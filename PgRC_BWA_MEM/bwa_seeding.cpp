//
// Created by ixiaohu on 2023/5/6.
//

#include "bwa_seeding.h"

#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>

#include "../Matching/SimplePgMatcher.h"
#include "../Utils/helper.h"
#include "../FM_index/bntseq.h"
#include "../bwalib/bwa.h"
#include "../cstl/kvec.h"
#include "../cstl/kthread.h"
#include "../bwalib/utils.h"

void BWA_seeding::load_all_PGs(std::ifstream &in) {
	fprintf(stderr, "1. Decompressing HQ reads...\n");
	PseudoGenomeHeader hq_pgh(in);
	auto *hq_prop = new ReadsSetProperties(in);
	assert(!confirm_text_read_mode(in));
	hq_reads_list = load_extend_reads_list(
			in,
			hq_pgh.get_max_read_length(),
			hq_pgh.get_reads_count(),
			preserve_order_mode,
			false,
			false);
	hq_reads_list->reads_count = hq_prop->reads_count;
	hq_reads_list->read_length = hq_prop->max_read_length;

	fprintf(stderr, "2. Decompressing LQ reads...\n");
	PseudoGenomeHeader lq_pgh(in);
	auto *lq_prop = new ReadsSetProperties(in);
	assert(!confirm_text_read_mode(in));
	lq_reads_list = load_extend_reads_list(
			in,
			lq_pgh.get_max_read_length(),
			lq_pgh.get_reads_count(),
			preserve_order_mode,
			true,
			true);
	lq_reads_list->reads_count = lq_prop->reads_count;
	lq_reads_list->read_length = lq_prop->max_read_length;

	fprintf(stderr, "3. Decompressing N reads...\n");
	PseudoGenomeHeader n_pgh;
	ReadsSetProperties *n_prop = nullptr;
	if (separate_N) {
		n_pgh = PseudoGenomeHeader(in);
		n_prop = new ReadsSetProperties(in);
		assert(!confirm_text_read_mode(in));
		n_reads_list = load_extend_reads_list(
				in,
				n_pgh.get_max_read_length(),
				n_pgh.get_reads_count(),
				preserve_order_mode,
				true,
				true);
		n_reads_list->reads_count = n_prop->reads_count;
		n_reads_list->read_length = n_prop->max_read_length;
	}

	read_length = hq_prop->max_read_length;
	hq_reads_count = hq_prop->reads_count;
	lq_reads_count = lq_prop->reads_count;
	non_reads_count = hq_reads_count + lq_reads_count;
	n_reads_count = n_prop ?n_prop->reads_count :0;
	total_reads_count = non_reads_count + n_reads_count;

	hq_pg_length = hq_pgh.get_pg_length();
	non_pg_length = hq_pg_length + lq_pgh.get_pg_length();
	if (preserve_order_mode) {
		joined_pg_len_std = non_pg_length + n_pgh.get_pg_length() <= UINT32_MAX;
		if (joined_pg_len_std) {
			decompress_pg_position(in, pos_pg_32, total_reads_count, single_end_mode);
		} else {
			decompress_pg_position(in, pos_pg_64, total_reads_count, single_end_mode);
		}
	} else if (not single_end_mode) {
		restore_paired_idx(in, paired_idx);
	}

	restore_all_matched_pg(in, hq_pg_length, hq_pg, lq_pg, n_pg);
}

template<typename uint_pg_len>
void BWA_seeding::apply_rc_pair_to_pg(std::vector<uint_pg_len> &pg_pos) {
	if (preserve_order_mode) {
		int hq_idx = 0;
		const int pairs_count = total_reads_count / 2;
		for (int i = 0; i < pairs_count; i++) {
			uint_pg_len pos = pg_pos[i];
			if (pos < hq_pg_length) hq_idx++;
		}
		for (int i = pairs_count; i < total_reads_count; i++) {
			uint_pg_len pos = pg_pos[i];
			if (pos < hq_pg_length) {
				hq_reads_list->rev_comp[hq_idx] = !hq_reads_list->rev_comp[hq_idx];
				hq_idx++;
			}
		}
	} else {
		// Reverse complement each read from file2
		for (int i = 1; i < total_reads_count; i += 2) {
			int idx = paired_idx[i];
			if (idx < hq_reads_count) {
				hq_reads_list->rev_comp[idx] = !hq_reads_list->rev_comp[idx];
			}
		}
	}
}

mem_opt_t *mem_opt_init() {
	auto *o = (mem_opt_t*) calloc(1, sizeof(mem_opt_t));
	o->flag = 0;
	o->a = 1; o->b = 4;
	o->o_del = o->o_ins = 6;
	o->e_del = o->e_ins = 1;
	o->w = 100;
	o->T = 30;
	o->zdrop = 100;
	o->pen_unpaired = 17;
	o->pen_clip5 = o->pen_clip3 = 5;

	o->max_mem_intv = 20;

	o->min_seed_len = 19;
	o->split_width = 10;
	o->max_occ = 500;
	o->max_chain_gap = 10000;
	o->max_ins = 10000;
	o->mask_level = 0.50;
	o->drop_ratio = 0.50;
	o->XA_drop_ratio = 0.80;
	o->split_factor = 1.5;
	o->chunk_size = 10000000;
	o->n_threads = 1;
	o->max_XA_hits = 5;
	o->max_XA_hits_alt = 200;
	o->max_matesw = 50;
	o->mask_level_redun = 0.95;
	o->min_chain_weight = 0;
	o->max_chain_extend = 1<<30;
	o->mapQ_coef_len = 50; o->mapQ_coef_fac = log(o->mapQ_coef_len);
	bwa_fill_scmat(o->a, o->b, o->mat);
	return o;
}

smem_aux_t *smem_aux_init() {
	smem_aux_t *a;
	a = (smem_aux_t*) calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = (bwtintv_v*) calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = (bwtintv_v*) calloc(1, sizeof(bwtintv_v));
	return a;
}

void smem_aux_destroy(smem_aux_t *a) {
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

static void bwt_reverse_intvs(bwtintv_v *p) {
	if (p->n > 1) {
		int j;
		for (j = 0; j < p->n>>1; ++j) {
			bwtintv_t tmp = p->a[p->n - 1 - j];
			p->a[p->n - 1 - j] = p->a[j];
			p->a[j] = tmp;
		}
	}
}

int BWA_seeding::super_mem1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2]) {
	int i, j, c, ret;
	bwtintv_t ik, ok[4];
	bwtintv_v a[2], *prev, *curr, *swap;

	mem->n = 0;
	if (q[x] > 3) return x + 1;
	if (min_intv < 1) min_intv = 1; // the interval size should be at least 1
	kv_init(a[0]); kv_init(a[1]);
	prev = tmpvec && tmpvec[0]? tmpvec[0] : &a[0]; // use the temporary vector if provided
	curr = tmpvec && tmpvec[1]? tmpvec[1] : &a[1];
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	ik.info = x + 1;

	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (ik.x[2] < max_intv) { // an interval small enough
			kv_push(bwtintv_t, *curr, ik);
			break;
		} else if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			uint64_t t_start = __rdtsc();
			bwt_extend(bwt, &ik, ok, 0);
			bwa_bwt_ticks += __rdtsc() - t_start;
			global_bwa_bwt++;
			if (ok[c].x[2] != ik.x[2]) { // change of the interval size
				kv_push(bwtintv_t, *curr, ik);
				if (ok[c].x[2] < min_intv) break; // the interval size is too small to be extended further
			}
			ik = ok[c]; ik.info = i + 1;
		} else { // an ambiguous base
			kv_push(bwtintv_t, *curr, ik);
			break; // always terminate extension at an ambiguous base; in this case, i<len always stands
		}
	}
	if (i == len) kv_push(bwtintv_t, *curr, ik); // push the last interval if we reach the end
	bwt_reverse_intvs(curr); // s.t. smaller intervals (i.e. longer matches) visited first
	ret = curr->a[0].info; // this will be the returned value
	swap = curr; curr = prev; prev = swap;

	for (i = x - 1; i >= -1; --i) { // backward search for MEMs
		c = i < 0? -1 : q[i] < 4? q[i] : -1; // c==-1 if i<0 or q[i] is an ambiguous base
		for (j = 0, curr->n = 0; j < prev->n; ++j) {
			bwtintv_t *p = &prev->a[j];
			if (c >= 0 && ik.x[2] >= max_intv) {
				uint64_t t_start = __rdtsc();
				bwt_extend(bwt, p, ok, 1);
				bwa_bwt_ticks += __rdtsc() - t_start;
				global_bwa_bwt++;
			}
			if (c < 0 || ik.x[2] < max_intv || ok[c].x[2] < min_intv) { // keep the hit if reaching the beginning or an ambiguous base or the intv is small enough
				if (curr->n == 0) { // test curr->n>0 to make sure there are no longer matches
					if (mem->n == 0 || i + 1 < mem->a[mem->n-1].info>>32) { // skip contained matches
						ik = *p; ik.info |= (uint64_t)(i + 1)<<32;
						kv_push(bwtintv_t, *mem, ik);
					}
				} // otherwise the match is contained in another longer match
			} else if (curr->n == 0 || ok[c].x[2] != curr->a[curr->n-1].x[2]) {
				ok[c].info = p->info;
				kv_push(bwtintv_t, *curr, ok[c]);
			}
		}
		if (curr->n == 0) break;
		swap = curr; curr = prev; prev = swap;
	}

	bwt_reverse_intvs(mem); // s.t. sorted by the start coordinate

	if (tmpvec == 0 || tmpvec[0] == 0) free(a[0].a);
	if (tmpvec == 0 || tmpvec[1] == 0) free(a[1].a);
	return ret;
}

int BWA_seeding::seed_strategy1(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem) {
	int i, c;
	bwtintv_t ik, ok[4];

	memset(mem, 0, sizeof(bwtintv_t));
	if (q[x] > 3) return x + 1;
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	for (i = x + 1; i < len; ++i) { // forward search
		if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			uint64_t t_start = __rdtsc();
			bwt_extend(bwt, &ik, ok, 0);
			bwa_bwt_ticks += __rdtsc() - t_start;
			global_bwa_bwt++;
			if (ok[c].x[2] < max_intv && i - x >= min_len) {
				*mem = ok[c];
				mem->info = (uint64_t)x<<32 | (i + 1);
				return i + 1;
			}
			ik = ok[c];
		} else return i + 1;
	}
	return len;
}

void BWA_seeding::mem_collect_intv(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a) {
	int i, k, x = 0, old_n;
	int start_width = 1;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	a->mem.n = 0;
	// first pass: find all SMEMs
	long bwt_start = global_bwa_bwt;
	uint64_t stamp = __rdtsc();
	while (x < len) {
		if (seq[x] < 4) {
			x = super_mem1(bwt, len, seq, x, start_width, 0, &a->mem1, a->tmpv);
			for (i = 0; i < a->mem1.n; ++i) {
				bwtintv_t *p = &a->mem1.a[i];
				int slen = (uint32_t)p->info - (p->info>>32); // seed length
				if (slen >= opt->min_seed_len)
					kv_push(bwtintv_t, a->mem, *p);
			}
		} else ++x;
	}
	bwa_time.seeding += __rdtsc() - stamp;
	bwa_bwt_calls[1] += global_bwa_bwt - bwt_start;

	// second pass: find MEMs inside a long SMEM
	bwt_start = global_bwa_bwt;
	stamp = __rdtsc();
	old_n = a->mem.n;
	for (k = 0; k < old_n; ++k) {
		bwtintv_t *p = &a->mem.a[k];
		int start = p->info>>32, end = (int32_t)p->info;
		if (end - start < split_len || p->x[2] > opt->split_width) continue;
		super_mem1(bwt, len, seq, (start + end)>>1, p->x[2]+1, 0, &a->mem1, a->tmpv);
		for (i = 0; i < a->mem1.n; ++i)
			if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len) {
				kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
			}
	}
	bwa_time.reseed += __rdtsc() - stamp;
	bwa_bwt_calls[2] += global_bwa_bwt - bwt_start;

	// third pass: LAST-like
	bwt_start = global_bwa_bwt;
	stamp = __rdtsc();
	if (opt->max_mem_intv > 0) {
		x = 0;
		while (x < len) {
			if (seq[x] < 4) {
				if (1) {
					bwtintv_t m;
					x = seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
					if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
				} else { // for now, we never come to this block which is slower
					x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
					for (i = 0; i < a->mem1.n; ++i)
						kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
				}
			} else ++x;
		}
	}
	bwa_time.third += __rdtsc() - stamp;
	bwa_bwt_calls[3] += global_bwa_bwt - bwt_start;
}

SST::SST(const bwt_t *b) {
	bwt = b;
	SST_Node_t root;
	for (int i = 0; i < 4; i++) root.children[i] = i + 1;
	nodes.push_back(root);

	SST_Node_t child;
	for (uint8_t c = 0; c < 4; c++) {
		bwt_set_intv(bwt, c, child.match);
		nodes.push_back(child);
	}
}

int SST::query_child(int parent, uint8_t base, bool is_back) {
	if (nodes[parent].children[base] == -1) {
		uint64_t start = __rdtsc();
		bwt_extend(bwt, &nodes[parent].match, next, is_back);
		bwt_ticks += __rdtsc() - start;
		bwt_calls++;

		for (uint8_t c = 0; c < 4; c++) { // Is it necessary?
			SST_Node_t child;
			child.match = next[c];
			nodes[parent].children[c] = nodes.size();
			nodes.push_back(child);
		}
	}
	auto &c = nodes[nodes[parent].children[base]];
	// If find an empty node, BWT query is required
	if (c.match.x[0] + c.match.x[1] + c.match.x[2] == 0) {
		uint64_t start = __rdtsc();
		bwt_extend(bwt, &nodes[parent].match, next, is_back);
		bwt_ticks += __rdtsc() - start;
		bwt_calls++;
		c.match = next[base];
	}
	return nodes[parent].children[base];
}

static bool check_mem(const bwt_t *bwt, const std::string &s, const bwtintv_t &ref) {
	bwtintv_t ik, ok[4];
	bwt_set_intv(bwt, nst_nt4_table[s[0]], ik);
	for (int i = 1; i < s.size(); i++) {
		int c = nst_nt4_table[s[i]];
		bwt_extend(bwt, &ik, ok, 1);
		ik = ok[c];
	}
	bool correct = ref.x[0] == ik.x[0] and ref.x[1] == ik.x[1] and ref.x[2] == ik.x[2];
	if (not correct) {
		fprintf(stderr, "%s:[%ld, %ld, %ld] <-> [%ld, %ld, %ld]\n",
		  s.c_str(),
		  ik.x[0], ik.x[1], ik.x[2],
		  ref.x[0], ref.x[1], ref.x[2]);
	}
	return correct;
}

static bool check_mem(const bwt_t *bwt, const uint8_t *seq, int start, int end, const bwtintv_t &ref) {
	bwtintv_t ik, ok[4];
	bwt_set_intv(bwt, seq[start], ik);
	for (int i = start + 1; i < end; i++) {
		bwt_extend(bwt, &ik, ok, 0);
		ik = ok[3 - seq[i]];
	}
	return ref.x[0] == ik.x[0] and ref.x[1] == ik.x[1] and ref.x[2] == ik.x[2];
}

static void dfs_print_tree(SST *tree, int node, std::string &prefix) {
	for (uint8_t c = 0; c < 4; c++) {
		int next = tree->get_child(node, c);
		if (next == -1) continue;
		prefix.push_back("ACGT"[c]);
		auto intv = tree->get_intv(next);
		if (intv.x[0] + intv.x[1] + intv.x[2] > 0) {
			fprintf(stderr, "%s\tSA:[%ld, %ld, %ld]\n",
		        prefix.c_str(), intv.x[0], intv.x[1], intv.x[2]);
		}
		dfs_print_tree(tree, next, prefix);
		prefix.pop_back();
	}
}

int BWA_seeding::collect_smem_with_sst(const uint8_t *seq, int len, int pivot, int min_hits, thread_aux_t &aux) {
	aux.mem.clear();
	if (seq[pivot] > 3) return pivot + 1;
	bwtintv_t ik, next[4];
	int node_id = forward_sst->query_forward_child(0, seq[pivot]);
	ik = forward_sst->get_intv(node_id);
	ik.info = pivot + 1;
	aux.prev_intv.clear();
	int ret_pivot = len;
	for (int i = pivot + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = forward_sst->query_forward_child(node_id, c);
			next[c] = forward_sst->get_intv(node_id);
			if (next[c].x[2] != ik.x[2]) {
				aux.prev_intv.push_back(ik);
				if (next[c].x[2] < min_hits) {
					ret_pivot = i;
					break;
				}
			}
			ik = next[c];
			ik.info = i + 1;
		} else {
			aux.prev_intv.push_back(ik);
			ret_pivot = i + 1;
			break;
		}
	}
	if (ret_pivot == len) aux.prev_intv.push_back(ik);
	if (pivot == 0) { // Special judge
		aux.mem.push_back(aux.prev_intv.back());
		return ret_pivot;
	}
	std::reverse(aux.prev_intv.begin(), aux.prev_intv.end());

	// Collect node indexes for LEPs in backward SST
	for (auto &p : aux.prev_intv) {
		node_id = 0;
		for (int j = (int)p.info - 1; j >= pivot + 1; j--) {
			node_id = backward_sst->add_empty_child(node_id, seq[j]);
		}
		node_id = backward_sst->add_lep_child(node_id, seq[pivot], p.x);
		p.info |= (1UL * node_id << 32);
	}

	for (int i = pivot - 1; i >= -1; i--) {
		int c = (i == -1) ?4 :seq[i];
		aux.curr_intv.clear();
		for (int j = 0; j < aux.prev_intv.size(); j++) {
			ik = aux.prev_intv[j];
			node_id = ik.info >> 32;
			if (c < 4) {
				node_id = backward_sst->query_backward_child(node_id, c);
				next[c] = backward_sst->get_intv(node_id);
			}
			if (c > 3 or next[c].x[2] < min_hits) {
				if (aux.mem.empty() or i + 1 < aux.mem.back().info >> 32) {
					ik.info = (1UL * (i + 1) << 32) | (int32_t)ik.info;
					aux.mem.push_back(ik);
				}
			} else if (aux.curr_intv.empty() or next[c].x[2] != aux.curr_intv.back().x[2]) {
				next[c].info = (1UL * node_id << 32) | (int32_t)ik.info;
				aux.curr_intv.push_back(next[c]);
			}
		}
		if (aux.curr_intv.empty()) break;
		std::swap(aux.prev_intv, aux.curr_intv);
	}
	return ret_pivot;
}

int BWA_seeding::tem_forward_sst(const uint8_t *seq, int len, int start, int min_len, int max_intv, bwtintv_t *mem) {
	if (seq[start] > 3) return start + 1;
	memset(mem, 0, sizeof(bwtintv_t));
	int node_id = forward_sst->query_forward_child(0, seq[start]);
	bwtintv_t ik = forward_sst->get_intv(node_id);
	for (int i = start + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = forward_sst->query_forward_child(node_id, c);
			ik = forward_sst->get_intv(node_id);
			if (ik.x[2] < max_intv && i - start >= min_len) {
				*mem = ik;
				mem->info = (1UL * start << 32) | (i + 1);
				return i + 1;
			}
		} else return i + 1;
	}
	return len;
}

static inline bool mem_cmp(const bwtintv_t &a, const bwtintv_t &b) {
	return a.info < b.info;
}

static inline bool mem_merge_cmp(const bwtintv_t &a, const bwtintv_t &b) {
	return a.x[0] != b.x[0] ?a.x[0] < b.x[0] :a.x[2] < b.x[2];
}

static inline bool mem_eq(const bwtintv_t &a, const bwtintv_t &b) {
//	return a.info == b.info and a.x[0] == b.x[0] and
//		a.x[1] == b.x[1] and a.x[2] == b.x[2];
	// I have manually modified the variable x[1]
	return a.info == b.info and a.x[0] == b.x[0] and a.x[2] == b.x[2];
}

static inline int mem_beg(const bwtintv_t &a) {
	return a.info >> 32;
}

static inline int mem_end(const bwtintv_t &a) {
	return (int)a.info;
}

static inline int mem_len(const bwtintv_t &a) {
	return mem_end(a) - mem_beg(a);
}

static std::string mem_str(const bwtintv_t &a) {
	char buf[256];
	sprintf(buf, "[%d, %d) [%ld, %ld, %ld]", mem_beg(a), mem_end(a), a.x[0], a.x[1], a.x[2]);
	return std::string(buf);
}

static inline bool merge_cmp0(const bwtintv_t &a, const bwtintv_t &b) {
	return a.x[0] != b.x[0] ?a.x[0] < b.x[0] :a.x[2] < b.x[2];
}

static inline bool merge_cmp1(const bwtintv_t &a, const bwtintv_t &b) {
	return a.x[1] != b.x[1] ?a.x[1] < b.x[1] :a.x[2] < b.x[2];
}

void BWA_seeding::test_a_batch(const std::vector<long> &offset, std::vector<std::string> &batch) {
	uint64_t tick_start = __rdtsc();
	forward_sst->clear(); backward_sst->clear();
	for (int i = batch.size() - 1; i >= 0; i--) {
		auto &read = batch[i];
		for (int j = 0; j < read.length(); j++) {
			read[j] = nst_nt4_table[(int) read[j]];
		}

		long bwt_beg, bwt_end; uint64_t stamp;
		bwt_beg = forward_sst->bwt_calls + backward_sst->bwt_calls;
		stamp = __rdtsc();
		std::vector<bwtintv_t> &smems = batch_mem[i]; smems.clear();
		for (int j = 0; j < read.length(); ) {
			j = collect_smem_with_sst((const uint8_t*)read.c_str(), read.length(), j, 1, thr_aux);
			for (const auto &m : thr_aux.mem) {
				if ((int)m.info - (m.info >> 32) >= mem_opt->min_seed_len) {
					smems.push_back(m);
				}
			}
		}
		comp_time.seeding += __rdtsc() - stamp;
		bwt_end = forward_sst->bwt_calls + backward_sst->bwt_calls;
		comp_bwt_calls[1] += bwt_end - bwt_beg;

		bwt_beg = forward_sst->bwt_calls + backward_sst->bwt_calls;
		stamp = __rdtsc();
		int old_n = (int)smems.size();
		for (int j = 0; j < old_n; j++) {
			const auto &p = smems[j] ;
			int start = p.info >> 32, end = (int)p.info;
			if (end - start < (int)(mem_opt->min_seed_len * mem_opt->split_factor + .499)
			    or p.x[2] > mem_opt->split_width) continue;
			collect_smem_with_sst((const uint8_t*) read.c_str(), read.length(), (start + end) / 2, p.x[2] + 1, thr_aux);
			for (const auto &m : thr_aux.mem) {
				if ((int)m.info - (m.info >> 32) >= mem_opt->min_seed_len) {
					smems.push_back(m);
				}
			}
		}
		comp_time.reseed += __rdtsc() - stamp;
		bwt_end = forward_sst->bwt_calls + backward_sst->bwt_calls;
		comp_bwt_calls[2] += bwt_end - bwt_beg;

		bwt_beg = forward_sst->bwt_calls + backward_sst->bwt_calls;
		stamp = __rdtsc();
		if (mem_opt->max_mem_intv > 0) {
			for (int j = 0; j < read.length(); ) {
				if (read[j] < 4) {
					bwtintv_t m;
					j = tem_forward_sst((const uint8_t*) read.c_str(), read.length(), j, mem_opt->min_seed_len, mem_opt->max_mem_intv, &m);
					if (m.x[2] > 0) smems.push_back(m);
				} else {
					j++;
				}
			}
		}
		comp_time.third += __rdtsc() - stamp;
		bwt_end = forward_sst->bwt_calls + backward_sst->bwt_calls;
		comp_bwt_calls[3] += bwt_end - bwt_beg;
		std::sort(smems.begin(), smems.end(), mem_cmp);
	}

	std::vector<bwtintv_t> collected_mem;
	for (uint64_t i = 0; i < batch.size(); i++) {
		const auto &mem = batch_mem[i];
		for (const auto &m : mem) {
			collected_mem.push_back(m);
		}
	}
	// X0 repeat
	std::sort(collected_mem.begin(), collected_mem.end(), merge_cmp0);
	for (int i = 1; i < collected_mem.size(); i++) {
		const auto &p = collected_mem[i-1];
		const auto &m = collected_mem[i];
		if (p.x[0] != m.x[0] or p.x[2] != m.x[2]) {
			x0_unique++;
			x0_sal += std::min((int)p.x[2], mem_opt->max_occ);
		}
	}
	// X1 repeat
	std::sort(collected_mem.begin(), collected_mem.end(), merge_cmp1);
	for (int i = 1; i < collected_mem.size(); i++) {
		const auto &p = collected_mem[i-1];
		const auto &m = collected_mem[i];
		if (p.x[1] != m.x[1] or p.x[2] != m.x[2]) {
			x1_unique++;
			x1_sal += std::min((int)p.x[2], mem_opt->max_occ);
		}
	}
	total_mem += collected_mem.size();
	for (const auto &p : collected_mem) {
		total_sal += std::min((int)p.x[2], mem_opt->max_occ);
	}
	// O(N^2) loop
	for (int i = 0; i < collected_mem.size(); i++) {
		bool is_unique = true;
		const auto &c = collected_mem[i];
		for (int j = 0; j < i; j++) {
			const auto &p = collected_mem[j];
			if (p.x[0] == c.x[0] and p.x[2] == c.x[2]) { is_unique = false; break; }
			if (p.x[1] == c.x[1] and p.x[2] == c.x[2]) { is_unique = false; break; }
		}
		if (is_unique) {
			x_unique++;
			x_sal += std::min((int)c.x[2], mem_opt->max_occ);
		}
	}
	fprintf(stderr, "Sort0 MEM(%.2f) SAL(%.2f); Sort1 MEM(%.2f) SAL(%.2f); MAX MEM(%.2f) SAL(%.2f)\n",
	        100.0 * x0_unique / total_mem, 100.0 * x0_sal / total_sal,
	        100.0 * x1_unique / total_mem, 100.0 * x1_sal / total_sal,
	        100.0 * x_unique / total_mem, 100.0 * x_sal / total_sal);

	uint64_t sal_tick = __rdtsc();
	// Find unique MEMs by merging and sorting
	merge_mem.clear();
	for (uint64_t i = 0; i < batch.size(); i++) {
		auto &mem = batch_mem[i];
		for (uint64_t j = 0; j < mem.size(); j++) {
			mem[j].x[1] = i << 32U | j;
			merge_mem.push_back(mem[j]);
		}
	}
	std::sort(merge_mem.begin(), merge_mem.end(), mem_merge_cmp);
	// Lookup suffix array for unique MEMs
	int inc_pointer = 0; unique_seed.clear();
	for (int i = 0; i < merge_mem.size(); i++) {
		const auto *c = &merge_mem[i]; // Current MEM
		const auto *p = i == 0 ? nullptr : &merge_mem[i-1]; // Previous MEM
		if (p == nullptr or p->x[0] != c->x[0] or p->x[2] != c->x[2]) { // Another unique MEM
			inc_pointer = unique_seed.size();
			uint64_t step = c->x[2] > mem_opt->max_occ ? c->x[2] / mem_opt->max_occ : 1;
			int slen = (int)c->info - (c->info>>32), count = 0;
			for (uint64_t k = 0; k < c->x[2] && count < mem_opt->max_occ; k += step, count++) {
				mem_seed_t s;
				s.rbeg = bwt_sa(bwa_idx->bwt, c->x[0] + k); // the base coordinate in forward-reverse reference
				comp_sal++;
				unique_seed.push_back(s);
			}
		}
		int id1 = c->x[1] >> 32, id2 = (int)c->x[1];
		batch_mem[id1][id2].x[1] = inc_pointer;
	}
	// Distribute seeds back to each read
	for (int i = 0; i < batch.size(); i++) {
		const auto &mem = batch_mem[i];
		auto &seed = batch_seed[i]; seed.clear();
		for (const auto &m : mem) {
			int size = std::min((int)m.x[2], mem_opt->max_occ);
			int qbeg = mem_beg(m), slen = mem_len(m);
			for (int j = 0; j < size; j++) {
				mem_seed_t s = unique_seed[m.x[1] + j];
				s.qbeg = qbeg; s.score = s.len = slen;
				s.rid = bns_intv2rid(bwa_idx->bns, s.rbeg, s.rbeg + s.len);
				seed.push_back(s);
			}
		}
	}
	comp_time.sal += __rdtsc() - sal_tick;
	comp_time.total += __rdtsc() - tick_start;


	tick_start = __rdtsc();
	for (int i = 0; i < batch.size(); i++) {
		const auto &read = batch[i];
		mem_collect_intv(mem_opt, bwa_idx->bwt, read.length(), (const uint8_t*) read.c_str(), mem_aux);

		std::vector<bwtintv_t> &truth_set = truth_mem[i]; truth_set.clear();
		for (int j = 0; j < mem_aux->mem.n; j++) {
			truth_set.push_back(mem_aux->mem.a[j]);
		}
		std::sort(truth_set.begin(), truth_set.end(), mem_cmp);

		sal_tick = __rdtsc();
		auto &seed = truth_seed[i]; seed.clear();
		for (const auto &m : truth_set) {
			const bwtintv_t *p = &m;
			int step, count, slen = (uint32_t)p->info - (p->info>>32); // seed length
			int64_t k;
			// if (slen < opt->min_seed_len) continue; // ignore if too short or too repetitive
			step = p->x[2] > mem_opt->max_occ? p->x[2] / mem_opt->max_occ : 1;
			for (k = count = 0; k < p->x[2] && count < mem_opt->max_occ; k += step, ++count) {
				mem_seed_t s;
				s.rbeg = bwt_sa(bwa_idx->bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
				bwa_sal++;
				s.qbeg = p->info >> 32;
				s.score = s.len = slen;
				s.rid = bns_intv2rid(bwa_idx->bns, s.rbeg, s.rbeg + s.len);
				seed.push_back(s);
			}
		}
		bwa_time.sal += __rdtsc() - sal_tick;
	}
	bwa_time.total += __rdtsc() - tick_start;

	// Verify correctness
	for (int i = 0; i < batch.size(); i++) {
		const auto &mem = batch_mem[i];
		const auto &bwa_set = truth_mem[i];
		assert(mem.size() == bwa_set.size());
		for (int j = 0; j < mem.size(); j++) assert(mem_eq(mem[j], bwa_set[j]));

		const auto &seed = batch_seed[i];
		const auto &bwa_seed = truth_seed[i];
		assert(seed.size() == bwa_seed.size());
		for (int j = 0; j < seed.size(); j++) {
			const auto &s = seed[j], &b = bwa_seed[j];
			assert(s.rbeg == b.rbeg);
			assert(s.rid == b.rid);
			assert(s.qbeg == b.qbeg);
			assert(s.len == b.len);
		}

		bool perfect_match = false;
		for (auto &m : mem) {
			perfect_match |= (mem_len(m) == read_length) ;
		}
		full_read_match += perfect_match;
	}

	fprintf(stderr, "Batch %d pass BWA(%.3f) COMP(%.3f) Gain(%.2f)\n",
		 ++batch_id,
		 1.0 * bwa_time.total / cpu_frequency,
		 1.0 * comp_time.total / cpu_frequency,
		 1.0 * comp_time.total / bwa_time.total);
}

void BWA_seeding::seeding_SE() {
	std::vector<std::string> read_batch;
	std::vector<long> offset;

	// Decompressing HQ Reads
	char *buffer = new char[read_length * 2]; memset(buffer, 0, read_length * 2);
	long curr_position = 0, curr_mis_cnt = 0;
	forward_sst = new SST(bwa_idx->bwt);
	backward_sst = new SST(bwa_idx->bwt);
	mem_aux = smem_aux_init();

	for (int i = 0; i < hq_reads_list->reads_count; i++) {
		curr_position += hq_reads_list->off[i];
		memcpy(buffer, (hq_pg.data() + curr_position), hq_reads_list->read_length);
		if (hq_reads_list->rev_comp[i]) {
			rev_comp_in_place(buffer, hq_reads_list->read_length);
		}
		for (uint8_t j = 0; j < hq_reads_list->mis_cnt[i]; j++) {
			uint8_t mis_pos = hq_reads_list->mis_off[curr_mis_cnt];
			uint8_t mis_code = hq_reads_list->mis_sym_code[curr_mis_cnt];
			buffer[mis_pos] = code_to_mismatch(buffer[mis_pos], mis_code);
			curr_mis_cnt++;
		}

		// Reverse-complement it temporarily
		if (hq_reads_list->rev_comp[i]) {
			rev_comp_in_place(buffer, hq_reads_list->read_length);
		}
		// Output strand-corrected reads and shifted distance between consecutive reads
		read_batch.emplace_back(std::string(buffer));
		offset.push_back(curr_position);
		if (read_batch.size() >= COMP_BATCH_SIZE or i == hq_reads_list->reads_count - 1) {
			test_a_batch(offset, read_batch);
			read_batch.clear();
			offset.clear();
		}
		if (batch_id > 3000) break;
	}
	delete forward_sst;
	delete backward_sst;
	smem_aux_destroy(mem_aux);

	// Restoring LQ reads
	curr_position = 0;
	for (int i = 0; i < lq_reads_list->reads_count; i++) {
		curr_position += lq_reads_list->off[i];
		memcpy(buffer, (lq_pg.data() + curr_position), lq_reads_list->read_length);
	}

	// Restoring N reads
	curr_position = 0;
	for (int i = 0; i < n_reads_list->reads_count; i++) {
		curr_position += n_reads_list->off[i];
		memcpy(buffer, (n_pg.data() + curr_position), n_reads_list->read_length);
	}

	delete [] buffer;
}

void BWA_seeding::compressive_seeding() {
	// Loading FM-index
	bwa_idx = bwa_idx_load_from_shm(index_name.c_str());
	if (bwa_idx == nullptr) {
		bwa_idx = bwa_idx_load(index_name.c_str(), BWA_IDX_ALL);
		if (bwa_idx == nullptr) {
			fprintf(stderr, "Load the FM-index failed\n");
			exit(EXIT_FAILURE);
		} else {
			fprintf(stderr, "Load the FM-index from disk\n");
		}
	} else {
		fprintf(stderr, "Load the FM-index from shared memory\n");
	}
	mem_opt = mem_opt_init();

	std::ifstream in(archive_name); assert(in.is_open());
	for (int i = 0; i < strlen(PGRC_HEADER); i++) assert(in.get() == PGRC_HEADER[i]);
	char version_mode = in.get(); assert(version_mode == '#');
	char version_major = in.get();
	char version_minor = in.get();
	char version_revision = in.get();
	if (version_major != PGRC_VERSION_MAJOR or version_minor != PGRC_VERSION_MINOR) {
		fprintf(stderr, "Archive is packed with a different version PgRC %d.%d\n",
		        version_major, version_minor);
		exit(EXIT_FAILURE);
	}
	compression_level = in.get();
	char pgrc_mode = in.get();
	if (pgrc_mode != PGRC_SE_MODE) {
		fprintf(stderr, "PgRC-BWA-MEM now supports SE mode only\n");
		exit(EXIT_FAILURE);
	}
	const std::string five_modes[] = {"SE", "PE", "SE_ORD", "PE_ORD"};
	fprintf(stderr, "Compression mode: %s\n", five_modes[pgrc_mode].c_str());
	separate_N = (bool) in.get();
	if (pgrc_mode == PGRC_PE_MODE or pgrc_mode == PGRC_ORD_PE_MODE) {
		rev_comp_pair = (bool) in.get();
	}
	std::string tmp_directory; in >> tmp_directory; tmp_directory += "/";
	in.get();

	preserve_order_mode = (pgrc_mode == PGRC_ORD_SE_MODE or pgrc_mode == PGRC_ORD_PE_MODE);
	single_end_mode = (pgrc_mode == PGRC_SE_MODE or pgrc_mode == PGRC_ORD_SE_MODE);

	load_all_PGs(in);

	if (rev_comp_pair) {
		fprintf(stderr, "Reverse complement second read from a pair\n");
		if (joined_pg_len_std) apply_rc_pair_to_pg(pos_pg_32);
		else apply_rc_pair_to_pg(pos_pg_64);
	}

	if (pgrc_mode == PGRC_SE_MODE) seeding_SE();

	fprintf(stderr, "Decompressed %d reads in total\n", total_reads_count);
	in.close();
}

void BWA_seeding::on_dec_reads(const char *fn) {
	// Loading FM-index
	bwa_idx = bwa_idx_load_from_shm(index_name.c_str());
	if (bwa_idx == nullptr) {
		bwa_idx = bwa_idx_load(index_name.c_str(), BWA_IDX_ALL);
		if (bwa_idx == nullptr) {
			fprintf(stderr, "Load the FM-index failed\n");
			exit(EXIT_FAILURE);
		} else {
			fprintf(stderr, "Load the FM-index from disk\n");
		}
	} else {
		fprintf(stderr, "Load the FM-index from shared memory\n");
	}
	mem_opt = mem_opt_init();

	uint64_t stamp = __rdtsc(); sleep(1); cpu_frequency = __rdtsc() - stamp;

	fprintf(stderr, "Input test data with strand-corrected reads and overlapping information\n");
	std::ifstream in(fn); assert(in.is_open());
	std::vector<std::string> read_batch;
	std::vector<long> offset;
	char *buffer = new char[1024]; memset(buffer, 0, 1024 * sizeof(char));
	long off, curr_position = 0, curr_mis_cnt = 0;
	forward_sst = new SST(bwa_idx->bwt);
	backward_sst = new SST(bwa_idx->bwt);
	mem_aux = smem_aux_init();
	while (in >> buffer >> off) {
		curr_position += off;
		read_batch.emplace_back(std::string(buffer));
		offset.push_back(curr_position);
		read_length = (int)strlen(buffer);
		hq_reads_count++;
		total_reads_count++;
		if (read_batch.size() >= COMP_BATCH_SIZE or in.eof()) {
			test_a_batch(offset, read_batch);
			read_batch.clear();
			offset.clear();
		}
		if (batch_id > 3000) break;
	}

	fprintf(stderr, "Input %d reads in total\n", total_reads_count);
	fprintf(stderr, "Overall BWA\tCOMP\tGAIN\n");
	fprintf(stderr, "  Total %.2f\t%.2f\t%.2f\n", __time(bwa_time.total), __time(comp_time.total), 1.0 * comp_time.total / bwa_time.total);
	fprintf(stderr, "  S1    %.2f\t%.2f\t%.2f\n", __time(bwa_time.seeding), __time(comp_time.seeding), 1.0 * comp_time.seeding / bwa_time.seeding);
	fprintf(stderr, "  S2    %.2f\t%.2f\t%.2f\n", __time(bwa_time.reseed), __time(comp_time.reseed), 1.0 * comp_time.reseed / bwa_time.reseed);
	fprintf(stderr, "  S3    %.2f\t%.2f\t%.2f\n", __time(bwa_time.third), __time(comp_time.third), 1.0 * comp_time.third / bwa_time.third);
	fprintf(stderr, "  S4    %.2f\t%.2f\t%.2f\n", __time(bwa_time.sal), __time(comp_time.sal), 1.0 * comp_time.sal / bwa_time.sal);
	fprintf(stderr, "Perfect Matched Reads: %d (%.2f %%)\n", full_read_match, 100.0 * full_read_match / hq_reads_count);
	fprintf(stderr, "BWT Gain: Seeding(%.2f) Reseed(%.2f) Third(%.2f) BWT_SUM(%.2f) SAL(%.2f)\n",
		 1.0 * comp_bwt_calls[1] / bwa_bwt_calls[1],
		 1.0 * comp_bwt_calls[2] / bwa_bwt_calls[2],
		 1.0 * comp_bwt_calls[3] / bwa_bwt_calls[3],
		 1.0 * (comp_bwt_calls[1] + comp_bwt_calls[2] + comp_bwt_calls[3]) / (bwa_bwt_calls[1] + bwa_bwt_calls[2] + bwa_bwt_calls[3]),
		 1.0 * comp_sal / bwa_sal);
	fprintf(stderr, "Comp BWT percentage %.2f\n", 100.0 * (forward_sst->bwt_ticks + backward_sst->bwt_ticks) / comp_time.total);
	fprintf(stderr, "BWA  BWT percentage %.2f\n", 100.0 * bwa_bwt_ticks / bwa_time.total);

	free(mem_opt);
	in.close();
	delete [] buffer;
	delete forward_sst;
	delete backward_sst;
	smem_aux_destroy(mem_aux);
}

int main(int argc, char *argv[]) {
	if (argc == 1) {
		fprintf(stderr, "Usage: PBM <bwa_index> <pgrc_archieve>\n");
		return 1;
	}
	BWA_seeding worker;
	worker.set_index_name(argv[1]);
	if (std::string(argv[2]) == "test") {
		worker.on_dec_reads(argv[3]);
	} else {
		worker.set_archive_name(argv[2]);
		worker.compressive_seeding();
	}
}
