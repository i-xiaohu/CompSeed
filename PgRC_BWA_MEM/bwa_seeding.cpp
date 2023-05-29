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
	int node_id = aux.forward_sst->query_forward_child(0, seq[pivot]);
	ik = aux.forward_sst->get_intv(node_id);
	ik.info = pivot + 1;
	aux.prev_intv.clear();
	int ret_pivot = len;
	for (int i = pivot + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = aux.forward_sst->query_forward_child(node_id, c);
			next[c] = aux.forward_sst->get_intv(node_id);
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
			node_id = aux.backward_sst->add_empty_child(node_id, seq[j]);
		}
		node_id = aux.backward_sst->add_lep_child(node_id, seq[pivot], p.x);
		p.info |= (1UL * node_id << 32);
	}

	for (int i = pivot - 1; i >= -1; i--) {
		int c = (i == -1) ?4 :seq[i];
		aux.curr_intv.clear();
		for (int j = 0; j < aux.prev_intv.size(); j++) {
			ik = aux.prev_intv[j];
			node_id = ik.info >> 32;
			if (c < 4) {
				node_id = aux.backward_sst->query_backward_child(node_id, c);
				next[c] = aux.backward_sst->get_intv(node_id);
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

int BWA_seeding::tem_forward_sst(const uint8_t *seq, int len, int start, int min_len, int max_intv, bwtintv_t *mem, thread_aux_t &aux) {
	if (seq[start] > 3) return start + 1;
	memset(mem, 0, sizeof(bwtintv_t));
	int node_id = aux.forward_sst->query_forward_child(0, seq[start]);
	bwtintv_t ik = aux.forward_sst->get_intv(node_id);
	for (int i = start + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = aux.forward_sst->query_forward_child(node_id, c);
			ik = aux.forward_sst->get_intv(node_id);
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

void BWA_seeding::test_a_batch(const std::vector<long> &offset, std::vector<std::string> &batch, thread_aux_t &aux) {
	// Agency of thread auxiliary
	auto *mem_aux = aux.mem_aux;
	auto &truth_mem = aux.truth_mem;
	auto &truth_seed = aux.truth_seed;
	auto &bwa_time = aux.bwa_time;
	auto &comp_time = aux.comp_time;
	auto *forward_sst = aux.forward_sst;
	auto *backward_sst = aux.backward_sst;
	auto &batch_mem = aux.batch_mem;
	auto &batch_seed = aux.batch_seed;
	auto &unique_sal = aux.unique_sal;

	forward_sst->clear(); backward_sst->clear();
	long ref_position = -1;
	for (int i = batch.size() - 1; i >= 0; i--) {
		auto &read = batch[i];
		for (auto &c : read) c = c < 5 ?c :nst_nt4_table[c];

		long bwt_beg, bwt_end; uint64_t stamp;
		bwt_beg = forward_sst->bwt_calls + backward_sst->bwt_calls;
		std::vector<bwtintv_t> &smems = batch_mem[i]; smems.clear();
		bool full_match = false;
		if (ref_position != -1 and ref_position - offset[i] <= read.length() / 3) {
			// Lookup forward SST to test if this read can be fully matched
			int node_id = 0;
			for (int j = ref_position - offset[i]; j < read.length(); j++) {
				if (read[j] < 4) {
					if (j == ref_position - offset[i]) node_id = forward_sst->get_child(node_id, read[j]);
					else node_id = forward_sst->get_child(node_id, 3 - read[j]);
				} else node_id = -1;
				if (node_id == -1) break;
			}
			if (node_id != -1) {
				// matched in forward SST till the end of read
				bwtintv_t ik = forward_sst->get_intv(node_id), ok[4];
				// It is a potential full-match read
				long temp_bwt_calls = 0;
				for (int j = ref_position - offset[i] - 1; j >= 0; j--) {
					if (read[j] < 4) {
						bwt_extend(bwa_idx->bwt, &ik, ok, 1);
						temp_bwt_calls++;
						ik = ok[read[j]];
					} else ik.x[2] = 0;
					if (ik.x[2] == 0) break;
				}
				aux.comp_bwt_calls[1] += temp_bwt_calls;
				if (ik.x[2] > 0) {
					full_match = true;
					aux.shortcut++;
					ik.info = read.length();
					if (mem_len(ik) >= mem_opt->min_seed_len) smems.push_back(ik);
				} else {
					aux.comp_bwt_calls[4] += temp_bwt_calls;
				}
			}
		}

		if (not full_match) {
			for (int j = 0; j < read.length(); ) {
				j = collect_smem_with_sst((const uint8_t*)read.c_str(), read.length(), j, 1, aux);
				for (const auto &m : aux.mem) {
					if (mem_len(m) >= mem_opt->min_seed_len) {
						smems.push_back(m);
					}
				}
			}
			ref_position = offset[i];
		}
		bwt_end = forward_sst->bwt_calls + backward_sst->bwt_calls;
		aux.comp_bwt_calls[1] += bwt_end - bwt_beg;

		bwt_beg = forward_sst->bwt_calls + backward_sst->bwt_calls;
		int old_n = (int)smems.size();
		for (int j = 0; j < old_n; j++) {
			const auto &p = smems[j] ;
			int start = p.info >> 32, end = (int)p.info;
			if (end - start < (int)(mem_opt->min_seed_len * mem_opt->split_factor + .499)
			    or p.x[2] > mem_opt->split_width) continue;
			collect_smem_with_sst((const uint8_t*) read.c_str(), read.length(), (start + end) / 2, p.x[2] + 1, aux);
			for (const auto &m : aux.mem) {
				if ((int)m.info - (m.info >> 32) >= mem_opt->min_seed_len) {
					smems.push_back(m);
				}
			}
		}
		bwt_end = forward_sst->bwt_calls + backward_sst->bwt_calls;
		aux.comp_bwt_calls[2] += bwt_end - bwt_beg;

		bwt_beg = forward_sst->bwt_calls + backward_sst->bwt_calls;
		if (mem_opt->max_mem_intv > 0) {
			for (int j = 0; j < read.length(); ) {
				if (read[j] < 4) {
					bwtintv_t m;
					j = tem_forward_sst((const uint8_t*) read.c_str(), read.length(), j, mem_opt->min_seed_len, mem_opt->max_mem_intv, &m, aux);
					if (m.x[2] > 0) smems.push_back(m);
				} else {
					j++;
				}
			}
		}
		bwt_end = forward_sst->bwt_calls + backward_sst->bwt_calls;
		aux.comp_bwt_calls[3] += bwt_end - bwt_beg;
		std::sort(smems.begin(), smems.end(), mem_cmp);
	}

	// Find unique hit locations by merging and sorting
	unique_sal.clear();
	for (int read_id = 0; read_id < batch.size(); read_id++) {
		const auto &mem = batch_mem[read_id];
		auto &seed = batch_seed[read_id]; seed.clear();
		int array_id = 0;
		for (const auto &m : mem) {
			uint64_t step = m.x[2] > mem_opt->max_occ ? m.x[2] / mem_opt->max_occ : 1;
			for (uint64_t k = 0, count = 0; k < m.x[2] && count < mem_opt->max_occ; k += step, count++) {
				mem_seed_t s;
				s.qbeg = mem_beg(m);
				s.score = s.len = mem_len(m);
				seed.push_back(s);
				unique_sal.emplace_back(SAL_Packed(m.x[0] + k, read_id, array_id));
				array_id++;
			}
		}
	}
	std::sort(unique_sal.begin(), unique_sal.end());
	// Lookup suffix array for hit locations
	uint64_t coordinate = 0;
	for (int i = 0; i < unique_sal.size(); i++) {
		const auto &p = unique_sal[i];
		if (i == 0 or unique_sal[i-1].hit_location != p.hit_location) {
			coordinate = bwt_sa(bwa_idx->bwt, p.hit_location);
			aux.comp_sal++;
		}
		auto &s = batch_seed[p.read_id][p.array_id];
		s.rbeg = coordinate;
		s.rid = bns_intv2rid(bwa_idx->bns, s.rbeg, s.rbeg + s.len);
	}

//	for (int i = 0; i < batch.size(); i++) {
//		for (char c : batch[i]) fprintf(stdout, "%c", "ACGTN"[c]); fprintf(stdout, "\n");
//		for (const auto &s : batch_seed[i]) {
//			fprintf(stdout, "[%d:%ld %d_%d]\n", s.rid, s.rbeg, s.qbeg, s.len);
//		}
//	}
}

void BWA_seeding::display_profile() {
	thread_aux_t total;
	for (int i = 0; i < threads_n; i++) total += thr_aux[i];
	const auto &bwa_time = total.bwa_time;
	const auto &comp_time = total.comp_time;
	const auto *bwa_bwt_calls = total.bwa_bwt_calls;
	const auto *comp_bwt_calls = total.comp_bwt_calls;
	auto full_read_match = total.full_read_match;
	auto shortcut = total.shortcut;
	auto comp_sal = total.comp_sal, bwa_sal = total.bwa_sal;
	auto bwa_bwt_ticks = total.bwa_bwt_ticks;
	fprintf(stderr, "Input %d reads in total\n", total_reads_count);
	fprintf(stderr, "Perfect Matched Reads: %d (%.2f %%)\n", full_read_match, 100.0 * full_read_match / total_reads_count);
	fprintf(stderr, "Reads go shortcut:     %d (%.2f %%)\n", shortcut, 100.0 * shortcut / total_reads_count);
	fprintf(stderr, "BWT Gain: Seeding(%ld) Reseed(%ld) Third(%ld) BWT_SUM(%ld) SAL(%ld)\n",
	        comp_bwt_calls[1],
	        comp_bwt_calls[2],
	        comp_bwt_calls[3],
	        (comp_bwt_calls[1] + comp_bwt_calls[2] + comp_bwt_calls[3]),
	        comp_sal);
	fprintf(stderr, "Wasted Comp BWT:    %.2f\n", 100.0 * comp_bwt_calls[4] / comp_bwt_calls[1]);
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
	uint64_t stamp = __rdtsc(); sleep(1); cpu_frequency = __rdtsc() - stamp;

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

//	if (pgrc_mode == PGRC_SE_MODE) seeding_SE();

	fprintf(stderr, "Decompressed %d reads in total\n", total_reads_count);
	in.close();
}

void BWA_seeding::seeding_with_thread(int batch_id, int tid) {
	std::vector<std::string> batch;
	std::vector<long> offset;
	int start = batch_id * BATCH_SIZE, end = std::min((batch_id + 1) * BATCH_SIZE, (int)all_batch.size());
	assert(start < all_batch.size() and end <= all_batch.size());
	for (int i = start; i < end; i++) {
		batch.push_back(all_batch[i]);
		offset.push_back(all_offset[i]);
	}
//	fprintf(stdout, "Thread %d processed reads from %d to %d\n", tid, start, end);
//	for (auto & i : batch) {
//		fprintf(stdout, "%s\n", i.c_str());
//	}
	test_a_batch(offset, batch, thr_aux[tid]);
}

void thread_worker(void *data, long batch_id, int tid) {
	auto *owner = (BWA_seeding*)data;
	owner->seeding_with_thread(batch_id, tid);
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

	// Prepare for thread auxiliary
	thr_aux = new thread_aux_t[threads_n];
	for (int i = 0; i < threads_n; i++) {
		auto &a = thr_aux[i];
		a.forward_sst = new SST(bwa_idx->bwt);
		a.backward_sst = new SST(bwa_idx->bwt);
		a.mem_aux = smem_aux_init();
	}

	fprintf(stderr, "Input test data with strand-corrected reads and overlapping information\n");
	std::ifstream in(fn); assert(in.is_open());
	char *buffer = new char[1024]; memset(buffer, 0, 1024 * sizeof(char));
	long off, curr_position = 0, bytes = 0;
	while (in >> buffer >> off) {
		curr_position += off;
		all_batch.emplace_back(std::string(buffer));
		all_offset.push_back(curr_position);
		read_length = (int)strlen(buffer);
		bytes += read_length;
		hq_reads_count++;
		total_reads_count++;
		if (bytes >= 10 * 1024 * 1024) {
			kt_for(threads_n, thread_worker, this, std::max(1UL, all_batch.size() / BATCH_SIZE));
			fprintf(stderr, "%d Reads processed\n", total_reads_count);
			bytes = 0;
			all_batch.clear();
			all_offset.clear();
		}
	}
	display_profile();

	free(mem_opt);
	for (int i = 0; i < threads_n; i++) {
		auto &a = thr_aux[i];
		delete a.forward_sst;
		delete a.backward_sst;
		smem_aux_destroy(a.mem_aux);
	}
	delete [] thr_aux;
	in.close();
	delete [] buffer;
}

int main(int argc, char *argv[]) {
	if (argc == 1) {
		fprintf(stderr, "Usage: PBM <bwa_index> <pgrc_archieve>\n");
		return 1;
	}
	BWA_seeding worker;
	int c; bool input_test = false;
	while ((c = getopt(argc, argv, "t:1")) >= 0) {
		if (c == 't') {
			worker.set_threads(atoi(optarg));
		} else if (c == '1') {
			input_test = true;
		} else {
			exit(EXIT_FAILURE);
		}
	}
	worker.set_index_name(argv[optind]);
	if (input_test) {
		worker.on_dec_reads(argv[optind + 1]);
	} else {
		worker.set_archive_name(argv[optind + 1]);
		worker.compressive_seeding();
	}
}
