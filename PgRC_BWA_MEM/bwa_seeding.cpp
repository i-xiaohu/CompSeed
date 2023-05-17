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

void mem_collect_intv(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a) {
	int i, k, x = 0, old_n;
	int start_width = 1;
	int split_len = (int)(opt->min_seed_len * opt->split_factor + .499);
	a->mem.n = 0;
	// first pass: find all SMEMs
	while (x < len) {
		if (seq[x] < 4) {
			x = bwt_smem1(bwt, len, seq, x, start_width, &a->mem1, a->tmpv);
			for (i = 0; i < a->mem1.n; ++i) {
				bwtintv_t *p = &a->mem1.a[i];
				int slen = (uint32_t)p->info - (p->info>>32); // seed length
				if (slen >= opt->min_seed_len)
					kv_push(bwtintv_t, a->mem, *p);
			}
		} else ++x;
	}
	// second pass: find MEMs inside a long SMEM
	old_n = a->mem.n;
	for (k = 0; k < old_n; ++k) {
		bwtintv_t *p = &a->mem.a[k];
		int start = p->info>>32, end = (int32_t)p->info;
		if (end - start < split_len || p->x[2] > opt->split_width) continue;
		bwt_smem1(bwt, len, seq, (start + end)>>1, p->x[2]+1, &a->mem1, a->tmpv);
		for (i = 0; i < a->mem1.n; ++i)
			if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len) {
				kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
			}
	}
	// third pass: LAST-like
	if (opt->max_mem_intv > 0) {
		x = 0;
		while (x < len) {
			if (seq[x] < 4) {
				if (1) {
					bwtintv_t m;
					x = bwt_seed_strategy1(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m);
					if (m.x[2] > 0) kv_push(bwtintv_t, a->mem, m);
				} else { // for now, we never come to this block which is slower
					x = bwt_smem1a(bwt, len, seq, x, start_width, opt->max_mem_intv, &a->mem1, a->tmpv);
					for (i = 0; i < a->mem1.n; ++i)
						kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
				}
			} else ++x;
		}
	}
}

SST::SST(bwt_t *b) {
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
	auto &p = nodes[parent];
	if (p.children[base] == -1) {
		bwt_extend(bwt, &p.match, next, is_back);
		for (uint8_t c = 0; c < 4; c++) {
			SST_Node_t child;
			child.match = next[c];
			p.children[c] = (int)nodes.size();
			nodes.push_back(child);
		}
	}
	auto &c = nodes[p.children[base]];
	// If find an empty node, BWT query is required
	if (c.match.x[0] + c.match.x[1] + c.match.x[2] == 0) {
		bwt_extend(bwt, &p.match, next, is_back);
		c.match = next[base];
	}
	return p.children[base];
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

int SST::add_lep_child(int parent, uint8_t base, const uint64_t *x) {
	monitor.push_back("ACGT"[base]);
	auto &p = nodes[parent];
	if (p.children[base] == -1) {
		p.children[base] = (int)nodes.size();
		SST_Node_t child;
		child.match.x[0] = x[0];
		child.match.x[1] = x[1];
		child.match.x[2] = x[2];
		nodes.push_back(child);
	} else {
		auto &c = nodes[p.children[base]];
		if (c.match.x[0] + c.match.x[1] + c.match.x[2] == 0) {
			c.match.x[0] = x[0];
			c.match.x[1] = x[1];
			c.match.x[2] = x[2];
			assert(check_mem(bwt, monitor, c.match));
		} else {
			assert(check_mem(bwt, monitor, c.match));
			for (int i = 0; i < 3; i++) {
				assert(c.match.x[i] == x[i]);
			}
		}
	}
	return p.children[base];
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

int SST::add_empty_child(int parent, uint8_t base) {
	monitor.push_back("ACGT"[base]);
	auto &p = nodes[parent];
	if (p.children[base] == -1) {
		p.children[base] = (int)nodes.size();
		SST_Node_t child;
		child.match.x[0] = 0;
		child.match.x[1] = 0;
		child.match.x[2] = 0;
		nodes.push_back(child);
	} // else do nothing whatever the child node is empty or not
	return p.children[base];
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
	if (seq[pivot] > 3) return pivot + 1;
	bwtintv_t ik, next[4];
	int node_id = forward_sst->query_child(0, seq[pivot], true);
	ik = forward_sst->get_intv(node_id);
	ik.info = pivot + 1;
	aux.prev_intv.clear();
	int ret_pivot = len;
	for (int i = pivot + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = forward_sst->query_child(node_id, c, false);
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
	std::reverse(aux.prev_intv.begin(), aux.prev_intv.end());
	for (const auto &m : aux.prev_intv) {
		assert(check_mem(bwa_idx->bwt, seq, pivot, m.info, m));
	}

	// Collect node indexes for LEPs in backward SST
//	for (int i = 0; i < len; i++) fprintf(stderr, "%c", "ACGTN"[seq[i]]); fprintf(stderr, "\n");
	aux.prev_node.clear();
	for (const auto &p : aux.prev_intv) {
		int start = pivot, end = (int)p.info;
//		fprintf(stderr, "LEP [%d, %d] [%ld, %ld, %ld]\n", start, end, p.x[0], p.x[1], p.x[2]);
//		for (int i = start; i < end; i++) fprintf(stderr, "%c", "ACGT"[seq[i]]); fprintf(stderr, "\n");

		node_id = 0;
		backward_sst->monitor.clear();
		for (int i = (int)p.info - 1; i >= pivot + 1; i--) {
			node_id = backward_sst->add_empty_child(node_id, seq[i]);
		}
		node_id = backward_sst->add_lep_child(node_id, seq[pivot], p.x);
		aux.prev_node.push_back(node_id);
	}

	aux.mem.clear();
	for (int i = pivot - 1; i >= -1; i--) {
		int c = (i == -1) ?4 :seq[i];
		aux.curr_intv.clear(); aux.curr_node.clear();
		for (int j = 0; j < aux.prev_node.size(); j++) {
			node_id = aux.prev_node[j];
			ik = aux.prev_intv[j];
			if (c < 4) {
				node_id = backward_sst->query_child(node_id, c, true);
				next[c] = backward_sst->get_intv(node_id);
			}
			if (c > 3 or next[c].x[2] < min_hits) {
				if (not aux.curr_intv.empty()) {
					fprintf(stderr, "Left extend from [%d, %d] [%ld, %ld, %ld] to [%d, %d] [%ld, %ld, %ld]\n",
						 i + 1, (int) ik.info, ik.x[0], ik.x[1], ik.x[2],
						 i, (int) ik.info, next[c].x[0], next[c].x[1], next[c].x[2]);
					for (const auto &a : aux.curr_intv) {
						fprintf(stderr, "[%d, %d] [%ld, %ld, %ld]\n", i, (int)a.info, a.x[0], a.x[1], a.x[2]);
					}
				}
				assert(aux.curr_intv.empty() and aux.curr_node.empty());
				if (aux.mem.empty() or i + 1 < aux.mem.back().info >> 32) {
					ik.info |= (1UL * (i + 1) << 32);
					aux.mem.push_back(ik);
				}
			} else {
				if (aux.curr_intv.empty() or next[c].x[2] != aux.curr_intv.back().x[2]) {
					next[c].info = ik.info;
					aux.curr_intv.push_back(next[c]);
					aux.curr_node.push_back(node_id);
				}
			}
		}
		if (aux.curr_intv.empty()) break;
		std::swap(aux.prev_intv, aux.curr_intv);
		std::swap(aux.prev_node, aux.curr_node);
	}

	for (const auto &p : aux.mem) {
		int start = (int)(p.info >> 32), end = (int)p.info;
//		fprintf(stderr, "[%d, %d) [%ld, %ld, %ld]\n", start, end, p.x[0], p.x[1], p.x[2]);
		assert(check_mem(bwa_idx->bwt, seq, start, end, p));
	}
	return ret_pivot;
}

void BWA_seeding::test_a_batch(const std::vector<long> &offset, std::vector<std::string> &batch) {
	forward_sst->clear(); backward_sst->clear();
	for (int i = 0; i < batch.size(); i++) {
		auto &read = batch[i];
		for (int j = 0; j < read.length(); j++) {
			read[j] = nst_nt4_table[(int) read[j]];
		}
		std::vector<bwtintv_t> smems;
		for (int j = 0; j < read.length(); ) {
			j = collect_smem_with_sst((const uint8_t*)read.c_str(), read.length(), j, 1, thr_aux);
			for (const auto &m : thr_aux.mem) {
				if ((int)m.info - (m.info >> 32) >= mem_opt->min_seed_len) {
					smems.push_back(m);
				}
			}
		}

		int old_n = (int)smems.size();
		for (int j = 0; j < old_n; j++) {
			const auto &m = smems[j] ;
			int start = m.info >> 32, end = (int)m.info;
			if (end - start < (int)(mem_opt->min_seed_len * mem_opt->split_factor + .499)
				or m.x[2] > mem_opt->split_width) continue;
			collect_smem_with_sst((const uint8_t*)read.c_str(), read.length(), (end + start) / 2, m.x[2] + 1, thr_aux);
			for (const auto &rem : thr_aux.mem) {
				if ((int)rem.info - (rem.info >> 32) >= mem_opt->min_seed_len) {
					smems.push_back(rem);
				}
			}
		}

		if (mem_opt->max_mem_intv > 0) {
			for (int j = 0; j < read.length(); ) {
				if (read[j] < 4) {
					bwtintv_t m;
					j = bwt_seed_strategy1(bwa_idx->bwt, read.length(), (const uint8_t*) read.c_str(), j, mem_opt->min_seed_len, mem_opt->max_mem_intv, &m);
					if (m.x[2] > 0) smems.push_back(m);
				} else {
					j++;
				}
			}
		}
	}
	fprintf(stderr, "Batch %d done\n", ++batch_id);
}

void BWA_seeding::seeding_SE() {
	std::vector<std::string> read_batch;
	std::vector<long> offset;

	// Decompressing HQ Reads
	char *buffer = new char[read_length * 2]; memset(buffer, 0, read_length * 2);
	long curr_position = 0, curr_mis_cnt = 0;
	forward_sst = new SST(bwa_idx->bwt);
	backward_sst = new SST(bwa_idx->bwt);
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
		read_batch.emplace_back(std::string(buffer));
		offset.push_back(curr_position);
		if (read_batch.size() >= COMP_BATCH_SIZE or i == hq_reads_list->reads_count - 1) {
			test_a_batch(offset, read_batch);
//			exit(EXIT_SUCCESS);
			read_batch.clear();
			offset.clear();
		}
	}
	delete forward_sst;
	delete backward_sst;

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
	mem_aux = smem_aux_init();


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

int main(int argc, char *argv[]) {
	if (argc == 1) {
		fprintf(stderr, "Usage: PBM <bwa_index> <pgrc_archieve>\n");
		return 1;
	}
	BWA_seeding worker;
	worker.set_index_name(argv[1]);
	worker.set_archive_name(argv[2]);
	worker.compressive_seeding();
}
