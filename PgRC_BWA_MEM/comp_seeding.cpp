//
// Created by ixiaohu on 2023/5/6.
//

#include "comp_seeding.h"

#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>

#include "../Matching/SimplePgMatcher.h"
#include "../Utils/helper.h"
#include "../FM_index/bntseq.h"
#include "../bwalib/bwa.h"
#include "../bwalib/ksw.h"
#include "../cstl/kvec.h"
#include "../cstl/kthread.h"
#include "../cstl/kbtree.h"
#include "../cstl/ksort.h"
#include "../cstl/kstring.h"
#include "../bwalib/utils.h"

/********************************
 * Bookmark 1: Seeding with SST *
 ********************************/
CompAligner::CompAligner() {
	uint64_t stamp = __rdtsc(); sleep(1); cpu_frequency = __rdtsc() - stamp;
	opt = mem_opt_init();
}

void CompAligner::load_index(const char *fn) {
	bwa_idx = bwa_idx_load_from_shm(fn);
	if (bwa_idx == nullptr) {
		bwa_idx = bwa_idx_load(fn, BWA_IDX_ALL);
		if (bwa_idx == nullptr) {
			fprintf(stderr, "Load the FM-index `%s` failed\n", fn);
			exit(EXIT_FAILURE);
		} else {
			fprintf(stderr, "Load the FM-index from disk\n");
		}
	} else {
		fprintf(stderr, "Load the FM-index from shared memory\n");
	}
	bns = bwa_idx->bns;
	bwt = bwa_idx->bwt;
	pac = bwa_idx->pac;
}

int CompAligner::collect_mem_with_sst(const uint8_t *seq, int len, int pivot, int min_hits, thread_aux &aux) {
	aux.super_mem.clear();
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
	if (pivot == 0) { // Quickly return full-length matched reads
		aux.super_mem.push_back(aux.prev_intv.back());
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
				if (aux.super_mem.empty() or i + 1 < aux.super_mem.back().info >> 32) {
					ik.info = (1UL * (i + 1) << 32) | (int32_t)ik.info;
					aux.super_mem.push_back(ik);
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

int CompAligner::tem_forward_sst(const uint8_t *seq, int len, int start, bwtintv_t *mem, thread_aux &aux) {
	if (seq[start] > 3) return start + 1;
	memset(mem, 0, sizeof(bwtintv_t));
	int node_id = aux.forward_sst->query_forward_child(0, seq[start]);
	bwtintv_t ik = aux.forward_sst->get_intv(node_id);
	for (int i = start + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = aux.forward_sst->query_forward_child(node_id, c);
			ik = aux.forward_sst->get_intv(node_id);
			if (ik.x[2] < opt->max_mem_intv && i - start >= opt->min_seed_len ) {
				*mem = ik;
				mem->info = (1UL * start << 32) | (i + 1);
				return i + 1;
			}
		} else return i + 1;
	}
	return len;
}

static inline bool mem_cmp(const bwtintv_t &a, const bwtintv_t &b) { return a.info < b.info; }

static inline int mem_beg(const bwtintv_t &a) { return a.info >> 32; }

static inline int mem_end(const bwtintv_t &a) { return (int)a.info; }

static inline int mem_len(const bwtintv_t &a) { return mem_end(a) - mem_beg(a); }

static std::string mem_str(const bwtintv_t &a) {
	char buf[256];
	sprintf(buf, "[%d, %d) [%ld, %ld, %ld]", mem_beg(a), mem_end(a), a.x[0], a.x[1], a.x[2]);
	return std::string(buf);
}

/***********************
 * Bookmark 2 Chaining *
 ***********************/

#define chain_cmp(a, b) (((b).anchor < (a).anchor) - ((a).anchor < (b).anchor))
KBTREE_INIT(chn, seed_chain, chain_cmp)

bool CompAligner::add_seed_to_chain(seed_chain *c, const seed_hit &s) {
	if (s.rid != c->rid) return false; // Seed is not on same chromosome with chain
	if (s.qbeg >= c->que_beg() and s.qbeg + s.len <= c->que_end()  and
		s.rbeg >= c->ref_beg() and s.rbeg + s.len <= c->ref_end()) {
		return true; // Seed in contained in chain
	}
	const auto &last = c->seeds[c->n - 1];
	if ((last.rbeg < bns->l_pac or c->ref_beg() < bns->l_pac) and s.rbeg >= bns->l_pac) {
		return false; // Seed is on different strand from chain
	}
	int64_t x = s.qbeg - last.qbeg; // Non-negative; guaranteed by sorting matches
	int64_t y = s.rbeg - last.rbeg;
	if (// Chained seeds should be co-linear on both read and reference
		y >= 0 and
		// Minimum insertion/deletions should be within DP matrix bandwidth
		std::abs(x - y) <= opt->w and
		// Gap between chain and seed on both read and reference is limited
		x - last.len < opt->max_chain_gap and y - last.len < opt->max_chain_gap
	) {
		c->push_back(s); // Grow the chain
		return true;
	} else {
		return false;
	}
}

void CompAligner::print_chains_to(const std::vector<seed_chain> &chains, kstring_t *s) {
	for (int i = 0; i < chains.size(); i++) {
		const auto &c = chains[i];
		ksprintf(s, "* Found CHAIN(%d): n=%d; weight=%d", i, c.n, c.calc_weight());
		for (int j = 0; j < c.n; j++) {
			bwtint_t pos; int is_rev;
			pos = bns_depos(bns, c.seeds[j].rbeg, &is_rev);
			if (is_rev) pos -= c.seeds[j].len - 1;
			ksprintf(s, "\t%d;%d;%d,%ld(%s:%c%ld)",
				c.seeds[j].score, c.seeds[j].len,
				c.seeds[j].qbeg, (long)c.seeds[j].rbeg,
				bns->anns[c.rid].name, "+-"[is_rev],
				(long)(pos - bns->anns[c.rid].offset) + 1);
		}
		ksprintf(s, "\n");
	}
}

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(flt, seed_chain, flt_lt)

std::vector<seed_chain> CompAligner::chaining(const std::vector<seed_hit> &seed) {
	kbtree_chn_t *tree = kb_init_chn(KB_DEFAULT_SIZE) ;
	for (const auto &s : seed) {
		if (s.rid < 0) continue;
		bool need_new_chain = true;
		seed_chain tmp(s.rbeg), *lower, *upper;
		if (kb_size(tree) > 0) {
			kb_intervalp_chn(tree, &tmp, &lower, &upper);
			if (lower) {
				// Check if the seed could merged into the closet chain
				need_new_chain = not add_seed_to_chain(lower, s);
			}
		}
		if (need_new_chain) {
			tmp.add_first_seed(s);
			tmp.is_alt = bns->anns[s.rid].is_alt;
			kb_putp_chn(tree, &tmp);
		}
	}

	std::vector<seed_chain> chain;
	chain.reserve(kb_size(tree));
	#define traverse_func(p) (chain.push_back(*(p)))
		__kb_traverse(seed_chain, tree, traverse_func);
	#undef traverse_func
	kb_destroy(chn, tree);

	if (chain.empty()) return chain;

	// Drop chains with the weight smaller than the specified minimum
	const int DOMINATE = 3; // Dominating chains with large weight
	const int SHADOWED = 2; // Shadowed by another chain
	const int DISCARD = 0; // Dropped chain
	int kept_chain_n = 0;
	for (auto &c : chain) {
		c.kept = DISCARD;
		c.first_cover = -1;
		c.w = c.calc_weight();
		if (c.w < opt->min_chain_weight) c.destroy();
		else chain[kept_chain_n++] = c;
	}
	chain.resize(kept_chain_n);
	// Sort chains by decreasing weight
	auto *chain_array = chain.data();
	ks_introsort_flt(kept_chain_n, chain_array);
	// Pairwise comparing chains
	std::vector<int> kept_index; // Storing the indexes of kept chains
	chain[0].kept = DOMINATE;
	kept_index.push_back(0);
	for (int i = 1; i < chain.size(); i++) {
		bool masked = false, dropped = false;
		for (auto j : kept_index) { // Check whether the current chain is overlapped by any kept chain
			int beg_max = std::max(chain[i].que_beg(), chain[j].que_beg());
			int end_min = std::min(chain[i].que_end(), chain[j].que_end());
			if (end_min > beg_max and not (chain[j].is_alt and not chain[i].is_alt)) { // Overlap found
				// When the kept chain is alternative and the current chain is primary, don't consider overlap
				int min_len = std::min(
					chain[i].que_end() - chain[i].que_beg(),
					chain[j].que_end() - chain[j].que_beg()
				);
				if ((float)(end_min - beg_max) >= (float)min_len * opt->mask_level
					and min_len < opt->max_chain_gap) {
					masked = true;
					// Record the first or heaviest chain that a kept chain shadows for more accurate MAPQ
					if (chain[j].first_cover < 0) chain[j].first_cover = i;
					if ((float)chain[i].w < (float) chain[j].w * opt->drop_ratio
						and chain[j].w - chain[i].w >= opt->min_seed_len * 2) {
						dropped = true;
						break;
					}
				}
			}
		}
		if (not dropped) {
			chain[i].kept = masked ?SHADOWED :DOMINATE;
			kept_index.push_back(i);
		}
	}
	for (auto i : kept_index) {
		const auto &c = chain[i];
		if (c.first_cover >= 0) {
			chain[c.first_cover].kept = 1; // the function of value 1?
		}
	}
	// Do not extend too many insignificant chains
	int to_extend = 0;
	for (auto &i : kept_index) {
		auto &c = chain[i];
		if (c.kept == DOMINATE) continue;
		to_extend++;
		if (to_extend >= opt->max_chain_extend) c.kept = DISCARD;
	}
	kept_chain_n = 0;
	for (const auto &c : chain) {
		if (c.kept == DISCARD) c.destroy();
		else chain[kept_chain_n++] = c;
	}
	chain.resize(kept_chain_n);
	return chain;
}

int CompAligner::seed_sw_score(const ngs_read &read, const seed_hit &s) {
	const int MEM_SHORT_LEN = 200;
	const int MEM_SHORT_EXT = 50;
	// The seed is longer than the max-extend; no need for SW
	if (s.len >= MEM_SHORT_LEN) return -1;

	int qb = std::max(0, s.qbeg - MEM_SHORT_EXT);
	int qe = std::min((int)read.len, s.qbeg + s.len + MEM_SHORT_EXT);
	int64_t rb = std::max(0L, s.rbeg - MEM_SHORT_EXT);
	int64_t re = std::min(bns->l_pac * 2, s.rbeg + s.len + MEM_SHORT_EXT);
	int64_t mid = (s.rbeg + s.rbeg + s.len) / 2;
	if (rb < bns->l_pac and bns->l_pac < re) {
		if (mid < bns->l_pac) re = bns->l_pac;
		else rb = bns->l_pac;
	}
	// The seed seems good enough; no need for SW
	if (qe - qb >= MEM_SHORT_LEN or re - rb >= MEM_SHORT_LEN) return -1;

	int rid;
	uint8_t *ref = bns_fetch_seq(bns, pac, &rb, mid, &re, &rid);
	kswr_t ans = ksw_align2(
			qe - qb, (uint8_t*)(read.bases + qb),
			re - rb, ref,
			5, opt->mat,
			opt->o_del, opt->e_del,
			opt->o_ins, opt->e_ins,
			KSW_XSTART,
			nullptr);
	free(ref);
	return ans.score;
}

void CompAligner::filter_seed_in_chain(const ngs_read &read, std::vector<seed_chain> &chain) {
	const float MEM_HSP_COEF = 1.1f;
	const float MEM_MINSC_COEF = 5.5f;
	const float MEM_SEEDSW_COEF = 0.05f;

	double min_l = opt->min_chain_weight ? MEM_HSP_COEF * (float)opt->min_chain_weight
										 : MEM_MINSC_COEF * log(read.len);
	if (min_l > MEM_SEEDSW_COEF * (float)read.len) return; // Don't suit to short reads
	int min_hsp_score = (int)(opt->a * min_l + 0.499);
	for (auto &c : chain) {
		int kept_seed_n = 0;
		for (int i = 0; i < c.n; i++) {
			auto &s = c.seeds[i];
			s.score = seed_sw_score(read, s);
			if (s.score < 0 or s.score >= min_hsp_score) {
				s.score = s.score < 0 ?read.len * opt->a :s.score;
				c.seeds[kept_seed_n++] = s;
			}
		}
		c.n = kept_seed_n;
	}
}

static inline int calc_max_gap(const mem_opt_t *opt, int len) {
	int l_del = (int) ((double)(len * opt->a - opt->o_del) / opt->e_del + 1.0);
	int l_ins = (int) ((double)(len * opt->a - opt->o_ins) / opt->e_ins + 1.0);
	int l = std::max(l_del, l_ins); // Maximum insertions or deletions
	l = std::max(l, 1); // At least one
	return std::min(l, opt->w * 2); // Should not exceed two times bandwidth
}

std::vector<align_region> CompAligner::
	extend_chain(const ngs_read &read, std::vector<seed_chain> &chain) {
	// Filter poor seeds in chain
	filter_seed_in_chain(read, chain);

	// Extend the chains of the read
	std::vector<align_region> regions;
	for (const auto &c : chain) {
		// All seeds are poor and discarded
		if (c.n == 0) { c.destroy(); continue; }

		// Get the max possible span
		int64_t span_l = bns->l_pac * 2, span_r = 0;
		for (int i = 0; i < c.n; i++) {
			const auto &s = c.seeds[i];
			// Margins on both sides of the seed
			int margin_l = s.qbeg, margin_r = (read.len - (s.qbeg + s.len));
			int64_t beg = s.rbeg - (margin_l + calc_max_gap(opt, margin_l));
			int64_t end = s.rbeg + s.len + (margin_r + calc_max_gap(opt, margin_r));
			span_l = std::min(span_l, beg);
			span_r = std::max(span_r, end);
		}
		span_l = std::max(span_l, 0L);
		span_r = std::min(span_r, bns->l_pac * 2);
		// Choose one side if crossing the forward-reverse boundary
		if (span_l < bns->l_pac and bns->l_pac < span_r) {
			// It is safe because all seeds on the chain are guaranteed to be on the same strand
			if (c.ref_beg() < bns->l_pac) span_r = bns->l_pac;
			else span_l = bns->l_pac;
		}
		// Prefetch the reference sequence for all possible seed extensions
		int rid; // On the same chromosome with the chain
		auto *ref_seq = bns_fetch_seq(bns, pac, &span_l, c.ref_beg(), &span_r, &rid);
		assert(rid == c.rid);

		// Sort seeds by score
		auto *sorted = (uint64_t*) malloc(c.n * sizeof(uint64_t));
		for (uint64_t i = 0; i < c.n; i++) {
			sorted[i] = (uint64_t)c.seeds[i].score << 32U | i;
		}
		ks_introsort_64(c.n, sorted);

		// Extend seeds from high score
		for (int k = c.n - 1; k >= 0; k--) {
			const auto &s = c.seeds[(uint32_t)sorted[k]];
			// Test whether extension has been made before
			bool contained = false;
			for (const auto &p : regions) {
				if (s.rbeg < p.rb or s.rbeg + s.len > p.re or  // The seed is not fully contained in the alignment
				    s.qbeg < p.qb or s.qbeg + s.len > p.qe) continue;
				if (s.len - p.seed_len0 > 0.1 * read.len) continue; // The seed might extend into a better alignment
				int que_dis = s.qbeg - p.qb; // The distance ahead of the seed on query
				int64_t ref_dis = s.rbeg - p.rb; // The distance ahead of the seed on reference
				// The maximal gap allowed in the region ahead of the seed
				int max_gap = calc_max_gap(opt, std::min((int64_t)que_dis, ref_dis));
				int width = std::min(max_gap, p.band_width); // Bounded by the band width
				// The seed is "around" the previous alignment
				if (abs(que_dis - ref_dis) < width) { contained = true; break; }

				// Check for the region behind the seed
				que_dis = p.qe - (s.qbeg + s.len); // The distance behind the seed on query
				ref_dis = p.re - (s.rbeg + s.len); // The distance behind the seed on reference
				max_gap = calc_max_gap(opt, std::min((int64_t)que_dis, ref_dis));
				width = std::min(max_gap, p.band_width);
				if (abs(que_dis - ref_dis) < width) { contained = true; break; }
			}
			// The seed is almost contained contained in an existing aligned region.
			// Further testing is needed to confirm it is not leading to a difference alignment.
			if (contained) {
				bool overlapped = false;
				for (int i = k + 1; i < c.n; i++) { // Loop for seeds have been extended
					if (sorted[i] == 0) continue; // The extension is skipped
					const auto &t = c.seeds[(uint32_t)sorted[i]];
					if (t.len < s.len * 0.95) continue; // Check overlap only if t is long enough
					if (s.qbeg <= t.qbeg and s.qbeg + s.len - t.qbeg >= s.len / 4 and
						t.qbeg - s.qbeg != t.rbeg - s.rbeg) { overlapped = true; break; }
					if (t.qbeg <= s.qbeg and t.qbeg + t.len - s.qbeg >= s.len / 4 and
					    s.qbeg - t.qbeg != s.rbeg - t.rbeg) { overlapped = true; break; }
				}
				if (not overlapped) { // No overlapping seeds, skip the extension
					sorted[k] = 0; // Mark this seed extension has not been performed
					continue;
				}
			}

			// Perform extension for the seed that could lead to a new aligned region
			align_region a = {0};
			a.rid = c.rid;
			a.frac_rep = c.frac_rep;
			a.seed_len0 = s.len;
			a.local_score = -1;
			a.true_score = -1;

			// Left Extension (target length is long enough for local alignment)
			int que_len = s.qbeg, tar_len = s.rbeg - span_l;
			int actual_bw_left = opt->w; // Actual bandwidth of the DP matrix
			const int MAX_BAND_TRY = 2;
			if (que_len > 0) {
				auto *qs = (uint8_t*) malloc(que_len * sizeof(uint8_t));
				for (int i = 0; i < que_len; i++) qs[i] = read.bases[s.qbeg - 1 - i];
				auto *ts = (uint8_t*) malloc(tar_len * sizeof(uint8_t));
				for (int i = 0; i < tar_len; i++) ts[i] = ref_seq[tar_len - 1 - i];
				// Terminals on query and target of the extension for local alignment
				int local_que_ext, local_tar_ext;
				// Semi-global alignment if the query sequence is entirely aligned.
				// The extension has reached the end of query, only recording terminal on target.
				int global_tar_ext, global_score;
				int max_off; // Maximum offset from the DP matrix diagonal
				for (int i = 0; i < MAX_BAND_TRY; i++) {
					int prev_score = a.local_score;
					actual_bw_left = opt->w << i;
					a.local_score = ksw_extend2(
							que_len, qs, tar_len, ts,
							5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
							actual_bw_left,
							opt->pen_clip5,
							opt->zdrop,
							s.len * opt->a,
							&local_que_ext, &local_tar_ext,
							&global_tar_ext, &global_score,
							&max_off);
					// Stop increasing the bandwidth if the DP score does not change or
					// maximum offsetting diagonal distance < 75% of the applied bandwidth
					if (a.local_score == prev_score or max_off < actual_bw_left / 2 + actual_bw_left / 4) break;
				}
				// Semi-global alignment is preferred if its score greater than the
				// best local score minus the penalty of left soft clipping.
				if (global_score > 0 and global_score > a.local_score - opt->pen_clip5) {
					a.qb = 0;
					a.rb = s.rbeg - global_tar_ext;
					a.true_score = global_score;
				} else { // Otherwise, choose local alignment
					a.qb = s.qbeg - local_que_ext;
					a.rb = s.rbeg - local_tar_ext;
					a.true_score = a.local_score;
				}
				free(qs); free(ts);
			} else {
				// The seed reaches the left end of the read
				a.local_score = a.true_score = s.len * opt->a;
				a.qb = 0; a.rb = s.rbeg;
			}

			// Right extension (similar to left extension)
			que_len = read.len - (s.qbeg + s.len);
			tar_len = span_r - (s.rbeg + s.len);
			int actual_bw_right = opt->w;
			if (que_len > 0) {
				int score0 = a.local_score; // The initial DP score
				int local_que_ext, local_tar_ext;
				int global_tar_ext, global_score;
				int max_off;
				for (int i = 0; i < MAX_BAND_TRY; i++) {
					int prev_score = a.local_score;
					actual_bw_right = opt->w << i;
					a.local_score = ksw_extend2(
							que_len, (uint8_t*)read.bases + (s.qbeg + s.len),
							tar_len, ref_seq + (s.rbeg + s.len) - span_l,
							5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
							actual_bw_right,
							opt->pen_clip3,
							opt->zdrop,
							score0,
							&local_que_ext, &local_tar_ext,
							&global_tar_ext, &global_score,
							&max_off);
					if (a.local_score == prev_score or max_off < actual_bw_right / 2 + actual_bw_right / 4) break;
				}
				if (global_score > 0 and global_score > a.local_score - opt->pen_clip3) {
					a.qe = read.len;
					a.re = s.rbeg + s.len + global_tar_ext;
					a.true_score += global_score - score0;
				} else {
					a.qe = s.qbeg + s.len + local_que_ext;
					a.re = s.rbeg + s.len + local_tar_ext;
					a.true_score += a.local_score - score0;
				}
			} else {
				a.qe = read.len;
				a.re = s.rbeg + s.len;
			}

			a.band_width = std::max(actual_bw_left, actual_bw_right);
			// Compute seed coverage on this alignment
			a.seed_cover = 0;
			for (int i = 0; i < c.n; i++) {
				const auto &t = c.seeds[i];
				if (t.qbeg >= a.qb and t.qbeg + t.len <= a.qe and
					t.rbeg >= a.rb and t.rbeg + t.len <= a.re) { // Seed fully contained
					// This is not very accurate, but good enough for approximate MAPQ
					a.seed_cover += t.len;
				}
			}
			regions.push_back(a);
		}

		free(sorted); free(ref_seq);
		c.destroy();
	}
	return regions;
}

void CompAligner::display_profile() {
	thread_aux total;
	for (int i = 0; i < threads_n; i++) total += thr_aux[i];
	auto full_read_match = total.full_read_match;
	auto shortcut = total.shortcut;
	fprintf(stderr, "Input %d reads in total\n", total_reads_count);
	fprintf(stderr, "Perfect Matched Reads: %d (%.2f %%)\n", full_read_match, 100.0 * full_read_match / total_reads_count);
	fprintf(stderr, "Reads go shortcut:     %d (%.2f %%)\n", shortcut, 100.0 * shortcut / total_reads_count);
}

void CompAligner::seed_and_extend(int _start, int _end, int tid) {
	auto &aux = thr_aux[tid];
	aux.forward_sst->clear(); aux.backward_sst->clear();
	_end = std::min(_end, (int) reads.size()); // Out of right boundary happens
	int n = _end - _start;
	long ref_position = -1; // Starting from this position might avoid many BWT-extension for fully matched reads
	for (int i = n - 1; i >= 0; i--) {
		auto &read = reads[_start + i];
		auto *bases = (uint8_t*) read.bases;
		for (int j = 0; j < read.len; j++) { // Convert ACGTN to 01234 if hasn't done so far
			if (bases[j] > 4) bases[j] = nst_nt4_table[bases[j]];
		}

		std::vector<bwtintv_t> &match = aux.match[i]; match.clear();
		bool full_match = false; // Whether this read could be full-length matched
		if (ref_position != -1 and ref_position - read.offset <= read.len / 3) {
			// Lookup forward SST to test if this read can be fully matched
			int node_id = 0;
			for (int j = ref_position - read.offset; j < read.len; j++) {
				if (bases[j] < 4) {
					if (j == ref_position - read.offset) node_id = aux.forward_sst->get_child(node_id, bases[j]);
					else node_id = aux.forward_sst->get_child(node_id, 3 - bases[j]);
				} else node_id = -1;
				if (node_id == -1) break;
			}
			if (node_id != -1) { // matched in forward SST till the end of read
				bwtintv_t ik = aux.forward_sst->get_intv(node_id), ok[4];
				// It is a potential full-match read
				for (int j = ref_position - read.offset - 1; j >= 0; j--) {
					if (bases[j] < 4) {
						bwt_extend(bwt, &ik, ok, 1);
						ik = ok[bases[j]];
					} else ik.x[2] = 0;
					if (ik.x[2] == 0) break;
				}
				if (ik.x[2] > 0) {
					full_match = true;
					aux.shortcut++;
					ik.info = read.len;
					if (mem_len(ik) >= opt->min_seed_len) match.push_back(ik);
				}
			}
		}

		if (not full_match) { // Go to the regular SMEM searching pass
			for (int j = 0; j < read.len; ) {
				j = collect_mem_with_sst(bases, read.len, j, 1, aux);
				for (const auto &m : aux.super_mem) {
					if (mem_len(m) >= opt->min_seed_len)
						match.push_back(m);
				}
			}
			ref_position = read.offset;
		}
		aux.full_read_match += (match.size() == 1 and mem_len(match[0]) == read.len);

		int old_n = (int)match.size();
		for (int j = 0; j < old_n; j++) {
			const auto &p = match[j] ;
			int beg = mem_beg(p), end = mem_end(p);
			if (end - beg < (int)(1.0 * opt->min_seed_len * opt->split_factor + .499) or p.x[2] > opt->split_width) continue;
			collect_mem_with_sst(bases, read.len, (beg + end) / 2, p.x[2] + 1, aux);
			for (const auto &m : aux.super_mem) {
				if (mem_len(m) >= opt->min_seed_len)
					match.push_back(m);
			}
		}

		if (opt->max_mem_intv > 0) {
			for (int j = 0; j < read.len; ) {
				if (bases[j] < 4) {
					bwtintv_t m;
					j = tem_forward_sst(bases, read.len, j, &m, aux);
					if (m.x[2] > 0) match.push_back(m);
				} else {
					j++;
				}
			}
		}
		std::sort(match.begin(), match.end(), mem_cmp);
	}

	// Find unique hit locations by merging and sorting
	auto unique_sal = aux.unique_sal; unique_sal.clear();
	for (int read_id = 0; read_id < n; read_id++) {
		const auto &mem = aux.match[read_id];
		auto &seed = aux.seed[read_id]; seed.clear();
		int array_id = 0;
		for (const auto &m : mem) {
			uint64_t step = m.x[2] > opt->max_occ ? m.x[2] / opt->max_occ : 1;
			for (uint64_t k = 0, count = 0; k < m.x[2] && count < opt->max_occ; k += step, count++) {
				seed_hit s = {0};
				s.qbeg = mem_beg(m);
				s.score = s.len = mem_len(m);
				seed.push_back(s);
				unique_sal.emplace_back(sal_request(m.x[0] + k, read_id, array_id));
				array_id++;
			}
		}
	}
	std::sort(unique_sal.begin(), unique_sal.end());
	// Lookup suffix array for hit locations
	uint64_t coordinate = 0;
	for (int i = 0; i < unique_sal.size(); i++) {
		const auto &p = unique_sal[i];
		if (i == 0 or unique_sal[i-1].que_location != p.que_location) {
			coordinate = bwt_sa(bwa_idx->bwt, p.que_location);
		}
		auto &s = aux.seed[p.read_id][p.array_id];
		s.rbeg = coordinate;
		s.rid = bns_intv2rid(bwa_idx->bns, s.rbeg, s.rbeg + s.len);
	}

	for (int i = 0; i < n; i++) {
		const auto &read = reads[_start + i];

		// Chaining co-linear seeds
		const auto &seed = aux.seed[i];
		auto chain = chaining(seed);
		// Calculate repetition fraction of chains
		const auto &mem = aux.match[i];
		int beg = 0, end = 0, l_repetition = 0;
		for (const auto &m : mem) { // Matches are sorted by interval already
			if (m.x[2] <= opt->max_occ) continue; // Not a repetitive match
			int b = mem_beg(m), e = mem_end(m);
			if (b > end) { l_repetition += end - beg; beg = b; end = e; }
			else end = std::max(end, e);
		}
		l_repetition += end - beg;
		for (auto &c : chain) {
			c.frac_rep = (float)l_repetition / read.len;
		}

		// Extend chains to alignments
		auto region = extend_chain(read, chain);

		// Output used variables of regions
		auto *out= &debug_out[_start + i];
		ksprintf(out, "Read %d has %ld chains, %ld aligned regions\n", _start + i, chain.size(), region.size());
		for (const auto &r : region) {
			ksprintf(out, "[%d,%d) => [%ld,%ld)\n", r.qb, r.qe, r.rb, r.re);
			ksprintf(out, "SW=%d, BW=%d, Final=%d, SeedCover=%d, Seed0=%d\n",
				r.local_score, r.band_width, r.true_score, r.seed_cover, r.seed_len0);
			ksprintf(out, "RID=%d, Freq=%.6f\n", r.rid, r.frac_rep);
			ksprintf(out, "\n");
		}
	}
}

void CompAligner::run(const char *fn) {
	// Prepare for thread auxiliary
	thr_aux = new thread_aux[threads_n];
	for (int i = 0; i < threads_n; i++) {
		auto &a = thr_aux[i];
		a.forward_sst = new SST(bwa_idx->bwt);
		a.backward_sst = new SST(bwa_idx->bwt);
	}

	fprintf(stderr, "Input test data with strand-corrected reads and overlapping information\n");
	std::ifstream in(fn); assert(in.is_open());
	char *buffer = new char[1024]; memset(buffer, 0, 1024 * sizeof(char));
	long off, global_offset = 0;
	long bytes = 0, round_n = 0;
	while (in >> buffer >> off) {
		global_offset += off;
		ngs_read read1;
		read1.bases = strdup(buffer);
		read1.len = strlen(buffer);
		read1.offset = global_offset;
		reads.push_back(read1);
		bytes += read1.len;
		total_reads_count++;
		if (bytes >= chunk_size * threads_n) {
			debug_out = (kstring_t*) calloc(reads.size(), sizeof(kstring_t));
			kt_for(
				threads_n,
				[](void *d, long i, int t) -> void {
					((CompAligner*)d)->seed_and_extend(i*BATCH_SIZE, (i+1)*BATCH_SIZE, t);
				},
				this,
				((int)reads.size() + BATCH_SIZE - 1) / BATCH_SIZE
			);
			for (auto &r : reads) { free(r.bases); free(r.sam); }
			fprintf(stderr, "Seed and Extend: %d reads processed\n", total_reads_count);
			for (int i = 0; i < reads.size(); i++) {
				fprintf(stdout, "%s\n", debug_out[i].s);
				free(debug_out[i].s);
			}
			free(debug_out);
			reads.clear(); bytes = 0; round_n++;
			if (round_n >= input_round_limit) break;
		}
	}
	display_profile();

	free(opt);
	bwa_idx_destroy(bwa_idx);
	for (int i = 0; i < threads_n; i++) {
		auto &a = thr_aux[i];
		delete a.forward_sst;
		delete a.backward_sst;
	}
	delete [] thr_aux;
	in.close();
	delete [] buffer;
}

void CompAligner::bwamem(const char *fn) {
	opt->n_threads = threads_n;

	fprintf(stderr, "Input test data with strand-corrected reads and overlapping information\n");
	std::ifstream in(fn); assert(in.is_open());
	char *buffer = new char[1024]; memset(buffer, 0, 1024 * sizeof(char));
	long off, bytes = 0, round_n = 0, n_processed = 0;
	std::vector<bseq1_t> bwa_seqs;
	while (in >> buffer >> off) {
		bseq1_t read1{0};
		read1.id = total_reads_count;
		read1.name = strdup(std::to_string(total_reads_count).c_str());
		read1.seq = strdup(buffer);
		read1.l_seq = (int) strlen(buffer);
		bwa_seqs.push_back(read1);
		total_reads_count++;
		bytes += read1.l_seq;
		if (bytes >= chunk_size * threads_n) {
			int n = bwa_seqs.size();
			mem_process_seqs(opt, bwt, bns, pac, n_processed, n, bwa_seqs.data(), nullptr);
			n_processed += n;
			for (int i = 0; i < n; i++) {
				const auto &s = bwa_seqs[i];
//				fprintf(stdout, "%s", s.sam);
				free(s.name); free(s.seq); free(s.sam);
			}

			fprintf(stderr, "BWA-MEM: %d reads processed\n", total_reads_count);
			bwa_seqs.clear(); bytes = 0; round_n++;
			if (round_n >= input_round_limit) break;
		}
	}

	free(opt);
	bwa_idx_destroy(bwa_idx);
	in.close();
	delete [] buffer;
}

int main(int argc, char *argv[]) {
	if (argc == 1) {
		fprintf(stderr, "Usage: PBM <bwa_index> <reads>\n");
		return 1;
	}
	CompAligner worker;
	int c; bool input_test = false;
	std::string mode = "test";
	while ((c = getopt(argc, argv, "t:m:r:")) >= 0) {
		if (c == 't') {
			worker.threads_n = strtol(optarg, nullptr, 10);
		} else if (c == 'r') {
			worker.input_round_limit = strtol(optarg, nullptr, 10);
		} else if (c == 'm') {
			mode = std::string(optarg);
		} else {
			fprintf(stderr, "Unrecognized option\n");
			exit(EXIT_FAILURE);
		}
	}
	worker.load_index(argv[optind]);
	if (mode == "comp") {
		worker.run(argv[optind + 1]);
	} else if (mode == "bwa") {
		worker.bwamem(argv[optind + 1]);
	}
	return 0;
}
