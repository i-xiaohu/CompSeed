//
// Created by ixiaohu on 2023/5/6.
//

#include "comp_seeding.h"

#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <fstream>

#include "../bwalib/ksw.h"
#include "../cstl/kthread.h"
#include "../cstl/kbtree.h"
#include "../cstl/ksort.h"
#include "../bwalib/utils.h"

//static time_profile t_prof[1024];

/********************************
 * Bookmark 1: Seeding with SST *
 ********************************/

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

int CompAligner::tem_forward_sst(const uint8_t *seq, int len, int start, bwtintv_t *mem, thread_aux &aux) const {
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

/*****************************
 * Bookmark 3 Smith-Waterman *
 *****************************/

static inline int calc_max_gap(const mem_opt_t *opt, int len) {
	int l_del = (int) ((double)(len * opt->a - opt->o_del) / opt->e_del + 1.0);
	int l_ins = (int) ((double)(len * opt->a - opt->o_ins) / opt->e_ins + 1.0);
	int l = std::max(l_del, l_ins); // Maximum insertions or deletions
	l = std::max(l, 1); // At least one
	return std::min(l, opt->w * 2); // Should not exceed two times bandwidth
}

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

int CompAligner::smith_waterman(int qlen, const uint8_t *query, int tlen, const uint8_t *target,
                                int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w,
                                int end_bonus,
                                int zdrop,
                                int h0,
                                int *_qle, int *_tle,
                                int *_gtle, int *_gscore,
                                int *_max_off,
                                thread_aux &aux) {

	struct eh_t { int32_t h, e; };
	eh_t *eh; // score array
	int8_t *qp; // query profile
	int i, j, k, oe_del = o_del + e_del, oe_ins = o_ins + e_ins, beg, end, max, max_i, max_j, max_ins, max_del, max_ie, gscore, max_off;
	assert(h0 > 0);
	// allocate memory
	qp = (int8_t*) malloc(qlen * m);
	eh = (eh_t*) calloc(qlen + 1, 8);
	// generate the query profile
	for (k = i = 0; k < m; ++k) {
		const int8_t *p = &mat[k * m];
		for (j = 0; j < qlen; ++j) qp[i++] = p[query[j]];
	}
	// fill the first row
	eh[0].h = h0; eh[1].h = h0 > oe_ins? h0 - oe_ins : 0;
	for (j = 2; j <= qlen && eh[j-1].h > e_ins; ++j)
		eh[j].h = eh[j-1].h - e_ins;
	// adjust $w if it is too large
	k = m * m;
	for (i = 0, max = 0; i < k; ++i) // get the max score
		max = max > mat[i]? max : mat[i];
	max_ins = (int)((double)(qlen * max + end_bonus - o_ins) / e_ins + 1.);
	max_ins = max_ins > 1? max_ins : 1;
	w = w < max_ins? w : max_ins;
	max_del = (int)((double)(qlen * max + end_bonus - o_del) / e_del + 1.);
	max_del = max_del > 1? max_del : 1;
	w = w < max_del? w : max_del;
	// DP loop
	max = h0, max_i = max_j = -1; max_ie = -1, gscore = -1;
	max_off = 0;
	beg = 0, end = qlen;
	for (i = 0; LIKELY(i < tlen); ++i) {
		int t, f = 0, h1, m = 0, mj = -1;
		int8_t *q = &qp[target[i] * qlen];
		// apply the band and the constraint (if provided)
		if (beg < i - w) beg = i - w;
		if (end > i + w + 1) end = i + w + 1;
		if (end > qlen) end = qlen;
		// compute the first column
		if (beg == 0) {
			h1 = h0 - (o_del + e_del * (i + 1));
			if (h1 < 0) h1 = 0;
		} else h1 = 0;
		for (j = beg; LIKELY(j < end); ++j) {
			// At the beginning of the loop: eh[j] = { H(i-1,j-1), E(i,j) }, f = F(i,j) and h1 = H(i,j-1)
			// Similar to SSE2-SW, cells are computed in the following order:
			//   H(i,j)   = max{H(i-1,j-1)+S(i,j), E(i,j), F(i,j)}
			//   E(i+1,j) = max{H(i,j)-gapo, E(i,j)} - gape
			//   F(i,j+1) = max{H(i,j)-gapo, F(i,j)} - gape
			eh_t *p = &eh[j];
			int h, M = p->h, e = p->e; // get H(i-1,j-1) and E(i-1,j)
			p->h = h1;          // set H(i,j-1) for the next row
			M = M? M + q[j] : 0;// separating H and M to disallow a cigar like "100M3I3D20M"
			h = M > e? M : e;   // e and f are guaranteed to be non-negative, so h>=0 even if M<0
			h = h > f? h : f;
			h1 = h;             // save H(i,j) to h1 for the next column
			mj = m > h? mj : j; // record the position where max score is achieved
			m = m > h? m : h;   // m is stored at eh[mj+1]
			t = M - oe_del;
			t = t > 0? t : 0;
			e -= e_del;
			e = e > t? e : t;   // computed E(i+1,j)
			p->e = e;           // save E(i+1,j) for the next row
			t = M - oe_ins;
			t = t > 0? t : 0;
			f -= e_ins;
			f = f > t? f : t;   // computed F(i,j+1)
		}
		eh[end].h = h1; eh[end].e = 0;
		if (j == qlen) {
			max_ie = gscore > h1? max_ie : i;
			gscore = gscore > h1? gscore : h1;
		}
		if (m > max) {
			max = m, max_i = i, max_j = mj;
			max_off = max_off > abs(mj - i)? max_off : abs(mj - i);
		} else if (zdrop > 0) {
			if (i - max_i > mj - max_j) {
				if (max - m - ((i - max_i) - (mj - max_j)) * e_del > zdrop) break;
			} else {
				if (max - m - ((mj - max_j) - (i - max_i)) * e_ins > zdrop) break;
			}
		}
		// update beg and end for the next round
		for (j = beg; LIKELY(j < end) && eh[j].h == 0 && eh[j].e == 0; ++j);
		beg = j;
		for (j = end; LIKELY(j >= beg) && eh[j].h == 0 && eh[j].e == 0; --j);
		end = j + 2 < qlen? j + 2 : qlen;
		//beg = 0; end = qlen; // uncomment this line for debugging
	}
	free(eh); free(qp);
	if (_qle) *_qle = max_j + 1;
	if (_tle) *_tle = max_i + 1;
	if (_gtle) *_gtle = max_ie + 1;
	if (_gscore) *_gscore = gscore;
	if (_max_off) *_max_off = max_off;
	return max;
}

std::vector<align_region> CompAligner::
	extend_chain(const ngs_read &read, std::vector<seed_chain> &chain, thread_aux &aux) {
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
					a.local_score = smith_waterman(
							que_len, qs, tar_len, ts,
							5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
							actual_bw_left,
							opt->pen_clip5,
							opt->zdrop,
							s.len * opt->a,
							&local_que_ext, &local_tar_ext,
							&global_tar_ext, &global_score,
							&max_off, aux);
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
					a.local_score = smith_waterman(
							que_len, (uint8_t*)read.bases + (s.qbeg + s.len),
							tar_len, ref_seq + (s.rbeg + s.len) - span_l,
							5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins,
							actual_bw_right,
							opt->pen_clip3,
							opt->zdrop,
							score0,
							&local_que_ext, &local_tar_ext,
							&global_tar_ext, &global_score,
							&max_off, aux);
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

void CompAligner::display_profile(const thread_aux &total) {
	auto full_read_match = total.full_read_match;
	auto shortcut = total.shortcut;
	fprintf(stderr, "Input reads in total:  %ld\n", processed_n);
	fprintf(stderr, "Read length:           %d\n", read_length);
	fprintf(stderr, "Perfect Matched Reads: %d (%.2f %%)\n", full_read_match, 100.0 * full_read_match / processed_n);
	fprintf(stderr, "Reads go shortcut:     %d (%.2f %%)\n", shortcut, 100.0 * shortcut / processed_n);
	fprintf(stderr, "Average calling of BWT-extend for each read: %.2f\n", 1.0 * total.bwt_call_times / processed_n);
	fprintf(stderr, "Average calling of SAL for each read:        %.2f\n", 1.0 * total.sal_call_times / processed_n);
	fprintf(stderr, "Seeding cost %.2f CPU seconds\n", total.seeding_cpu_sec);
	fprintf(stderr, "Real %.2f seconds in total threads\n", total.seeding_real_sec);
}

void CompAligner::seed_and_extend(int _start, int _end, int tid) {
	double real_start = realtime();
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
//		if (ref_position != -1 and ref_position - read.offset <= read.len / 3) {
//			// Lookup forward SST to test if this read can be fully matched
//			int node_id = 0;
//			for (int j = ref_position - read.offset; j < read.len; j++) {
//				if (bases[j] < 4) {
//					if (j == ref_position - read.offset) node_id = aux.forward_sst->get_child(node_id, bases[j]);
//					else node_id = aux.forward_sst->get_child(node_id, 3 - bases[j]);
//				} else node_id = -1;
//				if (node_id == -1) break;
//			}
//			if (node_id != -1) { // matched in forward SST till the end of read
//				bwtintv_t ik = aux.forward_sst->get_intv(node_id), ok[4];
//				// It is a potential full-match read
//				for (int j = ref_position - read.offset - 1; j >= 0; j--) {
//					if (bases[j] < 4) {
//						bwt_extend(bwt, &ik, ok, 1);
//						ik = ok[bases[j]];
//					} else ik.x[2] = 0;
//					if (ik.x[2] == 0) break;
//				}
//				if (ik.x[2] > 0) {
//					full_match = true;
//					aux.shortcut++;
//					ik.info = read.len;
//					if (mem_len(ik) >= opt->min_seed_len) match.push_back(ik);
//				}
//			}
//		}

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
		thr_aux[tid].full_read_match += (match.size() == 1 and mem_len(match[0]) == read.len);

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

	if (print_seed) { // Print all seeds to debug output buffer
		for (int i = 0; i < n; i++) {
			const auto &seed = aux.seed[i];
			kstring_t *o = &debug_out[_start + i];
			ksprintf(o, "Read %ld, Seed %ld\n", processed_n + _start + i, seed.size());
			for (int j = 0; j < seed.size(); j++) {
				ksprintf(o, "%d:{%d, %ld, %d}\n", j+1, seed[j].qbeg, seed[j].rbeg, seed[j].len);
			}
		}
	}
	// The code below that extends seeds to full alignments is deleted
	thr_aux[tid].seeding_real_sec += realtime() - real_start;
}

void CompAligner::run(const char *fn) {
	// Prepare for thread auxiliary
	threads_n = opt->n_threads;
	thread_aux total;

	fprintf(stderr, "Running compressive seeding...\n");
	gzFile in = gzopen(fn, "r"); assert(in != nullptr);
	char *buffer = new char[1024]; memset(buffer, 0, 1024 * sizeof(char));
	while (true) {
		// Input a batch of data
		long bytes = 0; reads.clear();
		while (gzgets(in, buffer, 1024)) {
		    int len = strlen(buffer); buffer[--len] = '\0';
			ngs_read read1;
			read1.bases = strdup(buffer);
			read_length = read1.len = len;
			reads.push_back(read1);
			bytes += read1.len;
			if (bytes >= actual_chunk_size) break;
		}
		if (reads.empty()) break; // End Of File

		// Normalize all reordered reads (restoring offset and correcting strand)

		// Processing (I/O thread not supported yet)
		if (print_seed) debug_out = (kstring_t*) calloc(reads.size(), sizeof(kstring_t));
		thr_aux = new thread_aux[threads_n];
		for (int i = 0; i < threads_n; i++) {
			auto &a = thr_aux[i];
			a.forward_sst = new SST(bwa_idx->bwt);
			a.backward_sst = new SST(bwa_idx->bwt);
		}

		double cpu_start = cputime();
		kt_for(
				threads_n,
				[](void *d, long i, int t) -> void {
					((CompAligner*)d)->seed_and_extend(i*BATCH_SIZE, (i+1)*BATCH_SIZE, t);
				},
				this,
				((int)reads.size() + BATCH_SIZE - 1) / BATCH_SIZE
		);
		thr_aux[0].seeding_cpu_sec += cputime() - cpu_start;

		for (int i = 0; i < threads_n; i++) total += thr_aux[i];
		for (int i = 0; i < threads_n; i++) {
			auto &a = thr_aux[i];
			delete a.forward_sst;
			delete a.backward_sst;
		}
		delete [] thr_aux;
		for (auto &r : reads) { free(r.bases); free(r.sam); }
		if (print_seed) {
			for (int i = 0; i < reads.size(); i++) {
				fprintf(stdout, "%s\n", debug_out[i].s);
				free(debug_out[i].s);
			}
			free(debug_out);
		}
		processed_n += reads.size();
		fprintf(stderr, "%ld reads processed\n", processed_n);
	}
	display_profile(total);

	bwa_idx_destroy(bwa_idx);
	gzclose(in);
	delete [] buffer;
}

/********************************
 * Bookmark 4: BWA-MEM Seeding  *
 ********************************/

mem_opt_t *mem_opt_init()
{
	mem_opt_t *o;
	o = (mem_opt_t*) calloc(1, sizeof(mem_opt_t));
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
	kv_resize(bwtintv_t, a->mem, 512);
	kv_resize(bwtintv_t, a->mem1, 512);
	kv_resize(bwtintv_t, *a->tmpv[0], 512);
	kv_resize(bwtintv_t, *a->tmpv[1], 512);
	return a;
}

void smem_aux_destroy(smem_aux_t *a) {
	free(a->tmpv[0]->a); free(a->tmpv[0]);
	free(a->tmpv[1]->a); free(a->tmpv[1]);
	free(a->mem.a); free(a->mem1.a);
	free(a);
}

#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)
static void mem_collect_intv(const mem_opt_t *opt, const bwt_t *bwt, int len, const uint8_t *seq, smem_aux_t *a)
{
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
			if ((uint32_t)a->mem1.a[i].info - (a->mem1.a[i].info>>32) >= opt->min_seed_len)
				kv_push(bwtintv_t, a->mem, a->mem1.a[i]);
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
	// sort
	ks_introsort(mem_intv, a->mem.n, a->mem.a);
}

void CompAligner::bwa_collect_seed(int seq_id, int tid) {
	double real_start = realtime();
	auto &read = reads[seq_id];
	auto *bases = (uint8_t*) read.bases;
	for (int j = 0; j < read.len; j++) { // Convert ACGTN to 01234 if hasn't done so far
		if (bases[j] > 4) bases[j] = nst_nt4_table[bases[j]];
	}
	auto *aux = smem_aux_init();
	mem_collect_intv(opt, bwt, read.len, bases, aux);
	for (int i = 0; i < aux->mem.n; i++) {
		if (mem_len(aux->mem.a[i]) == read.len) {
			thr_aux[tid].full_read_match++;
			break;
		}
	}
//	auto &seed = thr_aux[tid].seed[0]; seed.clear();
	for (int i = 0; i < aux->mem.n; i++) {
		const auto &p = aux->mem.a[i];
		int slen = mem_len(p); // seed length
		int64_t step = p.x[2] > opt->max_occ? p.x[2] / opt->max_occ : 1;
		for (int64_t k = 0, count = 0; k < p.x[2] && count < opt->max_occ; k += step, ++count) {
			seed_hit s = {0};
			s.rbeg = bwt_sa(bwt, p.x[0] + k); // this is the base coordinate in the forward-reverse reference
			s.qbeg = mem_beg(p);
			s.score= s.len = slen;
			s.rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
//			seed.push_back(s);
		}
	}
	smem_aux_destroy(aux);

	if (print_seed) {
		kstring_t *o = &debug_out[seq_id];
//		ksprintf(o, "Read %ld, Seed %ld\n", processed_n + seq_id, seed.size());
//		for (int j = 0; j < seed.size(); j++) {
//			ksprintf(o, "%d:{%d, %ld, %d}\n", j+1, seed[j].qbeg, seed[j].rbeg, seed[j].len);
//		}
	}
	thr_aux[tid].seeding_real_sec += realtime() - real_start;
}

void CompAligner::bwamem(const char *fn) {
	fprintf(stderr, "Running BWA-MEM seeding\n");
	// Prepare for thread auxiliary
	threads_n = opt->n_threads;
	thread_aux total;

	gzFile in = gzopen(fn, "r"); assert(in != nullptr);
	char *buffer = new char[1024]; memset(buffer, 0, 1024 * sizeof(char));
	while (true) {
		// Input a batch of data
		long bytes = 0; reads.clear();
		while (gzgets(in, buffer, 1024)) {
			int len = strlen(buffer); buffer[--len] = '\0';
			ngs_read read1;
			read1.bases = strdup(buffer);
			read_length = read1.len = len;
			reads.push_back(read1);
			bytes += read1.len;
			if (bytes >= actual_chunk_size) break;
		}
		if (reads.empty()) break; // End Of File

		// Normalize all reordered reads (restoring offset and correcting strand)

		// Processing (I/O thread not supported yet)
		if (print_seed) debug_out = (kstring_t*) calloc(reads.size(), sizeof(kstring_t));
		thr_aux = new thread_aux[threads_n];
		for (int i = 0; i < threads_n; i++) thr_aux[i].mem_aux = smem_aux_init();

		double cpu_start = cputime();
		kt_for(
				threads_n,
				[](void *d, long i, int t) -> void {
					((CompAligner*)d)->bwa_collect_seed(i, t);
				},
				this,
				(int)reads.size()
		);
		thr_aux[0].seeding_cpu_sec += cputime() - cpu_start;

		for (int i = 0; i < threads_n; i++) total += thr_aux[i];
		for (int i = 0; i < threads_n; i++) smem_aux_destroy(thr_aux[i].mem_aux);
		delete [] thr_aux;
		for (auto &r : reads) { free(r.bases); free(r.sam); }
		if (print_seed) {
			for (int i = 0; i < reads.size(); i++) {
				fprintf(stdout, "%s\n", debug_out[i].s);
				free(debug_out[i].s);
			}
			free(debug_out);
		}
		processed_n += reads.size();
		fprintf(stderr, "%ld reads processed\n", processed_n);
	}
	display_profile(total);

	bwa_idx_destroy(bwa_idx);
	gzclose(in);
	delete [] buffer;
}

struct worker {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	smem_aux_t **aux;
	ngs_read *seqs;
};

void bwa_worker(void *data, long seq_id, int tid) {
	auto *w = (worker*) data;
	const auto *opt = w->opt;
	const auto *bwt = w->bwt;
	const auto *bns = w->bns;
	const auto *pac = w->pac;
	const auto &read = w->seqs[seq_id];
	for (int i = 0; i < read.len; i++) {
		if (read.bases[i] > 4) {
			read.bases[i] = nst_nt4_table[read.bases[i]];
		}
	}
	auto *aux = w->aux[tid];
	mem_collect_intv(opt, bwt, read.len, (const uint8_t*) read.bases, aux);
	for (int i = 0; i < aux->mem.n; i++) {
		const auto &p = aux->mem.a[i];
		int slen = mem_len(p); // seed length
		int64_t step = p.x[2] > opt->max_occ? p.x[2] / opt->max_occ : 1;
		for (int64_t k = 0, count = 0; k < p.x[2] && count < opt->max_occ; k += step, ++count) {
			seed_hit s = {0};
			s.rbeg = bwt_sa(bwt, p.x[0] + k); // this is the base coordinate in the forward-reverse reference
			s.qbeg = mem_beg(p);
			s.score= s.len = slen;
			s.rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
		}
	}
}

void bwa_c_style(const char *index_fn, const char *read_fn, int actual_chunk_size, const mem_opt_t *opt) {
	fprintf(stderr, "Running BWA-MEM seeding in C style\n");
	bwaidx_t *idx = bwa_idx_load(index_fn, BWA_IDX_ALL);
	gzFile in = gzopen(read_fn, "r"); assert(in != nullptr);
	char *buffer = new char[1024]; memset(buffer, 0, 1024 * sizeof(char));
	double total_cpu_time = 0;
	long processed_n = 0;
	while (true) {
		// Input a batch of data
		long bytes = 0;
		std::vector<ngs_read> reads;
		while (gzgets(in, buffer, 1024)) {
			int len = strlen(buffer); buffer[--len] = '\0';
			ngs_read read1;
			read1.bases = strdup(buffer);
			read1.len = len;
			reads.push_back(read1);
			bytes += read1.len;
			if (bytes >= actual_chunk_size) break;
		}
		if (reads.empty()) break; // End Of File

		worker w = {};
		w.opt = opt;
		w.bwt = idx->bwt;
		w.bns = idx->bns;
		w.pac = idx->pac;
		w.seqs = reads.data();
		w.aux = (smem_aux_t**) malloc(opt->n_threads * sizeof(smem_aux_t*));
		for (int i = 0; i < opt->n_threads; i++) w.aux[i] = smem_aux_init();

		double cpu_start = cputime();
		kt_for(
				opt->n_threads,
				bwa_worker,
				&w,
				(int)reads.size()
		);
		total_cpu_time += cputime() - cpu_start;

		for (int i = 0; i < opt->n_threads; i++) smem_aux_destroy(w.aux[i]); free(w.aux);
		for (auto &r : reads) { free(r.bases); free(r.sam); }
		processed_n += reads.size();
		fprintf(stderr, "%ld reads processed\n", processed_n);
	}
	fprintf(stderr, "Seeding cost %.2f CPU seconds\n", total_cpu_time);

	bwa_idx_destroy(idx);
	gzclose(in);
	delete [] buffer;
}