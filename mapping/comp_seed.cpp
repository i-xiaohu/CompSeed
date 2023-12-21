//
// Created by ixiaohu on 2023/5/6.
//

#include "comp_seed.h"

#include <cstring>
#include <cassert>
#include <algorithm>
#include <cmath>
#include <climits>

#include "../cstl/kthread.h"
#include "../cstl/kbtree.h"
#include "../cstl/ksort.h"
#include "../bwalib/ksw.h"
#include "../bwalib/utils.h"
#include "macro.h"
#include "bandedSWA.h"
#include "memcpy_bwamem.h"

extern thread_aux_t tprof;

typedef kvec_t(int) int_v;

mem_opt_t *mem_opt_init() {
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

/********************
 * Seeding with SST *
 ********************/

int collect_mem_with_sst(const uint8_t *seq, int len, int pivot, int min_hits, thread_aux_t &aux) {
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
			aux.bwt_query_times++;
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
				aux.bwt_query_times++;
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

int tem_forward_sst(const mem_opt_t *opt, const uint8_t *seq, int len, int start, bwtintv_t *mem, thread_aux_t &aux) {
	memset(mem, 0, sizeof(bwtintv_t));
	if (seq[start] > 3) return start + 1;
	int node_id = aux.forward_sst->query_forward_child(0, seq[start]);
	bwtintv_t ik = aux.forward_sst->get_intv(node_id);
	for (int i = start + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = aux.forward_sst->query_forward_child(node_id, c);
			ik = aux.forward_sst->get_intv(node_id);
			aux.bwt_query_times++;
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

/************
 * Chaining *
 ************/

typedef struct { size_t n, m; mem_chain_t *a;  } mem_chain_v;

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, mem_chain_t, chain_cmp)

// return 1 if the seed is merged into the chain
static int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c, const mem_seed_t *p, int seed_rid) {
	int64_t qend, rend, x, y;
	const mem_seed_t *last = &c->seeds[c->n-1];
	qend = last->qbeg + last->len;
	rend = last->rbeg + last->len;
	if (seed_rid != c->rid) return 0; // different chr; request a new chain
	if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend && p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
		return 1; // contained seed; do nothing
	assert((last->rbeg < l_pac) == (c->seeds[0].rbeg < l_pac));
	if ((last->rbeg < l_pac || c->seeds[0].rbeg < l_pac) && p->rbeg >= l_pac) return 0; // don't chain if on different strand
	x = p->qbeg - last->qbeg; // always non-negtive
	y = p->rbeg - last->rbeg;
	if (y >= 0 && x - y <= opt->w && y - x <= opt->w && x - last->len < opt->max_chain_gap && y - last->len < opt->max_chain_gap) { // grow the chain
		if (c->n == c->m) {
			c->m <<= 1;
			c->seeds = (mem_seed_t*) realloc(c->seeds, c->m * sizeof(mem_seed_t));
		}
		c->seeds[c->n++] = *p;
		return 1;
	}
	return 0; // request to add a new chain
}

int mem_chain_weight(const mem_chain_t *c) {
	int64_t end;
	int j, w = 0, tmp;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->qbeg >= end) w += s->len;
		else if (s->qbeg + s->len > end) w += s->qbeg + s->len - end;
		end = end > s->qbeg + s->len? end : s->qbeg + s->len;
	}
	tmp = w; w = 0;
	for (j = 0, end = 0; j < c->n; ++j) {
		const mem_seed_t *s = &c->seeds[j];
		if (s->rbeg >= end) w += s->len;
		else if (s->rbeg + s->len > end) w += s->rbeg + s->len - end;
		end = end > s->rbeg + s->len? end : s->rbeg + s->len;
	}
	w = w < tmp? w : tmp;
	return w < 1<<30? w : (1<<30)-1;
}

void mem_print_chain(const bntseq_t *bns, mem_chain_v *chn) {
	int i, j;
	for (i = 0; i < chn->n; ++i) {
		mem_chain_t *p = &chn->a[i];
		err_printf("* Found CHAIN(%d): n=%d; weight=%d", i, p->n, mem_chain_weight(p));
		for (j = 0; j < p->n; ++j) {
			bwtint_t pos;
			int is_rev;
			pos = bns_depos(bns, p->seeds[j].rbeg, &is_rev);
			if (is_rev) pos -= p->seeds[j].len - 1;
			err_printf("\t%d;%d;%d,%ld(%s:%c%ld)", p->seeds[j].score, p->seeds[j].len, p->seeds[j].qbeg, (long)p->seeds[j].rbeg, bns->anns[p->rid].name, "+-"[is_rev], (long)(pos - bns->anns[p->rid].offset) + 1);
		}
		err_putchar('\n');
	}
}

mem_chain_v mem_chain(const mem_opt_t *opt, const bntseq_t *bns, int len, const std::vector<bwtintv_t> &mem, const std::vector<mem_seed_t> &seed) {
	mem_chain_v chain; kv_init(chain);
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match

	kbtree_t(chn) *tree = kb_init(chn, KB_DEFAULT_SIZE);
	for (const auto &s : seed) {
		mem_chain_t tmp, *lower, *upper;
		int rid, to_add = 0;
		tmp.pos = s.rbeg; // this is the base coordinate in the forward-reverse reference
		rid = bns_intv2rid(bns, s.rbeg, s.rbeg + s.len);
		if (rid < 0) continue; // bridging multiple reference sequences or the forward-reverse boundary; TODO: split the seed; don't discard it!!!
		if (kb_size(tree)) {
			kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
			if (!lower || !test_and_merge(opt, bns->l_pac, lower, &s, rid)) to_add = 1;
		} else to_add = 1;
		if (to_add) { // add the seed as a new chain
			tmp.n = 1; tmp.m = 4;
			tmp.seeds = (mem_seed_t*) calloc(tmp.m, sizeof(mem_seed_t));
			tmp.seeds[0] = s;
			tmp.rid = rid;
			tmp.is_alt = !!bns->anns[rid].is_alt;
			kb_putp(chn, tree, &tmp);
		}
	}

	kv_resize(mem_chain_t, chain, kb_size(tree));
	#define traverse_func(p_) (chain.a[chain.n++] = *(p_))
	__kb_traverse(mem_chain_t, tree, traverse_func);
	#undef traverse_func

	// Calculate repetition fraction
	int beg = 0, end = 0, l_rep = 0;
	for (const auto &m : mem) { // Matches are sorted by interval already
		if (m.x[2] <= opt->max_occ) continue; // Not a repetitive match
		int b = mem_beg(m), e = mem_end(m);
		if (b > end) { l_rep += end - beg; beg = b; end = e; }
		else end = std::max(end, e);
	}
	l_rep += end - beg;
	for (int i = 0; i < chain.n; i++) chain.a[i].frac_rep = (float)l_rep / len;
	if (bwa_verbose >= 4) printf("* fraction of repetitive seeds: %.3f\n", (float)l_rep / len);

	kb_destroy(chn, tree);
	return chain;
}

/********************
 * Filtering chains *
 ********************/

#define chn_beg(ch) ((ch).seeds->qbeg)
#define chn_end(ch) ((ch).seeds[(ch).n-1].qbeg + (ch).seeds[(ch).n-1].len)

#define flt_lt(a, b) ((a).w > (b).w)
KSORT_INIT(mem_flt, mem_chain_t, flt_lt)

int mem_chain_flt(const mem_opt_t *opt, int n_chn, mem_chain_t *a) {
	int i, k;
	int_v chains = {0,0,0}; // this keeps int indices of the non-overlapping chains
	if (n_chn == 0) return 0; // no need to filter
	// compute the weight of each chain and drop chains with small weight
	for (i = k = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		c->first = -1; c->kept = 0;
		c->w = mem_chain_weight(c);
		if (c->w < opt->min_chain_weight) free(c->seeds);
		else a[k++] = *c;
	}
	n_chn = k;
	ks_introsort(mem_flt, n_chn, a);
	// pairwise chain comparisons
	a[0].kept = 3;
	kv_push(int, chains, 0);
	for (i = 1; i < n_chn; ++i) {
		int large_ovlp = 0;
		for (k = 0; k < chains.n; ++k) {
			int j = chains.a[k];
			int b_max = chn_beg(a[j]) > chn_beg(a[i])? chn_beg(a[j]) : chn_beg(a[i]);
			int e_min = chn_end(a[j]) < chn_end(a[i])? chn_end(a[j]) : chn_end(a[i]);
			if (e_min > b_max && (!a[j].is_alt || a[i].is_alt)) { // have overlap; don't consider ovlp where the kept chain is ALT while the current chain is primary
				int li = chn_end(a[i]) - chn_beg(a[i]);
				int lj = chn_end(a[j]) - chn_beg(a[j]);
				int min_l = li < lj? li : lj;
				if (e_min - b_max >= min_l * opt->mask_level && min_l < opt->max_chain_gap) { // significant overlap
					large_ovlp = 1;
					if (a[j].first < 0) a[j].first = i; // keep the first shadowed hit s.t. mapq can be more accurate
					if (a[i].w < a[j].w * opt->drop_ratio && a[j].w - a[i].w >= opt->min_seed_len<<1)
						break;
				}
			}
		}
		if (k == chains.n) {
			kv_push(int, chains, i);
			a[i].kept = large_ovlp? 2 : 3;
		}
	}
	for (i = 0; i < chains.n; ++i) {
		mem_chain_t *c = &a[chains.a[i]];
		if (c->first >= 0) a[c->first].kept = 1;
	}
	free(chains.a);
	for (i = k = 0; i < n_chn; ++i) { // don't extend more than opt->max_chain_extend .kept=1/2 chains
		if (a[i].kept == 0 || a[i].kept == 3) continue;
		if (++k >= opt->max_chain_extend) break;
	}
	for (; i < n_chn; ++i)
		if (a[i].kept < 3) a[i].kept = 0;
	for (i = k = 0; i < n_chn; ++i) { // free discarded chains
		mem_chain_t *c = &a[i];
		if (c->kept == 0) free(c->seeds);
		else a[k++] = a[i];
	}
	return k;
}

/*********************************
 * Test if a seed is good enough *
 *********************************/

#define MEM_SHORT_EXT 50
#define MEM_SHORT_LEN 200

#define MEM_HSP_COEF 1.1f
#define MEM_MINSC_COEF 5.5f
#define MEM_SEEDSW_COEF 0.05f

int mem_seed_sw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_seed_t *s) {
	int qb, qe, rid;
	int64_t rb, re, mid, l_pac = bns->l_pac;
	uint8_t *rseq = 0;
	kswr_t x;

	if (s->len >= MEM_SHORT_LEN) return -1; // the seed is longer than the max-extend; no need to do SW
	qb = s->qbeg, qe = s->qbeg + s->len;
	rb = s->rbeg, re = s->rbeg + s->len;
	mid = (rb + re) >> 1;
	qb -= MEM_SHORT_EXT; qb = qb > 0? qb : 0;
	qe += MEM_SHORT_EXT; qe = qe < l_query? qe : l_query;
	rb -= MEM_SHORT_EXT; rb = rb > 0? rb : 0;
	re += MEM_SHORT_EXT; re = re < l_pac<<1? re : l_pac<<1;
	if (rb < l_pac && l_pac < re) {
		if (mid < l_pac) re = l_pac;
		else rb = l_pac;
	}
	if (qe - qb >= MEM_SHORT_LEN || re - rb >= MEM_SHORT_LEN) return -1; // the seed seems good enough; no need to do SW

	rseq = bns_fetch_seq(bns, pac, &rb, mid, &re, &rid);
	x = ksw_align2(qe - qb, (uint8_t*)query + qb, re - rb, rseq, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, KSW_XSTART, 0);
	free(rseq);
	return x.score;
}

void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, int n_chn, mem_chain_t *a) {
	double min_l = opt->min_chain_weight? MEM_HSP_COEF * opt->min_chain_weight : MEM_MINSC_COEF * log(l_query);
	int i, j, k, min_HSP_score = (int)(opt->a * min_l + .499);
	if (min_l > MEM_SEEDSW_COEF * l_query) return; // don't run the following for short reads
	for (i = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		for (j = k = 0; j < c->n; ++j) {
			mem_seed_t *s = &c->seeds[j];
			s->score = mem_seed_sw(opt, bns, pac, l_query, query, s);
			if (s->score < 0 || s->score >= min_HSP_score) {
				s->score = s->score < 0? s->len * opt->a : s->score;
				c->seeds[k++] = *s;
			}
		}
		c->n = k;
	}
}

/****************************************
 * Construct the alignment from a chain *
 ****************************************/

static inline int cal_max_gap(const mem_opt_t *opt, int qlen) {
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

#define MAX_BAND_TRY  2

void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av) {
	int i, k, rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
	int64_t l_pac = bns->l_pac, rmax[2], tmp, max = 0;
	const mem_seed_t *s;
	uint8_t *rseq = 0;
	uint64_t *srt;

	if (c->n == 0) return;
	// get the max possible span
	rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) {
		int64_t b, e;
		const mem_seed_t *t = &c->seeds[i];
		b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
		e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) + cal_max_gap(opt, l_query - t->qbeg - t->len));
		rmax[0] = rmax[0] < b? rmax[0] : b;
		rmax[1] = rmax[1] > e? rmax[1] : e;
		if (t->len > max) max = t->len;
	}
	rmax[0] = rmax[0] > 0? rmax[0] : 0;
	rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
	if (rmax[0] < l_pac && l_pac < rmax[1]) { // crossing the forward-reverse boundary; then choose one side
		if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac; // this works because all seeds are guaranteed to be on the same strand
		else rmax[0] = l_pac;
	}
	// retrieve the reference sequence
	rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
	assert(c->rid == rid);

	srt = (uint64_t*) malloc(c->n * sizeof(uint64_t));
	for (i = 0; i < c->n; ++i)
		srt[i] = (uint64_t)c->seeds[i].score<<32 | i;
	ks_introsort_64(c->n, srt);

	for (k = c->n - 1; k >= 0; --k) {
		mem_alnreg_t *a;
		s = &c->seeds[(uint32_t)srt[k]];

		for (i = 0; i < av->n; ++i) { // test whether extension has been made before
			mem_alnreg_t *p = &av->a[i];
			int64_t rd;
			int qd, w, max_gap;
			if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb || s->qbeg + s->len > p->qe) continue; // not fully contained
			if (s->len - p->seedlen0 > .1 * l_query) continue; // this seed may give a better alignment
			// qd: distance ahead of the seed on query; rd: on reference
			qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
			max_gap = cal_max_gap(opt, qd < rd? qd : rd); // the maximal gap allowed in regions ahead of the seed
			w = max_gap < p->w? max_gap : p->w; // bounded by the band width
			if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
			// similar to the previous four lines, but this time we look at the region behind
			qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
			max_gap = cal_max_gap(opt, qd < rd? qd : rd);
			w = max_gap < p->w? max_gap : p->w;
			if (qd - rd < w && rd - qd < w) break;
		}
		if (i < av->n) { // the seed is (almost) contained in an existing alignment; further testing is needed to confirm it is not leading to a different aln
			if (bwa_verbose >= 4)
				printf("** Seed(%d) [%ld;%ld,%ld] is almost contained in an existing alignment [%d,%d) <=> [%ld,%ld)\n",
				       k, (long)s->len, (long)s->qbeg, (long)s->rbeg, av->a[i].qb, av->a[i].qe, (long)av->a[i].rb, (long)av->a[i].re);
			for (i = k + 1; i < c->n; ++i) { // check overlapping seeds in the same chain
				const mem_seed_t *t;
				if (srt[i] == 0) continue;
				t = &c->seeds[(uint32_t)srt[i]];
				if (t->len < s->len * .95) continue; // only check overlapping if t is long enough; TODO: more efficient by early stopping
				if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 && t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
				if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 && s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
			}
			if (i == c->n) { // no overlapping seeds; then skip extension
				srt[k] = 0; // mark that seed extension has not been performed
				continue;
			}
			if (bwa_verbose >= 4)
				printf("** Seed(%d) might lead to a different alignment even though it is contained. Extension will be performed.\n", k);
		}

		a = kv_pushp(mem_alnreg_t, *av);
		memset(a, 0, sizeof(mem_alnreg_t));
		a->w = aw[0] = aw[1] = opt->w;
		a->score = a->truesc = -1;
		a->rid = c->rid;

		if (bwa_verbose >= 4) err_printf("** ---> Extending from seed(%d) [%ld;%ld,%ld] @ %s <---\n", k, (long)s->len, (long)s->qbeg, (long)s->rbeg, bns->anns[c->rid].name);
		if (s->qbeg) { // left extension
			uint8_t *rs, *qs;
			int qle, tle, gtle, gscore;
			qs = (uint8_t*) malloc(s->qbeg * sizeof(uint8_t));
			for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
			tmp = s->rbeg - rmax[0];
			rs = (uint8_t *) malloc(tmp * sizeof(uint8_t));
			for (i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i];
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[0] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Left ref:   "); for (j = 0; j < tmp; ++j) putchar("ACGTN"[(int)rs[j]]); putchar('\n');
					printf("*** Left query: "); for (j = 0; j < s->qbeg; ++j) putchar("ACGTN"[(int)qs[j]]); putchar('\n');
				}
				a->score = ksw_extend2(s->qbeg, qs, tmp, rs, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[0], opt->pen_clip5, opt->zdrop, s->len * opt->a, &qle, &tle, &gtle, &gscore, &max_off[0]);
				if (bwa_verbose >= 4) { printf("*** Left extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[0], max_off[0]); fflush(stdout); }
				if (a->score == prev || max_off[0] < (aw[0]>>1) + (aw[0]>>2)) break;
			}
			// check whether we prefer to reach the end of the query
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip5) { // local extension
				a->qb = s->qbeg - qle, a->rb = s->rbeg - tle;
				a->truesc = a->score;
			} else { // to-end extension
				a->qb = 0, a->rb = s->rbeg - gtle;
				a->truesc = gscore;
			}
			free(qs); free(rs);
		} else a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;

		if (s->qbeg + s->len != l_query) { // right extension
			int qle, tle, qe, re, gtle, gscore, sc0 = a->score;
			qe = s->qbeg + s->len;
			re = s->rbeg + s->len - rmax[0];
			assert(re >= 0);
			for (i = 0; i < MAX_BAND_TRY; ++i) {
				int prev = a->score;
				aw[1] = opt->w << i;
				if (bwa_verbose >= 4) {
					int j;
					printf("*** Right ref:   "); for (j = 0; j < rmax[1] - rmax[0] - re; ++j) putchar("ACGTN"[(int)rseq[re+j]]); putchar('\n');
					printf("*** Right query: "); for (j = 0; j < l_query - qe; ++j) putchar("ACGTN"[(int)query[qe+j]]); putchar('\n');
				}
				a->score = ksw_extend2(l_query - qe, query + qe, rmax[1] - rmax[0] - re, rseq + re, 5, opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, aw[1], opt->pen_clip3, opt->zdrop, sc0, &qle, &tle, &gtle, &gscore, &max_off[1]);
				if (bwa_verbose >= 4) { printf("*** Right extension: prev_score=%d; score=%d; bandwidth=%d; max_off_diagonal_dist=%d\n", prev, a->score, aw[1], max_off[1]); fflush(stdout); }
				if (a->score == prev || max_off[1] < (aw[1]>>1) + (aw[1]>>2)) break;
			}
			// similar to the above
			if (gscore <= 0 || gscore <= a->score - opt->pen_clip3) { // local extension
				a->qe = qe + qle, a->re = rmax[0] + re + tle;
				a->truesc += a->score - sc0;
			} else { // to-end extension
				a->qe = l_query, a->re = rmax[0] + re + gtle;
				a->truesc += gscore - sc0;
			}
		} else a->qe = l_query, a->re = s->rbeg + s->len;
		if (bwa_verbose >= 4) printf("*** Added alignment region: [%d,%d) <=> [%ld,%ld); score=%d; {left,right}_bandwidth={%d,%d}\n", a->qb, a->qe, (long)a->rb, (long)a->re, a->score, aw[0], aw[1]);

		// compute seedcov
		for (i = 0, a->seedcov = 0; i < c->n; ++i) {
			const mem_seed_t *t = &c->seeds[i];
			if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe && t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
				a->seedcov += t->len; // this is not very accurate, but for approx. mapQ, this is good enough
		}
		a->w = aw[0] > aw[1]? aw[0] : aw[1];
		a->seedlen0 = s->len;

		a->frac_rep = c->frac_rep;
	}
	free(srt); free(rseq);
}

/******************************
 * De-overlap single-end hits *
 ******************************/

#define alnreg_slt2(a, b) ((a).re < (b).re)
KSORT_INIT(mem_ars2, mem_alnreg_t, alnreg_slt2)

#define alnreg_slt(a, b) ((a).score > (b).score || ((a).score == (b).score && ((a).rb < (b).rb || ((a).rb == (b).rb && (a).qb < (b).qb))))
KSORT_INIT(mem_ars, mem_alnreg_t, alnreg_slt)

#define alnreg_hlt(a, b)  ((a).score > (b).score || ((a).score == (b).score && ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash, mem_alnreg_t, alnreg_hlt)

#define alnreg_hlt2(a, b) ((a).is_alt < (b).is_alt || ((a).is_alt == (b).is_alt && ((a).score > (b).score || ((a).score == (b).score && (a).hash < (b).hash))))
KSORT_INIT(mem_ars_hash2, mem_alnreg_t, alnreg_hlt2)

#define PATCH_MAX_R_BW 0.05f
#define PATCH_MIN_SC_RATIO 0.90f

int mem_patch_reg(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, const mem_alnreg_t *a, const mem_alnreg_t *b, int *_w) {
	int w, score, q_s, r_s;
	double r;
	if (bns == 0 || pac == 0 || query == 0) return 0;
	assert(a->rid == b->rid && a->rb <= b->rb);
	if (a->rb < bns->l_pac && b->rb >= bns->l_pac) return 0; // on different strands
	if (a->qb >= b->qb || a->qe >= b->qe || a->re >= b->re) return 0; // not colinear
	w = (a->re - b->rb) - (a->qe - b->qb); // required bandwidth
	w = w > 0? w : -w; // l = abs(l)
	r = (double)(a->re - b->rb) / (b->re - a->rb) - (double)(a->qe - b->qb) / (b->qe - a->qb); // relative bandwidth
	r = r > 0.? r : -r; // r = fabs(r)
	if (bwa_verbose >= 4)
		printf("* potential hit merge between [%d,%d)<=>[%ld,%ld) and [%d,%d)<=>[%ld,%ld), @ %s; w=%d, r=%.4g\n",
		       a->qb, a->qe, (long)a->rb, (long)a->re, b->qb, b->qe, (long)b->rb, (long)b->re, bns->anns[a->rid].name, w, r);
	if (a->re < b->rb || a->qe < b->qb) { // no overlap on query or on ref
		if (w > opt->w<<1 || r >= PATCH_MAX_R_BW) return 0; // the bandwidth or the relative bandwidth is too large
	} else if (w > opt->w<<2 || r >= PATCH_MAX_R_BW*2) return 0; // more permissive if overlapping on both ref and query
	// global alignment
	w += a->w + b->w;
	w = w < opt->w<<2? w : opt->w<<2;
	if (bwa_verbose >= 4) printf("* test potential hit merge with global alignment; w=%d\n", w);
	bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w, bns->l_pac, pac, b->qe - a->qb, query + a->qb, a->rb, b->re, &score, 0, 0);
	q_s = (int)((double)(b->qe - a->qb) / ((b->qe - b->qb) + (a->qe - a->qb)) * (b->score + a->score) + .499); // predicted score from query
	r_s = (int)((double)(b->re - a->rb) / ((b->re - b->rb) + (a->re - a->rb)) * (b->score + a->score) + .499); // predicted score from ref
	if (bwa_verbose >= 4) printf("* score=%d;(%d,%d)\n", score, q_s, r_s);
	if ((double)score / (q_s > r_s? q_s : r_s) < PATCH_MIN_SC_RATIO) return 0;
	*_w = w;
	return score;
}

int mem_sort_dedup_patch(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, int n, mem_alnreg_t *a) {
	int m, i, j;
	if (n <= 1) return n;
	ks_introsort(mem_ars2, n, a); // sort by the END position, not START!
	for (i = 0; i < n; ++i) a[i].n_comp = 1;
	for (i = 1; i < n; ++i) {
		mem_alnreg_t *p = &a[i];
		if (p->rid != a[i-1].rid || p->rb >= a[i-1].re + opt->max_chain_gap) continue; // then no need to go into the loop below
		for (j = i - 1; j >= 0 && p->rid == a[j].rid && p->rb < a[j].re + opt->max_chain_gap; --j) {
			mem_alnreg_t *q = &a[j];
			int64_t Or, oq, mr, mq;
			int score, w;
			if (q->qe == q->qb) continue; // a[j] has been excluded
			Or = q->re - p->rb; // overlap length on the reference
			oq = q->qb < p->qb? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
			mr = q->re - q->rb < p->re - p->rb? q->re - q->rb : p->re - p->rb; // min ref len in alignment
			mq = q->qe - q->qb < p->qe - p->qb? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
			if (Or > opt->mask_level_redun * mr && oq > opt->mask_level_redun * mq) { // one of the hits is redundant
				if (p->score < q->score) {
					p->qe = p->qb;
					break;
				} else q->qe = q->qb;
			} else if (q->rb < p->rb && (score = mem_patch_reg(opt, bns, pac, query, q, p, &w)) > 0) { // then merge q into p
				p->n_comp += q->n_comp + 1;
				p->seedcov = p->seedcov > q->seedcov? p->seedcov : q->seedcov;
				p->sub = p->sub > q->sub? p->sub : q->sub;
				p->csub = p->csub > q->csub? p->csub : q->csub;
				p->qb = q->qb, p->rb = q->rb;
				p->truesc = p->score = score;
				p->w = w;
				q->qb = q->qe;
			}
		}
	}
	for (i = 0, m = 0; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	n = m;
	ks_introsort(mem_ars, n, a);
	for (i = 1; i < n; ++i) { // mark identical hits
		if (a[i].score == a[i-1].score && a[i].rb == a[i-1].rb && a[i].qb == a[i-1].qb)
			a[i].qe = a[i].qb;
	}
	for (i = 1, m = 1; i < n; ++i) // exclude identical hits
		if (a[i].qe > a[i].qb) {
			if (m != i) a[m++] = a[i];
			else ++m;
		}
	return m;
}

/************************
 * Integrated interface *
 ************************/

int mem_approx_mapq_se(const mem_opt_t *opt, const mem_alnreg_t *a) {
	int mapq, l, sub = a->sub? a->sub : opt->min_seed_len * opt->a;
	double identity;
	sub = a->csub > sub? a->csub : sub;
	if (sub >= a->score) return 0;
	l = a->qe - a->qb > a->re - a->rb? a->qe - a->qb : a->re - a->rb;
	identity = 1. - (double)(l * opt->a - a->score) / (opt->a + opt->b) / l;
	if (a->score == 0) {
		mapq = 0;
	} else if (opt->mapQ_coef_len > 0) {
		double tmp;
		tmp = l < opt->mapQ_coef_len? 1. : opt->mapQ_coef_fac / log(l);
		tmp *= identity * identity;
		mapq = (int)(6.02 * (a->score - sub) / opt->a * tmp * tmp + .499);
	} else {
		mapq = (int)(MEM_MAPQ_COEF * (1. - (double)sub / a->score) * log(a->seedcov) + .499);
		mapq = identity < 0.95? (int)(mapq * identity * identity + .499) : mapq;
	}
	if (a->sub_n > 0) mapq -= (int)(4.343 * log(a->sub_n+1) + .499);
	if (mapq > 60) mapq = 60;
	if (mapq < 0) mapq = 0;
	mapq = (int)(mapq * (1. - a->frac_rep) + .499);
	return mapq;
}

void mem_mark_primary_se_core(const mem_opt_t *opt, int n, mem_alnreg_t *a, int_v *z) {
	int i, k, tmp;
	tmp = opt->a + opt->b;
	tmp = opt->o_del + opt->e_del > tmp? opt->o_del + opt->e_del : tmp;
	tmp = opt->o_ins + opt->e_ins > tmp? opt->o_ins + opt->e_ins : tmp;
	z->n = 0;
	kv_push(int, *z, 0);
	for (i = 1; i < n; ++i) {
		for (k = 0; k < z->n; ++k) {
			int j = z->a[k];
			int b_max = a[j].qb > a[i].qb? a[j].qb : a[i].qb;
			int e_min = a[j].qe < a[i].qe? a[j].qe : a[i].qe;
			if (e_min > b_max) { // have overlap
				int min_l = a[i].qe - a[i].qb < a[j].qe - a[j].qb? a[i].qe - a[i].qb : a[j].qe - a[j].qb;
				if (e_min - b_max >= min_l * opt->mask_level) { // significant overlap
					if (a[j].sub == 0) a[j].sub = a[i].score;
					if (a[j].score - a[i].score <= tmp && (a[j].is_alt || !a[i].is_alt))
						++a[j].sub_n;
					break;
				}
			}
		}
		if (k == z->n) kv_push(int, *z, i);
		else a[i].secondary = z->a[k];
	}
}

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id) {
	int i, n_pri;
	int_v z = {0,0,0};
	if (n == 0) return 0;
	for (i = n_pri = 0; i < n; ++i) {
		a[i].sub = a[i].alt_sc = 0, a[i].secondary = a[i].secondary_all = -1, a[i].hash = hash_64(id+i);
		if (!a[i].is_alt) ++n_pri;
	}
	ks_introsort(mem_ars_hash, n, a);
	mem_mark_primary_se_core(opt, n, a, &z);
	for (i = 0; i < n; ++i) {
		mem_alnreg_t *p = &a[i];
		p->secondary_all = i; // keep the rank in the first round
		if (!p->is_alt && p->secondary >= 0 && a[p->secondary].is_alt)
			p->alt_sc = a[p->secondary].score;
	}
	if (n_pri >= 0 && n_pri < n) {
		kv_resize(int, z, n);
		if (n_pri > 0) ks_introsort(mem_ars_hash2, n, a);
		for (i = 0; i < n; ++i) z.a[a[i].secondary_all] = i;
		for (i = 0; i < n; ++i) {
			if (a[i].secondary >= 0) {
				a[i].secondary_all = z.a[a[i].secondary];
				if (a[i].is_alt) a[i].secondary = INT_MAX;
			} else a[i].secondary_all = -1;
		}
		if (n_pri > 0) { // mark primary for hits to the primary assembly only
			for (i = 0; i < n_pri; ++i) a[i].sub = 0, a[i].secondary = -1;
			mem_mark_primary_se_core(opt, n_pri, a, &z);
		}
	} else {
		for (i = 0; i < n; ++i)
			a[i].secondary_all = a[i].secondary;
	}
	free(z.a);
	return n_pri;
}

void mem_reorder_primary5(int T, mem_alnreg_v *a) {
	int k, n_pri = 0, left_st = INT_MAX, left_k = -1;
	mem_alnreg_t t;
	for (k = 0; k < a->n; ++k)
		if (a->a[k].secondary < 0 && !a->a[k].is_alt && a->a[k].score >= T) ++n_pri;
	if (n_pri <= 1) return; // only one alignment
	for (k = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		if (p->secondary >= 0 || p->is_alt || p->score < T) continue;
		if (p->qb < left_st) left_st = p->qb, left_k = k;
	}
	assert(a->a[0].secondary < 0);
	if (left_k == 0) return; // no need to reorder
	t = a->a[0], a->a[0] = a->a[left_k], a->a[left_k] = t;
	for (k = 1; k < a->n; ++k) { // update secondary and secondary_all
		mem_alnreg_t *p = &a->a[k];
		if (p->secondary == 0) p->secondary = left_k;
		else if (p->secondary == left_k) p->secondary = 0;
		if (p->secondary_all == 0) p->secondary_all = left_k;
		else if (p->secondary_all == left_k) p->secondary_all = 0;
	}
}

/*****************************
 * Basic hit->SAM conversion *
 *****************************/

static inline int infer_bw(int l1, int l2, int score, int a, int q, int r) {
	int w;
	if (l1 == l2 && l1 * a - score < (q + r - a)<<1) return 0; // to get equal alignment length, we need at least two gaps
	w = ((double)((l1 < l2? l1 : l2) * a - score - q) / r + 2.);
	if (w < abs(l1 - l2)) w = abs(l1 - l2);
	return w;
}

mem_aln_t mem_reg2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const char *query_, const mem_alnreg_t *ar) {
	mem_aln_t a;
	int i, w2, tmp, qb, qe, NM, score, is_rev, last_sc = -(1<<30), l_MD;
	int64_t pos, rb, re;
	uint8_t *query;

	memset(&a, 0, sizeof(mem_aln_t));
	if (ar == 0 || ar->rb < 0 || ar->re < 0) { // generate an unmapped record
		a.rid = -1; a.pos = -1; a.flag |= 0x4;
		return a;
	}
	qb = ar->qb, qe = ar->qe;
	rb = ar->rb, re = ar->re;
	query = (uint8_t*)malloc(l_query * sizeof(uint8_t));
	for (i = 0; i < l_query; ++i) // convert to the nt4 encoding
		query[i] = query_[i] < 5? query_[i] : nst_nt4_table[(int)query_[i]];
	a.mapq = ar->secondary < 0? mem_approx_mapq_se(opt, ar) : 0;
	if (ar->secondary >= 0) a.flag |= 0x100; // secondary alignment
	tmp = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_del, opt->e_del);
	w2  = infer_bw(qe - qb, re - rb, ar->truesc, opt->a, opt->o_ins, opt->e_ins);
	w2 = w2 > tmp? w2 : tmp;
	if (bwa_verbose >= 4) printf("* Band width: inferred=%d, cmd_opt=%d, alnreg=%d\n", w2, opt->w, ar->w);
	if (w2 > opt->w) w2 = w2 < ar->w? w2 : ar->w;
	i = 0; a.cigar = 0;
	do {
		free(a.cigar);
		w2 = w2 < opt->w<<2? w2 : opt->w<<2;
		a.cigar = bwa_gen_cigar2(opt->mat, opt->o_del, opt->e_del, opt->o_ins, opt->e_ins, w2, bns->l_pac, pac, qe - qb, (uint8_t*)&query[qb], rb, re, &score, &a.n_cigar, &NM);
		if (bwa_verbose >= 4) printf("* Final alignment: w2=%d, global_sc=%d, local_sc=%d\n", w2, score, ar->truesc);
		if (score == last_sc || w2 == opt->w<<2) break; // it is possible that global alignment and local alignment give different scores
		last_sc = score;
		w2 <<= 1;
	} while (++i < 3 && score < ar->truesc - opt->a);
	l_MD = strlen((char*)(a.cigar + a.n_cigar)) + 1;
	a.NM = NM;
	pos = bns_depos(bns, rb < bns->l_pac? rb : re - 1, &is_rev);
	a.is_rev = is_rev;
	if (a.n_cigar > 0) { // squeeze out leading or trailing deletions
		if ((a.cigar[0]&0xf) == 2) {
			pos += a.cigar[0]>>4;
			--a.n_cigar;
			memmove(a.cigar, a.cigar + 1, a.n_cigar * 4 + l_MD);
		} else if ((a.cigar[a.n_cigar-1]&0xf) == 2) {
			--a.n_cigar;
			memmove(a.cigar + a.n_cigar, a.cigar + a.n_cigar + 1, l_MD); // MD needs to be moved accordingly
		}
	}
	if (qb != 0 || qe != l_query) { // add clipping to CIGAR
		int clip5, clip3;
		clip5 = is_rev? l_query - qe : qb;
		clip3 = is_rev? qb : l_query - qe;
		a.cigar = (uint32_t*) realloc(a.cigar, sizeof(uint32_t) * (a.n_cigar + 2) + l_MD);
		if (clip5) {
			memmove(a.cigar+1, a.cigar, a.n_cigar * 4 + l_MD); // make room for 5'-end clipping
			a.cigar[0] = clip5<<4 | 3;
			++a.n_cigar;
		}
		if (clip3) {
			memmove(a.cigar + a.n_cigar + 1, a.cigar + a.n_cigar, l_MD); // make room for 3'-end clipping
			a.cigar[a.n_cigar++] = clip3<<4 | 3;
		}
	}
	a.rid = bns_pos2rid(bns, pos);
	assert(a.rid == ar->rid);
	a.pos = pos - bns->anns[a.rid].offset;
	a.score = ar->score; a.sub = ar->sub > ar->csub? ar->sub : ar->csub;
	a.is_alt = ar->is_alt; a.alt_sc = ar->alt_sc;
	free(query);
	return a;
}

static inline int get_rlen(int n_cigar, const uint32_t *cigar) {
	int k, l;
	for (k = l = 0; k < n_cigar; ++k) {
		int op = cigar[k]&0xf;
		if (op == 0 || op == 2)
			l += cigar[k]>>4;
	}
	return l;
}

static inline void add_cigar(const mem_opt_t *opt, mem_aln_t *p, kstring_t *str, int which) {
	int i;
	if (p->n_cigar) { // aligned
		for (i = 0; i < p->n_cigar; ++i) {
			int c = p->cigar[i]&0xf;
			if (!(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt && (c == 3 || c == 4))
				c = which? 4 : 3; // use hard clipping for supplementary alignments
			kputw(p->cigar[i]>>4, str); kputc("MIDSH"[c], str);
		}
	} else kputc('*', str); // having a coordinate but unaligned (e.g. when copy_mate is true)
}

void mem_aln2sam(const mem_opt_t *opt, const bntseq_t *bns, kstring_t *str, bseq1_t *s, int n, const mem_aln_t *list, int which, const mem_aln_t *m_) {
	int i, l_name;
	mem_aln_t ptmp = list[which], *p = &ptmp, mtmp, *m = 0; // make a copy of the alignment to convert

	if (m_) mtmp = *m_, m = &mtmp;
	// set flag
	p->flag |= m? 0x1 : 0; // is paired in sequencing
	p->flag |= p->rid < 0? 0x4 : 0; // is mapped
	p->flag |= m && m->rid < 0? 0x8 : 0; // is mate mapped
	if (p->rid < 0 && m && m->rid >= 0) // copy mate to alignment
		p->rid = m->rid, p->pos = m->pos, p->is_rev = m->is_rev, p->n_cigar = 0;
	if (m && m->rid < 0 && p->rid >= 0) // copy alignment to mate
		m->rid = p->rid, m->pos = p->pos, m->is_rev = p->is_rev, m->n_cigar = 0;
	p->flag |= p->is_rev? 0x10 : 0; // is on the reverse strand
	p->flag |= m && m->is_rev? 0x20 : 0; // is mate on the reverse strand

	// print up to CIGAR
	l_name = strlen(s->name);
	ks_resize(str, str->l + s->l_seq + l_name + (s->qual? s->l_seq : 0) + 20);
	kputsn(s->name, l_name, str); kputc('\t', str); // QNAME
	kputw((p->flag&0xffff) | (p->flag&0x10000? 0x100 : 0), str); kputc('\t', str); // FLAG
	if (p->rid >= 0) { // with coordinate
		kputs(bns->anns[p->rid].name, str); kputc('\t', str); // RNAME
		kputl(p->pos + 1, str); kputc('\t', str); // POS
		kputw(p->mapq, str); kputc('\t', str); // MAPQ
		add_cigar(opt, p, str, which);
	} else kputsn("*\t0\t0\t*", 7, str); // without coordinte
	kputc('\t', str);

	// print the mate position if applicable
	if (m && m->rid >= 0) {
		if (p->rid == m->rid) kputc('=', str);
		else kputs(bns->anns[m->rid].name, str);
		kputc('\t', str);
		kputl(m->pos + 1, str); kputc('\t', str);
		if (p->rid == m->rid) {
			int64_t p0 = p->pos + (p->is_rev? get_rlen(p->n_cigar, p->cigar) - 1 : 0);
			int64_t p1 = m->pos + (m->is_rev? get_rlen(m->n_cigar, m->cigar) - 1 : 0);
			if (m->n_cigar == 0 || p->n_cigar == 0) kputc('0', str);
			else kputl(-(p0 - p1 + (p0 > p1? 1 : p0 < p1? -1 : 0)), str);
		} else kputc('0', str);
	} else kputsn("*\t0\t0", 5, str);
	kputc('\t', str);

	// print SEQ and QUAL
	if (p->flag & 0x100) { // for secondary alignments, don't write SEQ and QUAL
		kputsn("*\t*", 3, str);
	} else if (!p->is_rev) { // the forward strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) { // have cigar && not the primary alignment && not softclip all
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qb += p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qe -= p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qb; i < qe; ++i) str->s[str->l++] = "ACGTN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qb; i < qe; ++i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	} else { // the reverse strand
		int i, qb = 0, qe = s->l_seq;
		if (p->n_cigar && which && !(opt->flag&MEM_F_SOFTCLIP) && !p->is_alt) {
			if ((p->cigar[0]&0xf) == 4 || (p->cigar[0]&0xf) == 3) qe -= p->cigar[0]>>4;
			if ((p->cigar[p->n_cigar-1]&0xf) == 4 || (p->cigar[p->n_cigar-1]&0xf) == 3) qb += p->cigar[p->n_cigar-1]>>4;
		}
		ks_resize(str, str->l + (qe - qb) + 1);
		for (i = qe-1; i >= qb; --i) str->s[str->l++] = "TGCAN"[(int)s->seq[i]];
		kputc('\t', str);
		if (s->qual) { // printf qual
			ks_resize(str, str->l + (qe - qb) + 1);
			for (i = qe-1; i >= qb; --i) str->s[str->l++] = s->qual[i];
			str->s[str->l] = 0;
		} else kputc('*', str);
	}

	// print optional tags
	if (p->n_cigar) {
		kputsn("\tNM:i:", 6, str); kputw(p->NM, str);
		kputsn("\tMD:Z:", 6, str); kputs((char*)(p->cigar + p->n_cigar), str);
	}
	if (m && m->n_cigar) { kputsn("\tMC:Z:", 6, str); add_cigar(opt, m, str, which); }
	if (p->score >= 0) { kputsn("\tAS:i:", 6, str); kputw(p->score, str); }
	if (p->sub >= 0) { kputsn("\tXS:i:", 6, str); kputw(p->sub, str); }
	if (bwa_rg_id[0]) { kputsn("\tRG:Z:", 6, str); kputs(bwa_rg_id, str); }
	if (!(p->flag & 0x100)) { // not multi-hit
		for (i = 0; i < n; ++i)
			if (i != which && !(list[i].flag&0x100)) break;
		if (i < n) { // there are other primary hits; output them
			kputsn("\tSA:Z:", 6, str);
			for (i = 0; i < n; ++i) {
				const mem_aln_t *r = &list[i];
				int k;
				if (i == which || (r->flag&0x100)) continue; // proceed if: 1) different from the current; 2) not shadowed multi hit
				kputs(bns->anns[r->rid].name, str); kputc(',', str);
				kputl(r->pos+1, str); kputc(',', str);
				kputc("+-"[r->is_rev], str); kputc(',', str);
				for (k = 0; k < r->n_cigar; ++k) {
					kputw(r->cigar[k]>>4, str); kputc("MIDSH"[r->cigar[k]&0xf], str);
				}
				kputc(',', str); kputw(r->mapq, str);
				kputc(',', str); kputw(r->NM, str);
				kputc(';', str);
			}
		}
		if (p->alt_sc > 0)
			ksprintf(str, "\tpa:f:%.3f", (double)p->score / p->alt_sc);
	}
	if (p->XA) { kputsn("\tXA:Z:", 6, str); kputs(p->XA, str); }
	if (s->comment) { kputc('\t', str); kputs(s->comment, str); }
	if ((opt->flag&MEM_F_REF_HDR) && p->rid >= 0 && bns->anns[p->rid].anno != 0 && bns->anns[p->rid].anno[0] != 0) {
		int tmp;
		kputsn("\tXR:Z:", 6, str);
		tmp = str->l;
		kputs(bns->anns[p->rid].anno, str);
		for (i = tmp; i < str->l; ++i) // replace TAB in the comment to SPACE
			if (str->s[i] == '\t') str->s[i] = ' ';
	}
	kputc('\n', str);
}

static inline int get_pri_idx(double XA_drop_ratio, const mem_alnreg_t *a, int i)
{
	int k = a[i].secondary_all;
	if (k >= 0 && a[i].score >= a[k].score * XA_drop_ratio) return k;
	return -1;
}

// Okay, returning strings is bad, but this has happened a lot elsewhere. If I have time, I need serious code cleanup.
char **mem_gen_alt(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, const mem_alnreg_v *a, int l_query, const char *query) // ONLY work after mem_mark_primary_se()
{
	int i, k, r, *cnt, tot;
	kstring_t *aln = 0, str = {0,0,0};
	char **XA = 0, *has_alt;

	cnt = (int*) calloc(a->n, sizeof(int));
	has_alt = (char*) calloc(a->n, sizeof(char));
	for (i = 0, tot = 0; i < a->n; ++i) {
		r = get_pri_idx(opt->XA_drop_ratio, a->a, i);
		if (r >= 0) {
			++cnt[r], ++tot;
			if (a->a[i].is_alt) has_alt[r] = 1;
		}
	}
	if (tot == 0) goto end_gen_alt;
	aln = (kstring_t*) calloc(a->n, sizeof(kstring_t));
	for (i = 0; i < a->n; ++i) {
		mem_aln_t t;
		if ((r = get_pri_idx(opt->XA_drop_ratio, a->a, i)) < 0) continue;
		if (cnt[r] > opt->max_XA_hits_alt || (!has_alt[r] && cnt[r] > opt->max_XA_hits)) continue;
		t = mem_reg2aln(opt, bns, pac, l_query, query, &a->a[i]);
		str.l = 0;
		kputs(bns->anns[t.rid].name, &str);
		kputc(',', &str); kputc("+-"[t.is_rev], &str); kputl(t.pos + 1, &str);
		kputc(',', &str);
		for (k = 0; k < t.n_cigar; ++k) {
			kputw(t.cigar[k]>>4, &str);
			kputc("MIDSHN"[t.cigar[k]&0xf], &str);
		}
		kputc(',', &str); kputw(t.NM, &str);
		kputc(';', &str);
		free(t.cigar);
		kputsn(str.s, str.l, &aln[r]);
	}
	XA = (char**) calloc(a->n, sizeof(char*));
	for (k = 0; k < a->n; ++k)
		XA[k] = aln[k].s;

	end_gen_alt:
	free(has_alt); free(cnt); free(aln); free(str.s);
	return XA;
}

// TODO (future plan): group hits into a uint64_t[] array. This will be cleaner and more flexible
void mem_reg2sam(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, bseq1_t *s, mem_alnreg_v *a, int extra_flag, const mem_aln_t *m) {
	kstring_t str;
	kvec_t(mem_aln_t) aa;
	int k, l;
	char **XA = 0;

	if (!(opt->flag & MEM_F_ALL))
		XA = mem_gen_alt(opt, bns, pac, a, s->l_seq, s->seq);
	kv_init(aa);
	str.l = str.m = 0; str.s = 0;
	for (k = l = 0; k < a->n; ++k) {
		mem_alnreg_t *p = &a->a[k];
		mem_aln_t *q;
		if (p->score < opt->T) continue;
		if (p->secondary >= 0 && (p->is_alt || !(opt->flag&MEM_F_ALL))) continue;
		if (p->secondary >= 0 && p->secondary < INT_MAX && p->score < a->a[p->secondary].score * opt->drop_ratio) continue;
		q = kv_pushp(mem_aln_t, aa);
		*q = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, p);
		assert(q->rid >= 0); // this should not happen with the new code
		q->XA = XA? XA[k] : 0;
		q->flag |= extra_flag; // flag secondary
		if (p->secondary >= 0) q->sub = -1; // don't output sub-optimal score
		if (l && p->secondary < 0) // if supplementary
			q->flag |= (opt->flag&MEM_F_NO_MULTI)? 0x10000 : 0x800;
		if (!(opt->flag & MEM_F_KEEP_SUPP_MAPQ) && l && !p->is_alt && q->mapq > aa.a[0].mapq)
			q->mapq = aa.a[0].mapq; // lower mapq for supplementary mappings, unless -5 or -q is applied
		++l;
	}
	if (aa.n == 0) { // no alignments good enough; then write an unaligned record
		mem_aln_t t;
		t = mem_reg2aln(opt, bns, pac, s->l_seq, s->seq, 0);
		t.flag |= extra_flag;
		mem_aln2sam(opt, bns, &str, s, 1, &t, 0, m);
	} else {
		for (k = 0; k < aa.n; ++k)
			mem_aln2sam(opt, bns, &str, s, aa.n, aa.a, k, m);
		for (k = 0; k < aa.n; ++k) free(aa.a[k].cigar);
		free(aa.a);
	}
	s->sam = str.s;
	if (XA) {
		for (k = 0; k < a->n; ++k) free(XA[k]);
		free(XA);
	}
}

typedef struct {
	SeqPair *seqPairArrayAux[MAX_THREADS];
	SeqPair *seqPairArrayLeft128[MAX_THREADS];
	SeqPair *seqPairArrayRight128[MAX_THREADS];

	int64_t wsize[MAX_THREADS];

	int64_t wsize_buf_ref[MAX_THREADS*CACHE_LINE];
	int64_t wsize_buf_qer[MAX_THREADS*CACHE_LINE];

	uint8_t *seqBufLeftRef[MAX_THREADS*CACHE_LINE];
	uint8_t *seqBufRightRef[MAX_THREADS*CACHE_LINE];
	uint8_t *seqBufLeftQer[MAX_THREADS*CACHE_LINE];
	uint8_t *seqBufRightQer[MAX_THREADS*CACHE_LINE];

//	SMEM *matchArray[MAX_THREADS];
//	int32_t *min_intv_ar[MAX_THREADS];
//	int32_t *rid[MAX_THREADS];
	int32_t *lim[MAX_THREADS];
//	int16_t *query_pos_ar[MAX_THREADS];
//	uint8_t *enc_qdb[MAX_THREADS];

	int64_t wsize_mem[MAX_THREADS];
} mem_cache;

typedef struct {
	int n;
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	thread_aux_t *aux;
	bseq1_t *seqs;
	int64_t n_processed;
	mem_cache mmc;
} worker_t;

void* _mm_realloc(void *ptr, int64_t csize, int64_t nsize, int16_t dsize) {
	if (nsize <= csize)
	{
		fprintf(stderr, "Shringking not supported yet.\n");
		return ptr;
	}
	void *nptr = _mm_malloc(nsize * dsize, 64);
	assert(nptr != NULL);
	memcpy_bwamem(nptr, nsize * dsize, ptr, csize, (char*)__FILE__, __LINE__);
	_mm_free(ptr);

	return nptr;
}

inline void sortPairsLenExt(SeqPair *pairArray, int32_t count, SeqPair *tempArray,
                            int32_t *hist, int &numPairs128, int &numPairs16,
                            int &numPairs1, int score_a)
{
	int32_t i;
	numPairs128 = numPairs16 = numPairs1 = 0;

	int32_t *hist2 = hist + MAX_SEQ_LEN8;
	int32_t *hist3 = hist + MAX_SEQ_LEN8 + MAX_SEQ_LEN16;

	for(i = 0; i <= MAX_SEQ_LEN8 + MAX_SEQ_LEN16; i+=1)
		//_mm256_store_si256((__m256i *)(hist + i), zero256);
		hist[i] = 0;

	int *arr = (int*) calloc (count, sizeof(int));
	assert(arr != NULL);

	for(i = 0; i < count; i++)
	{
		SeqPair sp = pairArray[i];
		// int minval = sp.h0 + max_(sp.len1, sp.len2);
		int minval = sp.h0 + min_(sp.len1, sp.len2) * score_a;
		if (sp.len1 < MAX_SEQ_LEN8 && sp.len2 < MAX_SEQ_LEN8 && minval < MAX_SEQ_LEN8)
			hist[minval]++;
		else if(sp.len1 < MAX_SEQ_LEN16 && sp.len2 < MAX_SEQ_LEN16 && minval < MAX_SEQ_LEN16)
			hist2[minval] ++;
		else
			hist3[0] ++;

		arr[i] = 0;
	}

	int32_t cumulSum = 0;
	for(i = 0; i < MAX_SEQ_LEN8; i++)
	{
		int32_t cur = hist[i];
		hist[i] = cumulSum;
		cumulSum += cur;
	}
	for(i = 0; i < MAX_SEQ_LEN16; i++)
	{
		int32_t cur = hist2[i];
		hist2[i] = cumulSum;
		cumulSum += cur;
	}
	hist3[0] = cumulSum;

	for(i = 0; i < count; i++)
	{
		SeqPair sp = pairArray[i];
		// int minval = sp.h0 + max_(sp.len1, sp.len2);
		int minval = sp.h0 + min_(sp.len1, sp.len2) * score_a;

		if (sp.len1 < MAX_SEQ_LEN8 && sp.len2 < MAX_SEQ_LEN8 && minval < MAX_SEQ_LEN8)
		{
			int32_t pos = hist[minval];
			tempArray[pos] = sp;
			hist[minval]++;
			numPairs128 ++;
			if (arr[pos] != 0)
			{
				fprintf(stderr, "[%s] Error log: repeat, pos: %d, arr: %d, minval: %d, (%d %d)\n",
				        __func__, pos, arr[pos], minval, sp.len1, sp.len2);
				exit(EXIT_FAILURE);
			}
			arr[pos] = 1;
		}
		else if (sp.len1 < MAX_SEQ_LEN16 && sp.len2 < MAX_SEQ_LEN16 && minval < MAX_SEQ_LEN16) {
			int32_t pos = hist2[minval];
			tempArray[pos] = sp;
			hist2[minval]++;
			numPairs16 ++;
			if (arr[pos] != 0)
			{
				SeqPair spt = pairArray[arr[pos]-1];
				fprintf(stderr, "[%s] Error log: repeat, "
				                "i: %d, pos: %d, arr: %d, hist2: %d, minval: %d, (%d %d %d) (%d %d %d)\n",
				        __func__, i, pos, arr[pos], hist2[minval],  minval, sp.h0, sp.len1, sp.len2,
				        spt.h0, spt.len1, spt.len2);
				exit(EXIT_FAILURE);
			}
			arr[pos] = i + 1;
		}
		else {
			int32_t pos = hist3[0];
			tempArray[pos] = sp;
			hist3[0]++;
			arr[pos] = i + 1;
			numPairs1 ++;
		}
	}

	for(i = 0; i < count; i++) {
		pairArray[i] = tempArray[i];
	}

	free(arr);
}

inline void sortPairsLen(SeqPair *pairArray, int32_t count, SeqPair *tempArray, int32_t *hist)
{

	int32_t i;
#if ((!__AVX512BW__) & (__AVX2__ | __SSE2__))
	for(i = 0; i <= MAX_SEQ_LEN16; i++) hist[i] = 0;
#else
	__m512i zero512 = _mm512_setzero_si512();
    for(i = 0; i <= MAX_SEQ_LEN16; i+=16)
    {
        _mm512_store_si512((__m512i *)(hist + i), zero512);
    }
#endif

	for(i = 0; i < count; i++)
	{
		SeqPair sp = pairArray[i];
		hist[sp.len1]++;
	}
	int32_t cumulSum = 0;
	for(i = 0; i <= MAX_SEQ_LEN16; i++)
	{
		int32_t cur = hist[i];
		hist[i] = cumulSum;
		cumulSum += cur;
	}

	for(i = 0; i < count; i++)
	{
		SeqPair sp = pairArray[i];
		int32_t pos = hist[sp.len1];

		tempArray[pos] = sp;
		hist[sp.len1]++;
	}

	for(i = 0; i < count; i++) {
		pairArray[i] = tempArray[i];
	}
}

/* Restructured BSW parent function */
#define FAC 8
#define PFD 2
void mem_chain2aln_across_reads_V2(const mem_opt_t *opt, const bntseq_t *bns,
                                   const uint8_t *pac, bseq1_t *seq_, int nseq,
                                   mem_chain_v* chain_ar, mem_alnreg_v *av_v,
                                   mem_cache *mmc, int tid)
{
	SeqPair *seqPairArrayAux      = mmc->seqPairArrayAux[tid];
	SeqPair *seqPairArrayLeft128  = mmc->seqPairArrayLeft128[tid];
	SeqPair *seqPairArrayRight128 = mmc->seqPairArrayRight128[tid];
	int64_t *wsize_pair = &(mmc->wsize[tid]);

	uint8_t *seqBufLeftRef  = mmc->seqBufLeftRef[tid*CACHE_LINE];
	uint8_t *seqBufRightRef = mmc->seqBufRightRef[tid*CACHE_LINE];
	uint8_t *seqBufLeftQer  = mmc->seqBufLeftQer[tid*CACHE_LINE];
	uint8_t *seqBufRightQer = mmc->seqBufRightQer[tid*CACHE_LINE];
	int64_t *wsize_buf_ref = &(mmc->wsize_buf_ref[tid*CACHE_LINE]);
	int64_t *wsize_buf_qer = &(mmc->wsize_buf_qer[tid*CACHE_LINE]);

	// int32_t *lim_g = mmc->lim + (BATCH_SIZE + 32) * tid;
	int32_t *lim_g = mmc->lim[tid];

	mem_seed_t *s;
	int64_t l_pac = bns->l_pac, rmax[8] __attribute__((aligned(64)));
	// std::vector<int8_t> nexitv(nseq, 0);

	int numPairsLeft = 0, numPairsRight = 0;
	int numPairsLeft1 = 0, numPairsRight1 = 0;
	int numPairsLeft128 = 0, numPairsRight128 = 0;
	int numPairsLeft16 = 0, numPairsRight16 = 0;

	int64_t leftRefOffset = 0, rightRefOffset = 0;
	int64_t leftQerOffset = 0, rightQerOffset = 0;

	int srt_size = MAX_SEEDS_PER_READ, fac = FAC;
	uint64_t *srt = (uint64_t *) malloc(srt_size * 8);
	uint32_t *srtgg = (uint32_t*) malloc(nseq * SEEDS_PER_READ * fac * sizeof(uint32_t));

	int spos = 0;

	// uint64_t timUP = __rdtsc();
	for (int l=0; l<nseq; l++)
	{
		int max = 0;
		uint8_t *rseq = 0;

		uint32_t *srtg = srtgg;
		lim_g[l+1] = 0;

		const uint8_t *query = (uint8_t *) seq_[l].seq;
		int l_query = seq_[l].l_seq;

		mem_chain_v *chn = &chain_ar[l];
		mem_alnreg_v *av = &av_v[l];  // alignment
		mem_chain_t *c;

		_mm_prefetch((const char*) query, _MM_HINT_NTA);

		// aln mem allocation
		av->m = 0;
		for (int j=0; j<chn->n; j++) {
			c = &chn->a[j]; av->m += c->n;
		}
		av->a = (mem_alnreg_t*)calloc(av->m, sizeof(mem_alnreg_t));

		// aln mem allocation ends
		for (int j=0; j<chn->n; j++)
		{
			c = &chn->a[j];
//			assert(c->seqid == l);

			int64_t tmp = 0;
			if (c->n == 0) continue;

			_mm_prefetch((const char*) (srtg + spos + 64), _MM_HINT_NTA);
			_mm_prefetch((const char*) (lim_g), _MM_HINT_NTA);

			// get the max possible span
			rmax[0] = l_pac<<1; rmax[1] = 0;

			for (int i = 0; i < c->n; ++i) {
				int64_t b, e;
				const mem_seed_t *t = &c->seeds[i];
				b = t->rbeg - (t->qbeg + cal_max_gap(opt, t->qbeg));
				e = t->rbeg + t->len + ((l_query - t->qbeg - t->len) +
				                        cal_max_gap(opt, l_query - t->qbeg - t->len));

				tmp = rmax[0];
				rmax[0] = tmp < b? rmax[0] : b;
				rmax[1] = (rmax[1] > e)? rmax[1] : e;
				if (t->len > max) max = t->len;
			}

			rmax[0] = rmax[0] > 0? rmax[0] : 0;
			rmax[1] = rmax[1] < l_pac<<1? rmax[1] : l_pac<<1;
			if (rmax[0] < l_pac && l_pac < rmax[1])
			{
				if (c->seeds[0].rbeg < l_pac) rmax[1] = l_pac;
				else rmax[0] = l_pac;
			}

			/* retrieve the reference sequence */
			{
				int rid = 0;
				// free rseq
				rseq = bns_fetch_seq(bns, pac, &rmax[0], c->seeds[0].rbeg, &rmax[1], &rid);
//				rseq = bns_fetch_seq_v2(bns, pac, &rmax[0],
//				                        c->seeds[0].rbeg,
//				                        &rmax[1], &rid, ref_string,
//				                        (uint8_t*) seqPairArrayAux);
				assert(c->rid == rid);
			}

			_mm_prefetch((const char*) rseq, _MM_HINT_NTA);
			// _mm_prefetch((const char*) rseq + 64, _MM_HINT_NTA);

			// assert(c->n < MAX_SEEDS_PER_READ);  // temp
			if (c->n > srt_size) {
				srt_size = c->n + 10;
				srt = (uint64_t *) realloc(srt, srt_size * 8);
			}

			for (int i = 0; i < c->n; ++i)
				srt[i] = (uint64_t)c->seeds[i].score<<32 | i;

			if (c->n > 1)
				ks_introsort_64(c->n, srt);

			// assert((spos + c->n) < SEEDS_PER_READ * FAC * nseq);
			if ((spos + c->n) > SEEDS_PER_READ * fac * nseq) {
				fac <<= 1;
				srtgg = (uint32_t *) realloc(srtgg, nseq * SEEDS_PER_READ * fac * sizeof(uint32_t));
			}

			for (int i = 0; i < c->n; ++i)
				srtg[spos++] = srt[i];

			lim_g[l+1] += c->n;

			// uint64_t tim = __rdtsc();
			for (int k=c->n-1; k >= 0; k--)
			{
				s = &c->seeds[(uint32_t)srt[k]];

				mem_alnreg_t *a;
				// a = kv_pushp(mem_alnreg_t, *av);
				a = &av->a[av->n++];
				memset(a, 0, sizeof(mem_alnreg_t));

				s->aln = av->n-1;

				a->w = opt->w;
				a->score = a->truesc = -1;
				a->rid = c->rid;
				a->frac_rep = c->frac_rep;
				a->seedlen0 = s->len;
				a->c = c; //ptr
				a->rb = a->qb = a->re = a->qe = H0_;

//				tprof[PE19][tid] ++;

				int flag = 0;
				std::pair<int, int> pr;
				if (s->qbeg)  // left extension
				{
					SeqPair sp;
					sp.h0 = s->len * opt->a;
//					sp.seqid = c->seqid;
					sp.seqid = l;
					sp.regid = av->n - 1;

					if (numPairsLeft >= *wsize_pair) {
						fprintf(stderr, "[0000][%0.4d] Re-allocating seqPairArrays, in Left\n", tid);
						*wsize_pair += 1024;
						seqPairArrayAux = (SeqPair *) realloc(seqPairArrayAux,
						                                      (*wsize_pair + MAX_LINE_LEN)
						                                      * sizeof(SeqPair));
						mmc->seqPairArrayAux[tid] = seqPairArrayAux;
						seqPairArrayLeft128 = (SeqPair *) realloc(seqPairArrayLeft128,
						                                          (*wsize_pair + MAX_LINE_LEN)
						                                          * sizeof(SeqPair));
						mmc->seqPairArrayLeft128[tid] = seqPairArrayLeft128;
						seqPairArrayRight128 = (SeqPair *) realloc(seqPairArrayRight128,
						                                           (*wsize_pair + MAX_LINE_LEN)
						                                           * sizeof(SeqPair));
						mmc->seqPairArrayRight128[tid] = seqPairArrayRight128;
					}


					sp.idq = leftQerOffset;
					sp.idr = leftRefOffset;

					leftQerOffset += s->qbeg;
					if (leftQerOffset >= *wsize_buf_qer)
					{
						fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufQers in %s (left)\n",
						        tid, __func__);
						int64_t tmp = *wsize_buf_qer;
						*wsize_buf_qer *= 2;
						uint8_t *seqBufQer_ = (uint8_t*)
								_mm_realloc(seqBufLeftQer, tmp, *wsize_buf_qer, sizeof(uint8_t));
						mmc->seqBufLeftQer[tid*CACHE_LINE] = seqBufLeftQer = seqBufQer_;

						seqBufQer_ = (uint8_t*)
								_mm_realloc(seqBufRightQer, tmp, *wsize_buf_qer, sizeof(uint8_t));
						mmc->seqBufRightQer[tid*CACHE_LINE] = seqBufRightQer = seqBufQer_;
					}

					uint8_t *qs = seqBufLeftQer + sp.idq;
					for (int i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];

					tmp = s->rbeg - rmax[0];
					leftRefOffset += tmp;
					if (leftRefOffset >= *wsize_buf_ref)
					{
						fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufRefs in %s (left)\n",
						        tid, __func__);
						int64_t tmp = *wsize_buf_ref;
						*wsize_buf_ref *= 2;
						uint8_t *seqBufRef_ = (uint8_t*)
								_mm_realloc(seqBufLeftRef, tmp, *wsize_buf_ref, sizeof(uint8_t));
						mmc->seqBufLeftRef[tid*CACHE_LINE] = seqBufLeftRef = seqBufRef_;

						seqBufRef_ = (uint8_t*)
								_mm_realloc(seqBufRightRef, tmp, *wsize_buf_ref, sizeof(uint8_t));
						mmc->seqBufRightRef[tid*CACHE_LINE] = seqBufRightRef = seqBufRef_;
					}

					uint8_t *rs = seqBufLeftRef + sp.idr;
					for (int64_t i = 0; i < tmp; ++i) rs[i] = rseq[tmp - 1 - i]; //seq1

					sp.len2 = s->qbeg;
					sp.len1 = tmp;
					int minval = sp.h0 + min_(sp.len1, sp.len2) * opt->a;

					if (sp.len1 < MAX_SEQ_LEN8 && sp.len2 < MAX_SEQ_LEN8 && minval < MAX_SEQ_LEN8) {
						numPairsLeft128++;
					}
					else if (sp.len1 < MAX_SEQ_LEN16 && sp.len2 < MAX_SEQ_LEN16 && minval < MAX_SEQ_LEN16){
						numPairsLeft16++;
					}
					else {
						numPairsLeft1++;
					}

					seqPairArrayLeft128[numPairsLeft] = sp;
					numPairsLeft ++;
					a->qb = s->qbeg; a->rb = s->rbeg;
				}
				else
				{
					flag = 1;
					a->score = a->truesc = s->len * opt->a, a->qb = 0, a->rb = s->rbeg;
				}

				if (s->qbeg + s->len != l_query)  // right extension
				{
					int64_t qe = s->qbeg + s->len;
					int64_t re = s->rbeg + s->len - rmax[0];
					assert(re >= 0);
					SeqPair sp;

					sp.h0 = H0_; //random number
//					sp.seqid = c->seqid;
					sp.seqid = l;
					sp.regid = av->n - 1;

					if (numPairsRight >= *wsize_pair)
					{
						fprintf(stderr, "[0000] [%0.4d] Re-allocating seqPairArrays Right\n", tid);
						*wsize_pair += 1024;
						seqPairArrayAux = (SeqPair *) realloc(seqPairArrayAux,
						                                      (*wsize_pair + MAX_LINE_LEN)
						                                      * sizeof(SeqPair));
						mmc->seqPairArrayAux[tid] = seqPairArrayAux;
						seqPairArrayLeft128 = (SeqPair *) realloc(seqPairArrayLeft128,
						                                          (*wsize_pair + MAX_LINE_LEN)
						                                          * sizeof(SeqPair));
						mmc->seqPairArrayLeft128[tid] = seqPairArrayLeft128;
						seqPairArrayRight128 = (SeqPair *) realloc(seqPairArrayRight128,
						                                           (*wsize_pair + MAX_LINE_LEN)
						                                           * sizeof(SeqPair));
						mmc->seqPairArrayRight128[tid] = seqPairArrayRight128;
					}

					sp.len2 = l_query - qe;
					sp.len1 = rmax[1] - rmax[0] - re;

					sp.idq = rightQerOffset;
					sp.idr = rightRefOffset;

					rightQerOffset += sp.len2;
					if (rightQerOffset >= *wsize_buf_qer)
					{
						fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufQers in %s (right)\n",
						        tid, __func__);
						int64_t tmp = *wsize_buf_qer;
						*wsize_buf_qer *= 2;
						uint8_t *seqBufQer_ = (uint8_t*)
								_mm_realloc(seqBufLeftQer, tmp, *wsize_buf_qer, sizeof(uint8_t));
						mmc->seqBufLeftQer[tid*CACHE_LINE] = seqBufLeftQer = seqBufQer_;

						seqBufQer_ = (uint8_t*)
								_mm_realloc(seqBufRightQer, tmp, *wsize_buf_qer, sizeof(uint8_t));
						mmc->seqBufRightQer[tid*CACHE_LINE] = seqBufRightQer = seqBufQer_;
					}

					rightRefOffset += sp.len1;
					if (rightRefOffset >= *wsize_buf_ref)
					{
						fprintf(stderr, "[%0.4d] Re-allocating (doubling) seqBufRefs in %s (right)\n",
						        tid, __func__);
						int64_t tmp = *wsize_buf_ref;
						*wsize_buf_ref *= 2;
						uint8_t *seqBufRef_ = (uint8_t*)
								_mm_realloc(seqBufLeftRef, tmp, *wsize_buf_ref, sizeof(uint8_t));
						mmc->seqBufLeftRef[tid*CACHE_LINE] = seqBufLeftRef = seqBufRef_;

						seqBufRef_ = (uint8_t*)
								_mm_realloc(seqBufRightRef, tmp, *wsize_buf_ref, sizeof(uint8_t));
						mmc->seqBufRightRef[tid*CACHE_LINE] = seqBufRightRef = seqBufRef_;
					}

//					tprof[PE23][tid] += sp.len1 + sp.len2;

					uint8_t *qs = seqBufRightQer + sp.idq;
					uint8_t *rs = seqBufRightRef + sp.idr;

					for (int i = 0; i < sp.len2; ++i) qs[i] = query[qe + i];

					for (int i = 0; i < sp.len1; ++i) rs[i] = rseq[re + i]; //seq1

					int minval = sp.h0 + min_(sp.len1, sp.len2) * opt->a;

					if (sp.len1 < MAX_SEQ_LEN8 && sp.len2 < MAX_SEQ_LEN8 && minval < MAX_SEQ_LEN8) {
						numPairsRight128++;
					}
					else if(sp.len1 < MAX_SEQ_LEN16 && sp.len2 < MAX_SEQ_LEN16 && minval < MAX_SEQ_LEN16) {
						numPairsRight16++;
					}
					else {
						numPairsRight1++;
					}
					seqPairArrayRight128[numPairsRight] = sp;
					numPairsRight ++;
					a->qe = qe; a->re = rmax[0] + re;
				}
				else
				{
					a->qe = l_query, a->re = s->rbeg + s->len;
					// seedcov business, this "if" block should be redundant, check and remove.
					if (a->rb != H0_ && a->qb != H0_)
					{
						int i;
						for (i = 0, a->seedcov = 0; i < c->n; ++i)
						{
							const mem_seed_t *t = &c->seeds[i];
							if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
							    t->rbeg >= a->rb && t->rbeg + t->len <= a->re) // seed fully contained
								a->seedcov += t->len;
						}
					}
				}
			}
			free(rseq);
			// tprof[MEM_ALN2_DOWN1][tid] += __rdtsc() - tim;
		}
	}
	// tprof[MEM_ALN2_UP][tid] += __rdtsc() - timUP;


	int32_t *hist = (int32_t *)_mm_malloc((MAX_SEQ_LEN8 + MAX_SEQ_LEN16 + 32) *
	                                      sizeof(int32_t), 64);

	/* Sorting based score is required as that affects the use of SIMD lanes */
	sortPairsLenExt(seqPairArrayLeft128, numPairsLeft, seqPairArrayAux, hist,
	                numPairsLeft128, numPairsLeft16, numPairsLeft1, opt->a);
	assert(numPairsLeft == (numPairsLeft128 + numPairsLeft16 + numPairsLeft1));


	// SWA
	// uint64_t timL = __rdtsc();
	int nthreads = 1;

	// Now, process all the collected seq-pairs
	// First, left alignment, move out these calls
	BandedPairWiseSW bswLeft(opt->o_del, opt->e_del, opt->o_ins,
	                         opt->e_ins, opt->zdrop, opt->pen_clip5,
	                         opt->mat, opt->a, opt->b, nthreads);

	BandedPairWiseSW bswRight(opt->o_del, opt->e_del, opt->o_ins,
	                          opt->e_ins, opt->zdrop, opt->pen_clip3,
	                          opt->mat, opt->a, opt->b, nthreads);

	int i;
	// Left
	SeqPair *pair_ar = seqPairArrayLeft128 + numPairsLeft128 + numPairsLeft16;
	SeqPair *pair_ar_aux = seqPairArrayAux;
	int nump = numPairsLeft1;

	// scalar
	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int32_t w = opt->w << i;
		// uint64_t tim = __rdtsc();
		bswLeft.scalarBandedSWAWrapper(pair_ar,
		                               seqBufLeftRef,
		                               seqBufLeftQer,
		                               nump,
		                               nthreads,
		                               w);
		// tprof[PE5][0] += nump;
		// tprof[PE6][0] ++;
		// tprof[MEM_ALN2_B][tid] += __rdtsc() - tim;

		int num = 0;
		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev
			int prev = a->score;
			a->score = sp->score;

			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
			    i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip5) {
					a->qb -= sp->qle; a->rb -= sp->tle;
					a->truesc = a->score;
				} else {
					a->qb = 0; a->rb -= sp->gtle;
					a->truesc = sp->gscore;
				}

				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i){
						const mem_seed_t *t = &(a->c->seeds[i]);
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
						    t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}

			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}


	//****************** Left - vector int16 ***********************
	assert(numPairsLeft == (numPairsLeft128 + numPairsLeft16 + numPairsLeft1));

	pair_ar = seqPairArrayLeft128 + numPairsLeft128;
	pair_ar_aux = seqPairArrayAux;

	nump = numPairsLeft16;
	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int32_t w = opt->w << i;
		// int64_t tim = __rdtsc();
#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		bswLeft.scalarBandedSWAWrapper(pair_ar, seqBufLeftRef, seqBufLeftQer, nump, nthreads, w);
#else
		sortPairsLen(pair_ar, nump, seqPairArrayAux, hist);
		bswLeft.getScores16(pair_ar,
		                    seqBufLeftRef,
		                    seqBufLeftQer,
		                    nump,
		                    nthreads,
		                    w);
#endif

//		tprof[PE5][0] += nump;
//		tprof[PE6][0] ++;
		// tprof[MEM_ALN2_B][tid] += __rdtsc() - tim;

		int num = 0;
		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev

			int prev = a->score;
			a->score = sp->score;


			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
			    i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip5) {
					a->qb -= sp->qle; a->rb -= sp->tle;
					a->truesc = a->score;
				} else {
					a->qb = 0; a->rb -= sp->gtle;
					a->truesc = sp->gscore;
				}

				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i){
						const mem_seed_t *t = &(a->c->seeds[i]);
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
						    t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}

	//****************** Left - vector int8 ***********************
	pair_ar = seqPairArrayLeft128;
	pair_ar_aux = seqPairArrayAux;

	nump = numPairsLeft128;
	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int32_t w = opt->w << i;
		// int64_t tim = __rdtsc();

#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		bswLeft.scalarBandedSWAWrapper(pair_ar, seqBufLeftRef, seqBufLeftQer, nump, nthreads, w);
#else
		sortPairsLen(pair_ar, nump, seqPairArrayAux, hist);
		bswLeft.getScores8(pair_ar,
		                   seqBufLeftRef,
		                   seqBufLeftQer,
		                   nump,
		                   nthreads,
		                   w);
#endif

//		tprof[PE1][0] += nump;
//		tprof[PE2][0] ++;
		// tprof[MEM_ALN2_D][tid] += __rdtsc() - tim;

		int num = 0;
		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev

			int prev = a->score;
			a->score = sp->score;

			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
			    i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip5) {
					a->qb -= sp->qle; a->rb -= sp->tle;
					a->truesc = a->score;
				} else {
					a->qb = 0; a->rb -= sp->gtle;
					a->truesc = sp->gscore;
				}

				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i){
						const mem_seed_t *t = &(a->c->seeds[i]);
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
						    t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}

	// tprof[CLEFT][tid] += __rdtsc() - timL;

	// uint64_t timR = __rdtsc();
	// **********************************************************
	// Right, scalar
	for (int l=0; l<numPairsRight; l++) {
		mem_alnreg_t *a;
		SeqPair *sp = &seqPairArrayRight128[l];
		a = &(av_v[sp->seqid].a[sp->regid]); // prev
		sp->h0 = a->score;
	}

	sortPairsLenExt(seqPairArrayRight128, numPairsRight, seqPairArrayAux,
	                hist, numPairsRight128, numPairsRight16, numPairsRight1, opt->a);

	assert(numPairsRight == (numPairsRight128 + numPairsRight16 + numPairsRight1));

	pair_ar = seqPairArrayRight128 + numPairsRight128 + numPairsRight16;
	pair_ar_aux = seqPairArrayAux;
	nump = numPairsRight1;

	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int32_t w = opt->w << i;
		// tim = __rdtsc();
		bswRight.scalarBandedSWAWrapper(pair_ar,
		                                seqBufRightRef,
		                                seqBufRightQer,
		                                nump,
		                                nthreads,
		                                w);
		// tprof[PE7][0] += nump;
		// tprof[PE8][0] ++;
		// tprof[MEM_ALN2_C][tid] += __rdtsc() - tim;
		int num = 0;

		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev
			int prev = a->score;
			a->score = sp->score;

			// no further banding
			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
			    i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip3) {
					a->qe += sp->qle, a->re += sp->tle;
					a->truesc += a->score - sp->h0;
				} else {
					int l_query = seq_[sp->seqid].l_seq;
					a->qe = l_query, a->re += sp->gtle;
					a->truesc += sp->gscore - sp->h0;
				}
				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i) {
						const mem_seed_t *t = &a->c->seeds[i];
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
						    t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}

	// ************************* Right - vector int16 **********************
	pair_ar = seqPairArrayRight128 + numPairsRight128;
	pair_ar_aux = seqPairArrayAux;
	nump = numPairsRight16;

	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int32_t w = opt->w << i;
		// uint64_t tim = __rdtsc();
#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		bswRight.scalarBandedSWAWrapper(pair_ar, seqBufRightRef, seqBufRightQer, nump, nthreads, w);
#else
		sortPairsLen(pair_ar, nump, seqPairArrayAux, hist);
		bswRight.getScores16(pair_ar,
		                     seqBufRightRef,
		                     seqBufRightQer,
		                     nump,
		                     nthreads,
		                     w);
#endif

//		tprof[PE7][0] += nump;
//		tprof[PE8][0] ++;
		// tprof[MEM_ALN2_C][tid] += __rdtsc() - tim;

		int num = 0;

		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev
			//OutScore *o = outScoreArray + l;
			int prev = a->score;
			a->score = sp->score;

			// no further banding
			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
			    i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip3) {
					a->qe += sp->qle, a->re += sp->tle;
					a->truesc += a->score - sp->h0;
				} else {
					int l_query = seq_[sp->seqid].l_seq;
					a->qe = l_query, a->re += sp->gtle;
					a->truesc += sp->gscore - sp->h0;
				}
				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i) {
						const mem_seed_t *t = &a->c->seeds[i];
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
						    t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}


	// ************************* Right, vector int8 **********************
	pair_ar = seqPairArrayRight128;
	pair_ar_aux = seqPairArrayAux;
	nump = numPairsRight128;

	for ( i=0; i<MAX_BAND_TRY; i++)
	{
		int32_t w = opt->w << i;
		// uint64_t tim = __rdtsc();

#if ((!__AVX512BW__) && (!__AVX2__) && (!__SSE2__))
		bswRight.scalarBandedSWAWrapper(pair_ar, seqBufRightRef, seqBufRightQer, nump, nthreads, w);
#else
		sortPairsLen(pair_ar, nump, seqPairArrayAux, hist);
		bswRight.getScores8(pair_ar,
		                    seqBufRightRef,
		                    seqBufRightQer,
		                    nump,
		                    nthreads,
		                    w);
#endif

//		tprof[PE3][0] += nump;
//		tprof[PE4][0] ++;
		// tprof[MEM_ALN2_E][tid] += __rdtsc() - tim;
		int num = 0;

		for (int l=0; l<nump; l++)
		{
			mem_alnreg_t *a;
			SeqPair *sp = &pair_ar[l];
			a = &(av_v[sp->seqid].a[sp->regid]); // prev
			//OutScore *o = outScoreArray + l;
			int prev = a->score;
			a->score = sp->score;
			// no further banding
			if (a->score == prev || sp->max_off < (w >> 1) + (w >> 2) ||
			    i+1 == MAX_BAND_TRY)
			{
				if (sp->gscore <= 0 || sp->gscore <= a->score - opt->pen_clip3) {
					a->qe += sp->qle, a->re += sp->tle;
					a->truesc += a->score - sp->h0;
				} else {
					int l_query = seq_[sp->seqid].l_seq;
					a->qe = l_query, a->re += sp->gtle;
					a->truesc += sp->gscore - sp->h0;
				}
				a->w = max_(a->w, w);
				if (a->rb != H0_ && a->qb != H0_ && a->qe != H0_ && a->re != H0_)
				{
					int i = 0;
					for (i = 0, a->seedcov = 0; i < a->c->n; ++i) {
						const mem_seed_t *t = &a->c->seeds[i];
						if (t->qbeg >= a->qb && t->qbeg + t->len <= a->qe &&
						    t->rbeg >= a->rb && t->rbeg + t->len <= a->re)
							a->seedcov += t->len;
					}
				}
			} else {
				pair_ar_aux[num++] = *sp;
			}
		}
		nump = num;
		SeqPair *tmp = pair_ar;
		pair_ar = pair_ar_aux;
		pair_ar_aux = tmp;
	}

	_mm_free(hist);
	// tprof[CRIGHT][tid] += __rdtsc() - timR;

	if (numPairsLeft >= *wsize_pair || numPairsRight >= *wsize_pair)
	{   // refine it!
		fprintf(stderr, "Error: Unexpected behaviour!!!\n");
		fprintf(stderr, "Error: assert failed for seqPair size, "
		                "numPairsLeft: %d, numPairsRight %d\nExiting.\n",
		        numPairsLeft, numPairsRight);
		exit(EXIT_FAILURE);
	}
	/* Discard seeds and hence their alignemnts */

	lim_g[0] = 0;
	for (int l=1; l<nseq; l++)
		lim_g[l] += lim_g[l-1];

	// uint64_t tim = __rdtsc();
	int *lim = (int *) calloc(BATCH_SIZE, sizeof(int));
	assert(lim != NULL);

	for (int l=0; l<nseq; l++)
	{
		int s_start = 0, s_end = 0;
		uint32_t *srtg = srtgg + lim_g[l];

		int l_query = seq_[l].l_seq;
		mem_chain_v *chn = &chain_ar[l];
		mem_alnreg_v *av = &av_v[l];  // alignment
		mem_chain_t *c;

		for (int j=0; j<chn->n; j++)
		{
			c = &chn->a[j];
//			assert(c->seqid == l);

			s_end = s_start + c->n;

			uint32_t *srt2 = srtg + s_start;
			s_start += c->n;

			int k = 0;
			for (k = c->n-1; k >= 0; k--)
			{
				s = &c->seeds[srt2[k]];
				int i = 0;
				int v = 0;
				for (i = 0; i < av->n && v < lim[l]; ++i)  // test whether extension has been made before
				{
					mem_alnreg_t *p = &av->a[i];
					if (p->qb == -1 && p->qe == -1) {
						continue;
					}

					int64_t rd;
					int qd, w, max_gap;
					if (s->rbeg < p->rb || s->rbeg + s->len > p->re || s->qbeg < p->qb
					    || s->qbeg + s->len > p->qe) {
						v++; continue; // not fully contained
					}

					if (s->len - p->seedlen0 > .1 * l_query) { v++; continue;}
					// qd: distance ahead of the seed on query; rd: on reference
					qd = s->qbeg - p->qb; rd = s->rbeg - p->rb;
					// the maximal gap allowed in regions ahead of the seed
					max_gap = cal_max_gap(opt, qd < rd? qd : rd);
					w = max_gap < p->w? max_gap : p->w; // bounded by the band width
					if (qd - rd < w && rd - qd < w) break; // the seed is "around" a previous hit
					// similar to the previous four lines, but this time we look at the region behind
					qd = p->qe - (s->qbeg + s->len); rd = p->re - (s->rbeg + s->len);
					max_gap = cal_max_gap(opt, qd < rd? qd : rd);
					w = max_gap < p->w? max_gap : p->w;
					if (qd - rd < w && rd - qd < w) break;

					v++;
				}

				// the seed is (almost) contained in an existing alignment;
				// further testing is needed to confirm it is not leading
				// to a different aln
				if (v < lim[l])
				{
					for (v = k + 1; v < c->n; ++v)
					{
						const mem_seed_t *t;
						if (srt2[v] == UINT_MAX) continue;
						t = &c->seeds[srt2[v]];
						//if (t->done == H0_) continue;  //check for interferences!!!
						// only check overlapping if t is long enough;
						// TODO: more efficient by early stopping
						if (t->len < s->len * .95) continue;
						if (s->qbeg <= t->qbeg && s->qbeg + s->len - t->qbeg >= s->len>>2 &&
						    t->qbeg - s->qbeg != t->rbeg - s->rbeg) break;
						if (t->qbeg <= s->qbeg && t->qbeg + t->len - s->qbeg >= s->len>>2 &&
						    s->qbeg - t->qbeg != s->rbeg - t->rbeg) break;
					}
					if (v == c->n) {                  // no overlapping seeds; then skip extension
						mem_alnreg_t *ar = &(av_v[l].a[s->aln]);
						ar->qb = ar->qe = -1;         // purge the alingment
						srt2[k] = UINT_MAX;
//						tprof[PE18][tid]++;
						continue;
					}
				}
				lim[l]++;
			}
		}
	}
	free(srtgg);
	free(srt);
	free(lim);
	// tprof[MEM_ALN2_DOWN][tid] += __rdtsc() - tim;
}

static void seed_and_extend(worker_t *w, int _start, int _end, int tid) {
	_end = std::min(_end, w->n); // Out of right boundary happens
	int n = _end - _start;
	if (n == 0) return;

	// Collect exact matches
	uint64_t real_start = __rdtsc();
	auto *opt = w->opt;
	auto *bwt = w->bwt;
	auto *bns = w->bns;
	auto *pac = w->pac;
	auto &aux = w->aux[tid];
	aux.forward_sst->clear(); aux.backward_sst->clear();
	for (int r = 0; r < n; r++) {
		auto &read = w->seqs[_start + r];
		auto *seq = (uint8_t*) read.seq;
		for (int j = 0; j < read.l_seq; j++) { // Convert ACGTN to 01234 if hasn't done so far
			if (seq[j] > 4) seq[j] = nst_nt4_table[seq[j]];
		}

		std::vector<bwtintv_t> &match = aux.match[r]; match.clear();
		// Seeding
		uint64_t t_start = __rdtsc();
		for (int j = 0; j < read.l_seq; ) {
			j = collect_mem_with_sst(seq, read.l_seq, j, 1, aux);
			for (const auto &m : aux.super_mem) {
				if (mem_len(m) >= opt->min_seed_len)
					match.push_back(m);
			}
		}
		aux.first += __rdtsc() - t_start;
		// Re-seeding
		t_start = __rdtsc();
		int old_n = (int)match.size();
		for (int j = 0; j < old_n; j++) {
			const auto &p = match[j] ;
			int beg = mem_beg(p), end = mem_end(p);
			if (end - beg < (int)(1.0 * opt->min_seed_len * opt->split_factor + .499) or p.x[2] > opt->split_width) continue;
			collect_mem_with_sst(seq, read.l_seq, (beg + end) / 2, p.x[2] + 1, aux);
			for (const auto &m : aux.super_mem) {
				if (mem_len(m) >= opt->min_seed_len)
					match.push_back(m);
			}
		}
		aux.second += __rdtsc() - t_start;
		// 3rd round seeding
		t_start = __rdtsc();
		if (opt->max_mem_intv > 0) {
			for (int j = 0; j < read.l_seq; ) {
				if (seq[j] < 4) {
					bwtintv_t m;
					j = tem_forward_sst(opt, seq, read.l_seq, j, &m, aux);
					if (m.x[2] > 0) match.push_back(m);
				} else {
					j++;
				}
			}
		}
		aux.third += __rdtsc() - t_start;
		std::sort(match.begin(), match.end(), mem_cmp);
	}
	aux.bwt_call_times += aux.forward_sst->bwt_call + aux.backward_sst->bwt_call;
	aux.bwt_real += __rdtsc() - real_start;

	// SAL after merging identical seeds
	real_start = __rdtsc();
	auto unique_sal = aux.unique_sal; unique_sal.clear();
	for (int r = 0; r < n; r++) {
		const auto &mem = aux.match[r];
		auto &seed = aux.seed[r]; seed.clear();
		int array_id = 0;
		for (const auto &m : mem) {
			uint64_t step = m.x[2] > opt->max_occ ? m.x[2] / opt->max_occ : 1;
			for (uint64_t k = 0, count = 0; k < m.x[2] && count < opt->max_occ; k += step, count++) {
				mem_seed_t s = {0};
				s.qbeg = mem_beg(m);
				s.score = s.len = mem_len(m);
				s.rbeg = m.x[0] + k; // Temporarily store SAL request
				seed.push_back(s);
				unique_sal.emplace_back(sal_request_t(m.x[0] + k));
				aux.sal_query_times++;
				array_id++;
			}
		}
	}
	std::sort(unique_sal.begin(), unique_sal.end());
	int _size = 0;
	for (int i = 0; i < unique_sal.size(); i++) {
		if (i == 0 or unique_sal[i-1].que_location != unique_sal[i].que_location) {
			unique_sal[_size++] = unique_sal[i];
		}
	}
	unique_sal.resize(_size);
	// Perform SAL according to the original order
	for (int r = 0; r < n; r++) {
		for (auto &s : aux.seed[r]) {
			auto k = std::lower_bound(unique_sal.begin(), unique_sal.end(), sal_request_t(s.rbeg));
			assert(k != unique_sal.end() && k->que_location == s.rbeg);
			if (k->coordinate == -1) {
				k->coordinate = bwt_sa(bwt, s.rbeg);
				aux.sal_call_times++;
			}
			s.rbeg = k->coordinate;
		}
	}
	aux.sal_real += __rdtsc() - real_start;

	real_start = __rdtsc();
	mem_chain_v chain_ar[n];
	mem_alnreg_v reg_ar[n];
	for (int r = 0; r < n; r++) {
		kv_init(chain_ar[r]);
		kv_init(reg_ar[r]);
	}
	for (int r = 0; r < n; r++) {
		auto &read = w->seqs[_start + r];
		auto *seq = (uint8_t*) read.seq;
		read.sam = nullptr;
		// 1. Chain co-linear seeds
		mem_chain_v chn = mem_chain(opt, bns, read.l_seq, aux.match[r], aux.seed[r]);
		aux.match[r].clear();
		aux.seed[r].clear(); // Seeds are copied to chains, it is okay to release
		// 2. Filter out shadowed chains
		chn.n = mem_chain_flt(opt, chn.n, chn.a);
		// 3. Filter out poor seeds
		mem_flt_chained_seeds(opt, bns, pac, read.l_seq, seq, chn.n, chn.a);
		if (bwa_verbose >= 4) mem_print_chain(bns, &chn);
		chain_ar[r] = chn;
	}

	// 4. SIMD Banded Smith Waterman
	int64_t rst = __rdtsc();
	mem_chain2aln_across_reads_V2(opt, bns, pac, w->seqs + _start, n, chain_ar, reg_ar, &w->mmc, tid);
	aux.simd_real += __rdtsc() - rst;

	for (int r = 0; r < n; r++) {
		auto &read = w->seqs[_start + r];
		auto *seq = (uint8_t*) read.seq;
		mem_chain_v *chain = &chain_ar[r];
		for (int i = 0; i < chain->n; ++i) {
			free(chain->a[i].seeds);
		}
		free(chain->a);
		mem_alnreg_v regs = reg_ar[r];
		int m = 0;
		for (int i = 0; i < regs.n; ++i) { // exclude identical hits
			if (regs.a[i].qe > regs.a[i].qb) {
				if (m != i) regs.a[m++] = regs.a[i];
				else ++m;
			}
		}
		regs.n = m;
		// 5. Filter out shadowed hits
		regs.n = mem_sort_dedup_patch(opt, bns, pac, seq, regs.n, regs.a);
		if (bwa_verbose >= 4) {
			err_printf("* %ld chains remain after removing duplicated chains\n", regs.n);
			for (int i = 0; i < regs.n; ++i) {
				mem_alnreg_t *p = &regs.a[i];
				printf("** %d, [%d,%d) <=> [%ld,%ld)\n", p->score, p->qb, p->qe, (long)p->rb, (long)p->re);
			}
		}
		for (int i = 0; i < regs.n; ++i) {
			mem_alnreg_t *p = &regs.a[i];
			if (p->rid >= 0 && bns->anns[p->rid].is_alt)
				p->is_alt = 1;
		}
		// 6. Mark primary alignment
		if (bwa_verbose >= 4) printf("=====> Finalizing read '%s' <=====\n", read.name);
		mem_mark_primary_se(opt, regs.n, regs.a, w->n_processed + _start + r);
		if (opt->flag & MEM_F_PRIMARY5) mem_reorder_primary5(opt->T, &regs);
		// 7. Convert alignment to SAM format
		mem_reg2sam(opt, bns, pac, &read, &regs, 0, 0);
		free(regs.a);
	}
	aux.ext_real += __rdtsc() - real_start;
}

/*** Memory pre-allocations ***/
void memoryAlloc(const mem_opt_t *opt, worker_t &w, int32_t nreads, int32_t nthreads) {
	int32_t memSize = nreads;
	int32_t readLen = READ_LEN;

	/* Mem allocation section for core kernels */
//	w.regs = NULL; w.chain_ar = NULL; w.seedBuf = NULL;
//	w.regs = (mem_alnreg_v *) calloc(memSize, sizeof(mem_alnreg_v));
//	w.chain_ar = (mem_chain_v*) malloc (memSize * sizeof(mem_chain_v));
//	w.seedBuf = (mem_seed_t *) calloc(sizeof(mem_seed_t),  memSize * AVG_SEEDS_PER_READ);

//	assert(w.seedBuf  != NULL);
//	assert(w.regs     != NULL);
//	assert(w.chain_ar != NULL);

//	w.seedBufSize = BATCH_SIZE * AVG_SEEDS_PER_READ;

	/*** printing ***/
//	int64_t allocMem = memSize * sizeof(mem_alnreg_v) +
//	                   memSize * sizeof(mem_chain_v) +
//	                   sizeof(mem_seed_t) * memSize * AVG_SEEDS_PER_READ;
//	fprintf(stderr, "------------------------------------------\n");
//	fprintf(stderr, "1. Memory pre-allocation for Chaining: %0.4lf MB\n", allocMem/1e6);


	/* SWA mem allocation */
	int64_t wsize = BATCH_SIZE * SEEDS_PER_READ;
	for(int l=0; l<nthreads; l++)
	{
		w.mmc.seqBufLeftRef[l*CACHE_LINE]  = (uint8_t *)
				_mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
		w.mmc.seqBufLeftQer[l*CACHE_LINE]  = (uint8_t *)
				_mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);
		w.mmc.seqBufRightRef[l*CACHE_LINE] = (uint8_t *)
				_mm_malloc(wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN, 64);
		w.mmc.seqBufRightQer[l*CACHE_LINE] = (uint8_t *)
				_mm_malloc(wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN, 64);

		w.mmc.wsize_buf_ref[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_REF;
		w.mmc.wsize_buf_qer[l*CACHE_LINE] = wsize * MAX_SEQ_LEN_QER;

		assert(w.mmc.seqBufLeftRef[l*CACHE_LINE]  != NULL);
		assert(w.mmc.seqBufLeftQer[l*CACHE_LINE]  != NULL);
		assert(w.mmc.seqBufRightRef[l*CACHE_LINE] != NULL);
		assert(w.mmc.seqBufRightQer[l*CACHE_LINE] != NULL);
	}

	for(int l=0; l<nthreads; l++) {
		w.mmc.seqPairArrayAux[l]      = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
		w.mmc.seqPairArrayLeft128[l]  = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
		w.mmc.seqPairArrayRight128[l] = (SeqPair *) malloc((wsize + MAX_LINE_LEN)* sizeof(SeqPair));
		w.mmc.wsize[l] = wsize;

		assert(w.mmc.seqPairArrayAux[l] != NULL);
		assert(w.mmc.seqPairArrayLeft128[l] != NULL);
		assert(w.mmc.seqPairArrayRight128[l] != NULL);
	}


//	allocMem = (wsize * MAX_SEQ_LEN_REF * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads * 2+
//			   (wsize * MAX_SEQ_LEN_QER * sizeof(int8_t) + MAX_LINE_LEN) * opt->n_threads  * 2 +
//			    wsize * sizeof(SeqPair) * opt->n_threads * 3;
//	fprintf(stderr, "2. Memory pre-allocation for BSW: %0.4lf MB\n", allocMem/1e6);

	for (int l=0; l<nthreads; l++)
	{
		w.mmc.wsize_mem[l]     = BATCH_MUL * BATCH_SIZE *               readLen;
//		w.mmc.matchArray[l]    = (SMEM *) _mm_malloc(w.mmc.wsize_mem[l] * sizeof(SMEM), 64);
//		w.mmc.min_intv_ar[l]   = (int32_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int32_t));
//		w.mmc.query_pos_ar[l]  = (int16_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int16_t));
//		w.mmc.enc_qdb[l]       = (uint8_t *) malloc(w.mmc.wsize_mem[l] * sizeof(uint8_t));
//		w.mmc.rid[l]           = (int32_t *) malloc(w.mmc.wsize_mem[l] * sizeof(int32_t));
		w.mmc.lim[l]           = (int32_t *) _mm_malloc((BATCH_SIZE + 32) * sizeof(int32_t), 64); // candidate not for reallocation, deferred for next round of changes.
	}

//	allocMem = nthreads * BATCH_MUL * BATCH_SIZE * readLen * sizeof(SMEM) +
//	           nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t) +
//	           nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int16_t) +
//	           nthreads * BATCH_MUL * BATCH_SIZE * readLen *sizeof(int32_t) +
//	           nthreads * (BATCH_SIZE + 32) * sizeof(int32_t);
//	fprintf(stderr, "3. Memory pre-allocation for BWT: %0.4lf MB\n", allocMem/1e6);
//	fprintf(stderr, "------------------------------------------\n");
}

void memoryFree(worker_t &w, int nthreads) {
	for(int l=0; l<nthreads; l++) {
		_mm_free(w.mmc.seqBufLeftRef[l*CACHE_LINE]);
		_mm_free(w.mmc.seqBufRightRef[l*CACHE_LINE]);
		_mm_free(w.mmc.seqBufLeftQer[l*CACHE_LINE]);
		_mm_free(w.mmc.seqBufRightQer[l*CACHE_LINE]);
	}

	for(int l=0; l<nthreads; l++) {
		free(w.mmc.seqPairArrayAux[l]);
		free(w.mmc.seqPairArrayLeft128[l]);
		free(w.mmc.seqPairArrayRight128[l]);
	}

	for(int l=0; l<nthreads; l++) {
//		_mm_free(w.mmc.matchArray[l]);
//		free(w.mmc.min_intv_ar[l]);
//		free(w.mmc.query_pos_ar[l]);
//		free(w.mmc.enc_qdb[l]);
//		free(w.mmc.rid[l]);
		_mm_free(w.mmc.lim[l]);
	}
}

void mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0) {
	worker_t w;
	double ctime, rtime;

	ctime = cputime(); rtime = realtime();
	w.opt = opt; w.bwt = bwt; w.bns = bns; w.pac = pac;
	w.n = n; w.seqs = seqs; w.n_processed = n_processed;
	w.aux = new thread_aux_t[opt->n_threads];
	for (int i = 0; i < opt->n_threads; i++) {
		w.aux[i].forward_sst = new SST(bwt);
		w.aux[i].backward_sst = new SST(bwt);
	}
	memoryAlloc(opt, w, n, opt->n_threads);

	kt_for(
			opt->n_threads,
			[](void *d, long i, int t) -> void {
				seed_and_extend((worker_t *)d, i*BATCH_SIZE, (i+1)*BATCH_SIZE, t);
			},
			&w,
			(n + BATCH_SIZE - 1) / BATCH_SIZE
	);

	for (int i = 0; i < opt->n_threads; i++) {
		tprof += w.aux[i];
		delete w.aux[i].forward_sst;
		delete w.aux[i].backward_sst;
	}
	delete [] w.aux;
	memoryFree(w, opt->n_threads);

	if (bwa_verbose >= 3)
		fprintf(stderr, "[M::%s] Processed %d reads in %.3f CPU sec, %.3f real sec\n", __func__, n, cputime() - ctime, realtime() - rtime);
}