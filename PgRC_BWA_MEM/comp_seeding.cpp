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
typedef struct {
	int n, m, first, rid;
	uint32_t w:29, kept:2, is_alt:1;
	float frac_rep;
	int64_t pos;
	seed_hit *seeds;
} mem_chain_t;

typedef kvec_t(mem_chain_t) mem_chain_v;

#define chain_cmp(a, b) (((b).pos < (a).pos) - ((a).pos < (b).pos))
KBTREE_INIT(chn, mem_chain_t, chain_cmp)

// return 1 if the seed is merged into the chain
static int test_and_merge(const mem_opt_t *opt, int64_t l_pac, mem_chain_t *c, const seed_hit *p, int seed_rid) {
	int64_t qend, rend, x, y;
	const seed_hit *last = &c->seeds[c->n - 1];
	qend = last->qbeg + last->len;
	rend = last->rbeg + last->len;
	if (seed_rid != c->rid) return 0; // different chr; request a new chain
	if (p->qbeg >= c->seeds[0].qbeg && p->qbeg + p->len <= qend && p->rbeg >= c->seeds[0].rbeg && p->rbeg + p->len <= rend)
		return 1; // contained seed; do nothing
	if ((last->rbeg < l_pac || c->seeds[0].rbeg < l_pac) && p->rbeg >= l_pac) return 0; // don't chain if on different strand
	x = p->qbeg - last->qbeg; // always non-negtive
	y = p->rbeg - last->rbeg;
	if (y >= 0 && x - y <= opt->w && y - x <= opt->w && x - last->len < opt->max_chain_gap && y - last->len < opt->max_chain_gap) { // grow the chain
		if (c->n == c->m) {
			c->m <<= 1;
			c->seeds = (seed_hit*) realloc(c->seeds, c->m * sizeof(seed_hit));
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
		const seed_hit *s = &c->seeds[j];
		if (s->qbeg >= end) w += s->len;
		else if (s->qbeg + s->len > end) w += s->qbeg + s->len - end;
		end = end > s->qbeg + s->len? end : s->qbeg + s->len;
	}
	tmp = w; w = 0;
	for (j = 0, end = 0; j < c->n; ++j) {
		const seed_hit *s = &c->seeds[j];
		if (s->rbeg >= end) w += s->len;
		else if (s->rbeg + s->len > end) w += s->rbeg + s->len - end;
		end = end > s->rbeg + s->len? end : s->rbeg + s->len;
	}
	w = w < tmp? w : tmp;
	return w < 1<<30? w : (1<<30)-1;
}

mem_chain_v mem_chain(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, int len, const uint8_t *seq, const std::vector<bwtintv_t> &mem, const std::vector<seed_hit> &seeds) {
	mem_chain_v chain; kv_init(chain);
	if (len < opt->min_seed_len) return chain; // if the query is shorter than the seed length, no match
	kbtree_t(chn) *tree = kb_init(chn, KB_DEFAULT_SIZE);

	int b = 0, e = 0, l_rep = 0;
	for (const auto &p : mem) { // compute frac_rep
		int sb = mem_beg(p), se = mem_end(p);
		if (p.x[2] <= opt->max_occ) continue;
		if (sb > e) l_rep += e - b, b = sb, e = se;
		else e = e > se? e : se;
	}
	l_rep += e - b;

	for (const auto &s : seeds) {
		mem_chain_t tmp, *lower, *upper;
		int rid = s.rid, to_add = 0;
		if (rid < 0) continue; // bridging multiple reference sequences or the forward-reverse boundary;
		tmp.pos = s.rbeg;
		if (kb_size(tree)) {
			kb_intervalp(chn, tree, &tmp, &lower, &upper); // find the closest chain
			if (!lower || !test_and_merge(opt, bns->l_pac, lower, &s, rid)) to_add = 1;
		} else to_add = 1;
		if (to_add) { // add the seed as a new chain
			tmp.n = 1; tmp.m = 4;
			tmp.seeds = (seed_hit*) calloc(tmp.m, sizeof(seed_hit));
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
	kb_destroy(chn, tree);

	for (int i = 0; i < chain.n; ++i) chain.a[i].frac_rep = (float)l_rep / len;
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
	kvec_t(int) chains = {0,0,0}; // this keeps int indices of the non-overlapping chains
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

int mem_seed_sw(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const seed_hit *s)
{
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

void mem_flt_chained_seeds(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, int n_chn, mem_chain_t *a)
{
	double min_l = opt->min_chain_weight? MEM_HSP_COEF * opt->min_chain_weight : MEM_MINSC_COEF * log(l_query);
	int i, j, k, min_HSP_score = (int)(opt->a * min_l + .499);
	if (min_l > MEM_SEEDSW_COEF * l_query) return; // don't run the following for short reads
	for (i = 0; i < n_chn; ++i) {
		mem_chain_t *c = &a[i];
		for (j = k = 0; j < c->n; ++j) {
			seed_hit *s = &c->seeds[j];
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

static inline int cal_max_gap(const mem_opt_t *opt, int qlen)
{
	int l_del = (int)((double)(qlen * opt->a - opt->o_del) / opt->e_del + 1.);
	int l_ins = (int)((double)(qlen * opt->a - opt->o_ins) / opt->e_ins + 1.);
	int l = l_del > l_ins? l_del : l_ins;
	l = l > 1? l : 1;
	return l < opt->w<<1? l : opt->w<<1;
}

#define MAX_BAND_TRY  2

void mem_chain2aln(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, int l_query, const uint8_t *query, const mem_chain_t *c, mem_alnreg_v *av)
{
	int i, k, rid, max_off[2], aw[2]; // aw: actual bandwidth used in extension
	int64_t l_pac = bns->l_pac, rmax[2], tmp, max = 0;
	const seed_hit *s;
	uint8_t *rseq = 0;
	uint64_t *srt;

	if (c->n == 0) return;
	// get the max possible span
	rmax[0] = l_pac<<1; rmax[1] = 0;
	for (i = 0; i < c->n; ++i) {
		int64_t b, e;
		const seed_hit *t = &c->seeds[i];
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

	srt = (uint64_t*) malloc(c->n * 8);
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
				const seed_hit *t;
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
			qs = (uint8_t*) malloc(s->qbeg);
			for (i = 0; i < s->qbeg; ++i) qs[i] = query[s->qbeg - 1 - i];
			tmp = s->rbeg - rmax[0];
			rs = (uint8_t*) malloc(tmp);
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
			const seed_hit *t = &c->seeds[i];
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

int mem_patch_reg(const mem_opt_t *opt, const bntseq_t *bns, const uint8_t *pac, uint8_t *query, const mem_alnreg_t *a, const mem_alnreg_t *b, int *_w)
{
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
			int64_t _or, oq, mr, mq;
			int score, w;
			if (q->qe == q->qb) continue; // a[j] has been excluded
			_or = q->re - p->rb; // overlap length on the reference
			oq = q->qb < p->qb? q->qe - p->qb : p->qe - q->qb; // overlap length on the query
			mr = q->re - q->rb < p->re - p->rb? q->re - q->rb : p->re - p->rb; // min ref len in alignment
			mq = q->qe - q->qb < p->qe - p->qb? q->qe - q->qb : p->qe - p->qb; // min qry len in alignment
			if (_or > opt->mask_level_redun * mr && oq > opt->mask_level_redun * mq) { // one of the hits is redundant
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

typedef kvec_t(int) int_v;

static void mem_mark_primary_se_core(const mem_opt_t *opt, int n, mem_alnreg_t *a, int_v *z)
{ // similar to the loop in mem_chain_flt()
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

int mem_mark_primary_se(const mem_opt_t *opt, int n, mem_alnreg_t *a, int64_t id)
{
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
				seed_hit s;
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

	// Print seeds
	for (int i = 0; i < n; i++) {
		const auto &seed = aux.seed[i];
		kstring_t *d = &debug_out[_start + i];
		for (const auto &s : seed) {
			ksprintf(d, "qbeg=%d, rbeg=%d, len=%d\n", s.qbeg, s.rbeg, s.len);
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
