//
// Created by ixiaohu on 2023/6/27.
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "../bwalib/bwa.h"
#include "../bwalib/utils.h"
#include "../cstl/kvec.h"
#include "../cstl/kthread.h"
#include "../cstl/ksort.h"

#define min(a, b) ((a) < (b) ?(a) :(b))

typedef struct {
	int a, b;               // match score and mismatch penalty
	int o_del, e_del;
	int o_ins, e_ins;
	int pen_unpaired;       // phred-scaled penalty for unpaired reads
	int pen_clip5,pen_clip3;// clipping penalty. This score is not deducted from the DP score.
	int w;                  // band width
	int zdrop;              // Z-dropoff

	uint64_t max_mem_intv;

	int T;                  // output score threshold; only affecting output
	int flag;               // see MEM_F_* macros
	int min_seed_len;       // minimum seed length
	int min_chain_weight;
	int max_chain_extend;
	float split_factor;     // split into a seed if MEM is longer than min_seed_len*split_factor
	int split_width;        // split into a seed if its occurence is smaller than this value
	int max_occ;            // skip a seed if its occurence is larger than this value
	int max_chain_gap;      // do not chain seed if it is max_chain_gap-bp away from the closest seed
	int n_threads;          // number of threads
	int chunk_size;         // process chunk_size-bp sequences in a batch
	float mask_level;       // regard a hit as redundant if the overlap with another better hit is over mask_level times the min length of the two hits
	float drop_ratio;       // drop a chain if its seed coverage is below drop_ratio times the seed coverage of a better chain overlapping with the small chain
	float XA_drop_ratio;    // when counting hits for the XA tag, ignore alignments with score < XA_drop_ratio * max_score; only effective for the XA tag
	float mask_level_redun;
	float mapQ_coef_len;
	int mapQ_coef_fac;
	int max_ins;            // when estimating insert size distribution, skip pairs with insert longer than this value
	int max_matesw;         // perform maximally max_matesw rounds of mate-SW for each end
	int max_XA_hits, max_XA_hits_alt; // if there are max_hits or fewer, output them all
	int8_t mat[25];         // scoring matrix; mat[0] == 0 if unset
} mem_opt_t;

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

/** SMEM Search Tree in C implementation */
typedef struct {
	bwtintv_t match; // SA interval in FMD-index
	int children[4]; // Indexes of four children A,C,G,T
} SST_Node_t;
typedef kvec_t(SST_Node_t) SST_Node_v;

typedef struct {
	SST_Node_v nodes;
	const bwt_t *bwt;
	bwtintv_t next[4];
	int64_t bwt_calls;
} SST_t;

inline SST_t* sst_init(const bwt_t *bwt) {
	SST_t *sst = (SST_t*) calloc(1, sizeof(SST_t));
	kv_resize(SST_Node_t, sst->nodes, 1024);
	sst->bwt = bwt;

	SST_Node_t root;
	for (int i = 0; i < 4; i++) root.children[i] = i + 1;
	kv_push(SST_Node_t, sst->nodes, root);
	for (uint8_t c = 0; c < 4; c++) {
		SST_Node_t *p = kv_pushp(SST_Node_t, sst->nodes);
		p->children[0] = p->children[1] = p->children[2] = p->children[3] = -1;
		bwt_set_intv(bwt, c, p->match);
	}
	return sst;
}

inline int query_forward_child(SST_t *sst, int parent, uint8_t base) {
	if (sst->nodes.a[parent].children[base] == -1) {
		bwt_extend(sst->bwt, &sst->nodes.a[parent].match, sst->next, 0);
		sst->bwt_calls++;

		SST_Node_t *p = kv_pushp(SST_Node_t, sst->nodes);
		p->match = sst->next[base];
		p->children[0] = p->children[1] = p->children[2] = p->children[3] = -1;
		sst->nodes.a[parent].children[base] = sst->nodes.n - 1;
	}
	return sst->nodes.a[parent].children[base];
}

inline int query_backward_child(SST_t *sst, int parent, uint8_t base) {
	if (sst->nodes.a[parent].children[base] == -1) {
		bwt_extend(sst->bwt, &sst->nodes.a[parent].match, sst->next, 1);
		sst->bwt_calls++;

		SST_Node_t *p = kv_pushp(SST_Node_t, sst->nodes);
		p->match = sst->next[base];
		p->children[0] = p->children[1] = p->children[2] = p->children[3] = -1;
		sst->nodes.a[parent].children[base] = sst->nodes.n - 1;
	}
	SST_Node_t *c = &sst->nodes.a[sst->nodes.a[parent].children[base]];
	// If find an empty node, BWT query is required
	if (c->match.x[0] || c->match.x[1] || c->match.x[2]) {
		return sst->nodes.a[parent].children[base];
	} else {
		bwt_extend(sst->bwt, &sst->nodes.a[parent].match, sst->next, 1);
		sst->bwt_calls++;
		c->match = sst->next[base];
		return sst->nodes.a[parent].children[base];
	}
}

inline int add_lep_child(SST_t *sst,  int parent, uint8_t base, const uint64_t *x) {
	if (sst->nodes.a[parent].children[base] == -1) {
		SST_Node_t *p = kv_pushp(SST_Node_t, sst->nodes);
		p->match.x[0] = x[0];
		p->match.x[1] = x[1];
		p->match.x[2] = x[2];
		p->children[0] = p->children[1] = p->children[2] = p->children[3] = -1;
		sst->nodes.a[parent].children[base] = sst->nodes.n - 1;
	} else {
		SST_Node_t *c = &sst->nodes.a[sst->nodes.a[parent].children[base]];
		c->match.x[0] = x[0];
		c->match.x[1] = x[1];
		c->match.x[2] = x[2];
	}
	return sst->nodes.a[parent].children[base];
}

inline int add_empty_child(SST_t *sst, int parent, uint8_t base) {
	if (sst->nodes.a[parent].children[base] == -1) {
		SST_Node_t *p = kv_pushp(SST_Node_t, sst->nodes);
		p->match.x[0] = p->match.x[1] = p->match.x[2] = 0;
		p->children[0] = p->children[1] = p->children[2] = p->children[3] = -1;
		sst->nodes.a[parent].children[base] = sst->nodes.n - 1;
	} // else do nothing whatever the child node is empty or not
	return sst->nodes.a[parent].children[base];
}

inline void sst_clear(SST_t *sst) {
	sst->nodes.n = 5;
	for (int i = 1; i <= 4; i++) {
		for (int j = 0; j < 4; j++) {
			sst->nodes.a[i].children[j] = -1;
		}
	}
}

inline void sst_destroy(SST_t *sst) {
	free(sst->nodes.a);
	free(sst);
}

#define BATCH_SIZE 512

/** Seed hit is an exact match between target reference and query read */
typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
	int32_t rid, score;
} seed_hit_t;
typedef kvec_t(seed_hit_t) seed_hit_v;

/** A request for Suffix Array Lookup */
typedef struct {
	uint64_t que_location; // The request position in suffix array
	uint32_t read_id, array_id; // Which mem of which read made the request
//	sal_request(uint64_t h, uint32_t r, uint32_t a): que_location(h), read_id(r), array_id(a) {}
//	bool operator < (const sal_request &a) const { // Sort to merge duplicated SAL requests
//		return this->que_location < a.que_location;
//	}
} sal_request_t;
typedef kvec_t(sal_request_t) sal_request_v;

#define sal_lt(a, b) ((a).que_location < (b).que_location)
KSORT_INIT(sal, sal_request_t, sal_lt)

/** Auxiliary for each thread, maintaining essential buffers for alignment */
typedef struct {
	SST_t *forward_sst; // Forward SST caching BWT forward extension
	SST_t *backward_sst; // Backward SST caching BWT backward extension
	bwtintv_v prev_intv, curr_intv; // Buffer for forward search LEP and backward extension
	bwtintv_v super_mem; // SMEMs returned from the function collect-mem
	bwtintv_v match[BATCH_SIZE]; // Exact matches for each read (minimum seed length guaranteed)
	sal_request_v unique_sal; // Sorted suffix array hit locations for merging SAL operations
	seed_hit_v seed[BATCH_SIZE]; // Seed hits for each read

	// Profiling time cost and calls number for BWT extension and SAL
	long sal_call_times;
	long bwt_call_times;
	double seeding_real_sec;
	int full_read_match; // Number of full-length matched reads
	int shortcut; // Number of reads that are full-length matched and avoid regular SMEM search
} thread_aux_t;

void aux_add(thread_aux_t *a, const thread_aux_t *b) {
	a->full_read_match += b->full_read_match;
	a->shortcut += b->shortcut;
	a->sal_call_times += b->sal_call_times;
	a->bwt_call_times += b->bwt_call_times;
	a->seeding_real_sec += b->seeding_real_sec;
}

thread_aux_t *aux_init(const bwt_t *bwt) {
	thread_aux_t *aux = (thread_aux_t*) calloc(1, sizeof(thread_aux_t));
	aux->forward_sst = sst_init(bwt);
	aux->backward_sst = sst_init(bwt);
	for (int i = 0; i < BATCH_SIZE; i++) {
		kv_init(aux->match[i]);
		kv_init(aux->seed[i]);
	}
	return aux;
}

void aux_destroy(thread_aux_t *aux) {
	sst_destroy(aux->forward_sst);
	sst_destroy(aux->backward_sst);
	free(aux->prev_intv.a); free(aux->curr_intv.a);
	free(aux->super_mem.a);
	for (int i = 0; i < BATCH_SIZE; i++) free(aux->match[i].a);
	free(aux->unique_sal.a);
	for (int i = 0; i < BATCH_SIZE; i++) free(aux->seed[i].a);
	free(aux);
}

typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	thread_aux_t **aux;
	int n;
	bseq1_t *seqs;
} worker_t;

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

#define kv_empty(v) ((v).n == 0)
#define kv_back(v) ((v).a[(v).n-1])
#define kv_clear(v) ((v).n = 0)

static inline int mem_beg(const bwtintv_t *a) { return a->info >> 32; }

static inline int mem_end(const bwtintv_t *a) { return (int)a->info; }

static inline int mem_len(const bwtintv_t *a) { return mem_end(a) - mem_beg(a); }

int collect_mem_with_sst(const uint8_t *seq, int len, int pivot, int min_hits, thread_aux_t *aux) {
	kv_clear(aux->super_mem);
	if (seq[pivot] > 3) return pivot + 1;
	bwtintv_t ik, next[4];
	int node_id = query_forward_child(aux->forward_sst, 0, seq[pivot]);
	ik = aux->forward_sst->nodes.a[node_id].match;
	ik.info = pivot + 1;
	kv_clear(aux->prev_intv);
	int ret_pivot = len;
	for (int i = pivot + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = query_forward_child(aux->forward_sst, node_id, c);
			next[c] = aux->forward_sst->nodes.a[node_id].match;
			if (next[c].x[2] != ik.x[2]) {
				kv_push(bwtintv_t, aux->prev_intv, ik);
				if (next[c].x[2] < min_hits) {
					ret_pivot = i;
					break;
				}
			}
			ik = next[c];
			ik.info = i + 1;
		} else {
			kv_push(bwtintv_t, aux->prev_intv, ik);
			ret_pivot = i + 1;
			break;
		}
	}
	if (ret_pivot == len) kv_push(bwtintv_t, aux->prev_intv, ik);
	if (pivot == 0) { // Quickly return full-length matched reads
		kv_push(bwtintv_t, aux->super_mem, kv_back(aux->prev_intv));
		return ret_pivot;
	}
	bwt_reverse_intvs(&aux->prev_intv);

	// Collect node indexes for LEPs in backward SST
	for (int i = 0; i < aux->prev_intv.n; i++) {
		bwtintv_t *p = &aux->prev_intv.a[i];
		node_id = 0;
		for (int j = (int) p->info - 1; j >= pivot + 1; j--) {
			node_id = add_empty_child(aux->backward_sst, node_id, seq[j]);
		}
		node_id = add_lep_child(aux->backward_sst, node_id, seq[pivot], p->x);
		p->info |= (1UL * node_id << 32);
	}

	for (int i = pivot - 1; i >= -1; i--) {
		int c = (i == -1) ? 4 : seq[i];
		kv_clear(aux->curr_intv);
		for (int j = 0; j < aux->prev_intv.n; j++) {
			ik = aux->prev_intv.a[j];
			node_id = ik.info >> 32;
			if (c < 4) {
				node_id = query_backward_child(aux->backward_sst, node_id, c);
				next[c] = aux->backward_sst->nodes.a[node_id].match;
			}
			if (c > 3 || next[c].x[2] < min_hits) {
				if (kv_empty(aux->super_mem) || i + 1 < kv_back(aux->super_mem).info >> 32) {
					ik.info = (1UL * (i + 1) << 32) | (int32_t) ik.info;
					kv_push(bwtintv_t, aux->super_mem, ik);
				}
			} else if (kv_empty(aux->curr_intv) || next[c].x[2] != kv_back(aux->curr_intv).x[2]) {
				next[c].info = (1UL * node_id << 32) | (int32_t) ik.info;
				kv_push(bwtintv_t, aux->curr_intv, next[c]);
			}
		}
		if (kv_empty(aux->curr_intv)) break;
		bwtintv_v tmp = aux->prev_intv;
		aux->prev_intv = aux->curr_intv;
		aux->curr_intv = tmp;
	}
	return ret_pivot;
}

int tem_forward_sst(const mem_opt_t *opt, const uint8_t *seq, int len, int start, bwtintv_t *mem, thread_aux_t *aux) {
	if (seq[start] > 3) return start + 1;
	memset(mem, 0, sizeof(bwtintv_t));
	int node_id = query_forward_child(aux->forward_sst, 0, seq[start]);
	for (int i = start + 1; i < len; i++) {
		if (seq[i] < 4) {
			int c = 3 - seq[i];
			node_id = query_forward_child(aux->forward_sst, node_id, c);
			bwtintv_t ik = aux->forward_sst->nodes.a[node_id].match;
			if (ik.x[2] < opt->max_mem_intv && i - start >= opt->min_seed_len) {
				*mem = ik;
				mem->info = (1UL * start << 32) | (i + 1);
				return i + 1;
			}
		} else return i + 1;
	}
	return len;
}

#define intv_lt(a, b) ((a).info < (b).info)
KSORT_INIT(mem_intv, bwtintv_t, intv_lt)

void comp_worker(void *data, long seq_id, int tid) {
	double real_start = realtime();
	worker_t *w = (worker_t*)data;
	const mem_opt_t *opt = w->opt;
	thread_aux_t *aux = w->aux[tid];
	sst_clear(aux->forward_sst); sst_clear(aux->backward_sst);
	int _start = seq_id * BATCH_SIZE;
	int _end = min((seq_id + 1) * BATCH_SIZE, w->n); // Out of right boundary happens
	int n = _end - _start;
	long ref_position = -1; // Starting from this position might avoid many BWT-extension for fully matched reads
	for (int i = n - 1; i >= 0; i--) {
		bseq1_t *read = &w->seqs[_start + i];
		uint8_t *bases = (uint8_t*) read->seq;
		for (int j = 0; j < read->l_seq; j++) { // Convert ACGTN to 01234 if hasn't done so far
			if (bases[j] > 4) bases[j] = nst_nt4_table[bases[j]];
		}

		bwtintv_v *match = &aux->match[i]; match->n = 0;
		int full_match = 0; // Whether this read could be full-length matched
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

		if (full_match == 0) { // Go to the regular SMEM searching pass
			for (int j = 0; j < read->l_seq; ) {
				j = collect_mem_with_sst(bases, read->l_seq, j, 1, aux);
				for (int k = 0; k < aux->super_mem.n; k++) {
					const bwtintv_t *m = &aux->super_mem.a[k];
					if (mem_len(m) >= opt->min_seed_len) {
						kv_push(bwtintv_t, *match, *m);
					}
				}
			}
//			ref_position = read.offset;
		}
		w->aux[tid]->full_read_match += (kv_size(*match) == 1 && mem_len(&match->a[0]) == read->l_seq);

		int old_n = (int)(kv_size(*match));
		for (int j = 0; j < old_n; j++) {
			const bwtintv_t *p = &match->a[j] ;
			int beg = mem_beg(p), end = mem_end(p);
			if (end - beg < (int)(1.0 * opt->min_seed_len * opt->split_factor + .499) || p->x[2] > opt->split_width) continue;
			collect_mem_with_sst(bases, read->l_seq, (beg + end) / 2, p->x[2] + 1, aux);
			for (int k = 0; k < aux->super_mem.n; k++) {
				const bwtintv_t *m = &aux->super_mem.a[k];
				if (mem_len(m) >= opt->min_seed_len) {
					kv_push(bwtintv_t, *match, *m);
				}

			}
		}

		if (opt->max_mem_intv > 0) {
			for (int j = 0; j < read->l_seq; ) {
				if (bases[j] < 4) {
					bwtintv_t m;
					j = tem_forward_sst(opt, bases, read->l_seq, j, &m, aux);
					if (m.x[2] > 0) kv_push(bwtintv_t, *match, m) ;
				} else {
					j++;
				}
			}
		}
		ks_introsort(mem_intv, match->n, match->a);
	}

	// Find unique hit locations by merging and sorting
	sal_request_v *unique_sal = &aux->unique_sal; kv_clear(*unique_sal);
	for (int read_id = 0; read_id < n; read_id++) {
		const bwtintv_v *mem = &aux->match[read_id];
		seed_hit_v *seed = &aux->seed[read_id]; kv_clear(*seed);
		int array_id = 0;
		for (int i = 0; i < mem->n; i++) {
			const bwtintv_t *m = &mem->a[i];
			uint64_t step = m->x[2] > opt->max_occ ? m->x[2] / opt->max_occ : 1;
			for (uint64_t k = 0, count = 0; k < m->x[2] && count < opt->max_occ; k += step, count++) {
				seed_hit_t s;
				s.qbeg = mem_beg(m);
				s.score = s.len = mem_len(m);
				kv_push(seed_hit_t, *seed, s);

				sal_request_t t;
				t.que_location = m->x[0] + k;
				t.read_id = read_id;
				t.array_id = array_id;
				kv_push(sal_request_t, *unique_sal, t);
				array_id++;
			}
		}
	}
	ks_introsort(sal, unique_sal->n, unique_sal->a);
	// Lookup suffix array for hit locations
	uint64_t coordinate = 0;
	for (int i = 0; i < unique_sal->n; i++) {
		const sal_request_t *p = &unique_sal->a[i];
		if (i == 0 || unique_sal->a[i-1].que_location != p->que_location) {
			coordinate = bwt_sa(w->bwt, p->que_location);
		}
		seed_hit_t *s = &aux->seed[p->read_id].a[p->array_id];
		s->rbeg = coordinate;
		s->rid = bns_intv2rid(w->bns, s->rbeg, s->rbeg + s->len);
	}

	// The code below that extends seeds to full alignments is deleted
	aux->seeding_real_sec += realtime() - real_start;
}

typedef kvec_t(bseq1_t) bseq1_v;

static void run(const char *index_fn, const char *read_fn, int actual_chunk_size, const mem_opt_t *opt) {
	fprintf(stderr, "Running compressive seeding in C style...\n");

	// Prepare for thread auxiliary
	thread_aux_t total; memset(&total, 0, sizeof(total));
	bwaidx_t *idx = bwa_idx_load(index_fn, BWA_IDX_ALL); assert(idx != NULL);
	gzFile in = gzopen(read_fn, "r"); assert(in != NULL);
	char *buffer = (char*) calloc(1024, sizeof(char));
	double total_cpu_time = 0;
	long processed_n = 0;
	while (1) {
		// Input a batch of data
		long bytes = 0;
		bseq1_v seqs; kv_init(seqs);
		while (gzgets(in, buffer, 1024)) {
			int len = (int)strlen(buffer); buffer[--len] = '\0';
			bseq1_t b;
			b.seq = strdup(buffer);
			b.l_seq = len;
			kv_push(bseq1_t, seqs, b);
			bytes += b.l_seq;
			if (bytes >= actual_chunk_size) break;
		}
		if (seqs.n == 0) break; // End Of File

		// Normalize all reordered reads (restoring offset and correcting strand)

		// Processing (I/O thread not supported yet)
		worker_t w = {};
		w.opt = opt;
		w.bwt = idx->bwt;
		w.bns = idx->bns;
		w.pac = idx->pac;
		w.n = seqs.n;
		w.seqs = seqs.a;
		w.aux = (thread_aux_t **) malloc(opt->n_threads * sizeof(thread_aux_t *));
		for (int i = 0; i < opt->n_threads; i++) w.aux[i] = aux_init(idx->bwt);

		double cpu_start = cputime();
		kt_for(
				opt->n_threads,
				comp_worker,
				&w,
				((int)seqs.n + BATCH_SIZE - 1) / BATCH_SIZE
		);
		total_cpu_time += cputime() - cpu_start;

		for (int i = 0; i < seqs.n; i++) free(seqs.a[i].seq); free(seqs.a);
		for (int i = 0; i < opt->n_threads; i++) aux_add(&total, w.aux[i]);
		for (int i = 0; i < opt->n_threads; i++) aux_destroy(w.aux[i]); free(w.aux);
		processed_n += seqs.n;
		fprintf(stderr, "%ld reads processed\n", processed_n);
	}

	fprintf(stderr, "Seeding cost %.2f CPU seconds\n", total_cpu_time);
//	display_profile(total);

	bwa_idx_destroy(idx);
	gzclose(in);
	free(buffer);
}

int main(int argc, char *argv[]) {
	mem_opt_t *opt = mem_opt_init();
	if (argc == 1) { free(opt); return 1; }

	const char short_opts[] = "t:k:r:y:c:K:";
	int fixed_chunk_size = 0;
	while (1) {
		int c = getopt(argc, argv, short_opts);
		if (c < 0) break;
		else if (c == 't') opt->n_threads = strtol(optarg, NULL, 10);
		else if (c == 'k') opt->min_seed_len = strtol(optarg, NULL, 10);
		else if (c == 'r') opt->split_factor = strtof(optarg, NULL);
		else if (c == 'y') opt->max_mem_intv = strtol(optarg, NULL, 10);
		else if (c == 'c') opt->max_occ = strtol(optarg, NULL, 10);
		else if (c == 'K') fixed_chunk_size = strtol(optarg, NULL, 10);
		else {
			fprintf(stderr, "Unrecognized option\n");
			exit(EXIT_FAILURE);
		}
	}
	int actual_chunk_size = fixed_chunk_size == 0 ?opt->n_threads * opt->chunk_size :fixed_chunk_size;
	run(argv[optind], argv[optind + 1], actual_chunk_size, opt);
	free(opt);
	return 0;
}
