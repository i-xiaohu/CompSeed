//
// Created by ixiaohu on 2023/6/26.
//
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include "../bwalib/bwa.h"
#include "../bwalib/utils.h"
#include "../cstl/kvec.h"
#include "../cstl/ksort.h"
#include "../cstl/kthread.h"

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

typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
	long bwt_cnt, sal_cnt; // Profile
} smem_aux_t;

smem_aux_t *smem_aux_init() {
	smem_aux_t *a;
	a = (smem_aux_t*) calloc(1, sizeof(smem_aux_t));
	a->tmpv[0] = (bwtintv_v*) calloc(1, sizeof(bwtintv_v));
	a->tmpv[1] = (bwtintv_v*) calloc(1, sizeof(bwtintv_v));
	a->bwt_cnt = a->sal_cnt = 0;
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

int bwt_smem1_profile(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_intv, uint64_t max_intv, bwtintv_v *mem, bwtintv_v *tmpvec[2], long *bwt_cnt) {
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
	(*bwt_cnt)++;
	ik.info = x + 1;

	for (i = x + 1, curr->n = 0; i < len; ++i) { // forward search
		if (ik.x[2] < max_intv) { // an interval small enough
			kv_push(bwtintv_t, *curr, ik);
			break;
		} else if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			bwt_extend(bwt, &ik, ok, 0);
			(*bwt_cnt)++;
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
				bwt_extend(bwt, p, ok, 1);
				(*bwt_cnt)++;
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

int bwt_seed_profile(const bwt_t *bwt, int len, const uint8_t *q, int x, int min_len, int max_intv, bwtintv_t *mem, long *bwt_cnt) {
	int i, c;
	bwtintv_t ik, ok[4];

	memset(mem, 0, sizeof(bwtintv_t));
	if (q[x] > 3) return x + 1;
	bwt_set_intv(bwt, q[x], ik); // the initial interval of a single base
	for (i = x + 1; i < len; ++i) { // forward search
		if (q[i] < 4) { // an A/C/G/T base
			c = 3 - q[i]; // complement of q[i]
			bwt_extend(bwt, &ik, ok, 0); (*bwt_cnt)++;
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
			x = bwt_smem1_profile(bwt, len, seq, x, start_width, 0, &a->mem1, a->tmpv, &a->bwt_cnt);
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
		bwt_smem1_profile(bwt, len, seq, (start + end)>>1, p->x[2]+1, 0, &a->mem1, a->tmpv, &a->bwt_cnt);
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
					x = bwt_seed_profile(bwt, len, seq, x, opt->min_seed_len, opt->max_mem_intv, &m, &a->bwt_cnt);
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


typedef struct {
	const mem_opt_t *opt;
	const bwt_t *bwt;
	const bntseq_t *bns;
	const uint8_t *pac;
	smem_aux_t **aux;
	bseq1_t *seqs;
} worker;

typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
	int score;
} mem_seed_t; // unaligned memory

void bwa_worker(void *data, long seq_id, int tid) {
	worker *w = (worker*) data;
	const bseq1_t *b = &w->seqs[seq_id];
	for (int i = 0; i < b->l_seq; i++) {
		if (b->seq[i] > 4) {
			b->seq[i] = nst_nt4_table[b->seq[i]];
		}
	}
	smem_aux_t *aux = w->aux[tid];
	mem_collect_intv(w->opt, w->bwt, b->l_seq, (const uint8_t*) b->seq, aux);
	for (int i = 0; i < aux->mem.n; i++) {
		const bwtintv_t *p = &aux->mem.a[i];
		int slen = (int32_t)p->info - (p->info>>32);
		int64_t step = p->x[2] > w->opt->max_occ? p->x[2] / w->opt->max_occ : 1;
		for (int64_t k = 0, count = 0; k < p->x[2] && count < w->opt->max_occ; k += step, ++count) {
			mem_seed_t s;
			s.rbeg = bwt_sa(w->bwt, p->x[0] + k); // this is the base coordinate in the forward-reverse reference
			aux->sal_cnt++;
			s.qbeg = p->info >> 32;
			s.score= s.len = slen;
//			int rid = bns_intv2rid(w->bns, s.rbeg, s.rbeg + s.len);
		}
	}
}

typedef kvec_t(bseq1_t) bseq1_v;

void bwa_c_style(const char *index_fn, const char *read_fn, int actual_chunk_size, const mem_opt_t *opt) {
	fprintf(stderr, "Running BWA-MEM seeding implemented by C\n");
	bwaidx_t *idx = bwa_idx_load(index_fn, BWA_IDX_ALL); assert(idx != NULL);
	gzFile in = gzopen(read_fn, "r"); assert(in != NULL);
	char *buffer = (char*) calloc(1024, sizeof(char));
	double total_cpu_time = 0;
	long processed_n = 0;
	long bwt_cnt = 0, sal_cnt = 0;
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

		worker w = {};
		w.opt = opt;
		w.bwt = idx->bwt;
		w.bns = idx->bns;
		w.pac = idx->pac;
		w.seqs = seqs.a;
		w.aux = (smem_aux_t**) malloc(opt->n_threads * sizeof(smem_aux_t*));
		for (int i = 0; i < opt->n_threads; i++) w.aux[i] = smem_aux_init();

		double cpu_start = cputime();
		kt_for(
				opt->n_threads,
				bwa_worker,
				&w,
				(int)seqs.n
		);
		total_cpu_time += cputime() - cpu_start;

		for (int i = 0; i < opt->n_threads; i++) {
			bwt_cnt += w.aux[i]->bwt_cnt;
			sal_cnt += w.aux[i]->sal_cnt;
		}
		for (int i = 0; i < seqs.n; i++) free(seqs.a[i].seq); free(seqs.a);
		for (int i = 0; i < opt->n_threads; i++) smem_aux_destroy(w.aux[i]); free(w.aux);
		processed_n += seqs.n;
		fprintf(stderr, "%ld reads processed\n", processed_n);
	}
	fprintf(stderr, "BWT-extend: %ld calls\n", bwt_cnt);
	fprintf(stderr, "SA Lookup:  %ld calls\n", sal_cnt);
	fprintf(stderr, "Seeding cost %.2f CPU seconds\n", total_cpu_time);

	bwa_idx_destroy(idx);
	gzclose(in);
	free(buffer);
}

static void print_usage(const mem_opt_t *opt) {
	fprintf(stderr, "Usage: bwamem_seeding [options] <FM-index> <Reordered Reads>\n");
	fprintf(stderr, "    -t INT        number of threads [%d]\n", opt->n_threads);
	fprintf(stderr, "    -k INT        minimum seed length [%d]\n", opt->min_seed_len);
	fprintf(stderr, "    -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
	fprintf(stderr, "    -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
	fprintf(stderr, "    -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
	fprintf(stderr, "    -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility)\n");
}

int main(int argc, char *argv[]) {
	mem_opt_t *opt = mem_opt_init();
	if (argc == 1) { print_usage(opt); free(opt); return 1; }

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
	bwa_c_style(argv[optind], argv[optind + 1], actual_chunk_size, opt);
	free(opt);
	return 0;
}
