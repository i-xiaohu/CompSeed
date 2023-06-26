//
// Created by ixiaohu on 2023/5/6.
//

#ifndef COMP_SEED_COMP_SEEDING_H
#define COMP_SEED_COMP_SEEDING_H

#include <string>
#include <vector>
#include "../bwalib/bwa.h"
#include "../cstl/kvec.h"
#include "../cstl/kstring.h"
#include "SST.h"

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

mem_opt_t *mem_opt_init();

/** Storing decompressed NGS reads */
struct ngs_read {
	int16_t len; // Read length
	int16_t is_rc; // Reverse-complemented or not
	int32_t offset; // Offset in consensus built by compressor
	char *bases; // ACGTN (const char*) or their encodings 01234 (const uint8_t*)
	char *sam; // Alignment result in SAM format
	ngs_read() { len = is_rc = offset = 0; bases = sam = nullptr; }
};

/** Seed hit is an exact match between target reference and query read */
struct seed_hit {
	int64_t rbeg;
	int32_t qbeg, len;
	int32_t rid, score;
};

/** A request for Suffix Array Lookup */
struct sal_request {
	uint64_t que_location; // The request position in suffix array
	uint32_t read_id, array_id; // Which mem of which read made the request
	sal_request(uint64_t h, uint32_t r, uint32_t a): que_location(h), read_id(r), array_id(a) {}
	bool operator < (const sal_request &a) const { // Sort to merge duplicated SAL requests
		return this->que_location < a.que_location;
	}
};

/** Chain of co-linear seeds */
struct seed_chain { // Do not change this aligned struct; it affects the B-tree order
	int64_t anchor; // Anchor or the first seed location
	int rid; // Chain can't cross multiple chromosomes
	int first_cover; // Which chain that the current chain first shadows
	int n, m; seed_hit *seeds; // Storing chained seeds
	uint32_t w:29, kept:2, is_alt:1; // Chain weight; Kept level; alternative or not
	float frac_rep; // Repetitive segment fraction across the chain

	seed_chain() = default;

	seed_chain(int64_t anchor): anchor(anchor) {}

	int64_t ref_beg() const { return seeds[0].rbeg; }

	int64_t ref_end() const { return seeds[n-1].rbeg + seeds[n-1].len; }

	int que_beg() const { return seeds[0].qbeg; }

	int que_end() const { return seeds[n-1].qbeg + seeds[n-1].len; }

	void add_first_seed(const seed_hit &s) {
		n = 1; m = 16;
		seeds = (seed_hit*) malloc(m * sizeof(seed_hit));
		seeds[0] = s;
		rid = s.rid;
	}

	void push_back(const seed_hit &s) {
		if (n == m){
			m *= 2;
			seeds = (seed_hit*) realloc(seeds, m * sizeof(seed_hit));
		}
		seeds[n++] = s;
	}

	int calc_weight() const {
		// Seeds are non-decreasing on both read and reference
		int que_end = 0, weight_on_que = 0;
		for (int i = 0; i < n; i++) {
			const auto &s = seeds[i];
			if (s.qbeg >= que_end) weight_on_que += s.len;
			else if (s.qbeg + s.len > que_end) weight_on_que += s.qbeg + s.len - que_end;
			que_end = std::max(que_end, s.qbeg + s.len);
		}

		int64_t ref_end = 0, weight_on_ref = 0;
		for (int i = 0; i < n; i++) {
			const auto &s = seeds[i];
			if (s.rbeg >= ref_end) weight_on_ref += s.len;
			else if (s.rbeg + s.len > ref_end) weight_on_ref += s.rbeg + s.len - ref_end;
			ref_end = std::max(ref_end, s.rbeg + s.len);
		}

		return std::min((int)weight_on_ref, weight_on_que);
	}

	void destroy() const { free(seeds); }
};

struct align_region {
	int64_t rb, re; // [rb, re) reference sequence in the alignment
	int qb, qe;     // [qb, qe) query sequence in the alignment
	int rid;        // Reference sequence ID
	int local_score;// Best local Smith-Waterman score
	                // Actual score of the aligned region; could be the best local score
    int true_score;	// or the global alignment score which is possibly smaller the former
	int sub_score;  // Suboptimal SW score
	int alt_score;  //
	int tandem_sco; // SW score of a tandem hit
	int sub_n;      // approximate number of suboptimal hits
	int band_width; // Band width of SW matrix in seed extension
	int seed_cover; // Length of the aligned region covered by seeds
	int secondary;  // Point to the parent hit shadowing the current hit; < 0 if primary
	int second_all; //
	int seed_len0;  // Length of the staring seed that is extended to this alignment
	int32_t n_comp:30, is_alt:2; // Number of sub-alignments chained together
	float frac_rep; // Repetition fraction; equal to of the corresponding chain
	uint64_t hash;  //
};

/** How many reads to process at a time. Batch size can impact the performance of compressive aligner.
 * Increasing batch size could dig out more redundancy but lowering the load balancing of threads.
 * Lowering batch size could not make good use of benefits provided by compressors. */
#define BATCH_SIZE 512

/** Auxiliary of BWA-MEM seeding */
typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

/** Auxiliary for each thread, maintaining essential buffers for alignment */
struct thread_aux {
	SST *forward_sst = nullptr; // Forward SST caching BWT forward extension
	SST *backward_sst = nullptr; // Backward SST caching BWT backward extension
	std::vector<bwtintv_t> prev_intv, curr_intv; // Buffer for forward search LEP and backward extension
	std::vector<bwtintv_t> super_mem; // SMEMs returned from the function collect-mem
	std::vector<bwtintv_t> match[BATCH_SIZE]; // Exact matches for each read (minimum seed length guaranteed)
	std::vector<sal_request> unique_sal; // Sorted suffix array hit locations for merging SAL operations
	std::vector<seed_hit> seed[BATCH_SIZE]; // Seed hits for each read
	smem_aux_t *mem_aux;

	// Profiling time cost and calls number for BWT extension and SAL
	long sal_call_times = 0;
	long bwt_call_times = 0;
	double seeding_cpu_sec = 0;
	double seeding_real_sec = 0;
	int full_read_match = 0; // Number of full-length matched reads
	int shortcut = 0; // Number of reads that are full-length matched and avoid regular SMEM search
	void operator += (const thread_aux &a) {
		full_read_match += a.full_read_match;
		shortcut += a.shortcut;
		sal_call_times += a.sal_call_times;
		bwt_call_times += a.bwt_call_times;
		seeding_cpu_sec += a.seeding_cpu_sec;
		seeding_real_sec += a.seeding_real_sec;
	}
};

class CompAligner {
private:
	// FM-index
	bwaidx_t *bwa_idx = nullptr;
	bntseq_t *bns = nullptr;
	bwt_t *bwt = nullptr;
	uint8_t *pac = nullptr;

	// Working threads
	int threads_n = 1;
	thread_aux *thr_aux = nullptr;

	int read_length = 0;
	long processed_n = 0;
	std::vector<ngs_read> reads; // Reads input from compressed file

	kstring_t *debug_out = nullptr; // Debug information output to stdout

public:
	void load_index(const char *fn);

	mem_opt_t *opt = nullptr; // BWA-MEM built-in parameters
	bool print_seed = false; // Print all seeds to stdout for validation
	int actual_chunk_size = 0; // Size of each input

	/** Compressed super-mem1 algorithm with SST; used for seeding and re-seeding. */
	static int collect_mem_with_sst(const uint8_t *seq, int len, int pivot, int min_hits, thread_aux &aux);

	int tem_forward_sst(const uint8_t *seq, int len, int start, bwtintv_t *mem, thread_aux &aux) const;

	/** Chaining seeds and filtering out light-weight chains and seeds. */
	std::vector<seed_chain> chaining(const std::vector<seed_hit> &seed);

	bool add_seed_to_chain(seed_chain *c, const seed_hit &s);

	void print_chains_to(const std::vector<seed_chain> &chains, kstring_t *s);

	/** Filtering poor seed in chain (not for short reads) */
	void filter_seed_in_chain(const ngs_read &read, std::vector<seed_chain> &chain);

	/** Calculate Smith-Waterman score for a seed */
	int seed_sw_score(const ngs_read &read, const seed_hit &s);

	/** Extend chain to alignment with Smith-Waterman Algorithm */
	std::vector<align_region> extend_chain(const ngs_read &read,
			std::vector<seed_chain> &chain, thread_aux &aux);

	int smith_waterman(int qlen, const uint8_t *query, int tlen, const uint8_t *target,
					int m, const int8_t *mat, int o_del, int e_del, int o_ins, int e_ins, int w,
					int end_bonus,
					int zdrop,
					int h0,
					int *_qle, int *_tle,
					int *_gtle, int *_gscore,
					int *_max_off, thread_aux &aux);

	/** Restoring offset and corrected strand for reordered reads. */

	/** Align reads from start to end-1 with thread tid. (So far, only run seeding) */
	void seed_and_extend(int start, int end, int tid);

	void display_profile(const thread_aux &total);

	void run(const char *fn);

	void bwa_collect_seed(int seq_id, int tid);

	void bwamem(const char *fn);
};

void bwa_c_style(const char *index_fn, const char *read_fn, int actual_chunk_size, const mem_opt_t *opt);

#endif //COMP_SEED_COMP_SEEDING_H
