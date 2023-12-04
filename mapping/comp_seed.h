//
// Created by ixiaohu on 2023/5/6.
//

#ifndef COMP_SEED_COMP_SEED_H
#define COMP_SEED_COMP_SEED_H

#include <string>
#include <vector>
#include <x86intrin.h>
#include "../bwalib/bwa.h"
#include "../cstl/kvec.h"
#include "../cstl/kstring.h"
#include "SST.h"

#define MEM_MAPQ_COEF 30.0
#define MEM_MAPQ_MAX  60

#define MEM_F_PE        0x2
#define MEM_F_NOPAIRING 0x4
#define MEM_F_ALL       0x8
#define MEM_F_NO_MULTI  0x10
#define MEM_F_NO_RESCUE 0x20
#define MEM_F_REF_HDR	0x100
#define MEM_F_SOFTCLIP  0x200
#define MEM_F_SMARTPE   0x400
#define MEM_F_PRIMARY5  0x800
#define MEM_F_KEEP_SUPP_MAPQ 0x1000

/* How many reads to process at a time.
 * Batch size can impact the performance of compressive aligner.
 * Increasing batch size digs out more redundancy provided by compressors.
 * Decreasing batch size improves the load balance of multiple threads. */
#define BATCH_SIZE 512

/* CompSeed restricts maximum read length within 16-bit integer. */
#define MAX_READ_LEN 65535

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

/* I modified the structure and made the memory align given that CompSeed works for short reads compressed by tailored
 * NGS compressors, which are typically designed for fixed-length reads of hundreds of base pairs. */
struct mem_seed_t {
	int64_t rbeg;
	int16_t qbeg, len;
	int32_t score;
};

/* A request for Suffix Array Lookup(SAL) */
struct sal_request_t {
	uint64_t que_location; // The request position in suffix array
	uint32_t read_id, array_id; // Which mem of which read made the request
	sal_request_t(uint64_t h, uint32_t r, uint32_t a): que_location(h), read_id(r), array_id(a) {}
	bool operator < (const sal_request_t &a) const { // Sort to merge duplicated SAL requests
		return this->que_location < a.que_location;
	}
};

typedef struct {
	int64_t rb, re; // [rb,re): reference sequence in the alignment
	int qb, qe;     // [qb,qe): query sequence in the alignment
	int rid;        // reference seq ID
	int score;      // best local SW score
	int truesc;     // actual score corresponding to the aligned region; possibly smaller than $score
	int sub;        // 2nd best SW score
	int alt_sc;
	int csub;       // SW score of a tandem hit
	int sub_n;      // approximate number of suboptimal hits
	int w;          // actual band width used in extension
	int seedcov;    // length of regions coverged by seeds
	int secondary;  // index of the parent hit shadowing the current hit; <0 if primary
	int secondary_all;
	int seedlen0;   // length of the starting seed
	int n_comp:30, is_alt:2; // number of sub-alignments chained together
	float frac_rep;
	uint64_t hash;
} mem_alnreg_t;

typedef kvec_t(mem_alnreg_t) mem_alnreg_v;

typedef struct {
	int low, high;   // lower and upper bounds within which a read pair is considered to be properly paired
	int failed;      // non-zero if the orientation is not supported by sufficient data
	double avg, std; // mean and stddev of the insert size distribution
} mem_pestat_t;

// This struct is only used for the convenience of API.
typedef struct {
	int64_t pos;     // forward strand 5'-end mapping position
	int rid;         // reference sequence index in bntseq_t; <0 for unmapped
	int flag;        // extra flag
	uint32_t is_rev:1, is_alt:1, mapq:8, NM:22; // is_rev: whether on the reverse strand; mapq: mapping quality; NM: edit distance
	int n_cigar;     // number of CIGAR operations
	uint32_t *cigar; // CIGAR in the BAM encoding: opLen<<4|op; op to integer mapping: MIDSH=>01234
	char *XA;        // alternative mappings

	int score, sub, alt_sc;
} mem_aln_t;

/* Auxiliary for each thread, maintaining essential buffers for alignment */
struct thread_aux_t {
	SST *forward_sst = nullptr; // Forward SST caching BWT forward extension
	SST *backward_sst = nullptr; // Backward SST caching BWT backward extension
	std::vector<bwtintv_t> prev_intv, curr_intv; // Buffer for forward search LEP and backward extension
	std::vector<bwtintv_t> super_mem; // SMEMs returned from the function collect-mem
	std::vector<bwtintv_t> match[BATCH_SIZE]; // Exact matches for each read (minimum seed length guaranteed)
	std::vector<sal_request_t> unique_sal; // Sorted suffix array hit locations for merging SAL operations
	std::vector<mem_seed_t> seed[BATCH_SIZE]; // Seed hits for each read

	// Query and call number of BWT-extend and SAL
	long sal_call_times = 0, sal_query_times = 0;
	long bwt_call_times = 0, bwt_query_times = 0;
	// CPU cycles of each stage for each thread
	uint64_t bwt_real = 0, sal_real = 0, ext_real = 0;
	uint64_t first = 0, second = 0, third = 0;
	void operator += (const thread_aux_t &a) {
		bwt_call_times += a.bwt_call_times;
		bwt_query_times += a.bwt_query_times;
		sal_call_times += a.sal_call_times;
		sal_query_times += a.sal_query_times;
		bwt_real += a.bwt_real;
		sal_real += a.sal_real;
		ext_real += a.ext_real;
		first += a.first;
		second += a.second;
		third += a.third;
	}
};

mem_opt_t *mem_opt_init();

void mem_process_seqs(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0);

#endif //COMP_SEED_COMP_SEED_H
