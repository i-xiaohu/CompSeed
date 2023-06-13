//
// Created by ixiaohu on 2023/5/6.
//

#ifndef PGRC_LEARN_COMP_SEEDING_H
#define PGRC_LEARN_COMP_SEEDING_H

#include <string>
#include <vector>
#include "../PseudoGenome/SeparatedPG.h"
#include "../bwalib/bwa.h"
#include "../cstl/kvec.h"
#include "../cstl/kstring.h"
#include "bwamem.h"
#include "SST.h"

/** For recording time */
#if defined(__GNUC__) && !defined(__clang__) && defined(__x86_64__)
static inline unsigned long long __rdtsc(void) {
	unsigned hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#endif

struct time_record {
	uint64_t seeding, reseed, third, sal, total;
	time_record(): seeding(0), reseed(0), third(0), sal(0), total(0) {}
	void operator += (const time_record &t) {
		seeding += t.seeding;
		reseed += t.reseed;
		third += t.third;
		sal += t.sal;
		total += t.total;
	}
};

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
struct seed_chain { // Do not change this unaligned struct; it affects the B-tree order
	int64_t anchor; // Anchor or the first seed location
	int rid; // Chain can't cross multiple chromosomes
	int first; //
	int n, m; seed_hit *seeds; // Storing chained seeds
	uint32_t w:29, kept:2, is_alt:1; // Chain weight; Kept level; alternative or not
	float frac_rep; // Repetitive segment fraction across the chain

	seed_chain() = default;

	seed_chain(int64_t anchor): anchor(anchor) {}

	int64_t ref_beg() const { return seeds[0].rbeg; }

	int64_t ref_end() const { return seeds[n-1].rbeg + seeds[n-1].len; }

	int que_beg() const { return seeds[0].qbeg; }

	int que_end() const { return seeds[n-1].qbeg + seeds[n-1].len; }

	void first_seed(const seed_hit &s) {
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

	int weight() const {
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
};

struct PPPP {
	int n, m, first, rid;
	uint32_t w:29, kept:2, is_alt:1;
	float frac_rep;
	int64_t pos;
	seed_hit *seeds;
};

/** How many reads to process at a time. Batch size can impact the performance of compressive aligner.
 * Increasing batch size could dig out more redundancy but lowering the load balancing of threads.
 * Lowering batch size could not make good use of benefits provided by compressors. */
#define BATCH_SIZE 512

/** Auxiliary for each thread, maintaining essential buffers for alingment */
struct thread_aux {
	SST *forward_sst = nullptr; // Forward SST caching BWT forward extension
	SST *backward_sst = nullptr; // Backward SST caching BWT backward extension
	std::vector<bwtintv_t> prev_intv, curr_intv; // Buffer for forward search LEP and backward extension
	std::vector<bwtintv_t> super_mem; // SMEMs returned from the function collect-mem
	std::vector<bwtintv_t> match[BATCH_SIZE]; // Exact matches for each read (minimum seed length guaranteed)
	std::vector<sal_request> unique_sal; // Sorted suffix array hit locations for merging SAL operations
	std::vector<seed_hit> seed[BATCH_SIZE]; // Seed hits for each read

	// Profiling time cost and calls number for BWT extension and SAL
	long sal_times = 0;
	int full_read_match = 0; // Number of full-length matched reads
	int shortcut = 0; // Number of reads that are full-length matched and avoid regular SMEM search
	void operator += (const thread_aux &a) {
		full_read_match += a.full_read_match;
		shortcut += a.shortcut;
	}
};

class CompAligner {
private:
	bwaidx_t *bwa_idx = nullptr;
	bntseq_t *bns = nullptr;
	bwt_t *bwt = nullptr;
	uint8_t *pac = nullptr;
	mem_opt_t *opt = nullptr;
	thread_aux *thr_aux = nullptr;

	int total_reads_count = 0;
	int chunk_size = 5 * 1000 * 1000; // Chunk size for each thread; total input size = #theads * chunk size
	std::vector<ngs_read> reads; // Reads input from compressed file

	uint64_t cpu_frequency = 1;
	inline double _time(uint64_t x) { return 1.0 * x / cpu_frequency; }

public:
	CompAligner();

	void load_index(const char *fn);

	int threads_n = 1; // Working thread number
	int input_round_limit = 1; // Input round limit for testing
	kstring_t *debug_out = nullptr; // Debug information output to stdout

	/** Compressed super-mem1 algorithm with SST; used for seeding and re-seeding. */
	static int collect_mem_with_sst(const uint8_t *seq, int len, int pivot, int min_hits, thread_aux &aux);

	int tem_forward_sst(const uint8_t *seq, int len, int start, bwtintv_t *mem, thread_aux &aux);

	/** chaining seeds and filtering out light-weight chains and seeds. */
	std::vector<seed_chain> chaining(const std::vector<seed_hit> &seed);

	bool add_seed_to_chain(seed_chain *c, const seed_hit &s);

	void print_chains_to(const std::vector<seed_chain> &chains, kstring_t *s);

	/** Align reads from start to end-1 with thread tid. */
	void seed_and_extend(int start, int end, int tid);

	void display_profile();

	void run(const char *fn);

	void bwamem(const char *fn);
};

#endif //PGRC_LEARN_COMP_SEEDING_H
