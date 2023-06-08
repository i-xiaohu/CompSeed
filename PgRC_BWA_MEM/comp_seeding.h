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
#include "bwamem.h"
#include "SST.h"

static const char PGRC_SE_MODE = 0; // reordered Singled End
static const char PGRC_PE_MODE = 1; // reordered Paired End
static const char PGRC_ORD_SE_MODE = 2; // preserve Order Singled End
static const char PGRC_ORD_PE_MODE = 3; // preserve Order Paired End

static const char *const PGRC_HEADER = "PgRC";
static const char PGRC_VERSION_MAJOR = 1;
static const char PGRC_VERSION_MINOR = 2;
static const char PGRC_VERSION_REVISION = 2;

#if defined(__GNUC__) && !defined(__clang__)
#if defined(__i386__)
static inline unsigned long long __rdtsc(void)
{
    unsigned long long int x;
    __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
    return x;
}
#elif defined(__x86_64__)
static inline unsigned long long __rdtsc(void)
{
	unsigned hi, lo;
	__asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
	return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}
#endif
#endif

typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

struct time_rec_t {
	uint64_t seeding, reseed, third, sal, total;
	time_rec_t():seeding(0), reseed(0), third(0), sal(0), total(0) {}
	void operator += (const time_rec_t &t) {
		seeding += t.seeding;
		reseed += t.reseed;
		third += t.third;
		sal += t.sal;
		total += t.total;
	}
};

typedef struct {
	int64_t rbeg;
	int32_t qbeg, len;
	int rid, score;
} mem_seed_t;

struct SAL_Packed {
	uint64_t hit_location;
	uint32_t read_id, array_id;
	SAL_Packed(uint64_t h, uint32_t r, uint32_t a): hit_location(h), read_id(r), array_id(a) {}
	bool operator < (const SAL_Packed &a) const {
		return this->hit_location < a.hit_location;
	}
};

#define BATCH_SIZE 512

struct thread_aux_t {
	SST *forward_sst = nullptr; // Forward SST caching BWT forward extension
	SST *backward_sst = nullptr; // Backward SST caching BWT backward extension
	std::vector<bwtintv_t> prev_intv, curr_intv; // Buffer for forward search LEP and backward extension
	std::vector<int> prev_node, curr_node; // Buffer for SST node indexes
	std::vector<bwtintv_t> mem; // Storing all SMEMs but minimal seed length not considered
	std::vector<bwtintv_t> batch_mem[BATCH_SIZE]; // Exact matches for each read
	std::vector<mem_seed_t> batch_seed[BATCH_SIZE]; // Seeds for each read
	std::vector<SAL_Packed> unique_sal; // Sorted suffix array hit locations for merging SAL operations

	std::vector<bwtintv_t> truth_mem[BATCH_SIZE]; // Ground truth provided by BWA-MEM
	std::vector<mem_seed_t> truth_seed[BATCH_SIZE];
	smem_aux_t *mem_aux = nullptr;

	// Profiling time cost and calls number for BWT extension and SAL
	time_rec_t bwa_time, comp_time;
	long comp_bwt_calls[16] = {0};
	long bwa_bwt_calls[16] = {0};
	long comp_sal = 0, bwa_sal = 0;
	uint64_t bwa_bwt_ticks = 0;
	long global_bwa_bwt = 0;
	int full_read_match = 0; // #Reads are full-length matched
	int shortcut = 0; // #Reads that are full-length matched and avoid regular SMEM search

	void operator += (const thread_aux_t &a) {
		bwa_time += a.bwa_time;
		comp_time += a.comp_time;
		for (int i = 0; i < 16; i++) {
			bwa_bwt_calls[i] += a.bwa_bwt_calls[i];
			comp_bwt_calls[i] += a.comp_bwt_calls[i];
		}
		comp_sal += a.comp_sal; bwa_sal += a.bwa_sal;
		bwa_bwt_ticks += a.bwa_bwt_ticks;
		global_bwa_bwt += a.global_bwa_bwt;
		full_read_match += a.full_read_match;
		shortcut += a.shortcut;
	}
};

class CompSeeding {
private:
	bwaidx_t *bwa_idx = nullptr;
	bntseq_t *bns = nullptr;
	bwt_t *bwt = nullptr;
	uint8_t *pac = nullptr;
	mem_opt_t *opt = nullptr;
	thread_aux_t *thr_aux = nullptr;

	int total_reads_count = 0;
	int big_data_round_limit = 0;
	int chunk_size = 5 * 1000 * 1000;
	int threads_n = 1;
	std::vector<std::string> all_batch; // Storing 10*threads batches
	std::vector<long> all_offset; // Storing all offset values
	std::vector<mem_alnreg_v> all_regs;
	std::vector<std::string> all_sam; // Storing alignment results in SAM format
	int n_processed = 0;

	uint64_t cpu_frequency = 1;
	inline double __time(uint64_t x) { return 1.0 * x / cpu_frequency; }

public:
	CompSeeding();

	void load_index(const char *fn);

	void set_threads(int n) { this->threads_n = n; }

	int collect_smem_with_sst(const uint8_t *seq, int len, int pivot, int min_hits, thread_aux_t &aux);

	int tem_forward_sst(const uint8_t *seq, int len, int start, int min_len, int max_intv, bwtintv_t *mem, thread_aux_t &aux);

	void test_a_batch(int base, const std::vector<long> &offset, std::vector<std::string> &batch, thread_aux_t &aux);

	void seed_and_extend(int batch_id, int tid);

	void generate_sam(int seq_id);

	void display_profile();

	void on_dec_reads(const char *fn);

	void bwamem(const char *fn);
};

#endif //PGRC_LEARN_COMP_SEEDING_H
