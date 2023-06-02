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

/** SMEM Search Tree (SST) Node */
struct SST_Node_t {
	bwtintv_t match; // SA interval in FMD-index
	int children[4]; // Indexes of four children A,C,G,T
	SST_Node_t() { children[0] = children[1] = children[2] = children[3] = -1; }
};

class SST {
private:
	std::vector<SST_Node_t> nodes;
	const bwt_t *bwt;
	bwtintv_t next[4];

public:
	explicit SST(const bwt_t *b);

	long bwt_calls = 0;
	uint64_t bwt_ticks = 0;

	/** At parent node with prefix p, query child node for p + base.
	 * If the child node does not exist, query it from FMD-index. */
	int query_child(int parent, uint8_t base, bool is_back);

	inline int query_forward_child(int parent, uint8_t base);

	inline int query_backward_child(int parent, uint8_t base);

	/** Add [pivot, LEP] to backward SST. */
	inline int add_lep_child(int parent, uint8_t base, const uint64_t *x);

	/** Add suffixes of [pivot + 1, LEP] to backward SST. But SA interval of
	 * suffixes is unknown so {0,0,0} is used to mark such nodes in SST. */
	inline int add_empty_child(int parent, uint8_t base);

	inline bwtintv_t get_intv(int id) { return nodes[id].match; }

	inline int get_child(int parent, uint8_t base) { return nodes[parent].children[base]; }

	/** Only keep root and its four children */
	inline void clear() {
		nodes.resize(5);
		for (int i = 1; i <= 4; i++) {
			for (int &c : nodes[i].children) {
				c = -1;
			}
		}
	}
};

inline int SST::query_forward_child(int parent, uint8_t base) {
	if (nodes[parent].children[base] == -1) {
		uint64_t start = __rdtsc();
		bwt_extend(bwt, &nodes[parent].match, next, 0);
		bwt_ticks += __rdtsc() - start;
		bwt_calls++;

		SST_Node_t child;
		child.match = next[base];
		nodes[parent].children[base] = nodes.size();
		nodes.push_back(child);
	}
	return nodes[parent].children[base];
}

inline int SST::query_backward_child(int parent, uint8_t base) {
	if (nodes[parent].children[base] == -1) {
		uint64_t start = __rdtsc();
		bwt_extend(bwt, &nodes[parent].match, next, 1);
		bwt_ticks += __rdtsc() - start;
		bwt_calls++;

		SST_Node_t child;
		child.match = next[base];
		nodes[parent].children[base] = nodes.size();
		nodes.push_back(child);
	}
	auto &c = nodes[nodes[parent].children[base]];
	// If find an empty node, BWT query is required
	if (c.match.x[0] or c.match.x[1] or c.match.x[2]) {
		return nodes[parent].children[base];
	} else {
		uint64_t start = __rdtsc();
		bwt_extend(bwt, &nodes[parent].match, next, 1);
		bwt_ticks += __rdtsc() - start;
		bwt_calls++;
		c.match = next[base];
		return nodes[parent].children[base];
	}
}

inline int SST::add_lep_child(int parent, uint8_t base, const uint64_t *x) {
	if (nodes[parent].children[base] == -1) {
		nodes[parent].children[base] = nodes.size();
		SST_Node_t child;
		child.match.x[0] = x[0];
		child.match.x[1] = x[1];
		child.match.x[2] = x[2];
		nodes.push_back(child);
	} else {
		auto &c = nodes[nodes[parent].children[base]];
		c.match.x[0] = x[0];
		c.match.x[1] = x[1];
		c.match.x[2] = x[2];
	}
	return nodes[parent].children[base];
}

inline int SST::add_empty_child(int parent, uint8_t base) {
	if (nodes[parent].children[base] == -1) {
		nodes[parent].children[base] = nodes.size();
		SST_Node_t child;
		child.match.x[0] = child.match.x[1] = child.match.x[2] = 0;
		nodes.push_back(child);
	} // else do nothing whatever the child node is empty or not
	return nodes[parent].children[base];
}

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
	std::string archive_name;
	std::string index_name;
	uint8_t compression_level = 2;
	bool separate_N = true; // separate reads with N
	bool single_end_mode = false;
	bool preserve_order_mode = false;
	bool rev_comp_pair = false; // The second read from a pair is reverse-complemented

	ExtendReadsList *hq_reads_list = nullptr;
	ExtendReadsList *lq_reads_list = nullptr;
	ExtendReadsList *n_reads_list = nullptr;

	int read_length;
	int hq_reads_count;
	int lq_reads_count;
	int non_reads_count; // Non-N reads count
	int n_reads_count;
	int total_reads_count;

	long hq_pg_length;
	long non_pg_length; // Non-N PG length
	std::string hq_pg;
	std::string lq_pg;
	std::string n_pg;

	std::vector<int> paired_idx; // paired_idx[2i] and paired_idx[2i+1] store the indexes in PG for a pair of reads
	bool joined_pg_len_std;
	std::vector<uint32_t> pos_pg_32; // The position on PG of reads in original FASTQ order (32-bit)
	std::vector<uint64_t> pos_pg_64; // The position on PG of reads in original FASTQ order (64-bit)
	// org_idx[0] is PG position of the first read in FASTQ
	// ...
	// org_idx[n] is PG position of the n read in FASTQ

	bwaidx_t *bwa_idx = nullptr;
	bntseq_t *bns = nullptr;
	bwt_t *bwt = nullptr;
	uint8_t *pac = nullptr;
	mem_opt_t *opt = nullptr;
	thread_aux_t *thr_aux = nullptr;


	int big_data_round_limit = 20;
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
	void set_archive_name(const char *fn) { this->archive_name = fn; }

	void set_index_name(const char *fn) { this->index_name = fn; }

	void set_threads(int n) { this->threads_n = n; }

	void load_all_PGs(std::ifstream &in);

	template<typename uint_pg_len>
	void apply_rc_pair_to_pg(std::vector<uint_pg_len> &pg_pos);

	void compressive_seeding();

	int collect_smem_with_sst(const uint8_t *seq, int len, int pivot, int min_hits, thread_aux_t &aux);

	int tem_forward_sst(const uint8_t *seq, int len, int start, int min_len, int max_intv, bwtintv_t *mem, thread_aux_t &aux);

	void test_a_batch(int base, const std::vector<long> &offset, std::vector<std::string> &batch, thread_aux_t &aux);

	void seed_and_extend(int batch_id, int tid);

	void generate_sam(int seq_id);

	void display_profile();

//	void seeding_SE();

	void on_dec_reads(const char *fn);

	// BWA-MEM seeding
	void bwamem(const char *fn);
};

#endif //PGRC_LEARN_COMP_SEEDING_H
