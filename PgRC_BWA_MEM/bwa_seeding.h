//
// Created by ixiaohu on 2023/5/6.
//

#ifndef PGRC_LEARN_BWA_SEEDING_H
#define PGRC_LEARN_BWA_SEEDING_H

#include <string>
#include <vector>
#include "../PseudoGenome/SeparatedPG.h"
#include "../bwalib/bwa.h"

static const char PGRC_SE_MODE = 0; // reordered Singled End
static const char PGRC_PE_MODE = 1; // reordered Paired End
static const char PGRC_ORD_SE_MODE = 2; // preserve Order Singled End
static const char PGRC_ORD_PE_MODE = 3; // preserve Order Paired End

static const char *const PGRC_HEADER = "PgRC";
static const char PGRC_VERSION_MAJOR = 1;
static const char PGRC_VERSION_MINOR = 2;
static const char PGRC_VERSION_REVISION = 2;

typedef struct {
	bwtintv_v mem, mem1, *tmpv[2];
} smem_aux_t;

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

class BWA_seeding {
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

	const size_t CHUNK_SIZE_IN_BYTES = 100000; // 100KB

	bwaidx_t *bwa_idx = nullptr;
	mem_opt_t *mem_opt = nullptr;
	smem_aux_t *mem_aux = nullptr;

	int COMP_BATCH_SIZE = 1024; // Process 1024 sorted reads at a time
	int batch_id = 0;
	long divergent_column = 0, total_cs_length = 0;
	int read_with_n = 0, overflow_read = 0;
	long reseed_count = 0, identical_reseed = 0;
	long reseed_bin[512] = {0};
	long original_sal = 0, compressed_sal = 0;

public:
	void set_archive_name(const char *fn) { this->archive_name = fn; }

	void load_all_PGs(std::ifstream &in);

	template<typename uint_pg_len>
	void apply_rc_pair_to_pg(std::vector<uint_pg_len> &pg_pos);

	void compressive_seeding();

	void set_index_name(const char *fn) { this->index_name = fn; }

	void traditional_seeding(const std::string &cs, std::vector<long> &offset, std::vector<std::string> &batch);



	void seeding_SE();
};

#endif //PGRC_LEARN_BWA_SEEDING_H
