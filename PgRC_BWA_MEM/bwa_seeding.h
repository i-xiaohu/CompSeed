//
// Created by ixiaohu on 2023/5/6.
//

#ifndef PGRC_LEARN_BWA_SEEDING_H
#define PGRC_LEARN_BWA_SEEDING_H

#include <string>
#include <vector>
#include "../PseudoGenome/SeparatedPG.h"

static const char PGRC_SE_MODE = 0; // reordered Singled End
static const char PGRC_PE_MODE = 1; // reordered Paired End
static const char PGRC_ORD_SE_MODE = 2; // preserve Order Singled End
static const char PGRC_ORD_PE_MODE = 3; // preserve Order Paired End

static const char *const PGRC_HEADER = "PgRC";
static const char PGRC_VERSION_MAJOR = 1;
static const char PGRC_VERSION_MINOR = 2;
static const char PGRC_VERSION_REVISION = 2;


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

public:
	void set_archive_name(const char *fn) { this->archive_name = fn; }

	void load_all_PGs(std::ifstream &in);

	template<typename uint_pg_len>
	void apply_rc_pair_to_pg(std::vector<uint_pg_len> &pg_pos);

	void compressive_seeding();

	void set_index_name(const char *fn) { this->index_name = fn; }

	void seeding_SE();
};

#endif //PGRC_LEARN_BWA_SEEDING_H
