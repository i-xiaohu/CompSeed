//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_MANAGER_H
#define PGRC_LEARN_MANAGER_H

#include <string>
#include <vector>
#include "PseudoGenome/SeparatedPG.h"
using namespace std;

static const char PGRC_SE_MODE = 0; // reordered Singled End
static const char PGRC_PE_MODE = 1; // reordered Paired End
static const char PGRC_ORD_SE_MODE = 2; // preserve Order Singled End
static const char PGRC_ORD_PE_MODE = 3; // preserve Order Paired End
static const char PGRC_MIN_PE_MODE = 4; // Taking PE as SE

static const char *const PGRC_HEADER = "PgRC";
static const char PGRC_VERSION_MAJOR = 1;
static const char PGRC_VERSION_MINOR = 2;
static const char PGRC_VERSION_REVISION = 2;



class Manager {
private:
	string archive_name;
	uint8_t compression_level = 2;
	bool separate_N = true; // separate reads with N
	bool keep_pairing = true; // Keep pairing information
	bool preserve_order_mode = false;
	bool take_pe_as_se = false;
	bool single_end_mode = false;

	string pg_mapped_HQ_prefix;
	string pg_mapped_LQ_prefix;
	string pg_seq_final_HQ_prefix;
	string pg_seq_final_LQ_prefix;
	string pg_N_prefix;

	int read_length;
	int hq_reads_count;
	int lq_reads_count;
	int non_reads_count; // Non-N reads count
	int n_reads_count;
	int total_reads_count;

	long hq_pg_length;
	long non_pg_length; // Non-N PG length
	SeparatedPseudoGenome *hq_pg = nullptr;
	SeparatedPseudoGenome *lq_pg = nullptr;
	SeparatedPseudoGenome *n_pg = nullptr;

	bool joined_pg_len_std;
	vector<int> idx_order;
	const size_t CHUNK_SIZE_IN_BYTES = 100000; // 100KB
	vector<uint32_t> org_idx_32; // The position on PG of reads in original FASTQ order (32-bit)
	vector<uint64_t> org_idx_64; // The position on PG of reads in original FASTQ order (64-bit)
	// org_idx[0] is PG position of the first read in FASTQ
	// ...
	// org_idx[n] is PG position of the n read in FASTQ

public:
	void set_archive_name(const char *fn) { this->archive_name = fn; }

	void load_all_PGs(ifstream &in);

	void write_all_reads_SE(const std::string &out_fn) const;

	/** This function works for both SE_ORD and PE_ORD mode */
	template<typename uint_pg_len>
	void write_all_reads_ORD(const std::string &out_fn, const vector<uint_pg_len> &org_idx) const;

	void decompress(const std::string &out_fn);
};

#endif //PGRC_LEARN_MANAGER_H
