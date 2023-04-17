//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_MANAGER_H
#define PGRC_LEARN_MANAGER_H

#include <string>
using namespace std;

static const char PGRC_SE_MODE = 0; // reordered Singled End
static const char PGRC_PE_MODE = 1; // reordered Paired End
static const char PGRC_ORD_SE_MODE = 2; // preserve Order Singled End
static const char PGRC_ORD_PE_MODE = 3; // preserve Order Paired End
static const char PGRC_MIN_PE_MODE = 4; // Compatible Paired End

static const char *const PGRC_HEADER = "PgRC";
static const char PGRC_VERSION_MAJOR = 1;
static const char PGRC_VERSION_MINOR = 2;
static const char PGRC_VERSION_REVISION = 2;



class Manager {
private:
	string archive_name;
	uint8_t compression_level = 2;
	bool separate_N = true; // separate reads with N
	bool keep_pairing = false; // Keep pairing information
	bool preserve_order_mode = false;
	bool single_end_mode = false;

	string pg_mapped_HQ_prefix;
	string pg_mapped_LQ_prefix;
	string pg_seq_final_HQ_prefix;
	string pg_seq_final_LQ_prefix;
	string pg_N_prefix;
public:
	void set_archive_name(const char *fn) { this->archive_name = fn; }

	void load_all_PGs(ifstream &in);

	void decompress();
};

#endif //PGRC_LEARN_MANAGER_H
