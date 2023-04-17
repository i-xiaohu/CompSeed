//
// Created by ixiaohu on 2023/4/17.
//

#include "manager.h"

#include <fstream>
#include <cstring>
#include <cassert>

static const char *const GOOD_INFIX = "good";
static const char *const BAD_INFIX = "bad";
static const char *const N_INFIX = "N";

void Manager::load_all_PGs(ifstream &in) {
	fprintf(stderr, "1. Loading HQ pseudo-genome\n");

}

void Manager::decompress() {
	ifstream in(archive_name); assert(in.is_open());
	for (int i = 0; i < strlen(PGRC_HEADER); i++) {
		assert(in.get() == PGRC_HEADER[i]);
	}
	char pgrc_mode = in.get(); assert(pgrc_mode == '#');
	char version_major = in.get();
	char version_minor = in.get();
	char version_revision = in.get();
	if (version_major != PGRC_VERSION_MAJOR or version_minor != PGRC_VERSION_MINOR) {
		fprintf(stderr, "Archive is packed with a different version PgRC %d.%d\n",
		  version_major, version_minor);
		exit(EXIT_FAILURE);
	}
	pgrc_mode = in.get(); assert(0 <= pgrc_mode and pgrc_mode <= 4); // from 0 to 4
	separate_N = (bool) in.get();
	keep_pairing = (pgrc_mode == PGRC_PE_MODE or pgrc_mode == PGRC_ORD_PE_MODE);
	string tmp_directory; in >> tmp_directory; tmp_directory += "/";
    in.get();

    pg_mapped_HQ_prefix = tmp_directory + GOOD_INFIX;
    pg_mapped_LQ_prefix = tmp_directory + BAD_INFIX;
    pg_seq_final_HQ_prefix = tmp_directory + GOOD_INFIX;
    pg_seq_final_LQ_prefix = tmp_directory + BAD_INFIX;
    pg_N_prefix = tmp_directory + N_INFIX;

	preserve_order_mode = (pgrc_mode == PGRC_ORD_SE_MODE or pgrc_mode == PGRC_ORD_PE_MODE);
	single_end_mode = (pgrc_mode == PGRC_SE_MODE or pgrc_mode == PGRC_ORD_SE_MODE);



	in.close();
}