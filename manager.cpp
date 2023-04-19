//
// Created by ixiaohu on 2023/4/17.
//

#include "manager.h"

#include <fstream>
#include <cstring>
#include <cassert>

#include "PseudoGenome/PseudoGenomeBase.h"
#include "PseudoGenome/ReadList/ExtendReadsList.h"
#include "ReadsSet/ReadsSetBase.h"
#include "Utils/helper.h"

static const char *const GOOD_INFIX = "good";
static const char *const BAD_INFIX = "bad";
static const char *const N_INFIX = "N";

void Manager::load_all_PGs(ifstream &in) {
	fprintf(stderr, "1. Decompressing HQ reads...\n");
	PseudoGenomeHeader hq_pgh(in);
	ReadsSetProperties hq_prop(in);
	assert(!confirm_text_read_mode(in));
	assert(preserve_order_mode == false);
	auto *hq_reads = load_extend_reads_list(
			in,
			hq_pgh.get_max_read_length(),
			hq_pgh.get_reads_count(),
			preserve_order_mode,
			false,
			false);

	fprintf(stderr, "2. Decompressing LQ reads...\n");
	PseudoGenomeHeader lq_pgh(in);
	ReadsSetProperties lq_prop(in);
	assert(!confirm_text_read_mode(in));
	auto *lq_reads = load_extend_reads_list(
			in,
			lq_pgh.get_max_read_length(),
			lq_pgh.get_reads_count(),
			preserve_order_mode,
			true,
			true);

	fprintf(stderr, "3. Decompressing N reads...\n");
	PseudoGenomeHeader n_pgh;
	ReadsSetProperties n_prop;
	ExtendReadsList *n_reads = nullptr;
	if (separate_N) {
		n_pgh = PseudoGenomeHeader(in);
		n_prop = ReadsSetProperties(in);
		assert(!confirm_text_read_mode(in));
		n_reads = load_extend_reads_list(
				in,
				n_pgh.get_max_read_length(),
				n_pgh.get_reads_count(),
				preserve_order_mode,
				true,
				true);
	}
}

void Manager::decompress() {
	ifstream in(archive_name); assert(in.is_open());
	for (int i = 0; i < strlen(PGRC_HEADER); i++) {
		assert(in.get() == PGRC_HEADER[i]);
	}
	char version_mode = in.get(); assert(version_mode == '#');
	char version_major = in.get();
	char version_minor = in.get();
	char version_revision = in.get();
	if (version_major != PGRC_VERSION_MAJOR or version_minor != PGRC_VERSION_MINOR) {
		fprintf(stderr, "Archive is packed with a different version PgRC %d.%d\n",
		  version_major, version_minor);
		exit(EXIT_FAILURE);
	}
	compression_level = in.get();
	char pgrc_mode = in.get(); assert(0 <= pgrc_mode and pgrc_mode <= 4); // from 0 to 4
	separate_N = (bool) in.get();
	if (pgrc_mode == PGRC_PE_MODE or pgrc_mode == PGRC_ORD_PE_MODE) {
		keep_pairing = (bool) in.get();
	}
	string tmp_directory; in >> tmp_directory; tmp_directory += "/";
    in.get();

    pg_mapped_HQ_prefix = tmp_directory + GOOD_INFIX;
    pg_mapped_LQ_prefix = tmp_directory + BAD_INFIX;
    pg_seq_final_HQ_prefix = tmp_directory + GOOD_INFIX;
    pg_seq_final_LQ_prefix = tmp_directory + BAD_INFIX;
    pg_N_prefix = tmp_directory + N_INFIX;

	preserve_order_mode = (pgrc_mode == PGRC_ORD_SE_MODE or pgrc_mode == PGRC_ORD_PE_MODE);
	single_end_mode = (pgrc_mode == PGRC_SE_MODE or pgrc_mode == PGRC_ORD_SE_MODE);

	load_all_PGs(in);

	in.close();
}