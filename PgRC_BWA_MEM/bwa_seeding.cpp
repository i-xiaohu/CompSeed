//
// Created by ixiaohu on 2023/5/6.
//

#include "bwa_seeding.h"

#include <cstring>
#include <cassert>

#include "../Matching/SimplePgMatcher.h"
#include "../Utils/helper.h"

void BWA_seeding::load_all_PGs(std::ifstream &in) {
	fprintf(stderr, "1. Decompressing HQ reads...\n");
	PseudoGenomeHeader hq_pgh(in);
	auto *hq_prop = new ReadsSetProperties(in);
	assert(!confirm_text_read_mode(in));
	hq_reads_list = load_extend_reads_list(
			in,
			hq_pgh.get_max_read_length(),
			hq_pgh.get_reads_count(),
			preserve_order_mode,
			false,
			false);
	hq_reads_list->reads_count = hq_prop->reads_count;
	hq_reads_list->read_length = hq_prop->max_read_length;

	fprintf(stderr, "2. Decompressing LQ reads...\n");
	PseudoGenomeHeader lq_pgh(in);
	auto *lq_prop = new ReadsSetProperties(in);
	assert(!confirm_text_read_mode(in));
	lq_reads_list = load_extend_reads_list(
			in,
			lq_pgh.get_max_read_length(),
			lq_pgh.get_reads_count(),
			preserve_order_mode,
			true,
			true);
	lq_reads_list->reads_count = lq_prop->reads_count;
	lq_reads_list->read_length = lq_prop->max_read_length;

	fprintf(stderr, "3. Decompressing N reads...\n");
	PseudoGenomeHeader n_pgh;
	ReadsSetProperties *n_prop = nullptr;
	if (separate_N) {
		n_pgh = PseudoGenomeHeader(in);
		n_prop = new ReadsSetProperties(in);
		assert(!confirm_text_read_mode(in));
		n_reads_list = load_extend_reads_list(
				in,
				n_pgh.get_max_read_length(),
				n_pgh.get_reads_count(),
				preserve_order_mode,
				true,
				true);
		n_reads_list->reads_count = n_prop->reads_count;
		n_reads_list->read_length = n_prop->max_read_length;
	}

	read_length = hq_prop->max_read_length;
	hq_reads_count = hq_prop->reads_count;
	lq_reads_count = lq_prop->reads_count;
	non_reads_count = hq_reads_count + lq_reads_count;
	n_reads_count = n_prop ?n_prop->reads_count :0;
	total_reads_count = non_reads_count + n_reads_count;

	hq_pg_length = hq_pgh.get_pg_length();
	non_pg_length = hq_pg_length + lq_pgh.get_pg_length();
	if (preserve_order_mode) {
		joined_pg_len_std = non_pg_length + n_pgh.get_pg_length() <= UINT32_MAX;
		if (joined_pg_len_std) {
			decompress_pg_position(in, pos_pg_32, total_reads_count, single_end_mode);
		} else {
			decompress_pg_position(in, pos_pg_64, total_reads_count, single_end_mode);
		}
	} else if (not single_end_mode) {
		restore_paired_idx(in, paired_idx);
	}

	restore_all_matched_pg(in, hq_pg_length, hq_pg, lq_pg, n_pg);
}

template<typename uint_pg_len>
void BWA_seeding::apply_rc_pair_to_pg(std::vector<uint_pg_len> &pg_pos) {
	if (preserve_order_mode) {
		int hq_idx = 0;
		const int pairs_count = total_reads_count / 2;
		for (int i = 0; i < pairs_count; i++) {
			uint_pg_len pos = pg_pos[i];
			if (pos < hq_pg_length) hq_idx++;
		}
		for (int i = pairs_count; i < total_reads_count; i++) {
			uint_pg_len pos = pg_pos[i];
			if (pos < hq_pg_length) {
				hq_reads_list->rev_comp[hq_idx] = !hq_reads_list->rev_comp[hq_idx];
				hq_idx++;
			}
		}
	} else {
		// Reverse complement each read from file2
		for (int i = 1; i < total_reads_count; i += 2) {
			int idx = paired_idx[i];
			if (idx < hq_reads_count) {
				hq_reads_list->rev_comp[idx] = !hq_reads_list->rev_comp[idx];
			}
		}
	}
}

void BWA_seeding::seeding_SE() {
	char *buffer = new char[read_length * 2]; memset(buffer, 0, read_length * 2);
	long curr_position = 0, curr_mis_cnt = 0;
	for (int i = 0; i < hq_reads_list->reads_count; i++) {
		curr_position += hq_reads_list->off[i];
		memcpy(buffer, (hq_pg.data() + curr_position), hq_reads_list->read_length);
		if (hq_reads_list->rev_comp[i]) {
			rev_comp_in_place(buffer, hq_reads_list->read_length);
		}
		for (uint8_t j = 0; j < hq_reads_list->mis_cnt[i]; j++) {
			uint8_t mis_pos = hq_reads_list->mis_off[curr_mis_cnt];
			uint8_t mis_code = hq_reads_list->mis_sym_code[curr_mis_cnt];
			buffer[mis_pos] = code_to_mismatch(buffer[mis_pos], mis_code);
			curr_mis_cnt++;
		}
		fprintf(stdout, "%s\n", buffer);
	}

	curr_position = 0;
	for (int i = 0; i < lq_reads_list->reads_count; i++) {
		curr_position += lq_reads_list->off[i];
		memcpy(buffer, (lq_pg.data() + curr_position), lq_reads_list->read_length);
		fprintf(stdout, "%s\n", buffer);
	}

	curr_position = 0;
	for (int i = 0; i < n_reads_list->reads_count; i++) {
		curr_position += n_reads_list->off[i];
		memcpy(buffer, (n_pg.data() + curr_position), n_reads_list->read_length);
		fprintf(stdout, "%s\n", buffer);
	}

	delete [] buffer;
}

void BWA_seeding::compressive_seeding() {
	std::ifstream in(archive_name); assert(in.is_open());
	for (int i = 0; i < strlen(PGRC_HEADER); i++) assert(in.get() == PGRC_HEADER[i]);
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
	char pgrc_mode = in.get();
	if (pgrc_mode != PGRC_SE_MODE) {
		fprintf(stderr, "PgRC-BWA-MEM now supports SE mode only\n");
		exit(EXIT_FAILURE);
	}
	const std::string five_modes[] = {"SE", "PE", "SE_ORD", "PE_ORD"};
	fprintf(stderr, "Compression mode: %s\n", five_modes[pgrc_mode].c_str());
	separate_N = (bool) in.get();
	if (pgrc_mode == PGRC_PE_MODE or pgrc_mode == PGRC_ORD_PE_MODE) {
		rev_comp_pair = (bool) in.get();
	}
	std::string tmp_directory; in >> tmp_directory; tmp_directory += "/";
	in.get();

	preserve_order_mode = (pgrc_mode == PGRC_ORD_SE_MODE or pgrc_mode == PGRC_ORD_PE_MODE);
	single_end_mode = (pgrc_mode == PGRC_SE_MODE or pgrc_mode == PGRC_ORD_SE_MODE);

	load_all_PGs(in);

	if (rev_comp_pair) {
		fprintf(stderr, "Reverse complement second read from a pair\n");
		if (joined_pg_len_std) apply_rc_pair_to_pg(pos_pg_32);
		else apply_rc_pair_to_pg(pos_pg_64);
	}

	if (pgrc_mode == PGRC_SE_MODE) seeding_SE();

	fprintf(stderr, "Decompressed %d reads in total\n", total_reads_count);
	in.close();
}

int main(int argc, char *argv[]) {
	if (argc == 1) {
		fprintf(stderr, "Usage: PBM <bwa_index> <pgrc_archieve>\n");
		return 1;
	}
	BWA_seeding worker;
	worker.set_index_name(argv[1]);
	worker.set_archive_name(argv[2]);
	worker.compressive_seeding();
}
