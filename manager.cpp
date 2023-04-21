//
// Created by ixiaohu on 2023/4/17.
//

#include "manager.h"

#include <cstring>
#include <cassert>

#include "PseudoGenome/PseudoGenomeBase.h"
#include "PseudoGenome/ReadList/ExtendReadsList.h"
#include "PseudoGenome/SeparatedPG.h"
#include "Matching/SimplePgMatcher.h"
#include "Utils/helper.h"

static const char *const GOOD_INFIX = "good";
static const char *const BAD_INFIX = "bad";
static const char *const N_INFIX = "N";

void Manager::load_all_PGs(ifstream &in) {
	fprintf(stderr, "1. Decompressing HQ reads...\n");
	PseudoGenomeHeader hq_pgh(in);
	auto *hq_prop = new ReadsSetProperties(in);
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
	auto *lq_prop = new ReadsSetProperties(in);
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
	ReadsSetProperties *n_prop = nullptr;
	ExtendReadsList *n_reads = nullptr;
	if (separate_N) {
		n_pgh = PseudoGenomeHeader(in);
		n_prop = new ReadsSetProperties(in);
		assert(!confirm_text_read_mode(in));
		n_reads = load_extend_reads_list(
				in,
				n_pgh.get_max_read_length(),
				n_pgh.get_reads_count(),
				preserve_order_mode,
				true,
				true);
	}

	read_length = hq_prop->max_read_length;
	hq_reads_count = hq_prop->reads_count;
	lq_reads_count = lq_prop->reads_count;
	non_reads_count = hq_reads_count + lq_reads_count;
	n_reads_count = separate_N ?n_prop->reads_count :0;
	total_reads_count = non_reads_count + n_reads_count;

	hq_pg_length = hq_pgh.get_pg_length();
	non_pg_length = hq_pg_length + lq_pgh.get_pg_length();
	if (preserve_order_mode) {
		joined_pg_len_std = non_pg_length + n_pgh.get_pg_length() <= UINT32_MAX;
		// wait to implement
	} else {
		spg_decompress_reads_order(
				in,
				idx_order,
				preserve_order_mode,
				take_pe_as_se,
				single_end_mode);
	}

	string hq_pg_seq, lq_pg_seq, n_pg_seq;
	restore_all_matched_pg(in, hq_pg_length, hq_pg_seq, lq_pg_seq, n_pg_seq);
//	std::ofstream os("./pg.txt");
//	os << hq_pg_seq << endl;
//	os << lq_pg_seq << endl;
//	os << n_pg_seq << endl;
//	os.close();
	hq_pg = new SeparatedPseudoGenome(std::move(hq_pg_seq), hq_reads, hq_prop);
	lq_pg = new SeparatedPseudoGenome(std::move(lq_pg_seq), lq_reads, lq_prop);
	n_pg = new SeparatedPseudoGenome(std::move(n_pg_seq), n_reads, n_prop);
}

void Manager::write_all_reads_SE(const std::string &out_fn) const {
	fstream out(archive_name + "_out", ios_base::out | ios_base::binary | ios::trunc);
	string buffer, read1; read1.resize(read_length);
	const uint64_t buf_size_guard = CHUNK_SIZE_IN_BYTES;
	uint64_t total_size = total_reads_count * (read_length + 1); // Including '\n' in bases line
	buffer.reserve(total_size < buf_size_guard ?total_size :buf_size_guard + (read_length + 1));
	for (int i = 0; i < hq_reads_count; i++) {
		if (buffer.size() > buf_size_guard) {
			out << buffer;
			buffer.resize(0);
		}
		hq_pg->get_next_mis_read((char *) read1.data());
		buffer.append(read1);
		buffer.push_back('\n');
	}
	for (int i = 0; i < lq_reads_count; i++) {
		if (buffer.size() > buf_size_guard) {
			out << buffer;
			buffer.resize(0);
		}
		lq_pg->get_next_raw_read((char *) read1.data());
		buffer.append(read1);
		buffer.push_back('\n');
	}
	for (int i = 0; i < n_reads_count; i++) {
		if (buffer.size() > buf_size_guard) {
			out << buffer;
			buffer.resize(0);
		}
		n_pg->get_next_raw_read((char*) read1.data());
		buffer.append(read1);
		buffer.push_back('\n');
	}
	out << buffer;
	out.close();
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
	take_pe_as_se = (pgrc_mode == PGRC_MIN_PE_MODE);
	single_end_mode = (pgrc_mode == PGRC_SE_MODE or pgrc_mode == PGRC_ORD_SE_MODE);

	load_all_PGs(in);

	if (keep_pairing) {
		if (joined_pg_len_std) {

		} else {

		}
	}

	if (not preserve_order_mode) {
		if (single_end_mode) write_all_reads_SE("");
//		else write_all_reads_PE();
	} else {
//		if (joined_pg_len_std) write_all_reads_ORD<int>();
//		else write_all_reads_ORD<long>();
	}

	fprintf(stderr, "Decompressed %d reads in total\n", total_reads_count);

	if (hq_pg) { delete hq_pg; hq_pg = nullptr; }
	if (lq_pg) { delete lq_pg; lq_pg = nullptr; }
	if (n_pg) { delete n_pg; n_pg = nullptr; }
	in.close();
}