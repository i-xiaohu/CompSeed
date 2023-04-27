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

	string hq_pg_seq, lq_pg_seq, n_pg_seq;
	restore_all_matched_pg(in, hq_pg_length, hq_pg_seq, lq_pg_seq, n_pg_seq);
	hq_pg = new SeparatedPseudoGenome(std::move(hq_pg_seq), hq_reads, hq_prop);
	lq_pg = new SeparatedPseudoGenome(std::move(lq_pg_seq), lq_reads, lq_prop);
	n_pg = new SeparatedPseudoGenome(std::move(n_pg_seq), n_reads, n_prop);
}

void Manager::write_all_reads_SE(const std::string &out_fn) const {
	fstream out(out_fn, ios_base::out | ios_base::binary | ios::trunc);
	string buffer, read1; read1.resize(read_length);
	const uint64_t buf_size_guard = CHUNK_SIZE_IN_BYTES;
	uint64_t total_size = total_reads_count * (read_length + 1); // Including '\n' in bases line
	buffer.reserve(total_size < buf_size_guard ?total_size :buf_size_guard + (read_length + 1));
	for (int i = 0; i < hq_reads_count; i++) {
		if (buffer.size() > buf_size_guard) {
			out << buffer;
			buffer.resize(0);
		}
		hq_pg->get_next_mis_read((char*) read1.data());
		buffer.append(read1);
		buffer.push_back('\n');
	}
	for (int i = 0; i < lq_reads_count; i++) {
		if (buffer.size() > buf_size_guard) {
			out << buffer;
			buffer.resize(0);
		}
		lq_pg->get_next_raw_read((char*) read1.data());
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

void Manager::write_all_reads_PE(const std::string &out_fn) const {
	hq_pg->get_reads_list()->enable_constant_access(true);
	lq_pg->get_reads_list()->enable_constant_access(true);
	if (n_pg) n_pg->get_reads_list()->enable_constant_access(true);

	for (int p = 1; p <= 2; p++) {
		fstream os(out_fn + "_" + to_string(p), ios::out | ios::binary | ios::trunc);
		string buffer, read1; read1.resize(read_length);
		const uint64_t buf_size_guard = CHUNK_SIZE_IN_BYTES;
		uint64_t total_size = total_reads_count * (read_length + 1);
		buffer.reserve(total_size < buf_size_guard ?total_size :buf_size_guard + (read_length + 1));
		for (int i = p-1; i < total_reads_count; i += 2) {
			if (buffer.size() > buf_size_guard) {
				os << buffer;
				buffer.resize(0);
			}
			int idx = paired_idx[i];
			if (idx < hq_reads_count) hq_pg->get_mis_read(idx, (char*) read1.data());
			else {
				if (idx < non_reads_count) lq_pg->get_raw_read(idx - hq_reads_count, (char*) read1.data());
				else n_pg->get_raw_read(idx - non_reads_count, (char*) read1.data());
				if (p == 2) rev_comp_in_place((char*) read1.data(), read_length);
			}
			buffer.append(read1);
			buffer.push_back('\n');
		}
		os << buffer;
		os.close();
	}
}

template<typename uint_pg_len>
void Manager::write_all_reads_ORD(const std::string &out_fn, const vector<uint_pg_len> &pos_pg) const {
	int parts = single_end_mode ?1 :2;
	for (int p = 1; p <= parts; p++) {
		string out_filename = out_fn + (single_end_mode ?"" :"_" + to_string(p));
		fstream os(out_filename, ios::out | ios::binary | ios::trunc);
		string buffer, read1; read1.resize(read_length);
		const uint64_t buf_size_guard = CHUNK_SIZE_IN_BYTES;
		uint64_t total_size = total_reads_count * (read_length + 1);
		buffer.reserve(total_size < buf_size_guard ?total_size :buf_size_guard + (read_length + 1));
		for (int i = (total_reads_count / parts) * (p-1); i < (total_reads_count / parts) * p; i++) {
			if (buffer.size() > buf_size_guard) {
				os << buffer;
				buffer.resize(0);
			}
			uint_pg_len pos = pos_pg[i];
			if (pos < hq_pg_length) hq_pg->get_next_mis_read((char*) read1.data(), pos);
			else {
				if (pos < non_pg_length) lq_pg->get_next_raw_read((char*) read1.data(), pos - hq_pg_length);
				else n_pg->get_next_raw_read((char*) read1.data(), pos - non_pg_length);
				if (p == 2) rev_comp_in_place((char*) read1.data(), read_length);
			}
			buffer.append(read1);
			buffer.push_back('\n');
		}
		os << buffer;
		os.close();
	}
}

template<typename uint_pg_len>
void Manager::apply_rc_pair_to_pg(std::vector<uint_pg_len> &pg_pos) {
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
				hq_pg->get_reads_list()->rev_comp[hq_idx] = !hq_pg->get_reads_list()->rev_comp[hq_idx];
				hq_idx++;
			}
		}
	} else {
		// Reverse complement each read from file2
		for (int i = 1; i < total_reads_count; i += 2) {
			int idx = paired_idx[i];
			if (idx < hq_reads_count) {
				hq_pg->get_reads_list()->rev_comp[idx] = !hq_pg->get_reads_list()->rev_comp[idx];
			}
		}
	}
}

void Manager::decompress(const std::string &out_fn) {
	ifstream in(archive_name); assert(in.is_open());
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
	char pgrc_mode = in.get(); assert(0 <= pgrc_mode and pgrc_mode <= 4); // from 0 to 4
	const std::string five_modes[] = {"SE", "PE", "SE_ORD", "PE_ORD", "PE_MIN"};
	fprintf(stderr, "Compression mode: %s\n", five_modes[pgrc_mode].c_str());
	separate_N = (bool) in.get();
	if (pgrc_mode == PGRC_PE_MODE or pgrc_mode == PGRC_ORD_PE_MODE) {
		rev_comp_pair = (bool) in.get();
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
	ignore_pair_order = (pgrc_mode == PGRC_MIN_PE_MODE);

	load_all_PGs(in);

	if (rev_comp_pair) {
		if (joined_pg_len_std) apply_rc_pair_to_pg(pos_pg_32);
		else apply_rc_pair_to_pg(pos_pg_64);
	}

	if (not preserve_order_mode) {
		if (single_end_mode) write_all_reads_SE(out_fn);
		else write_all_reads_PE(out_fn);
	} else {
		if (joined_pg_len_std) write_all_reads_ORD(out_fn, pos_pg_32);
		else write_all_reads_ORD(out_fn, pos_pg_64);
	}

	fprintf(stderr, "Decompressed %d reads in total\n", total_reads_count);

	if (hq_pg) { delete hq_pg; hq_pg = nullptr; }
	if (lq_pg) { delete lq_pg; lq_pg = nullptr; }
	if (n_pg) { delete n_pg; n_pg = nullptr; }
	in.close();
}