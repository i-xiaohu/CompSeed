//
// Created by ixiaohu on 2023/4/17.
//

#include "ExtendReadsList.h"
#include "../Utils/LzmaLib.h"

void ExtendReadsList::enable_constant_access(bool disable_iteration_mode) {
	if (pos.empty()) {
		pos.reserve(reads_count + 1);
		long curr_pos = 0;
		for (int i = 0; i < reads_count; i++) {
			curr_pos += off[i];
			pos.push_back(curr_pos);
		}
		pos.push_back(pos.back() + read_length);
	}

	if (disable_iteration_mode) off.clear();
	if (not mis_cnt.empty()) {
		mis_cum_count.reserve(reads_count + 1);
		mis_cum_count.push_back(0);
		int cum_count = 0;
		for (int i = 0; i < reads_count; i++) {
			cum_count += mis_cnt[i];
			mis_cum_count.push_back(cum_count);
		}
		if (disable_iteration_mode) mis_cnt.clear();
	}
}

ExtendReadsList* load_extend_reads_list(
		std::istream &in,
		int max_read_length,
		int reads_count,
		bool preserve_order_mode,
		bool disable_rev_comp,
		bool disable_mismatch) {

	auto *res = new ExtendReadsList(max_read_length);
	if (not preserve_order_mode) {
		read_compressed(in, res->off);
	}
	if (not disable_rev_comp) {
		read_compressed(in, res->rev_comp);
	}
	if (not disable_mismatch) {
		// Input number of mismatches for all reads
		read_compressed(in, res->mis_cnt);
		// Input types of all mismatches
		read_compressed(in, res->mis_sym_code);
		// The maximum mismatches allowed
		uint8_t mismatches_count_limit = 0;
		read_value(in, mismatches_count_limit, false);
		std::vector<uint8_t> misCnt2srcIdx(UINT8_MAX, mismatches_count_limit);
		for (uint8_t m = 1; m < mismatches_count_limit; m++) {
			read_value(in, misCnt2srcIdx[m], false);
		}
		// Mismatch buckets, srcs[i] contains mismatch offset for all reads with i mismatches
		std::vector<uint8_t> srcs[UINT8_MAX];
		std::vector<int> src_counter(UINT8_MAX, 0);
		for (uint8_t m = 1; m <= mismatches_count_limit; m++) {
			read_compressed(in, srcs[m]);
		}
		// Restore mismatch offset for all reads
		res->mis_off.reserve(res->mis_sym_code.size());
		for (int i = 0; i < reads_count; i++) { // Lookup mismatch offsets by mismatch number
			uint8_t mis_cnt = res->mis_cnt[i];
			uint8_t src_idx = misCnt2srcIdx[mis_cnt];
			uint64_t mis_off_start_idx = res->mis_off.size();
			for (uint8_t m = 0; m < mis_cnt; m++) {
				res->mis_off.push_back(srcs[src_idx][src_counter[src_idx]++]);
			}
			// Convert the differential offset to 1-based position from 3' to 5'end
			convert_mis_rev_offset<uint8_t>(
					res->mis_off.data() + mis_off_start_idx,
					res->mis_cnt[i],
					res->read_length);
		}
	}
	fprintf(stderr, "Load PG reads list containing %d reads.\n", reads_count);
	return res;
}
