//
// Created by ixiaohu on 2023/4/17.
//

#include "ExtendReadsList.h"
#include "../../Utils/LzmaLib.h"

bool ExtendReadsList::move_next() {
	return false;
}


bool ExtendReadsList::rewind() {
	return false;
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
		// TODO: I can not understand the jump array
		std::vector<uint8_t> misCnt2srcIdx(UINT8_MAX, mismatches_count_limit);
		for (uint8_t m = 1; m < mismatches_count_limit; m++) {
			read_value(in, misCnt2srcIdx[m], false);
		}
		// Mismatch buckets, srcs[i] contains mismatch offset for all reads with i mismatches
		std::vector<uint8_t> srcs[UINT8_MAX]; // FIXME: is 8-bit enough for offset?
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
