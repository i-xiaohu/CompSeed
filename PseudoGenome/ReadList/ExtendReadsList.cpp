//
// Created by ixiaohu on 2023/4/17.
//

#include "ExtendReadsList.h"

bool ExtendReadsList::move_next() {
	return false;
}


bool ExtendReadsList::rewind() {
	return false;
}

ExtendReadsList* load_extend_reads_list(
		const std::istream &in,
		int max_read_length,
		int reads_count,
		bool preserve_order_mode,
		bool disable_rev_comp,
		bool disable_mismatch) {

	auto *res = new ExtendReadsList(max_read_length);
}
