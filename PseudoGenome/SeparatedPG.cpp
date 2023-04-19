//
// Created by ixiaohu on 2023/4/19.
//

#include "SeparatedPG.h"
#include "../Utils/LzmaLib.h"

void spg_decompress_reads_order(
		std::istream &in,
		std::vector<int> idx_order,
		bool complete_order_info,
		bool take_pe_as_se,
		bool single_end_mode) {

	if (single_end_mode) {
		if (not complete_order_info) return;
		read_compressed(in, idx_order);
		fprintf(stderr, "idx_order: %ld\n", idx_order.size());
	}
}
