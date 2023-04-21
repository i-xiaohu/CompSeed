//
// Created by ixiaohu on 2023/4/19.
//

#include "SeparatedPG.h"
#include "../Utils/LzmaLib.h"

SeparatedPseudoGenome::SeparatedPseudoGenome(
		std::string &&pg_seq, ExtendReadsList *rl,ReadsSetProperties *prop) {
	this->length = pg_seq.length();
	this->properties = prop;
	this->pg_sequence = std::move(pg_seq);
	this->reads_list = rl;
	if (this->reads_list) {
		this->reads_list->reads_count = prop->reads_count;
	}
}

SeparatedPseudoGenome::~SeparatedPseudoGenome() {
	delete reads_list;
	delete properties;
}

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
