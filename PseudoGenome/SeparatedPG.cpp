//
// Created by ixiaohu on 2023/4/19.
//

#include <cstring>
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

void SeparatedPseudoGenome::get_next_mis_read(char *ptr) {
	cur_pos += reads_list->off[next_idx];
	memcpy(ptr, (pg_sequence.data() + cur_pos), reads_list->read_length);
	if (reads_list->rev_comp[next_idx]) {
		rev_comp_in_place(ptr, reads_list->read_length);
	}
	uint8_t mis_cnt = reads_list->mis_cnt[next_idx];
	for (uint8_t i = 0; i < mis_cnt; i++) {
		uint8_t mis_pos = reads_list->mis_off[cur_mis_cnt];
		ptr[mis_pos] = code_to_mismatch(ptr[mis_pos], reads_list->mis_sym_code[cur_mis_cnt]);
		cur_mis_cnt++;
	}
	next_idx++;
}

void SeparatedPseudoGenome::get_next_raw_read(char *ptr) {
	cur_pos += reads_list->off[next_idx];
	memcpy(ptr, (pg_sequence.data() + cur_pos), reads_list->read_length);
	next_idx++;
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
