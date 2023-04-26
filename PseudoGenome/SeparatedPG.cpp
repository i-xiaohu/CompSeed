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

void SeparatedPseudoGenome::get_next_mis_read(char *ptr, uint64_t pos) {
	memcpy(ptr, (pg_sequence.data() + pos), reads_list->read_length);
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

void SeparatedPseudoGenome::get_next_mis_read(char *ptr) {
	cur_pos += reads_list->off[next_idx];
	get_next_mis_read(ptr, cur_pos);
}

void SeparatedPseudoGenome::get_next_raw_read(char *ptr, uint64_t pos) {
	memcpy(ptr, (pg_sequence.data() + pos), reads_list->read_length);
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

template<typename uint_pg_len>
void decompress_pg_position(std::istream &in, std::vector<uint_pg_len> &pg_pos,
                            int reads_count, bool se_mode) {
	if (se_mode) {
		read_compressed(in, pg_pos);
		return;
	}

	// Paired-end mode
	bool delta_enabled = (bool) in.get(); // enable delta pair encoding
	std::vector<uint8_t> offset_in_uint16_flag;
	std::vector<uint8_t> offset_is_base_first_flag;
	std::vector<uint16_t> offset_in_uint16_value;
	pg_pos.reserve(reads_count);
	read_compressed(in, pg_pos);
	read_compressed(in, offset_in_uint16_flag);
	read_compressed(in, offset_is_base_first_flag);
	read_compressed(in, offset_in_uint16_value);

	std::vector<uint8_t> delta_in_int16_flag;
	std::vector<uint8_t> delta_is_base_first_flag;
	std::vector<int16_t> delta_in_int16_value;
	if (delta_enabled) {
		read_compressed(in, delta_in_int16_flag);
		read_compressed(in, delta_is_base_first_flag);
		read_compressed(in, delta_in_int16_value);
	}

	std::vector<uint_pg_len> not_base_pair_pos;
	read_compressed(in, not_base_pair_pos);
	int pairs_count = reads_count / 2;
	std::vector<int> bpp_rank; bpp_rank.reserve(pairs_count);
	for (int p = 0; p < pairs_count; p++) bpp_rank.push_back(p);
}

template void decompress_pg_position(std::istream &in, std::vector<uint32_t> &pg_pos,
		int reads_count, bool se_mode);
template void decompress_pg_position(std::istream &in, std::vector<uint64_t> &pg_pos,
		int reads_count, bool se_mode);
