//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_HELPER_H
#define PGRC_LEARN_HELPER_H

#include <fstream>

bool confirm_text_read_mode(std::ifstream &in);

template<typename t_val>
void read_value(std::istream &in, t_val &value, bool plain_text_read_mode) {
	if (plain_text_read_mode) in >> value;
	else in.read((char*) &value, sizeof(t_val));
}

template<typename uint_read_len>
void convert_mis_rev_offset(uint_read_len *offsets, uint8_t mis_cnt, uint_read_len read_length) {
	// For example, three mismatches encoded as 0 0 2, and the read length is 101.
	// the real mismatched positions are 96, 99, 100
	// 101-(0+1) = 100
	// 100-(0+1) = 99
	// 99-(2+1)  = 96
	// This encoding strategy consider mismatches toward 5'end
	for (int i = 0; i < mis_cnt / 2; i++) {
		std::swap(offsets[i], offsets[mis_cnt - 1 - i]);
	}
	uint_read_len pos = read_length;
	for (int i = (int)mis_cnt - 1; i >= 0; i--) {
		pos -= offsets[i] + 1;
		offsets[i] = pos;
	}
}

void read_array(std::istream &in, void *dest_array, size_t array_size_in_bytes);

#endif //PGRC_LEARN_HELPER_H
