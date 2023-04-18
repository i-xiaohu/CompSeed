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

void read_array(std::istream &in, void *dest_array, size_t array_size_in_bytes);

#endif //PGRC_LEARN_HELPER_H
