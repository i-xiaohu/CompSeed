//
// Created by ixiaohu on 2023/4/17.
//

#include <cstring>

#include "ReadsSetBase.h"

ReadsSetProperties::ReadsSetProperties(std::ifstream &in) {
	in >> reads_count;
	in >> all_read_length;
	in >> constant_read_length;
	in >> min_read_length;
	in >> max_read_length;
	in >> symbols_count;
	in >> symbols_list;
	generate_symbol_order();
}

void ReadsSetProperties::generate_symbol_order() {
	memset(symbols_order, -1, sizeof(symbols_order));
	for (int i = 0; i < symbols_count; i++) {
		symbols_order[symbols_list[i]] = i;
	}
}