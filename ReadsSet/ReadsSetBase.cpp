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

	fprintf(stderr, "reads_count = %d\n", reads_count);
	fprintf(stderr, "all_read_length = %ld\n", all_read_length);
	fprintf(stderr, "constant_read_length = %d\n", constant_read_length);
	fprintf(stderr, "min_read_length = %d\n", min_read_length);
	fprintf(stderr, "max_read_length = %d\n", max_read_length);
	fprintf(stderr, "symbols_count = %d\n", symbols_count);
	fprintf(stderr, "symbols_list = %s\n", symbols_list);
}

void ReadsSetProperties::generate_symbol_order() {
	memset(symbols_order, -1, sizeof(symbols_order));
	for (int i = 0; i < symbols_count; i++) {
		symbols_order[symbols_list[i]] = i;
	}
}