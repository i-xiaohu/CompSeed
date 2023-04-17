//
// Created by ixiaohu on 2023/4/17.
//

#include "PseudoGenomeBase.h"
#include <cassert>

PseudoGenomeHeader::PseudoGenomeHeader(std::ifstream &in) {
	std::string header; in >> header; assert(header == PG_HEADER);
	in >> type;
	in >> constant_read_length;
	in >> max_read_length;
	in >> reads_count;
	in >> pg_length;
	in.get();

//	fprintf(stderr, "type: %s", type.c_str());
//	fprintf(stderr, "constant length: %d\n", constant_read_length);
//	fprintf(stderr, "max length: %u\n", max_read_length);
//	fprintf(stderr, "read count: %u\n", reads_count);
//	fprintf(stderr, "PG length: %lu\n", pg_length);
}