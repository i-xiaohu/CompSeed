//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_PSEUDOGENOMEBASE_H
#define PGRC_LEARN_PSEUDOGENOMEBASE_H

#include <fstream>
#include <cassert>

class PseudoGenomeHeader {
private:
	std::string PG_HEADER = "PGEN";

	std::string type;
	bool constant_read_length = true;
	int max_read_length = 0;
	int reads_count = 0;
	long pg_length = 0;


public:
	PseudoGenomeHeader() = default;

	explicit PseudoGenomeHeader(std::ifstream &in) {
		std::string header; in >> header; assert(header == PG_HEADER);
		in >> type;
		in >> constant_read_length;
		in >> max_read_length;
		in >> reads_count;
		in >> pg_length;
		in.get();
	}

	int get_max_read_length() const { return max_read_length; }

	int get_reads_count() const { return reads_count; }

	long get_pg_length() const { return pg_length; }
};

#endif //PGRC_LEARN_PSEUDOGENOMEBASE_H
