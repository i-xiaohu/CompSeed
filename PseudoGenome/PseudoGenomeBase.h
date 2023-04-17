//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_PSEUDOGENOMEBASE_H
#define PGRC_LEARN_PSEUDOGENOMEBASE_H

#include <fstream>

class PseudoGenomeHeader {
private:
	const std::string PG_HEADER = "PGEN";

	std::string type;
	bool constant_read_length = true;
	uint16_t max_read_length;
	uint32_t reads_count;
	uint64_t pg_length;


public:
	PseudoGenomeHeader(std::ifstream &in);
};

#endif //PGRC_LEARN_PSEUDOGENOMEBASE_H
