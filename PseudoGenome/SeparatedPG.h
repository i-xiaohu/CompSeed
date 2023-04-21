//
// Created by ixiaohu on 2023/4/19.
//

#ifndef PGRC_LEARN_SEPARATEDPG_H
#define PGRC_LEARN_SEPARATEDPG_H

#include <iostream>
#include <vector>
#include "ReadList/ExtendReadsList.h"

class SeparatedPseudoGenome {
private:
	long length;
	std::string pg_sequence;
	ExtendReadsList *reads_list = nullptr;
	ReadsSetProperties *properties = nullptr;

public:
	SeparatedPseudoGenome(std::string &&pg_seq, ExtendReadsList *rl, ReadsSetProperties *prop);

	~SeparatedPseudoGenome();
};

void spg_decompress_reads_order(
		std::istream &in,
		std::vector<int> idx_order,
		bool complete_order_info,
		bool take_pe_as_se,
		bool single_end_mode);

#endif //PGRC_LEARN_SEPARATEDPG_H
