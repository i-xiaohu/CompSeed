//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_READSSETBASE_H
#define PGRC_LEARN_READSSETBASE_H

#include <iostream>
#include <climits>
#include <fstream>

class ReadsSetProperties {
public:
	int reads_count = 0;
	long all_read_length = 0;
	bool constant_read_length = true;
	int min_read_length = -1;
	int max_read_length = 0;
	int symbols_count = 0;
	char symbols_list[UCHAR_MAX] = {0};
	int symbols_order[UCHAR_MAX] = {-1};

	ReadsSetProperties(std::ifstream &in);

	void generate_symbol_order();
};

#endif //PGRC_LEARN_READSSETBASE_H
