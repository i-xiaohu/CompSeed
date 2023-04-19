//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_EXTENDREADSLIST_H
#define PGRC_LEARN_EXTENDREADSLIST_H

#include <string>
#include <vector>
#include "ExtendIteratorInterface.h"
#include "../PseudoGenomeBase.h"
#include "../../ReadsSet/ReadsSetBase.h"

class ExtendReadsList : public ExtendIteratorInterface {
private:
	ReadListEntry entry;
	long current = -1;
	long cur_mis_cum_count = 0;

public:
	int read_length;
	int reads_count;
	std::vector<uint8_t> off;
	std::vector<uint8_t> rev_comp;
	std::vector<uint8_t> mis_cnt; // #mismatches in each read
	std::vector<uint8_t> mis_sym_code;
	std::vector<uint8_t> mis_off; // mismatch offset for all reads

	explicit ExtendReadsList(int len): read_length(len) {}

	bool move_next() override;
	bool rewind() override;
	ReadListEntry &peek_read_entry() override { return entry; }
};

// Reads List Factories
ExtendReadsList* load_extend_reads_list(
		std::istream &in,
		int max_read_length,
		int reads_count,
		bool preserve_order_mode = false,
		bool disable_rev_comp = false,
		bool disable_mismatch = false);

#endif //PGRC_LEARN_EXTENDREADSLIST_H
