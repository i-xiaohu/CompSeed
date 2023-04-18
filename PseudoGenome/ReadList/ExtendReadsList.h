//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_EXTENDREADSLIST_H
#define PGRC_LEARN_EXTENDREADSLIST_H

#include <string>
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

	explicit ExtendReadsList(int len): read_length(len) {}

	bool move_next() override;
	bool rewind() override;
	ReadListEntry &peek_read_entry() override { return entry; }
};

// Reads List Factories
ExtendReadsList* load_extend_reads_list(
		const std::istream &in,
		int max_read_length,
		int reads_count,
		bool preserve_order_mode = false,
		bool disable_rev_comp = false,
		bool disable_mismatch = false);

#endif //PGRC_LEARN_EXTENDREADSLIST_H
