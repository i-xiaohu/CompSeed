//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_EXTENDREADSLIST_H
#define PGRC_LEARN_EXTENDREADSLIST_H

#include <string>
#include <vector>
#include "PseudoGenomeHeader.h"
#include "../ReadsSet/ReadsSetBase.h"

class ReadListEntry {
public:
	long pos;
	int offset;
	int idx;
	bool rev_comp = false;
	int mismatches_count = 0;
	int mismatch_code[UCHAR_MAX];
	int mismatch_offset[UCHAR_MAX];

	explicit ReadListEntry(long p = 0, bool rc = false) {
		this->pos = p;
		this->rev_comp = rc;
		this->offset = 0;
	}

	void add_mismatch(int code, int os) {
		this->mismatch_code[this->mismatches_count] = code;
		this->mismatch_offset[this->mismatches_count] = os;
		this->mismatches_count++;
	}

	void advance_entry_by_position(long p, int i, bool rc = false) {
		this->offset = p - this->pos;
		this->pos = p;
		this->idx = i;
		this->rev_comp = rc;
		this->mismatches_count = 0;
	}

	void advance_entry_by_offset(int os, int i, bool rc = false) {
		advance_entry_by_position(this->pos + os, i, rc);
	}
};

class ExtendReadsList {
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

	// Constant access features
	std::vector<long> pos;
	std::vector<int> mis_cum_count;

	/**
	 * Load the exact (not cumulative) location offset and mismatch offset for all reads,
	 * allowing that constant access for every read. It works for paired-end reads where
	 * reads can not fetched by sequential iteration.
	 * Calling this function will increase memory cost.
	 *
	 * @param disable_iteration_mode if true, clear cumulative array.
	 */
	void enable_constant_access(bool disable_iteration_mode);

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
