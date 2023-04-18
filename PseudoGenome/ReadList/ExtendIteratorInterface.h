//
// Created by ixiaohu on 2023/4/17.
//

#ifndef PGRC_LEARN_EXTENDITERATORINTERFACE_H
#define PGRC_LEARN_EXTENDITERATORINTERFACE_H

#include <climits>

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


class ExtendIteratorInterface {
public:
	virtual ~ExtendIteratorInterface() {};
	virtual bool move_next() = 0;
	virtual bool rewind() = 0;

	virtual ReadListEntry& peek_read_entry() = 0;
};

#endif //PGRC_LEARN_EXTENDITERATORINTERFACE_H
