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

	// Iteration variables
	uint64_t cur_pos = 0;
	uint64_t next_idx = 0;
	uint64_t cur_mis_cnt = 0;

public:
	SeparatedPseudoGenome(std::string &&pg_seq, ExtendReadsList *rl, ReadsSetProperties *prop);

	~SeparatedPseudoGenome();

	/** Fetch next read from PG_hq[pos, pos+L) with mismatches */
	void get_next_mis_read(char *ptr, uint64_t pos);

	/** Fetch next read from PG_hq with mismatches in sorted order */
	void get_next_mis_read(char *ptr);

	/** Fetch next read from PG_lq[pos, pos+L) or Pg_n[pos, pos+L) with no mismatch */
	void get_next_raw_read(char *ptr, uint64_t pos);

	/** Fetch next read from PG_lq or PG_n without mismatch in sorted order */
	void get_next_raw_read(char *ptr);
};

void spg_decompress_reads_order(
		std::istream &in,
		std::vector<int> idx_order,
		bool complete_order_info,
		bool take_pe_as_se,
		bool single_end_mode);

template<typename uint_pg_len>
void decompress_pg_position(std::istream &in, std::vector<uint_pg_len> &pg_pos,
							int reads_count, bool se_mode);

#endif //PGRC_LEARN_SEPARATEDPG_H
