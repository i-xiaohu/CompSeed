//
// Created by ixiaohu on 2023/4/19.
//

#include "SimplePgMatcher.h"
#include "../Utils/helper.h"
#include "../Utils/LzmaLib.h"

void restore_matched_pg(std::istream &in, long hq_pg_len, std::string &hq_pg_seq,
                        std::string &lq_pg_seq, std::string &n_pg_seq) {
	long hq_mapped_len, lq_mapped_len, n_mapped_len;
	read_value(in, hq_mapped_len, false);
	read_value(in, lq_mapped_len, false);
	read_value(in, n_mapped_len, false);
	std::string combo_pg_mapped;
	read_compressed(in, combo_pg_mapped, 0);

}
