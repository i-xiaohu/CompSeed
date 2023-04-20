//
// Created by ixiaohu on 2023/4/19.
//

#include <sstream>
#include "SimplePgMatcher.h"
#include "../Utils/helper.h"
#include "../Utils/LzmaLib.h"

std::string restore_one_matched_pg(std::string &hq_pg_seq, size_t org_pg_len, const std::string &mapped_pg,
                                   std::istream &off_src, std::istream &len_src,
                                   bool rev_comp, bool text_mode, const std::string &name_pg) {
	std::string new_space; // Create new space for resulted PG sequence
	auto &res_pg = name_pg == "HQ" ? hq_pg_seq : new_space; // If PG_hq, gradually append mapped substrings to itself
	uint32_t min_match_length = 0;
	read_uint_byte_frugal(len_src, min_match_length);
	uint64_t start_pos = 0, mark_pos;
	while ((mark_pos = mapped_pg.find(MATCH_MARK, start_pos)) != std::string::npos) {
		res_pg.append(mapped_pg, start_pos, mark_pos - start_pos);
		start_pos = mark_pos + 1;
		uint64_t match_src_pos = 0;
		if (org_pg_len <= UINT32_MAX) {
			uint32_t tmp;
			read_value(off_src, tmp, text_mode);
			match_src_pos = tmp;
		} else {
			read_value(off_src, match_src_pos, text_mode);
		}

		uint64_t match_length = 0;
		read_uint_byte_frugal(len_src, match_length);
		match_length += min_match_length;
		if (rev_comp) {
			res_pg.append(reverse_complement(hq_pg_seq.substr(match_src_pos, match_length)));
		} else {
			res_pg.append(hq_pg_seq.substr(match_src_pos, match_length));
		}
	}
	res_pg.append(mapped_pg, start_pos, mapped_pg.length() - start_pos);
	fprintf(stderr, "Restore %s PG sequence of length: %ld\n", name_pg.c_str(), res_pg.length());
	return res_pg;
}

void restore_all_matched_pg(std::istream &in, long org_hq_len, std::string &hq_pg_seq,
                            std::string &lq_pg_seq, std::string &n_pg_seq) {
	uint64_t hq_mapped_len, lq_mapped_len, n_mapped_len;
	read_value(in, hq_mapped_len, false);
	read_value(in, lq_mapped_len, false);
	read_value(in, n_mapped_len, false);
	std::string combo_pg_mapped; // Including HQ, LQ, N mapped PG
	read_compressed(in, combo_pg_mapped, 0);
	std::string n_pg_mapped(combo_pg_mapped, hq_mapped_len + lq_mapped_len);
	combo_pg_mapped.resize(hq_mapped_len + lq_mapped_len);
	std::string lq_pg_mapped(combo_pg_mapped, hq_mapped_len);
	combo_pg_mapped.resize(hq_mapped_len);
	std::string hq_pg_mapped(std::move(combo_pg_mapped));

	// Restore the three mapped PGs using the mapping offset and length
	std::string pg_map_off, pg_map_len;
	std::istringstream map_off_src, map_len_src;
	read_compressed(in, pg_map_off, 0);
	read_compressed(in, pg_map_len, 0);
	map_off_src.str(pg_map_off); map_len_src.str(pg_map_len);
	hq_pg_seq = restore_one_matched_pg(hq_pg_seq, org_hq_len,
	                                   hq_pg_mapped, map_off_src, map_len_src,
	                                   true, false, "HQ");

	read_compressed(in, pg_map_off, 0);
	read_compressed(in, pg_map_len, 0);
	map_off_src.str(pg_map_off); map_len_src.str(pg_map_len);
	lq_pg_seq = restore_one_matched_pg(hq_pg_seq, org_hq_len,
									lq_pg_mapped, map_off_src, map_len_src,
									true, false, "LQ");

	if (n_mapped_len == 0) return;
	read_compressed(in, pg_map_off, 0);
	read_compressed(in, pg_map_len, 0);
	map_off_src.str(pg_map_off); map_len_src.str(pg_map_len);
	n_pg_seq = restore_one_matched_pg(hq_pg_seq, org_hq_len,
								   n_pg_mapped, map_off_src, map_len_src,
								   true, false, "N ");
}
