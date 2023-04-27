//
// Created by ixiaohu on 2023/4/19.
//

#include <cstring>
#include <algorithm>
#include "SeparatedPG.h"
#include "../Utils/LzmaLib.h"

SeparatedPseudoGenome::SeparatedPseudoGenome(
		std::string &&pg_seq, ExtendReadsList *rl,ReadsSetProperties *prop) {
	this->length = pg_seq.length();
	this->properties = prop;
	this->pg_sequence = std::move(pg_seq);
	this->reads_list = rl;
	if (this->reads_list) {
		this->reads_list->reads_count = prop->reads_count;
	}
}

SeparatedPseudoGenome::~SeparatedPseudoGenome() {
	delete reads_list;
	delete properties;
}

void SeparatedPseudoGenome::get_next_mis_read(char *ptr, uint64_t pos) {
	memcpy(ptr, (pg_sequence.data() + pos), reads_list->read_length);
	if (reads_list->rev_comp[next_idx]) {
		rev_comp_in_place(ptr, reads_list->read_length);
	}
	uint8_t mis_cnt = reads_list->mis_cnt[next_idx];
	for (uint8_t i = 0; i < mis_cnt; i++) {
		uint8_t mis_pos = reads_list->mis_off[cur_mis_cnt];
		ptr[mis_pos] = code_to_mismatch(ptr[mis_pos], reads_list->mis_sym_code[cur_mis_cnt]);
		cur_mis_cnt++;
	}
	next_idx++;
}

void SeparatedPseudoGenome::get_next_mis_read(char *ptr) {
	cur_pos += reads_list->off[next_idx];
	get_next_mis_read(ptr, cur_pos);
}

void SeparatedPseudoGenome::get_next_raw_read(char *ptr, uint64_t pos) {
	memcpy(ptr, (pg_sequence.data() + pos), reads_list->read_length);
}

void SeparatedPseudoGenome::get_next_raw_read(char *ptr) {
	cur_pos += reads_list->off[next_idx];
	memcpy(ptr, (pg_sequence.data() + cur_pos), reads_list->read_length);
	next_idx++;
}


void SeparatedPseudoGenome::get_mis_read(int idx, char *ptr) {
	memcpy(ptr, pg_sequence.data() + reads_list->pos[idx], reads_list->read_length);
	if (not reads_list->rev_comp.empty() and reads_list->rev_comp[idx]) {
		rev_comp_in_place(ptr, reads_list->read_length);
	}

	if (not reads_list->mis_cum_count.empty() or not reads_list->mis_off.empty()) {
		for (int i = reads_list->mis_cum_count[idx]; i < reads_list->mis_cum_count[idx+1]; i++) {
			uint8_t mis_pos = reads_list->mis_off[i];
			uint8_t mis_sym = reads_list->mis_sym_code[i];
			ptr[mis_pos] = code_to_mismatch(ptr[mis_pos], mis_sym);
		}
	}
}

void SeparatedPseudoGenome::get_raw_read(int idx, char *ptr) {
	memcpy(ptr, pg_sequence.data() + reads_list->pos[idx], reads_list->read_length);
}

void restore_paired_idx(std::istream &in, std::vector<int> &paired_idx) {
	// 1. offset_flag is a binary flag telling if the distance between two reads from the pair
	// is small enough to fit a byte; the distance is expressed in the number of aligned reads
	// between them (not in bases)
	std::vector<uint8_t> offset_flag;
	// 2. If offset_flag is positive, the value [0,255] goes to the offset_value stream.
	std::vector<uint8_t> offset_value;
	// 3. Otherwise, a binary delta_flag to tell if the difference between the current pair distance
	// and the "referential pair" distance fits a byte [-128,127]. If the distance between two reads of
	// a pair is large, then it happens quite often that such pairs of reads occur in ADJACENT positions,
	// considering their relative order along the pseudo-genome. Such a phenomenon usually corresponds to
	// the case of having relatively long matching strings in another region of the pseudo-genome. This is
	// why a referential pair is required. The distance between a paired reads does not change frequently.
	std::vector<uint8_t> delta_flag;
	// 4. If delta_flag is positive, the mentioned difference goes to the delta_value stream
	std::vector<int8_t>  delta_value;
	// 5. Otherwise, the current pair distance, represented in a 32-bit integer, goes to the full_offset stream
	std::vector<int32_t> full_offset;
	// 6. offset_pair_base is a binary flag for a pair of reads close enough (flag "true" in Stream #1), telling
	// if the first read of the pair, in the order of their location, comes from the file _1; usually, it is true.
	std::vector<uint8_t> close_flag;
	// 7. a binary flag for a pair of reads not close enough (flag "false" in Stream #1), telling if the first
	// read of the pair, in the order of their location, comes from the file _1. Stream #6 and #7 have different
	// distribution, so the classification can improve encoding.
	std::vector<uint8_t> far_flag;

	read_compressed(in, offset_flag);
	read_compressed(in, offset_value);
	read_compressed(in, delta_flag);
	read_compressed(in, delta_value);
	read_compressed(in, full_offset);
	read_compressed(in, close_flag);
	read_compressed(in, far_flag);

	int reads_count = (int)offset_flag.size() * 2;
	paired_idx.resize(reads_count);
	std::vector<bool> is_read_done(reads_count, false);

	long off_flag_idx = -1, del_flag_idx = -1;
	long off_idx = -1, del_idx = -1;
	long full_idx = -1;
	// "referential pair", the term defined above. It denotes either the currently last pair written to
	// Stream #4, or the currently last pair written to Stream #5 provided that there has been at least
	// two pairs written to Stream #5 after the last pair written to Stream #4.
	long referential_pair = 0;
	// Whether the current pair matches the referential pair
	bool match = false;
	long prev = 0;
	for (int i = 0; i < reads_count; i++) { // Enumerating reads in order of their PG locations.
		if (is_read_done[i]) continue;
		long pair_offset; // The distance (in number of reads) between i and its paired read
		if (offset_flag[++off_flag_idx]) { // Stream #1
			pair_offset = offset_value[++off_idx]; // Stream #2
		} else if (delta_flag[++del_flag_idx]) { // Stream #3
			pair_offset = referential_pair + delta_value[++del_idx]; // Stream #4
			referential_pair = pair_offset; // Updated in Stream #4
			match = true;
			prev = pair_offset;
		} else {
			pair_offset = full_offset[++full_idx]; // Stream #5
			if (not match or referential_pair != prev) referential_pair = pair_offset;
			match = false;
			prev = pair_offset;
		}
		paired_idx[off_flag_idx * 2] = i; // off_flag_idx = counter of pairs
		paired_idx[off_flag_idx * 2 + 1] = i + pair_offset;
		is_read_done[i + pair_offset] = true;
	}

	long close_idx = -1, far_idx = -1;
	for (int p = 0; p < reads_count / 2; p++) {
		bool swap_pair = offset_flag[p] ?close_flag[++close_idx] :far_flag[++far_idx];
		if (swap_pair) std::swap(paired_idx[p * 2], paired_idx[p * 2 + 1]);
	}
}

template<typename uint_pg_len>
void decompress_pg_position(std::istream &in, std::vector<uint_pg_len> &pg_pos,
                            int reads_count, bool se_mode) {
	if (se_mode) {
		read_compressed(in, pg_pos);
		return;
	}

	// Paired-end mode
	bool delta_enabled = (bool) in.get(); // enable delta pair encoding
	std::vector<uint8_t> offset_in_uint16_flag;
	std::vector<uint8_t> offset_is_base_first_flag;
	std::vector<uint16_t> offset_in_uint16_value;
	pg_pos.reserve(reads_count);
	read_compressed(in, pg_pos);
	read_compressed(in, offset_in_uint16_flag);
	read_compressed(in, offset_is_base_first_flag);
	read_compressed(in, offset_in_uint16_value);

	std::vector<uint8_t> delta_in_int16_flag;
	std::vector<uint8_t> delta_is_base_first_flag;
	std::vector<int16_t> delta_in_int16_value;
	if (delta_enabled) {
		read_compressed(in, delta_in_int16_flag);
		read_compressed(in, delta_is_base_first_flag);
		read_compressed(in, delta_in_int16_value);
	}
	std::vector<uint_pg_len> not_base_pair_pos;
	read_compressed(in, not_base_pair_pos);

	// Sort read indexes by PG position
	int pairs_count = reads_count / 2;
	std::vector<int> bpp_rank; bpp_rank.reserve(pairs_count);
	for (int p = 0; p < pairs_count; p++) bpp_rank.push_back(p);
	std::sort(bpp_rank.begin(), bpp_rank.end(),
		   [&](const int &l, const int &r) -> bool { return pg_pos[l] < pg_pos[r]; });

	pg_pos.resize(reads_count);
	for (int i = 0; i < pairs_count; i++) {
		int p = bpp_rank[i];

	}
}

template void decompress_pg_position(std::istream &in, std::vector<uint32_t> &pg_pos,
		int reads_count, bool se_mode);
template void decompress_pg_position(std::istream &in, std::vector<uint64_t> &pg_pos,
		int reads_count, bool se_mode);
