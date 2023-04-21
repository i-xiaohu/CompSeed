//
// Created by ixiaohu on 2023/4/17.
//

#include "helper.h"
#include <cassert>

bool confirm_text_read_mode(std::ifstream &in) {
	std::string read_mode; in >> read_mode; in.get();
	assert(read_mode == "BIN" or read_mode == "TXT");
	return read_mode == "TXT";
}

const size_t CHUNK_SIZE = 10000000;

void read_array(std::istream &in, void *dest_array, size_t array_size_in_bytes) {
	char *ptr = (char*) dest_array;
	size_t bytes_left = array_size_in_bytes;
	while (bytes_left > CHUNK_SIZE) {
		in.read(ptr, CHUNK_SIZE);
		ptr += CHUNK_SIZE;
		bytes_left -= CHUNK_SIZE;
	}
	in.read(ptr, bytes_left);
}


bool com_tab_init = false;
char complement_table[256] = {0};
static void init_com_tab() {
	complement_table['A'] = 'T'; complement_table['a'] = 'T';
	complement_table['C'] = 'G'; complement_table['c'] = 'G';
	complement_table['G'] = 'C'; complement_table['g'] = 'C';
	complement_table['T'] = 'A'; complement_table['t'] = 'A';
	complement_table['N'] = 'N'; complement_table['n'] = 'N';
	complement_table['U'] = 'A'; complement_table['u'] = 'A';
	complement_table['Y'] = 'R'; complement_table['y'] = 'R';
	complement_table['R'] = 'Y'; complement_table['r'] = 'Y';
	complement_table['K'] = 'M'; complement_table['k'] = 'M';
	complement_table['M'] = 'K'; complement_table['m'] = 'K';
	complement_table['B'] = 'V'; complement_table['b'] = 'V';
	complement_table['D'] = 'H'; complement_table['d'] = 'H';
	complement_table['H'] = 'D'; complement_table['h'] = 'D';
	complement_table['V'] = 'B'; complement_table['v'] = 'B';
	com_tab_init = true;
}

std::string reverse_complement(const std::string &seq) {
	if (not com_tab_init) init_com_tab();
	std::string res; res.resize(seq.length());
	size_t i = seq.length();
	for (auto c : seq) {
		res[--i] = complement_table[c];
	}
	return res;
}

void rev_comp_in_place(char *ptr, int len) {
	if (not com_tab_init) init_com_tab();
	for (int i = 0; i < len; i++) {
		ptr[i] = complement_table[ptr[i]];
	}
	for (int i = 0; i < len / 2; i++) {
		std::swap(ptr[i], ptr[len-1-i]);
	}
}

bool sym_init = false;
uint8_t sym2val[256] = {255};
char val2sym[256] = {' '};

char code_to_mismatch(char actual, uint8_t code) {
	if (not sym_init) {
		sym2val['A'] = 0; val2sym[0] = 'A';
		sym2val['C'] = 1; val2sym[1] = 'C';
		sym2val['G'] = 2; val2sym[2] = 'G';
		sym2val['T'] = 3; val2sym[3] = 'T';
		sym2val['N'] = 4; val2sym[4] = 'N';
		sym_init = true;
	}
	return val2sym[code < sym2val[actual] ?code :(code + 1)];
}