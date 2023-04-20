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

std::string reverse_complement(const std::string &seq) {
	if (not com_tab_init) {
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
	std::string res; res.resize(seq.length());
	size_t i = seq.length();
	for (auto c : seq) {
		res[--i] = complement_table[c];
	}
	return res;
}
