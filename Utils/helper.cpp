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
