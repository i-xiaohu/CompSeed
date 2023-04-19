//
// Created by ixiaohu on 2023/4/18.
//

#ifndef PGRC_LEARN_LZMALIB_H
#define PGRC_LEARN_LZMALIB_H

#include <vector>
#include <cassert>
#include "helper.h"

#include "../lzma/7zTypes.h"

#define LZMA_PROPS_SIZE 5

const static uint8_t LZMA_CODER  = 1;
const static uint8_t LZMA2_CODER = 2;
const static uint8_t PPMD7_CODER = 3;
const static uint8_t VARLEN_DNA_CODER = 11;
const static uint8_t COMPOUND_CODER_TYPE = 77;


void uncompress(char *dest, size_t dest_len, std::istream &in, size_t src_len, uint8_t coder_type);

void uncompress(char *dest, size_t dest_len, const char *src, size_t src_len, uint8_t coder_type);

template<typename T>
void read_compressed(std::istream &in, std::vector<T> &dest) {
	size_t dest_len = 0;
	read_value<uint64_t>(in, dest_len, false);
	if (dest_len == 0) return ;
	assert(dest_len % sizeof(T) == 0);
	dest.resize(dest_len / sizeof(T));
	size_t  src_len = 0;
	read_value<uint64_t>(in, src_len, false);
	uint8_t coder_type = 0;
	read_value<uint8_t>(in, coder_type, false);
	uncompress((char*)dest.data(), dest_len, in, src_len, coder_type);
}

void read_compressed(std::istream &in, std::string &dest, int level);

#endif //PGRC_LEARN_LZMALIB_H
