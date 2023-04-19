//
// Created by ixiaohu on 2023/4/19.
//

#include "VarLenDNACoder.h"
#include <cstring>

void VarLenDNACoder::init_using(const std::string &codes) {
	// Split input string by '\n' and put them into the code book
	int next_pos = 0, pos = 0, counter = 0;
	while ((next_pos = codes.find('\n', pos)) != std::string::npos and counter <= UINT8_MAX) {
		strncpy(code_book[counter++], codes.data() + pos, next_pos - pos);
		pos = next_pos + 1;
	}
	if (counter > UINT8_MAX) {
		fprintf(stderr, "Too many codes in the var-len coder book\n");
		exit(EXIT_FAILURE);
	}
	strncpy(code_book[counter++], codes.data() + pos, codes.length() - pos);
	if (strlen(code_book[NOT_FOUND_CODE]) != 1) {
		fprintf(stderr, "First code in the var-len coder book should be a single symbol\n");
		exit(EXIT_FAILURE);
	}

	code_lut = new uint8_t[CODE_LUT_SIZE];
	memset(code_lut, NOT_FOUND_CODE, CODE_LUT_SIZE);
	for (int i = 0; i < counter; i++) {
		uint32_t temp = 0;
		memcpy(&temp, code_book[i], MAX_CODE_LENGTH);
		code_lut[temp & CODE_LUT_MASK] = i;
	}
}

VarLenDNACoder::VarLenDNACoder(const std::string &codes) {
	init_using(codes);
}

VarLenDNACoder::~VarLenDNACoder() {
	delete [] code_lut;
}

int VarLenDNACoder::decode(unsigned char *dest, size_t exp_dest_len,
						   const unsigned char *src, size_t src_len) {
	uint8_t *dest_ptr = dest;
	for (size_t i = 0; i < src_len; i++) {
		char *ptr = code_book[src[i]]; // Each byte in the source stream indicates several symbols,
		while (*ptr) *(dest_ptr++) = *(ptr++); // which are encoded in the code book.
	}

	size_t dest_len = dest_ptr - dest;
	if (exp_dest_len != dest_len) {
		fprintf(stderr, "Unexpected decoded length: %ld (expected %ld)\n", dest_len, exp_dest_len);
		exit(EXIT_FAILURE);
	}
	return 0;
}

int VarLenDNACoder::uncompress(unsigned char *dest, size_t dest_len,
							   unsigned char *src, size_t src_len) {
	int coder_param = src[0];
	size_t header_size = VAR_LEN_PROPS_SIZE;
	VarLenDNACoder *coder;
	if (coder_param == STATIC_CODES_CODER_PARAM or
	    coder_param == DYNAMIC_CODES_CODER_PARAM) {
		std::string codes = std::string((char*)(src + header_size));
		coder = new VarLenDNACoder(codes);
		header_size += codes.length() + 1;
	} else {
		fprintf(stderr, "Unsupported var-len DNA coder parameter: %d\n", coder_param);
		exit(EXIT_FAILURE);
	}
	int res = coder->decode(dest, dest_len, src + header_size, src_len - header_size);
	delete (coder);
	return res;
}