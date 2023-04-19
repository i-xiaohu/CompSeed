//
// Created by ixiaohu on 2023/4/19.
//

#ifndef PGRC_LEARN_VARLENDNACODER_H
#define PGRC_LEARN_VARLENDNACODER_H

#include "helper.h"

class VarLenDNACoder {
private:
	const static uint32_t MAX_NUMBER_OF_CODES = UINT8_MAX + 1; // code from 0 to 255
	const static uint8_t MAX_CODE_LENGTH = 4;
	const static uint8_t NOT_FOUND_CODE = 0;
	char code_book[MAX_NUMBER_OF_CODES][MAX_CODE_LENGTH + 1] = {};

	const static size_t VAR_LEN_PROPS_SIZE = 2;

	const static uint32_t CODE_LUT_SIZE = 1U << 27;
	const static uint32_t CODE_LUT_MASK = CODE_LUT_SIZE - 1;
	uint8_t *code_lut = nullptr;

	void init_using(const std::string &codes);

public:
	explicit VarLenDNACoder(const std::string &codes);

	~VarLenDNACoder();

	int decode(unsigned char *dest, size_t exp_dest_len,
			   const unsigned char *src, size_t src_len);

	const static int STATIC_CODES_CODER_PARAM = 0;
	const static int DYNAMIC_CODES_CODER_PARAM = 1;
	static int uncompress(unsigned char *dest, size_t dest_len,
					      unsigned char *src, size_t src_len);
};



#endif //PGRC_LEARN_VARLENDNACODER_H
