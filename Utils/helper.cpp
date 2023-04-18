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
