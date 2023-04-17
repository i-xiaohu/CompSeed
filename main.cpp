#include <iostream>
#include <cstring>
#include "manager.h"

int main(int argc, char *argv[]) {
	Manager worker;
	if (!strcmp(argv[1], "dec")) worker.set_archive_name(argv[2]);
	worker.decompress();
	return 0;
}
