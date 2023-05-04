#include <iostream>
#include "manager.h"

const char *PROG_NAME = "PgRC-BWA-MEM";

static int usage() {
	fprintf(stderr, "Usage: PBM(%s) <archive> <decompression file>\n", PROG_NAME);
	return 1;
}

int main(int argc, char *argv[]) {
	if (argc == 1) return usage();
	Manager worker;
	worker.set_archive_name(argv[1]);
	worker.decompress(argv[2]);
	return 0;
}
