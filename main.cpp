//
// Created by ixiaohu on 2023/6/22.
//

#include <cstdio>
#include <getopt.h>
#include "seeding/comp_seeding.h"

static void print_usage(const mem_opt_t *opt) {
	fprintf(stderr, "This program implements the compressive version of BWA-MEM seeding, and benchmarks the speedup over\n");
	fprintf(stderr, "BWA-MEM seeding under different compressors, including SPRING, Minicom and PgRC.\n");
	fprintf(stderr, "Usage: comp_seed [options] <FM-index> <Reordered Reads>\n");
	fprintf(stderr, "comp_seed supports all built-in parameters of BWA-MEM seeding and generates identical seeds with BWA-MEM.\n");
	fprintf(stderr, "    -t INT        number of threads [%d]\n", opt->n_threads);
	fprintf(stderr, "    -k INT        minimum seed length [%d]\n", opt->min_seed_len);
	fprintf(stderr, "    -r FLOAT      look for internal seeds inside a seed longer than {-k} * FLOAT [%g]\n", opt->split_factor);
	fprintf(stderr, "    -y INT        seed occurrence for the 3rd round seeding [%ld]\n", (long)opt->max_mem_intv);
	fprintf(stderr, "    -c INT        skip seeds with more than INT occurrences [%d]\n", opt->max_occ);
	fprintf(stderr, "    -K INT        process INT input bases in each batch regardless of nThreads (for reproducibility)\n");
	fprintf(stderr, "Use option --bwa or --comp to specify the original or compressive seeding algorithm.\n");
	fprintf(stderr, "Use option --print to print all seeds to stdout. The output is in text mode, where a seed structure\n"
				    "contains start positions on query and reference sequence and the seed length. It allows the validation\n"
					"whether comp_seed produce identical seeds with BWA-MEM.\n");
	fprintf(stderr, "comp_seed takes the decompressed/reordered reads as input, not directly supporting compression format yet.\n");
}

int main(int argc, char *argv[]) {
	mem_opt_t *opt = mem_opt_init();
	if (argc == 1) { print_usage(opt); free(opt); return 1; }

	struct option long_opts[]  = {
		{"print", no_argument, nullptr, 3}
	};
	const char short_opts[] = "t:k:r:y:c:K";
	CompAligner worker;
	int fixed_chunk_size = 0;
	while (true) {
		int index_ptr, c;
		c = getopt_long(argc, argv, short_opts, long_opts, &index_ptr);
		if (c < 0) break;
		else if (c == 3) worker.print_seed = true;
		else if (c == 't') opt->n_threads = strtol(optarg, nullptr, 10);
		else if (c == 'k') opt->min_seed_len = strtol(optarg, nullptr, 10);
		else if (c == 'r') opt->split_factor = strtof(optarg, nullptr);
		else if (c == 'y') opt->max_mem_intv = strtol(optarg, nullptr, 10);
		else if (c == 'c') opt->max_occ = strtol(optarg, nullptr, 10);
		else if (c == 'K') fixed_chunk_size = strtol(optarg, nullptr, 10);
		else {
			fprintf(stderr, "Unrecognized option\n");
			exit(EXIT_FAILURE);
		}
	}
	worker.opt = opt;
	worker.actual_chunk_size = fixed_chunk_size == 0 ?opt->n_threads * opt->chunk_size :fixed_chunk_size;
	worker.load_index(argv[optind]);
	worker.run(argv[optind + 1]);
	free(opt);
	return 0;
}