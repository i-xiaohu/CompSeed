//
// Created by ixiaohu on 2023/11/27.
//

#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <cstdio>
#include <unistd.h>
#include <sys/time.h>
#include <x86intrin.h>

using namespace std::chrono;

static inline long gettimeofday_wrapper() {
	struct timeval tp;
	gettimeofday(&tp, nullptr);
	return tp.tv_sec  * 1000000 + tp.tv_usec;
}

static inline long clock_gettime_wrapper() {
	struct timespec tp;
	clock_gettime(CLOCK_MONOTONIC, &tp);
	return tp.tv_sec * 1000000000 + tp.tv_nsec;
}

int main() {
	long n = 1e8;
	time_point<system_clock, nanoseconds> t_start, t_end;
	nanoseconds elapsed;

	// std::chrono
	t_start = std::chrono::high_resolution_clock::now();
	long num = 0; // For validate
	for (long i = 0; i < n; i++) { // Run time of empty loop is 266ns
		auto x = std::chrono::high_resolution_clock::now();
		auto y = std::chrono::high_resolution_clock::now();
		num += (y - x).count();
	}
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Overhead of std::chrono %.2f ns (verify=%.2f ns)\n", (double)elapsed.count() / n / 2, (double)num / n);

	// gettimeofday
	t_start = std::chrono::high_resolution_clock::now();
	num = 0; // For validate
	for (long i = 0; i < n; i++) {
		auto x = gettimeofday_wrapper();
		auto y = gettimeofday_wrapper();
		num += y - x;
	}
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Overhead of gettimeofday %.2f ns (verify=%.2f ns)\n", (double)elapsed.count() / n / 2, (double)num / n);

	// clock_gettime
	t_start = std::chrono::high_resolution_clock::now();
	num = 0; // For validate
	for (long i = 0; i < n; i++) {
		auto x = clock_gettime_wrapper();
		auto y = clock_gettime_wrapper();
		num += y - x;
	}
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Overhead of clock_gettime %.2f ns (verify=%.2f ns)\n", (double)elapsed.count() / n / 2, (double)num / n);

	// __rdtsc
	uint64_t tim = __rdtsc();
	sleep(1);
	uint64_t proc_freq = __rdtsc() - tim;
	fprintf(stderr, "CPU frequency: %ld\n", proc_freq);
	t_start = std::chrono::high_resolution_clock::now();
	num = 0; // For validate
	for (long i = 0; i < n; i++) {
		uint64_t x = __rdtsc();
		uint64_t y = __rdtsc();
		num += y - x;
	}
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Overhead of __rdtsc %.2f ns (verify=%.2f cycles or %.2f ns)\n", (double)elapsed.count() / n / 2, (double)num / n, (double)num / n / proc_freq * 1e9);

	return 0;
}