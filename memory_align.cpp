//
// Created by ixiaohu on 2023/11/27.
//

#include <cstdint>
#include <cstdlib>
#include <chrono>
#include <cstdio>

#include "mm_malloc.h"
using namespace std::chrono;


void granularity_1(uint32_t n, void *a) {
	auto data8 = (uint8_t *)a;
	auto dataEnd = data8 + n;
	while (data8 != dataEnd) {
		*data8 = ~(*data8);
		data8++;
	}
}

void granularity_2(uint32_t n, void *a) {
	auto data16 = (uint16_t *)a;
	auto data16End = data16 + (n >> 1U);
	auto data8 = (uint8_t*)data16End;
	auto data8End = data8 + (n & 1U);
	while (data16 != data16End) {
		*data16 = ~(*data16);
		data16++;
	}
	while (data8 != data8End) {
		*data8 = ~(*data8);
		data8++;
	}
}

void granularity_q(uint32_t n, void *a) {
	auto data32 = (uint32_t *)a;
	auto data32End = data32 + (n >> 2U);
	auto data8 = (uint8_t*)data32End;
	auto data8End = data8 + (n & 3U);
	while (data32 != data32End) {
		*data32 = ~(*data32);
		data32++;
	}
	while (data8 != data8End) {
		*data8 = ~(*data8);
		data8++;
	}
}

void granularity_8(uint32_t n, void *a) {
	auto data64 = (uint64_t *)a;
	auto data64End = data64 + (n >> 3U);
	auto data8 = (uint8_t*)data64End;
	auto data8End = data8 + (n & 7U);
	while (data64 != data64End) {
		*data64 = ~(*data64);
		data64++;
	}
	while (data8 != data8End) {
		*data8 = ~(*data8);
		data8++;
	}
}

int main() {
	void *a;
	uint32_t n = 800 * 1024 * 1024;
	time_point<system_clock, nanoseconds> t_start, t_end;
	nanoseconds elapsed;

	// Unaligned memory block
	a = (uint8_t*) malloc((n + 1) * sizeof(uint8_t));
	a = (uint8_t*)a + 1; // 1 byte bias
	fprintf(stderr, "Pointer %% 8 = %lu\n", (uint64_t)a % 8);
	// Memory granularity: 1 byte
	t_start = std::chrono::high_resolution_clock::now();
	granularity_1(n, a);
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Granularity 1: %.0f milliseconds\n", elapsed.count() / 1000.0);
	// Memory granularity: 2 byte
	t_start = std::chrono::high_resolution_clock::now();
	granularity_2(n, a);
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Granularity 2: %.0f milliseconds\n", elapsed.count() / 1000.0);
	// Memory granularity: 4 byte
	t_start = std::chrono::high_resolution_clock::now();
	granularity_q(n, a);
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Granularity 4: %.0f milliseconds\n", elapsed.count() / 1000.0);
	// Memory granularity: 8 byte
	t_start = std::chrono::high_resolution_clock::now();
	granularity_8(n, a);
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Granularity 8: %.0f milliseconds\n", elapsed.count() / 1000.0);
	free((uint8_t*)a - 1);

	// Aligned memory block
	a = (uint8_t*) _mm_malloc(n * sizeof(uint8_t), 64);
	fprintf(stderr, "Pointer %% 8 = %lu\n", (uint64_t)a % 8);
	// Memory granularity: 1 byte
	t_start = std::chrono::high_resolution_clock::now();
	granularity_1(n, a);
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Granularity 1: %.0f milliseconds\n", elapsed.count() / 1000.0);
	// Memory granularity: 2 byte
	t_start = std::chrono::high_resolution_clock::now();
	granularity_2(n, a);
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Granularity 2: %.0f milliseconds\n", elapsed.count() / 1000.0);
	// Memory granularity: 4 byte
	t_start = std::chrono::high_resolution_clock::now();
	granularity_q(n, a);
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Granularity 4: %.0f milliseconds\n", elapsed.count() / 1000.0);
	// Memory granularity: 8 byte
	t_start = std::chrono::high_resolution_clock::now();
	granularity_8(n, a);
	t_end = std::chrono::high_resolution_clock::now();
	elapsed = t_end - t_start;
	fprintf(stderr, "Granularity 8: %.0f milliseconds\n", elapsed.count() / 1000.0);
	_mm_free(a);

	return 0;
}