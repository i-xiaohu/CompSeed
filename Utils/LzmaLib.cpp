//
// Created by ixiaohu on 2023/4/18.
//


#include "LzmaLib.h"
#include "../lzma/7zTypes.h"
#include "../lzma/LzmaDec.h"
#include "../lzma/Alloc.h"
#include "../lzma/Ppmd7.h"
#include "../Utils/VarLenDNACoder.h"

#define RC_INIT_SIZE 5
#define UNCOMPRESS_BUFFER_SIZE 8192

int lzma_uncompress(unsigned char* dest, size_t *dest_len, std::istream &in, size_t *src_len) {
	size_t props_size = LZMA_PROPS_SIZE;
	unsigned char props_buf[LZMA_PROPS_SIZE];
	read_array(in, (void*) props_buf, props_size);

	SRes res;
	SizeT outSize = *dest_len, inSize = *src_len - props_size;
	*dest_len = *src_len = 0;
	ELzmaStatus status = LZMA_STATUS_NOT_SPECIFIED;
	if (inSize < RC_INIT_SIZE) return SZ_ERROR_INPUT_EOF;

	CLzmaDec p; LzmaDec_Construct(&p);
	RINOK(LzmaDec_AllocateProbs(&p, props_buf, props_size, &g_Alloc));
	p.dic = dest;
	p.dicBufSize = outSize;
	LzmaDec_Init(&p);
	auto *srcBuf = new unsigned char[UNCOMPRESS_BUFFER_SIZE];
	Int64 srcLeftCount = inSize;
	do {
		size_t srcBufSize = UNCOMPRESS_BUFFER_SIZE > srcLeftCount ?srcLeftCount :UNCOMPRESS_BUFFER_SIZE;
		*src_len = srcBufSize;
		read_array(in, srcBuf, srcBufSize);
		srcLeftCount -= srcBufSize;
		res = LzmaDec_DecodeToDic(&p, outSize, srcBuf, src_len, LZMA_FINISH_ANY, &status);
	} while (srcLeftCount > 0 and status == LZMA_STATUS_NEEDS_MORE_INPUT);
	delete [] srcBuf;
	*dest_len = p.dicPos;
	if (res == SZ_OK and status == LZMA_STATUS_NEEDS_MORE_INPUT)
		res = SZ_ERROR_INPUT_EOF;
	LzmaDec_FreeProbs(&p, &g_Alloc);
	*src_len = inSize + props_size;
	return res;
}

class CByteInBufWrap {
public:
	IByteIn vt;
	const Byte *cur = nullptr;
	const Byte *lim = nullptr;
	Byte *buf;
	UInt32 size;
	bool extra;
	int res;

	std::istream *src;
	unsigned char *const src_buf = nullptr;
	Int64 src_left_count;

	bool reload_src_buf();

	CByteInBufWrap(std::istream &src, size_t buf_len);

	~CByteInBufWrap() { delete [] src_buf; }

	UInt64 get_processed() const { return src_buf ?(size - src_left_count) :(cur - buf); }
};

bool CByteInBufWrap::reload_src_buf() {
	if (src_left_count == 0) return false;
	size_t src_buf_size = UNCOMPRESS_BUFFER_SIZE > src_left_count ?src_left_count :UNCOMPRESS_BUFFER_SIZE;
	read_array(*src, src_buf, src_buf_size);
	cur = src_buf;
	lim = src_buf + src_buf_size;
	src_left_count -= src_buf_size;
	return true;
}

static Byte wrap_stream_read_byte(const IByteIn *pp) noexcept {
	CByteInBufWrap *p = CONTAINER_FROM_VTBL_CLS(pp, CByteInBufWrap, vt);
	if (p->cur != p->lim) return *p->cur++;
	if (p->reload_src_buf()) return *p->cur++; // reach to end
	p->res = SZ_ERROR_INPUT_EOF;
	return 0;
}

CByteInBufWrap::CByteInBufWrap(std::istream &src, size_t buf_len):
	src(&src), src_buf(new unsigned char[UNCOMPRESS_BUFFER_SIZE]), size(buf_len)
{
	src_left_count = buf_len;
	reload_src_buf();
	extra = false;
	res = SZ_OK;
	vt.Read = wrap_stream_read_byte; // read one byte
}

int ppmd_uncompress(unsigned char *dest, size_t *dest_len, std::istream &src, size_t *src_len) {
	CPpmd7 ppmd; Ppmd7_Construct(&ppmd);
	size_t props_size = LZMA_PROPS_SIZE;
	unsigned char props_buf[LZMA_PROPS_SIZE];
	read_array(src, (void*) props_buf, props_size);
	uint32_t mem_size = GetUi32(props_buf + 1);
	if (!Ppmd7_Alloc(&ppmd, mem_size, &g_Alloc)) return SZ_ERROR_MEM;

	unsigned int order = props_buf[0];
	Ppmd7_Init(&ppmd, order);
	CPpmd7z_RangeDec rDec;
	Ppmd7z_RangeDec_CreateVTable(&rDec);
	CByteInBufWrap _istream(src, *src_len - props_size);
	rDec.Stream = &_istream.vt;

	int res = SZ_OK;
	if (!Ppmd7z_RangeDec_Init(&rDec)) res = SZ_ERROR_DATA;
	else if (_istream.extra) res = (_istream.res != SZ_OK ?_istream.res :SZ_ERROR_DATA);
	else {
		size_t i;
		for (i = 0; i < *dest_len; i++) {
			int sym = Ppmd7_DecodeSymbol(&ppmd, &rDec.vt);
			if (_istream.extra or sym < 0) break;
			dest[i] = sym;
		}
		if (i != *dest_len) {
			res = (_istream.res != SZ_OK ?_istream.res :SZ_ERROR_DATA);
		} else if (_istream.get_processed() != *src_len - props_size or !Ppmd7z_RangeDec_IsFinishedOK(&rDec)) {
			res = SZ_ERROR_DATA;
		}
	}
	Ppmd7_Free(&ppmd, &g_Alloc);
	return res;
}

void uncompress(char *dest, size_t dest_len, std::istream &in, size_t src_len, uint8_t coder_type) {
	int res = 0;
	size_t out_len = dest_len;
	switch (coder_type) {
		case LZMA_CODER:
			res = lzma_uncompress((unsigned char*)dest, &out_len, in, &src_len);
			break;
		case PPMD7_CODER:
			res = ppmd_uncompress((unsigned char*)dest, &out_len, in, &src_len);
			break;
		default:
			fprintf(stderr, "Unsupported coder type: %d\n", coder_type);
			exit(EXIT_FAILURE);
	}
	assert(out_len == dest_len);
	if (res != SZ_OK) {
		fprintf(stderr, "Error code %d during decompression\n", res);
		exit(EXIT_FAILURE);
	}
}

void uncompress(char *dest, size_t dest_len, const char *src, size_t src_len, uint8_t coder_type) {
	int res;
	size_t out_len = dest_len;
	if (coder_type == VARLEN_DNA_CODER) {
		res = VarLenDNACoder::uncompress((unsigned char*)dest, dest_len,
								         (unsigned char*)src, src_len);
	} else {
		fprintf(stderr, "Unsupported coder type: %d\n", coder_type);
		exit(EXIT_FAILURE);
	}
	assert(out_len == dest_len);
	if (res != SZ_OK) {
		fprintf(stderr, "Error code %d during decompression\n", res);
		exit(EXIT_FAILURE);
	}
}

void read_compressed(std::istream &in, std::string &dest, int level) {
	size_t dest_len = 0; read_value(in, dest_len);
	dest.resize(dest_len);
	if (dest_len == 0) return;
	size_t src_len = 0; read_value(in, src_len);
	uint8_t coder_type = 0; read_value(in, coder_type);
//	fprintf(stderr, "level=%d, dest_len=%ld, src_len=%ld, coder_type=%d\n",
//		 level, dest_len, src_len, coder_type);
	if (coder_type == COMPOUND_CODER_TYPE) {
		read_value(in, coder_type);
		std::string component;
		read_compressed(in, component, level+1);
		assert(src_len == component.length());
		uncompress((char*)dest.data(), dest_len, component.data(), src_len, coder_type);
	} else {
		uncompress((char*)dest.data(), dest_len, in, src_len, coder_type);
	}
}
