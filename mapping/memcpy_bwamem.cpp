/*************************************************************************************
                           The MIT License

   BWA-MEM2  (Sequence alignment using Burrows-Wheeler Transform),
   Copyright (C) 2019  Intel Corporation, Heng Li.

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.

Authors: Sanchit Misra <sanchit.misra@intel.com>; Vasimuddin Md <vasimuddin.md@intel.com>;
*****************************************************************************************/

#include "memcpy_bwamem.h"
#include <cstdint>
#include <cstdio>
#include <cstring>

errno_t memcpy_bwamem(void *dest, rsize_t dmax, const void *src, rsize_t smax, char *file_name, int line_num)
{
    errno_t ret;
    int64_t bytes_copied;
    for(bytes_copied = 0; bytes_copied < smax; bytes_copied += RSIZE_MAX_MEM)
    {
        int64_t bytes_remaining = smax - bytes_copied;
        int64_t bytes_to_copy = (bytes_remaining > RSIZE_MAX_MEM) ? RSIZE_MAX_MEM : bytes_remaining;
        int64_t dest_bytes_remaining = dmax - bytes_copied;
        int64_t dest_bytes = (dest_bytes_remaining < bytes_to_copy) ? dest_bytes_remaining : bytes_to_copy;
        if((memcpy((char *)dest + bytes_copied, (const char *)src + bytes_copied, bytes_to_copy)) == 0)
        {
            fprintf(stderr, "[%s: %d] memcpy_s returned %d\n", file_name, line_num, ret);
            exit(EXIT_FAILURE);
        }
    }
    return 0;
}
