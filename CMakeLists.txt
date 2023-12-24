cmake_minimum_required(VERSION 3.10)
project(comp_seed)

include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-mavx512bw" AVX512_SUPPORTED)
check_cxx_compiler_flag("-mavx2" AVX2_SUPPORTED)
check_cxx_compiler_flag("-mavx" AVX_SUPPORTED)
check_cxx_compiler_flag("-msse4.2" SSSE_SUPPORTED)

if (AVX512_SUPPORTED)
    message("AVX512 instructions are supported on this machine.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx512bw")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx512bw")
elseif(AVX2_SUPPORTED)
    message("AVX2 instructions are supported on this machine.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx2")
elseif (AVX_SUPPORTED)
    message("AVX instructions are supported on this machine.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mavx")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mavx")
elseif (SSSE_SUPPORTED)
    message("AVX instructions are supported on this machine.")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -msse -msse2 -msse3 -mssse3 -msse4.1 -msse4.2")
endif()

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")

set(cstl
        cstl/kvec.h
        cstl/khash.h
        cstl/ksort.h
        cstl/kbtree.h
        cstl/kstring.c cstl/kstring.h
        cstl/kthread.c cstl/kthread.h)

set(FM_index
        FM_index/bntseq.c FM_index/bntseq.h
        FM_index/bwt.c FM_index/bwt.h
        FM_index/bwt_gen.c
        FM_index/is.c
        FM_index/QSufSort.c FM_index/QSufSort.h
        FM_index/rle.c FM_index/rle.h
        FM_index/rope.c FM_index/rope.h)

set(bwalib
        bwalib/kseq.h
        bwalib/bwashm.c
        bwalib/bwa.c bwalib/bwa.h
        bwalib/kopen.c bwalib/kopen.h
        bwalib/ksw.c bwalib/ksw.h
        bwalib/utils.c bwalib/utils.h)

add_executable(bwaidx ${FM_index} ${bwalib} FM_index/index_main.c)
target_link_libraries(bwaidx -lz -lrt)

add_executable(bwamem ${cstl} ${FM_index} ${bwalib}
        mapping/bwamem.c mapping/bwamem.h
        mapping/bwamem_pair.c
        mapping/bwamem_extra.c
        mapping/fastmap.c)
target_link_libraries(bwamem -lm -lz -lpthread -lrt)

add_executable(CompSeed ${cstl} ${FM_index} ${bwalib}
        mapping/comp_seed.cpp mapping/comp_seed.h
        mapping/SST.cpp mapping/SST.h
        mapping/bandedSWA.cpp mapping/bandedSWA.h
        mapping/memcpy_bwamem.cpp mapping/memcpy_bwamem.h
        mapping/macro.h
        main.cpp)
target_link_libraries(CompSeed -lm -lz -lpthread -lrt)