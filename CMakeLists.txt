cmake_minimum_required(VERSION 3.10)
project(comp_seed)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")
set(EXE "comp_seed")

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

set(BWA_Package ${cstl} ${FM_index} ${bwalib})
add_executable(bwaidx ${BWA_Package} FM_index/index_main.c)
target_link_libraries(bwaidx -lz -lrt -lpthread)

set(seeding
        seeding/comp_seeding.cpp seeding/comp_seeding.h
        seeding/SST.cpp seeding/SST.h)

add_executable(${EXE} ${BWA_Package} ${seeding} main.cpp)
target_link_libraries(${EXE} -lm -lz -lrt -lpthread)