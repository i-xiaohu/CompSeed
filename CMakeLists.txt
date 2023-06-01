cmake_minimum_required(VERSION 3.10)
project(pgrc_learn)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3")
set(EXE "PBM")

set(LZMA
        Utils/LzmaLib.cpp Utils/LzmaLib.h
        Utils/VarLenDNACoder.cpp Utils/VarLenDNACoder.h
        lzma/7zTypes.h
        lzma/Alloc.c lzma/Alloc.h
        lzma/CpuArch.c lzma/CpuArch.h
        lzma/LzFind.c lzma/LzFind.h
        lzma/LzFindMt.c lzma/LzFindMt.h
        lzma/LzHash.h
        lzma/LzmaDec.c lzma/LzmaDec.h
        lzma/LzmaEnc.c lzma/LzmaEnc.h
        lzma/Ppmd.h
        lzma/Ppmd7.c lzma/Ppmd7.h
        lzma/Ppmd7Dec.c
        lzma/Ppmd7Enc.c
        lzma/Threads.c lzma/Threads.h)

set(Matcher
        Matching/SimplePgMatcher.cpp Matching/SimplePgMatcher.h)

set(PseudoGenome
        PseudoGenome/PseudoGenomeHeader.h
        PseudoGenome/SeparatedPG.cpp PseudoGenome/SeparatedPG.h
        PseudoGenome/ExtendReadsList.cpp PseudoGenome/ExtendReadsList.h)

set(ReadsSet
        ReadsSet/ReadsSetBase.cpp ReadsSet/ReadsSetBase.h)

set(Utils
        Utils/helper.cpp Utils/helper.h)

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

set(PgRC_Package ${LZMA} ${Matcher} ${PseudoGenome} ${ReadsSet} ${Utils})
add_executable(pgrc_dec ${PgRC_Package} manager.cpp manager.h main.cpp)
target_link_libraries(pgrc_dec -lpthread)

add_executable(${EXE} ${BWA_Package} ${PgRC_Package}
        PgRC_BWA_MEM/comp_seeding.cpp PgRC_BWA_MEM/comp_seeding.h
        PgRC_BWA_MEM/bwamem.h PgRC_BWA_MEM/bwamem.c
        PgRC_BWA_MEM/bwamem_extra.c PgRC_BWA_MEM/bwamem_pair.c)
target_link_libraries(${EXE} -lm -lz -lrt -lpthread)