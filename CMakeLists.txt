cmake_minimum_required(VERSION 3.10)
project(pgrc_learn)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
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

add_executable(bwa_index
        ${FM_index}
        ${bwalib}
        FM_index/index_main.c)
target_link_libraries(bwa_index -lz -lrt)

add_executable(${EXE}
        ${LZMA}
        ${Matcher}
        ${PseudoGenome}
        ${ReadsSet}
        ${Utils}
        manager.cpp manager.h
        main.cpp)
target_link_libraries(${EXE} -lpthread)