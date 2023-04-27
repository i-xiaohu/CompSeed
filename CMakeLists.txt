cmake_minimum_required(VERSION 3.10)
project(pgrc_learn)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

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
        PseudoGenome/PseudoGenomeBase.h
        PseudoGenome/SeparatedPG.cpp PseudoGenome/SeparatedPG.h
        PseudoGenome/ReadList/ExtendIteratorInterface.h
        PseudoGenome/ReadList/ExtendReadsList.cpp PseudoGenome/ReadList/ExtendReadsList.h)

set(ReadsSet
        ReadsSet/ReadsSetBase.cpp ReadsSet/ReadsSetBase.h)

set(Utils
        Utils/helper.cpp Utils/helper.h
)

add_executable(pgrc_learn
        main.cpp
        ${LZMA}
        ${Matcher}
        ${PseudoGenome}
        ${ReadsSet}
        ${Utils}
        manager.cpp manager.h
)

target_link_libraries(pgrc_learn pthread)