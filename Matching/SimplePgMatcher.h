//
// Created by ixiaohu on 2023/4/19.
//

#ifndef PGRC_LEARN_SIMPLEPGMATCHER_H
#define PGRC_LEARN_SIMPLEPGMATCHER_H

#include <iostream>

static char MATCH_MARK = '%';

std::string restore_one_matched_pg(std::string &hq_pg_seq, size_t org_pg_len, const std::string &mapped_pg,
                                   std::istream &off_src, std::istream &len_src,
                                   bool rev_comp, bool text_mode, const std::string &name_pg);

void restore_all_matched_pg(std::istream &in, long org_hq_len, std::string &hq_pg_seq,
                            std::string &lq_pg_seq, std::string &n_pg_seq);

#endif //PGRC_LEARN_SIMPLEPGMATCHER_H
