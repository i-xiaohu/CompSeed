//
// Created by ixiaohu on 2021/9/22.
//

#ifndef ZIP_SEEDING_KOPEN_H
#define ZIP_SEEDING_KOPEN_H

void *kopen(const char *fn, int *_fd);
int kclose(void *a);

#endif //ZIP_SEEDING_KOPEN_H
