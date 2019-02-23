#ifndef RAND_H
#define RAND_H

#define _GNU_SOURCE

#include <unistd.h>

void random_bytes(unsigned char *x, size_t xlen);
int pseudo_random_bytes(unsigned char *r, unsigned int len, const unsigned char *seed);
int hash(const unsigned char *in, unsigned int len_in, unsigned char * out);
int gen_seed(unsigned char *in, unsigned int len_in, unsigned char * out);

#endif
