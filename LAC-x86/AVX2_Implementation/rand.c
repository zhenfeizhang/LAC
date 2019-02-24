#include "rand.h"
#include "fips202.h"
#include "lac_param.h"
#include <string.h>
#include <stdio.h>
#include <sys/types.h>
#include <errno.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <sys/syscall.h>

#define _GNU_SOURCE

static int fd = -1;
static void randombytes_fallback(unsigned char *x, size_t xlen)
{
	int i;
	
	if (fd == -1) {
		for (;;) {
			fd = open("/dev/urandom",O_RDONLY);
			if (fd != -1) break;
			sleep(1);
		}
	}

	while (xlen > 0) {
		if (xlen < 1048576) i = xlen; else i = 1048576;

		i = read(fd,x,i);
		if (i < 1) {
			sleep(1);
			continue;
		}

		x += i;
		xlen -= i;
	}
}

#ifdef SYS_getrandom
void random_bytes(unsigned char *buf,size_t buflen)
{
	size_t d = 0;
	int r;

	while(d<buflen)
	{
		errno = 0;
		r = syscall(SYS_getrandom, buf, buflen - d, 0); 
		if(r < 0) 
		{
			if (errno == EINTR) continue;
			randombytes_fallback(buf, buflen);
			return;
		}
		buf += r;
		d += r;
	}
}
#else
void random_bytes(unsigned char *buf,size_t buflen)
{
	randombytes_fallback(buf,buflen);
}
#endif

int pseudo_random_bytes(unsigned char *r, unsigned int len, const unsigned char *seed)
{
	//check  parameter
	if(r==NULL || seed==NULL)
	{
		return 1;
	}
	shake256(r, len, seed, SEED_LEN);
	return 0;
}

int hash(const unsigned char *in, unsigned int len_in, unsigned char *out)
{
	//check  parameter
	if(in==NULL || out==NULL)
	{
		return 1;
	}
	
	#if defined LAC128
	sha3_256(out,in,len_in);
	#endif
	
	#if defined LAC192
	sha3_256(out,in,len_in);
	#endif

	#if defined LAC256
	sha3_256(out,in,len_in);
	#endif

	
	return 0;
}

int gen_seed(unsigned char *in, unsigned int len_in, unsigned char *out)
{
	//check  parameter
	if(in==NULL || out==NULL)
	{
		return 1;
	}
	
	sha3_256(out,in,len_in);
	return 0;
}
