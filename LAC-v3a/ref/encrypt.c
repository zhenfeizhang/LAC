#include <string.h>
#include "api.h"
#include "rand.h"
#include "bch.h"
#include "ecc.h"
#include "bin-lwe.h"

#define RATIO 125
//compute parity
static inline int parity(const uint16_t *data, uint16_t *parity_c, uint8_t *parity_r)
{
	//8*16 two dimention parity check
	int i;
	uint16_t x[8];
	memcpy(x,data,16);
	//compute colume parity
	parity_c[0]=0;
	for(i=0;i<8;i++)
	{
		parity_c[0]^=x[i];
	}
	
	//compute row parity
	parity_r[0]=0;
	for(i=0;i<8;i++)
	{
		x[i]^=x[i]>>1;
		x[i]^=x[i]>>2;
		x[i]=(x[i] & 0x1111U)*(0x1111U);
		parity_r[0]=parity_r[0]^(((x[i]>>12)&1)<<i);
		
	}
	
	return 0;
}
//error correct using parity
static inline int parity_ecc(uint16_t *data, uint16_t parity_c, uint8_t parity_r)
{
	//16*16 two dimention parity check
	int i;
	uint16_t tmp_c;
	uint8_t  tmp_r;
	
	//compute  parity
	parity(data,&tmp_c,&tmp_r);
	
	//check
	tmp_c^=parity_c;
	tmp_r^=parity_r;
	for(i=0;i<8;i++)
	{
		data[i]^=(((tmp_r>>i)&1)*tmp_c);
	}
	
	return 0;
}

static int d2_decode(unsigned char *v1, unsigned char *v2, unsigned char *data, int len)
{
	int i;
	int temp1,temp2;
	int d2_bound=Q/2;
	int center_point=Q/2;
	int vec_bound=len/2;
	
	for(i=0;i<vec_bound;i++)
	{
		//D2 decoding:compute m*q/2+e1 + m*q/2+e2 in [0,2*Q]
		temp1=(v1[i]-v2[i]+Q)%Q;
		temp2=(v1[i+vec_bound]-v2[i+vec_bound]+Q)%Q;
		
		//shift
		if(temp1<center_point)
		{
			temp1=center_point-temp1+center_point;//mirror around Q/2
		}
		if(temp2<center_point)
		{
			temp2=center_point-temp2+center_point;//mirror around Q/2
		}
		//merge erors
		temp1+=temp2-Q;

		//recover m from m*q/2+e1 + m*q/2+e2, RATIO=q/2
		if(temp1<d2_bound)
		{
			data[i/8]=data[i/8]^(1<<(i%8));
		}
	}
	
	return 0;
}

//key generation
int crypto_encrypt_keypair( unsigned char *pk, unsigned char *sk)
{
	//check parameter
	if(pk==NULL || sk==NULL)
	{
		return -1;
	}
	kg(pk,sk);
	
	return 0;
}

//encryption
int crypto_encrypt( unsigned char *c, unsigned long long *clen, const unsigned char *m, unsigned long long mlen, const unsigned char *pk)
{
	//check parameter
	if(c==NULL || m==NULL || pk==NULL)
	{
		return -1;
	}
	if(mlen>MESSAGE_LEN)
	{
		return -1;
	}
	
	//call pke encryption function
	pke_enc(pk,m, mlen,c,clen);

	return 0;
}
//decryption
int crypto_encrypt_open(unsigned char *m, unsigned long long *mlen,const unsigned char *c, unsigned long long clen,const unsigned char *sk)
{
	//check parameter
	if(sk==NULL || m==NULL || c==NULL || mlen==NULL)
	{
		return -1;
	}
	
	//call pke decryption function
	pke_dec(sk,c,clen,m,mlen);

	return 0;
}

//key generation with seed
int kg_seed(unsigned char *pk, unsigned char *sk, unsigned char *seed)
{
	unsigned char seeds[3*SEED_LEN];
	unsigned char a[DIM_N];
	char e[DIM_N];
	//check pointer
	if(pk==NULL || sk==NULL)
	{
		return -1;
	}
	
	//generate three seeds for a,sk,e
	pseudo_random_bytes(seeds,3*SEED_LEN,seed);	
	
	//copy the seed of a to pk
	memcpy(pk,seeds,SEED_LEN);
	
	//generate a
	gen_a(a,pk);//print_bytes(sk,CRYPTO_SECRETKEYBYTES);
	
	//generate  sk,e
	gen_r((char*)sk,seeds+SEED_LEN);
	gen_e((char*)e,seeds+2*SEED_LEN);
	//compute pk=a*sk+e
	poly_aff(a,(char *)sk,e,pk+SEED_LEN,DIM_N);
	//copy pk=as+e to the second part of sk, now sk=s|pk
	memcpy(sk+DIM_N,pk,PK_LEN);
	
	return 0;
}

//key generation
int kg(unsigned char *pk, unsigned char *sk)
{
	unsigned char seed[SEED_LEN];
	
	//generate seed
	random_bytes(seed,SEED_LEN);		
	//key generation with seed 
	kg_seed(pk,sk,seed);	
	
	return 0;
}
// encryption
int pke_enc(const unsigned char *pk, const unsigned char *m, unsigned long long mlen, unsigned char *c, unsigned long long *clen)
{
	unsigned char seed[SEED_LEN];
	
	//generate seed
	random_bytes(seed,SEED_LEN);
	//encrypt with seed 
	pke_enc_seed(pk,m,mlen,c,clen,seed);	

	return 0;
}


// decrypt
int pke_dec(const unsigned char *sk, const unsigned char *c,unsigned long long clen, unsigned char *m, unsigned long long *mlen)
{
	//check parameter
	if(sk==NULL || m==NULL || c==NULL)
	{
		return -1;
	}
	
	unsigned char out[DIM_N];
	unsigned char c2[C2_VEC_NUM];
	int c2_len=(clen-DIM_N)*2;

	//c2 decompress
	poly_decompress(c+DIM_N,c2,c2_len);
	//c1*sk
	poly_mul(c,(char *)sk,out,c2_len);
	//compute c2-c1*sk and recover data from m*q/2+e
	
	#ifdef LAC_LIGHT
	//D2+ parity check decode
	uint16_t tmp_c,tmp_r;
	unsigned char code[C2_VEC_NUM];
	
	//compute mlen
	*mlen=c2_len/16-3;
	if(*mlen!=16)//only for 128 bit message
	{
		return 0;
	}
	memset(m,0,*mlen);
	memset(code,0,C2_VEC_NUM);
	//D2 decode
	d2_decode(c2,out,code,c2_len);
	//parity error correction
	memcpy(&tmp_c,code+16,2);
	memcpy(&tmp_r,code+18,1);
	parity_ecc((uint16_t *)code,tmp_c,tmp_r);
	memcpy(m,code,*mlen);
	
	#else
	//D2+BCH decode
		
	unsigned char code[CODE_LEN],*p_code;
	unsigned char m_buf[MESSAGE_LEN];

	//compute mlen
	*mlen=c2_len/16-ECC_LEN;
	//shif the pointer of ecc data
	p_code=code+(DATA_LEN-(*mlen));
	//init code
	memset(code,0,CODE_LEN);
	//D2 decode
	d2_decode(c2,out,p_code,c2_len);
	//bch decode to recover m
	ecc_dec(m_buf,code);
	//get plaintext
	memcpy(m,m_buf+(DATA_LEN-(*mlen)),*mlen);
	
	#endif
	
	
	return 0;
}

static int encode_to_e2(char *e2, const unsigned char *m, unsigned long long mlen, int *c2_len)
{
	int i;
	int8_t message;
	int vec_bound;
	
	#ifdef LAC_LIGHT
	if(mlen!=16)//only for 128 bit message
	{
		return 0;
	}
	//D2+ two dimention parity
	uint16_t parity_c[1];
	uint8_t  parity_r[1];
	unsigned char code[C2_VEC_NUM];
	//package  message + parity
	memcpy(code,m,mlen);
	parity((uint16_t*)code,parity_c,parity_r);
	memcpy(code+mlen,(unsigned char*)(parity_c),2);
	memcpy(code+mlen+2,parity_r,1);
	//compute the length of c2
	*c2_len=C2_VEC_NUM;
	vec_bound=*c2_len/2;
	//compute  code*q/2+e2, 
	for(i=0;i<vec_bound;i++)
	{
		//RATIO=q/2. add code*q/2 to e2
		message=RATIO*((code[i/8]>>(i%8))&1);
		e2[i]=e2[i]+message;
		//D2 encode, repeat at i+vec_bound
		e2[i+vec_bound]=e2[i+vec_bound]+message;
	}
	#else
	
	unsigned char m_buf[MESSAGE_LEN];
	unsigned char code[CODE_LEN],*p_code;
	//BCH+D2 encoding
	//package m_buf
	memset(m_buf,0,MESSAGE_LEN);
	//set data
	memcpy(m_buf+(DATA_LEN-mlen),m,mlen);
	//encode m with ecc code
	ecc_enc(m_buf,code);
	//set p_code
	p_code=code+(DATA_LEN-mlen);
	//compute the length of c2
	*c2_len=(mlen+ECC_LEN)*8*2;
	vec_bound=*c2_len/2;
	//compute  code*q/2+e2, 
	for(i=0;i<vec_bound;i++)
	{
		//RATIO=q/2. add code*q/2 to e2
		message=RATIO*((p_code[i/8]>>(i%8))&1);
		e2[i]=e2[i]+message;
		//D2 encode, repeat at i+vec_bound
		e2[i+vec_bound]=e2[i+vec_bound]+message;
	}
	#endif
	
	return 0;
}

// encryption with seed
int pke_enc_seed(const unsigned char *pk, const unsigned char *m, unsigned long long mlen, unsigned char *c, unsigned long long *clen, unsigned char *seed)
{
	unsigned char seeds[3*SEED_LEN];
	char r[DIM_N];
	char e1[DIM_N],e2[DIM_N];
	unsigned char c2[C2_VEC_NUM];
	unsigned char a[DIM_N];
	
	int c2_len;
	
	//check parameter
	if(pk==NULL || m==NULL || c==NULL )
	{
		return -1;
	}
	if(mlen>MESSAGE_LEN)
	{
		return -1;
	}
	
	//generate  a from seed in the first part of pk
	gen_a(a,pk);
	//generate three seeds for r,e1,e2
	pseudo_random_bytes(seeds,3*SEED_LEN,seed);
	//generate random vector r
	gen_r(r,seeds);
	//generate error vector e1
	gen_e(e1,seeds+SEED_LEN);
	//generate c1: c1=a*r+e1
	poly_aff(a,r,e1,c,DIM_N);
	//generate error vector e2
	gen_e(e2,seeds+2*SEED_LEN);
	//encode message to e2
	encode_to_e2(e2,m,mlen, &c2_len);
	//c2=b*r+e2+m*[q/2]
	poly_aff(pk+SEED_LEN,r,e2,c2,c2_len);
	//compress c2
	poly_compress(c2,c+DIM_N,c2_len);
	*clen=DIM_N+c2_len/2;

	return 0;

}

