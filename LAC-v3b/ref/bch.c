/*
 * This library provides encoding/decoding of binary
 * Bose-Chaudhuri-Hocquenghem (BCH) codes.
 *
 * Call init_bch to get a pointer to a newly allocated bch_control structure for
 * the given m (Galois field order), t (error correction capability) and
 * (optional) primitive polynomial parameters.
 *
 * Call encode_bch to compute and store ecc parity bytes to a given buffer.
 * Call decode_bch to detect and locate errors in received data.
 * Algorithmic details:
 *
 * Encoding is performed by processing 32 input bits in parallel, using 4
 * remainder lookup tables.
 *
 * The final stage of decoding involves the following internal steps:
 * a. Syndrome computation
 * b. Error locator polynomial computation using Berlekamp-Massey algorithm
 * c. Error locator root finding (by far the most expensive step)
 *
 
 */

# include <endian.h>
# include <errno.h>
# include <stdint.h>
# include <stdio.h>
# include <stdlib.h>
# include <string.h>
# include "bch.h"
# define MODULE_LICENSE(s)
# define MODULE_AUTHOR(s)
# define MODULE_DESCRIPTION(s)
# define EXPORT_SYMBOL_GPL(x)
# define GFP_KERNEL
# define kmalloc(size, flags) malloc(size)
# define kzalloc(size, flags) memset(malloc(size), 0, size)
# define kfree(ptr) free(ptr)
# define DIV_ROUND_UP(a, b) ((a + b - 1) / b)
# define ARRAY_SIZE(arr) (sizeof(arr) / sizeof(*(arr)))
# define cpu_to_be32(x) htobe32(x)
# define fls(x) ({ \
	unsigned int __tmp = x; \
	unsigned int __count = 0; \
	while (__tmp >>= 1) \
		__count++; \
	__count + 1; \
})

#define GF_M(_p)               ((_p)->m)
#define GF_T(_p)               ((_p)->t)
#define GF_N(_p)               ((_p)->n)

#define BCH_ECC_WORDS(_p)      DIV_ROUND_UP(GF_M(_p)*GF_T(_p), 32)
#define BCH_ECC_BYTES(_p)      DIV_ROUND_UP(GF_M(_p)*GF_T(_p), 8)

/*
 * represent a polynomial over GF(2^m)
 */
struct gf_poly {
	unsigned int deg;    /* polynomial degree */
	unsigned int c[0];   /* polynomial terms */
};

/* given its degree, compute a polynomial size in bytes */
#define GF_POLY_SZ(_d) (sizeof(struct gf_poly)+((_d)+1)*sizeof(unsigned int))

/* polynomial of degree 1 */
struct gf_poly_deg1 {
	struct gf_poly poly;
	unsigned int   c[2];
};

/*
 * same as encode_bch(), but process input data one byte at a time
 */
static void encode_bch_unaligned(struct bch_control *bch,
				 const unsigned char *data, unsigned int len,
				 uint32_t *ecc)
{
	int i;
	const uint32_t *p;
	const int l = BCH_ECC_WORDS(bch)-1;

	while (len--) {
		p = bch->mod8_tab + (l+1)*(((ecc[0] >> 24)^(*data++)) & 0xff);

		for (i = 0; i < l; i++)
			ecc[i] = ((ecc[i] << 8)|(ecc[i+1] >> 24))^(*p++);

		ecc[l] = (ecc[l] << 8)^(*p);
	}
}

/*
 * convert ecc bytes to aligned, zero-padded 32-bit ecc words
 */
static void load_ecc8(struct bch_control *bch, uint32_t *dst,
		      const uint8_t *src)
{
	uint8_t pad[4] = {0, 0, 0, 0};
	unsigned int i, nwords = BCH_ECC_WORDS(bch)-1;

	for (i = 0; i < nwords; i++, src += 4)
		dst[i] = (src[0] << 24)|(src[1] << 16)|(src[2] << 8)|src[3];

	memcpy(pad, src, BCH_ECC_BYTES(bch)-4*nwords);
	dst[nwords] = (pad[0] << 24)|(pad[1] << 16)|(pad[2] << 8)|pad[3];
}

/*
 * convert 32-bit ecc words to ecc bytes
 */
static void store_ecc8(struct bch_control *bch, uint8_t *dst,
		       const uint32_t *src)
{
	uint8_t pad[4];
	unsigned int i, nwords = BCH_ECC_WORDS(bch)-1;

	for (i = 0; i < nwords; i++) {
		*dst++ = (src[i] >> 24);
		*dst++ = (src[i] >> 16) & 0xff;
		*dst++ = (src[i] >>  8) & 0xff;
		*dst++ = (src[i] >>  0) & 0xff;
	}
	pad[0] = (src[nwords] >> 24);
	pad[1] = (src[nwords] >> 16) & 0xff;
	pad[2] = (src[nwords] >>  8) & 0xff;
	pad[3] = (src[nwords] >>  0) & 0xff;
	memcpy(dst, pad, BCH_ECC_BYTES(bch)-4*nwords);
}

/**
 * encode_bch - calculate BCH ecc parity of data
 * @bch:   BCH control structure
 * @data:  data to encode
 * @len:   data length in bytes
 * @ecc:   ecc parity data, must be initialized by caller
 *
 * The @ecc parity array is used both as input and output parameter, in order to
 * allow incremental computations. It should be of the size indicated by member
 * @ecc_bytes of @bch, and should be initialized to 0 before the first call.
 *
 * The exact number of computed ecc parity bits is given by member @ecc_bits of
 * @bch; it may be less than m*t for large values of t.
 */
void encode_bch(struct bch_control *bch, const uint8_t *data,
		unsigned int len, uint8_t *ecc)
{
	const unsigned int l = BCH_ECC_WORDS(bch)-1;
	unsigned int i, mlen;
	unsigned long m;
	uint32_t w, r[l+1];
	const uint32_t * const tab0 = bch->mod8_tab;
	const uint32_t * const tab1 = tab0 + 256*(l+1);
	const uint32_t * const tab2 = tab1 + 256*(l+1);
	const uint32_t * const tab3 = tab2 + 256*(l+1);
	const uint32_t *pdata, *p0, *p1, *p2, *p3;

	if (ecc) {
		/* load ecc parity bytes into internal 32-bit buffer */
		load_ecc8(bch, bch->ecc_buf, ecc);
	} else {
		memset(bch->ecc_buf, 0, sizeof(r));
	}

	/* process first unaligned data bytes */
	m = ((unsigned long)data) & 3;
	if (m) {
		mlen = (len < (4-m)) ? len : 4-m;
		encode_bch_unaligned(bch, data, mlen, bch->ecc_buf);
		data += mlen;
		len  -= mlen;
	}

	/* process 32-bit aligned data words */
	pdata = (uint32_t *)data;
	mlen  = len/4;
	data += 4*mlen;
	len  -= 4*mlen;
	memcpy(r, bch->ecc_buf, sizeof(r));

	/*
	 * split each 32-bit word into 4 polynomials of weight 8 as follows:
	 *
	 * 31 ...24  23 ...16  15 ... 8  7 ... 0
	 * xxxxxxxx  yyyyyyyy  zzzzzzzz  tttttttt
	 *                               tttttttt  mod g = r0 (precomputed)
	 *                     zzzzzzzz  00000000  mod g = r1 (precomputed)
	 *           yyyyyyyy  00000000  00000000  mod g = r2 (precomputed)
	 * xxxxxxxx  00000000  00000000  00000000  mod g = r3 (precomputed)
	 * xxxxxxxx  yyyyyyyy  zzzzzzzz  tttttttt  mod g = r0^r1^r2^r3
	 */
	while (mlen--) {
		/* input data is read in big-endian format */
		w = r[0]^cpu_to_be32(*pdata++);
		p0 = tab0 + (l+1)*((w >>  0) & 0xff);
		p1 = tab1 + (l+1)*((w >>  8) & 0xff);
		p2 = tab2 + (l+1)*((w >> 16) & 0xff);
		p3 = tab3 + (l+1)*((w >> 24) & 0xff);

		for (i = 0; i < l; i++)
			r[i] = r[i+1]^p0[i]^p1[i]^p2[i]^p3[i];

		r[l] = p0[l]^p1[l]^p2[l]^p3[l];
	}
	memcpy(bch->ecc_buf, r, sizeof(r));

	/* process last unaligned bytes */
	if (len)
		encode_bch_unaligned(bch, data, len, bch->ecc_buf);

	/* store ecc parity bytes into original parity buffer */
	if (ecc)
		store_ecc8(bch, ecc, bch->ecc_buf);
}

static inline int modulo(struct bch_control *bch, unsigned int v)
{
	const unsigned int n = GF_N(bch);
	while (v >= n) {
		v -= n;
		v = (v & n) + (v >> GF_M(bch));
	}
	return v;
}

/*
 * shorter and faster modulo function, only works when v < 2N.
 */
static inline int mod_s(struct bch_control *bch, unsigned int v)
{
	
	const unsigned int n = GF_N(bch);
	unsigned int tmp;
	tmp=(v<n)?0x0:0xffffffff;
	
	return v-(n&tmp);
}

static inline int deg(unsigned int poly)
{
	/* polynomial degree is the most-significant bit index */
	return fls(poly)-1;
}

/* Galois field basic operations: multiply, divide, inverse, etc. */

static inline unsigned int gf_mul(struct bch_control *bch, unsigned int a,
				  unsigned int b)
{
	unsigned int tmp,mask;
	
	tmp=bch->a_pow_tab[mod_s(bch, bch->a_log_tab[a]+bch->a_log_tab[b])];
	mask=(a && b)? 0xffffffff:0x0;
	
	return  tmp&mask;
}

static inline unsigned int gf_sqr(struct bch_control *bch, unsigned int a)
{
	unsigned int tmp,mask;
	
	tmp=bch->a_pow_tab[mod_s(bch, 2*bch->a_log_tab[a])];
	mask=(a)? 0xffffffff:0x0;
	
	return tmp&mask;
}

static inline unsigned int gf_div(struct bch_control *bch, unsigned int a,
				  unsigned int b)
{
	unsigned int tmp,mask;
	
	tmp=bch->a_pow_tab[mod_s(bch, bch->a_log_tab[a]+GF_N(bch)-bch->a_log_tab[b])];
	mask=(a && b)? 0xffffffff:0x0;
	
	return tmp&mask;
}

static inline unsigned int gf_inv(struct bch_control *bch, unsigned int a)
{
	return bch->a_pow_tab[GF_N(bch)-bch->a_log_tab[a]];
}

static inline unsigned int a_pow(struct bch_control *bch, int i)
{
	return bch->a_pow_tab[modulo(bch, i)];
}

static inline int a_log(struct bch_control *bch, unsigned int x)
{
	return bch->a_log_tab[x];
}

static inline int a_ilog(struct bch_control *bch, unsigned int x)
{
	return mod_s(bch, GF_N(bch)-bch->a_log_tab[x]);
}

/*
 * compute 2t syndromes of ecc polynomial, i.e. ecc(a^j) for j=1..2t
 */
static void compute_syndromes(struct bch_control *bch, uint32_t *ecc,
			      unsigned int *syn)
{
	int i, j, s;
	unsigned int m;
	uint32_t poly,mask_syn,syn_tmp;
	const int t = GF_T(bch);
	unsigned int w,w2,mask_w;

	s = bch->ecc_bits;
	/* make sure extra bits in last ecc word are cleared */
	m = ((unsigned int)s) & 31;
	if (m)
		ecc[s/32] &= ~((1u << (32-m))-1);
	memset(syn, 0, 2*t*sizeof(*syn));

	/* compute v(a^j) for j=1 .. 2t-1 */
	do {
		poly = *ecc++;
		s -= 32;
		
		for(i=31;i>=0 && i+s>=0;i--)
		{
			mask_syn=((poly>>i)&0x1)?0xffffffff:0x0;
			w=i+s;
			w2=w*2;
			for (j = 0; j < 2*t; j += 2)
			{
				syn_tmp=bch->a_pow_tab[w];
				syn[j] ^= (syn_tmp&mask_syn);
				w+=w2;
				mask_w=(w<bch->n)?0x0 : 0xffffffff;
				w-=(bch->n & mask_w);
			}
		}
	} while (s > 0);

	/* v(a^(2j)) = v(a^j)^2 */
	for (j = 0; j < t; j++)
		syn[2*j+1] = gf_sqr(bch, syn[j]);
}

static void gf_poly_copy(struct gf_poly *dst, struct gf_poly *src, int t)
{
	memcpy(dst, src, GF_POLY_SZ(t));
}

static int compute_error_locator_polynomial(struct bch_control *bch,
					    const unsigned int *syn)
{
	const unsigned int t = GF_T(bch);
	const unsigned int n = GF_N(bch);
	unsigned int i, j, tmp, l, pd = 1, d = syn[0];
	struct gf_poly *elp = bch->elp;
	struct gf_poly *pelp = bch->poly_2t[0];
	struct gf_poly *elp_copy = bch->poly_2t[1];
	struct gf_poly *elp_copy2 = bch->poly_2t[2];
	int k, pp = -1;
	uint16_t mask_d,mask_pelp;
	unsigned int mask_tmp,mask_deg,tmp_c,mask_tmp2;

	memset(pelp, 0, GF_POLY_SZ(2*t));
	memset(elp, 0, GF_POLY_SZ(2*t));

	pelp->deg = 0;
	pelp->c[0] = 1;
	elp->deg = 0;
	elp->c[0] = 1;

	/* use simplified binary Berlekamp-Massey algorithm */
	for (i = 0; i < t ; i++) 
	{
		mask_d=(d?0xffff:0x0);
		k = 2*i-pp;
			
		gf_poly_copy(elp_copy, elp,i);
			// e[i+1](X) = e[i](X)+di*dp^-1*X^2(i-p)*e[p](X) 
		tmp = a_log(bch, d)+n-a_log(bch, pd);

		for (j = 0; j <=i ; j++) 
		{
			mask_deg=((j <= pelp->deg)?0xffffffff:0x0);
			mask_pelp=(pelp->c[j]? 0xffff:0x0);
			l = a_log(bch, pelp->c[j]);
			tmp_c=a_pow(bch, tmp+l);
			elp->c[j+k] ^= (tmp_c&mask_d&mask_pelp&mask_deg);
		}
		// compute l[i+1] = max(l[i]->c[l[p]+2*(i-p]) 
		tmp = pelp->deg+(k&mask_d);
		mask_tmp =((tmp > elp->deg)?0xffffffff:0x0);
		mask_tmp2=((tmp > elp->deg)?0x0:0xffffffff);
		//memcpy cause 70 cycles gap
		if (tmp > elp->deg) 
		{
			gf_poly_copy(pelp, elp_copy,i);
			elp->deg = tmp;
		
		}
		else //just for constant time
		{
			gf_poly_copy(elp_copy2, pelp,i);
			tmp =  elp->deg;
		}
		//update according to mask_tmp	
		pd = (d & mask_tmp)^(pd &mask_tmp2);
		pp = ((2*i) & mask_tmp)^(pp &mask_tmp2);
			
		// di+1 = S(2i+3)+elp[i+1].1*S(2i+2)+...+elp[i+1].lS(2i+3-l) 
		unsigned int mask_j,tmp_d;
		if (i < t-1) 
		{
			d = syn[2*i+2];
			
			for (j = 1; j <=i+1 ; j++)
			{
				mask_j=(j <= elp->deg)?0xffffffff:0x00;
				tmp_d=gf_mul(bch, elp->c[j], syn[2*i+2-j]);
				d ^= (tmp_d&mask_j);
			}
		}
		
		
	}

	return (elp->deg > t) ? -1 : (int)elp->deg;
}

/*
 * build monic, log-based representation of a polynomial
 */
static void init_rep(struct bch_control *bch,
			   const struct gf_poly *a, uint16_t *rep, uint16_t *syn_mask, unsigned int pow_start)
{
	
	int i, d = a->deg, l = GF_N(bch)-a_log(bch, a->c[a->deg]);
	uint16_t mask_d;
    int t= bch->t;
	
	for (i = 0; i <= t; i++)
	{
		mask_d=(i<d) ? 0xFFFF : 0x0;
		rep[i] =  modulo(bch, a_log(bch, a->c[i])+l + pow_start*i)&mask_d;
	    syn_mask[i]=(a->c[i] ? 0xFFFF : 0x0)&mask_d;
	}
	
	rep[d]=modulo(bch, pow_start*d);
	syn_mask[d]=0xFFFF;
	
}


/*
 * exhaustive root search (Chien) implementation - not used, included only for
 * reference/comparison tests
 */
static int chien_search(struct bch_control *bch, unsigned int len,
			struct gf_poly *p, unsigned int *roots)
{
	
	unsigned int i, j, syn, syn0, count = 0;
	const unsigned int k = 8*len+bch->ecc_bits,bound=GF_N(bch)-bch->ecc_bits;
	uint16_t syn_mask[bch->t+1],syn_rep[bch->t+1];
	unsigned int mask_rep,syn_tmp;
	
	//use a log-based representation of polynomial 
	init_rep(bch,p,syn_rep,syn_mask,GF_N(bch)-k);
	
	syn0 = gf_div(bch, p->c[0], p->c[p->deg]);
	
	for (i = GF_N(bch)-k+1; i <= bound; i++) 
	{
		// compute elp(a^i) 
		for (j = 1, syn = syn0; j <= bch->t; j++) 
		{
			syn_rep[j]+=j;
			mask_rep=(syn_rep[j]<bch->n)?0x0:0xffffffff;
		    syn_rep[j]-=(bch->n & mask_rep);
			syn_tmp=bch->a_pow_tab[syn_rep[j]];
			syn ^= (syn_tmp&syn_mask[j]);
		}
		roots[count] = GF_N(bch)-i;
		count+=(syn==0);
		
	}
	
	return count;
}

/**
 * decode_bch - decode received codeword and find bit error locations
 * @bch:      BCH control structure
 * @data:     received data, ignored if @calc_ecc is provided
 * @len:      data length in bytes, must always be provided
 * @recv_ecc: received ecc, if NULL then assume it was XORed in @calc_ecc
 *
 * Returns:
 *  The number of errors found, or -EBADMSG if decoding failed, or -EINVAL if
 *  invalid parameters were provided
 */
int decode_bch(struct bch_control *bch, const uint8_t *data, unsigned int len,
	       const uint8_t *recv_ecc, const uint8_t *calc_ecc,
	       const unsigned int *syn, unsigned int *errloc)
{
	const unsigned int ecc_words = BCH_ECC_WORDS(bch);
	unsigned int nbits;
	int i, err ,t=bch->t;

	/* sanity check: make sure data length can be handled */
	if (8*len > (bch->n-bch->ecc_bits))
		return -EINVAL;
	//check data and ecc pointer
    if (!data || !recv_ecc)
		return -EINVAL;
		
	/* compute received data ecc into an internal buffer */
	encode_bch(bch, data, len, NULL);
	/* load received ecc  */
	load_ecc8(bch, bch->ecc_buf2, recv_ecc);
	/* XOR received and calculated ecc */
	for (i = 0; i < (int)ecc_words; i++) 
	{
		bch->ecc_buf[i] ^= bch->ecc_buf2[i];
	}
	
	//compute syndromes
	compute_syndromes(bch, bch->ecc_buf, bch->syn);
    //compute error locator polynomial
	compute_error_locator_polynomial(bch, bch->syn);
	//find roots
	err=chien_search(bch, len, bch->elp, errloc);

	//post process error location
	nbits = (len*8)+bch->ecc_bits;
	
	unsigned char mask_err;
	
	for (i = 0; i < t; i++) 
	{
		mask_err  = (i<err)? 0xff: 0x00;
		errloc[i] = (nbits-1-errloc[i])&mask_err;
		errloc[i] = ((errloc[i] & ~7)|(7-(errloc[i] & 7)))&mask_err;
	}

	return err;
}

/*
 * generate Galois field lookup tables
 */
static int build_gf_tables(struct bch_control *bch, unsigned int poly)
{
	unsigned int i, x = 1;
	const unsigned int k = 1 << deg(poly);

	/* primitive polynomial must be of degree m */
	if (k != (1u << GF_M(bch)))
		return -1;

	for (i = 0; i < GF_N(bch); i++) {
		bch->a_pow_tab[i] = x;
		bch->a_log_tab[x] = i;
		if (i && (x == 1))
			/* polynomial is not primitive (a^i=1 with 0<i<2^m-1) */
			return -1;
		x <<= 1;
		if (x & k)
			x ^= poly;
	}
	bch->a_pow_tab[GF_N(bch)] = 1;
	bch->a_log_tab[0] = 0;

	return 0;
}

/*
 * compute generator polynomial remainder tables for fast encoding
 */
static void build_mod8_tables(struct bch_control *bch, const uint32_t *g)
{
	int i, j, b, d;
	uint32_t data, hi, lo, *tab;
	const int l = BCH_ECC_WORDS(bch);
	const int plen = DIV_ROUND_UP(bch->ecc_bits+1, 32);
	const int ecclen = DIV_ROUND_UP(bch->ecc_bits, 32);

	memset(bch->mod8_tab, 0, 4*256*l*sizeof(*bch->mod8_tab));

	for (i = 0; i < 256; i++) {
		/* p(X)=i is a small polynomial of weight <= 8 */
		for (b = 0; b < 4; b++) {
			/* we want to compute (p(X).X^(8*b+deg(g))) mod g(X) */
			tab = bch->mod8_tab + (b*256+i)*l;
			data = i << (8*b);
			while (data) {
				d = deg(data);
				/* subtract X^d.g(X) from p(X).X^(8*b+deg(g)) */
				data ^= g[0] >> (31-d);
				for (j = 0; j < ecclen; j++) {
					hi = (d < 31) ? g[j] << (d+1) : 0;
					lo = (j+1 < plen) ?
						g[j+1] >> (31-d) : 0;
					tab[j] ^= hi|lo;
				}
			}
		}
	}
}

/*
 * build a base for factoring degree 2 polynomials
 */
static int build_deg2_base(struct bch_control *bch)
{
	const int m = GF_M(bch);
	int i, j, r;
	unsigned int sum, x, y, remaining, ak = 0, xi[m];

	/* find k s.t. Tr(a^k) = 1 and 0 <= k < m */
	for (i = 0; i < m; i++) {
		for (j = 0, sum = 0; j < m; j++)
			sum ^= a_pow(bch, i*(1 << j));

		if (sum) {
			ak = bch->a_pow_tab[i];
			break;
		}
	}
	/* find xi, i=0..m-1 such that xi^2+xi = a^i+Tr(a^i).a^k */
	remaining = m;
	memset(xi, 0, sizeof(xi));

	for (x = 0; (x <= GF_N(bch)) && remaining; x++) {
		y = gf_sqr(bch, x)^x;
		for (i = 0; i < 2; i++) {
			r = a_log(bch, y);
			if (y && (r < m) && !xi[r]) {
				bch->xi_tab[r] = x;
				xi[r] = 1;
				remaining--;
				
				break;
			}
			y ^= ak;
		}
	}
	/* should not happen but check anyway */
	return remaining ? -1 : 0;
}

static void *bch_alloc(size_t size, int *err)
{
	void *ptr;

	ptr = kmalloc(size, GFP_KERNEL);
	if (ptr == NULL)
		*err = 1;
	return ptr;
}

/*
 * compute generator polynomial for given (m,t) parameters.
 */
static uint32_t *compute_generator_polynomial(struct bch_control *bch)
{
	const unsigned int m = GF_M(bch);
	const unsigned int t = GF_T(bch);
	int n, err = 0;
	unsigned int i, j, nbits, r, word, *roots;
	struct gf_poly *g;
	uint32_t *genpoly;

	g = bch_alloc(GF_POLY_SZ(m*t), &err);
	roots = bch_alloc((bch->n+1)*sizeof(*roots), &err);
	genpoly = bch_alloc(DIV_ROUND_UP(m*t+1, 32)*sizeof(*genpoly), &err);

	if (err) {
		kfree(genpoly);
		genpoly = NULL;
		goto finish;
	}

	/* enumerate all roots of g(X) */
	memset(roots , 0, (bch->n+1)*sizeof(*roots));
	for (i = 0; i < t; i++) {
		for (j = 0, r = 2*i+1; j < m; j++) {
			roots[r] = 1;
			r = mod_s(bch, 2*r);
		}
	}
	/* build generator polynomial g(X) */
	g->deg = 0;
	g->c[0] = 1;
	for (i = 0; i < GF_N(bch); i++) {
		if (roots[i]) {
			/* multiply g(X) by (X+root) */
			r = bch->a_pow_tab[i];
			g->c[g->deg+1] = 1;
			for (j = g->deg; j > 0; j--)
				g->c[j] = gf_mul(bch, g->c[j], r)^g->c[j-1];

			g->c[0] = gf_mul(bch, g->c[0], r);
			g->deg++;
		}
	}
	/* store left-justified binary representation of g(X) */
	n = g->deg+1;
	i = 0;

	while (n > 0) {
		nbits = (n > 32) ? 32 : n;
		for (j = 0, word = 0; j < nbits; j++) {
			if (g->c[n-1-j])
				word |= 1u << (31-j);
		}
		genpoly[i++] = word;
		n -= nbits;
	}
	bch->ecc_bits = g->deg;

finish:
	kfree(g);
	kfree(roots);

	return genpoly;
}

/**
 * init_bch - initialize a BCH encoder/decoder
 * @m:          Galois field order, should be in the range 5-15
 * @t:          maximum error correction capability, in bits
 * @prim_poly:  user-provided primitive polynomial (or 0 to use default)
 *
 * Returns:
 *  a newly allocated BCH control structure if successful, NULL otherwise
 *
 * This initialization can take some time, as lookup tables are built for fast
 * encoding/decoding; make sure not to call this function from a time critical
 * path. Usually, init_bch() should be called on module/driver init and
 * free_bch() should be called to release memory on exit.
 *
 * You may provide your own primitive polynomial of degree @m in argument
 * @prim_poly, or let init_bch() use its default polynomial.
 *
 * Once init_bch() has successfully returned a pointer to a newly allocated
 * BCH control structure, ecc length in bytes is given by member @ecc_bytes of
 * the structure.
 */
struct bch_control *init_bch(int m, int t, unsigned int prim_poly)
{
	int err = 0;
	unsigned int i, words;
	uint32_t *genpoly;
	struct bch_control *bch = NULL;

	const int min_m = 5;
	const int max_m = 15;

	/* default primitive polynomials */
	static const unsigned int prim_poly_tab[] = {
		0x25, 0x43, 0x83, 0x11d, 0x211, 0x409, 0x805, 0x1053, 0x201b,
		0x402b, 0x8003,
	};

	if ((m < min_m) || (m > max_m))
		/*
		 * values of m greater than 15 are not currently supported;
		 * supporting m > 15 would require changing table base type
		 * (uint16_t) and a small patch in matrix transposition
		 */
		goto fail;

	/* sanity checks */
	if ((t < 1) || (m*t >= ((1 << m)-1)))
		/* invalid t value */
		goto fail;

	/* select a primitive polynomial for generating GF(2^m) */
	if (prim_poly == 0)
		prim_poly = prim_poly_tab[m-min_m];

	bch = kzalloc(sizeof(*bch), GFP_KERNEL);
	if (bch == NULL)
		goto fail;

	bch->m = m;
	bch->t = t;
	bch->n = (1 << m)-1;
	words  = DIV_ROUND_UP(m*t, 32);
	bch->ecc_bytes = DIV_ROUND_UP(m*t, 8);
	bch->a_pow_tab = bch_alloc((1+bch->n)*sizeof(*bch->a_pow_tab), &err);
	bch->a_log_tab = bch_alloc((1+bch->n)*sizeof(*bch->a_log_tab), &err);
	bch->mod8_tab  = bch_alloc(words*1024*sizeof(*bch->mod8_tab), &err);
	bch->ecc_buf   = bch_alloc(words*sizeof(*bch->ecc_buf), &err);
	bch->ecc_buf2  = bch_alloc(words*sizeof(*bch->ecc_buf2), &err);
	bch->xi_tab    = bch_alloc(m*sizeof(*bch->xi_tab), &err);
	bch->syn       = bch_alloc(2*t*sizeof(*bch->syn), &err);
	bch->cache     = bch_alloc(2*t*sizeof(*bch->cache), &err);
	bch->elp       = bch_alloc((t+1)*sizeof(struct gf_poly_deg1), &err);

	for (i = 0; i < ARRAY_SIZE(bch->poly_2t); i++)
		bch->poly_2t[i] = bch_alloc(GF_POLY_SZ(2*t), &err);

	if (err)
		goto fail;

	err = build_gf_tables(bch, prim_poly);
	if (err)
		goto fail;

	/* use generator polynomial for computing encoding tables */
	genpoly = compute_generator_polynomial(bch);
	if (genpoly == NULL)
		goto fail;

	build_mod8_tables(bch, genpoly);
	kfree(genpoly);

	err = build_deg2_base(bch);
	if (err)
		goto fail;

	return bch;

fail:
	free_bch(bch);
	return NULL;
}


/**
 *  free_bch - free the BCH control structure
 *  @bch:    BCH control structure to release
 */
void free_bch(struct bch_control *bch)
{
	unsigned int i;

	if (bch) {
		kfree(bch->a_pow_tab);
		kfree(bch->a_log_tab);
		kfree(bch->mod8_tab);
		kfree(bch->ecc_buf);
		kfree(bch->ecc_buf2);
		kfree(bch->xi_tab);
		kfree(bch->syn);
		kfree(bch->cache);
		kfree(bch->elp);

		for (i = 0; i < ARRAY_SIZE(bch->poly_2t); i++)
			kfree(bch->poly_2t[i]);

		kfree(bch);
	}
}