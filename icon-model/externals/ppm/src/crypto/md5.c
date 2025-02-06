/**
 * @file md5.c
 * @brief Public-Domain MD5 implementation massaged for POSIX-
 * and ScalES-PPM-conformance
 *
 * @copyright Copyright  (C)  2012  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords: MD5 Message-Digest Algorithm
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

/*	$OpenBSD: md5.c,v 1.2 2011/01/11 15:42:05 deraadt Exp $	*/

/*
 * This code implements the MD5 message-digest algorithm.
 * The algorithm is due to Ron Rivest.	This code was
 * written by Colin Plumb in 1993, no copyright is claimed.
 * This code is in the public domain; do with it what you wish.
 *
 * Equivalent code is available from RSA Data Security, Inc.
 * This code has been tested against that, and is equivalent,
 * except that you don't need to include two pages of legalese
 * with every copy.
 *
 * To compute the message digest of a chunk of bytes, declare an
 * MD5Context structure, pass it to MD5Init, call MD5Update as
 * needed on buffers full of bytes, and then call MD5Final, which
 * will fill a supplied 16-byte array with the digest.
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdlib.h>
#include <string.h>

#include "core/ppm_visibility.h"
#include "crypto/md5.h"

static void
PPM_MD5Transform(uint32_t state[4],
                 const uint8_t block[PPM_MD5_BLOCK_LENGTH]);

#define PUT_64BIT_LE(cp, value) do {					\
        (cp)[7] = (value) >> 56;					\
        (cp)[6] = (value) >> 48;					\
        (cp)[5] = (value) >> 40;					\
        (cp)[4] = (value) >> 32;					\
        (cp)[3] = (value) >> 24;					\
        (cp)[2] = (value) >> 16;					\
        (cp)[1] = (value) >> 8;						\
        (cp)[0] = (value); } while (0)

#define PUT_32BIT_LE(cp, value) do {					\
        (cp)[3] = (value) >> 24;					\
        (cp)[2] = (value) >> 16;					\
        (cp)[1] = (value) >> 8;						\
        (cp)[0] = (value); } while (0)

static uint8_t PADDING[PPM_MD5_BLOCK_LENGTH] = {
        0x80, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
};

/*
 * Start MD5 accumulation.  Set bit count to 0 and buffer to mysterious
 * initialization constants.
 */
void
PPM_MD5Init(PPM_MD5_CTX *ctx)
{
        ctx->count = 0;
        ctx->state[0] = UINT32_C(0x67452301);
        ctx->state[1] = UINT32_C(0xefcdab89);
        ctx->state[2] = UINT32_C(0x98badcfe);
        ctx->state[3] = UINT32_C(0x10325476);
}

/*
 * Update context to reflect the concatenation of another buffer full
 * of bytes.
 */
void
PPM_MD5Update(PPM_MD5_CTX *ctx, const uint8_t *input, size_t len)
{
  size_t have, need;

  /* Check how many bytes we already have and how many more we need. */
  have = (size_t)((ctx->count >> 3) & (PPM_MD5_BLOCK_LENGTH - 1));
  need = PPM_MD5_BLOCK_LENGTH - have;

  /* Update bitcount */
  ctx->count += (uint64_t)len << 3;

  if (len >= need) {
    if (have != 0) {
      memcpy(ctx->buffer + have, input, need);
      PPM_MD5Transform(ctx->state, ctx->buffer);
      input += need;
      len -= need;
      have = 0;
    }

    /* Process data in PPM_MD5_BLOCK_LENGTH-byte chunks. */
    while (len >= PPM_MD5_BLOCK_LENGTH) {
      PPM_MD5Transform(ctx->state, input);
      input += PPM_MD5_BLOCK_LENGTH;
      len -= PPM_MD5_BLOCK_LENGTH;
    }
  }

  /* Handle any remaining bytes of data. */
  if (len != 0)
    memcpy(ctx->buffer + have, input, len);
}

/*
 * Final wrapup - pad to 64-byte boundary with the bit pattern
 * 1 0* (64-bit count of bits processed, MSB-first)
 */
void
PPM_MD5Final(unsigned char digest[PPM_MD5_DIGEST_LENGTH], PPM_MD5_CTX *ctx)
{
        uint8_t count[8];
        size_t padlen;
        int i;

        /* Convert count to 8 bytes in little endian order. */
        PUT_64BIT_LE(count, ctx->count);

        /* Pad out to 56 mod 64. */
        padlen = PPM_MD5_BLOCK_LENGTH -
            ((ctx->count >> 3) & (PPM_MD5_BLOCK_LENGTH - 1));
        if (padlen < 1 + 8)
                padlen += PPM_MD5_BLOCK_LENGTH;
        PPM_MD5Update(ctx, PADDING, padlen - 8);		/* padlen - 8 <= 64 */
        PPM_MD5Update(ctx, count, 8);

        if (digest != NULL) {
                for (i = 0; i < 4; i++)
                        PUT_32BIT_LE(digest + i * 4, ctx->state[i]);
        }
        memset(ctx, 0, sizeof(*ctx));	/* in case it's sensitive */
}


/* The four core functions - F1 is optimized somewhat */

/* #define F1(x, y, z) (x & y | ~x & z) */
#define F1(x, y, z) (z ^ (x & (y ^ z)))
#define F2(x, y, z) F1(z, x, y)
#define F3(x, y, z) (x ^ y ^ z)
#define F4(x, y, z) (y ^ (x | ~z))

/* This is the central step in the MD5 algorithm. */
#define MD5STEP(f, w, x, y, z, data, s) \
        ( w += f(x, y, z) + data,  w = w<<s | w>>(32-s),  w += x )

/*
 * The core of the MD5 algorithm, this alters an existing MD5 hash to
 * reflect the addition of 16 longwords of new data.  PPM_MD5Update blocks
 * the data and converts bytes into longwords for this routine.
 */
static void
PPM_MD5Transform(uint32_t state[4], const uint8_t block[PPM_MD5_BLOCK_LENGTH])
{
        uint32_t a, b, c, d, in[PPM_MD5_BLOCK_LENGTH / 4];

#ifdef WORDS_BIGENDIAN
        for (a = 0; a < PPM_MD5_BLOCK_LENGTH / 4; a++) {
                in[a] = (uint32_t)(
                    (uint32_t)(block[a * 4 + 0]) |
                    (uint32_t)(block[a * 4 + 1]) <<  8 |
                    (uint32_t)(block[a * 4 + 2]) << 16 |
                    (uint32_t)(block[a * 4 + 3]) << 24);
        }
#else
        memcpy(in, block, sizeof(in));
#endif

        a = state[0];
        b = state[1];
        c = state[2];
        d = state[3];

        MD5STEP(F1, a, b, c, d, in[ 0] + UINT32_C(0xd76aa478),  7);
        MD5STEP(F1, d, a, b, c, in[ 1] + UINT32_C(0xe8c7b756), 12);
        MD5STEP(F1, c, d, a, b, in[ 2] + UINT32_C(0x242070db), 17);
        MD5STEP(F1, b, c, d, a, in[ 3] + UINT32_C(0xc1bdceee), 22);
        MD5STEP(F1, a, b, c, d, in[ 4] + UINT32_C(0xf57c0faf),  7);
        MD5STEP(F1, d, a, b, c, in[ 5] + UINT32_C(0x4787c62a), 12);
        MD5STEP(F1, c, d, a, b, in[ 6] + UINT32_C(0xa8304613), 17);
        MD5STEP(F1, b, c, d, a, in[ 7] + UINT32_C(0xfd469501), 22);
        MD5STEP(F1, a, b, c, d, in[ 8] + UINT32_C(0x698098d8),  7);
        MD5STEP(F1, d, a, b, c, in[ 9] + UINT32_C(0x8b44f7af), 12);
        MD5STEP(F1, c, d, a, b, in[10] + UINT32_C(0xffff5bb1), 17);
        MD5STEP(F1, b, c, d, a, in[11] + UINT32_C(0x895cd7be), 22);
        MD5STEP(F1, a, b, c, d, in[12] + UINT32_C(0x6b901122),  7);
        MD5STEP(F1, d, a, b, c, in[13] + UINT32_C(0xfd987193), 12);
        MD5STEP(F1, c, d, a, b, in[14] + UINT32_C(0xa679438e), 17);
        MD5STEP(F1, b, c, d, a, in[15] + UINT32_C(0x49b40821), 22);

        MD5STEP(F2, a, b, c, d, in[ 1] + UINT32_C(0xf61e2562),  5);
        MD5STEP(F2, d, a, b, c, in[ 6] + UINT32_C(0xc040b340),  9);
        MD5STEP(F2, c, d, a, b, in[11] + UINT32_C(0x265e5a51), 14);
        MD5STEP(F2, b, c, d, a, in[ 0] + UINT32_C(0xe9b6c7aa), 20);
        MD5STEP(F2, a, b, c, d, in[ 5] + UINT32_C(0xd62f105d),  5);
        MD5STEP(F2, d, a, b, c, in[10] + UINT32_C(0x02441453),  9);
        MD5STEP(F2, c, d, a, b, in[15] + UINT32_C(0xd8a1e681), 14);
        MD5STEP(F2, b, c, d, a, in[ 4] + UINT32_C(0xe7d3fbc8), 20);
        MD5STEP(F2, a, b, c, d, in[ 9] + UINT32_C(0x21e1cde6),  5);
        MD5STEP(F2, d, a, b, c, in[14] + UINT32_C(0xc33707d6),  9);
        MD5STEP(F2, c, d, a, b, in[ 3] + UINT32_C(0xf4d50d87), 14);
        MD5STEP(F2, b, c, d, a, in[ 8] + UINT32_C(0x455a14ed), 20);
        MD5STEP(F2, a, b, c, d, in[13] + UINT32_C(0xa9e3e905),  5);
        MD5STEP(F2, d, a, b, c, in[ 2] + UINT32_C(0xfcefa3f8),  9);
        MD5STEP(F2, c, d, a, b, in[ 7] + UINT32_C(0x676f02d9), 14);
        MD5STEP(F2, b, c, d, a, in[12] + UINT32_C(0x8d2a4c8a), 20);

        MD5STEP(F3, a, b, c, d, in[ 5] + UINT32_C(0xfffa3942),  4);
        MD5STEP(F3, d, a, b, c, in[ 8] + UINT32_C(0x8771f681), 11);
        MD5STEP(F3, c, d, a, b, in[11] + UINT32_C(0x6d9d6122), 16);
        MD5STEP(F3, b, c, d, a, in[14] + UINT32_C(0xfde5380c), 23);
        MD5STEP(F3, a, b, c, d, in[ 1] + UINT32_C(0xa4beea44),  4);
        MD5STEP(F3, d, a, b, c, in[ 4] + UINT32_C(0x4bdecfa9), 11);
        MD5STEP(F3, c, d, a, b, in[ 7] + UINT32_C(0xf6bb4b60), 16);
        MD5STEP(F3, b, c, d, a, in[10] + UINT32_C(0xbebfbc70), 23);
        MD5STEP(F3, a, b, c, d, in[13] + UINT32_C(0x289b7ec6),  4);
        MD5STEP(F3, d, a, b, c, in[ 0] + UINT32_C(0xeaa127fa), 11);
        MD5STEP(F3, c, d, a, b, in[ 3] + UINT32_C(0xd4ef3085), 16);
        MD5STEP(F3, b, c, d, a, in[ 6] + UINT32_C(0x04881d05), 23);
        MD5STEP(F3, a, b, c, d, in[ 9] + UINT32_C(0xd9d4d039),  4);
        MD5STEP(F3, d, a, b, c, in[12] + UINT32_C(0xe6db99e5), 11);
        MD5STEP(F3, c, d, a, b, in[15] + UINT32_C(0x1fa27cf8), 16);
        MD5STEP(F3, b, c, d, a, in[2 ] + UINT32_C(0xc4ac5665), 23);

        MD5STEP(F4, a, b, c, d, in[ 0] + UINT32_C(0xf4292244),  6);
        MD5STEP(F4, d, a, b, c, in[7 ] + UINT32_C(0x432aff97), 10);
        MD5STEP(F4, c, d, a, b, in[14] + UINT32_C(0xab9423a7), 15);
        MD5STEP(F4, b, c, d, a, in[5 ] + UINT32_C(0xfc93a039), 21);
        MD5STEP(F4, a, b, c, d, in[12] + UINT32_C(0x655b59c3),  6);
        MD5STEP(F4, d, a, b, c, in[3 ] + UINT32_C(0x8f0ccc92), 10);
        MD5STEP(F4, c, d, a, b, in[10] + UINT32_C(0xffeff47d), 15);
        MD5STEP(F4, b, c, d, a, in[1 ] + UINT32_C(0x85845dd1), 21);
        MD5STEP(F4, a, b, c, d, in[8 ] + UINT32_C(0x6fa87e4f),  6);
        MD5STEP(F4, d, a, b, c, in[15] + UINT32_C(0xfe2ce6e0), 10);
        MD5STEP(F4, c, d, a, b, in[6 ] + UINT32_C(0xa3014314), 15);
        MD5STEP(F4, b, c, d, a, in[13] + UINT32_C(0x4e0811a1), 21);
        MD5STEP(F4, a, b, c, d, in[4 ] + UINT32_C(0xf7537e82),  6);
        MD5STEP(F4, d, a, b, c, in[11] + UINT32_C(0xbd3af235), 10);
        MD5STEP(F4, c, d, a, b, in[2 ] + UINT32_C(0x2ad7d2bb), 15);
        MD5STEP(F4, b, c, d, a, in[9 ] + UINT32_C(0xeb86d391), 21);

        state[0] += a;
        state[1] += b;
        state[2] += c;
        state[3] += d;
}

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
