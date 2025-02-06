/**
 * @file xt_arithmetic_long.h
 *
 * @copyright Copyright  (C)  2016 Jörg Behrens <behrens@dkrz.de>
 *                                 Moritz Hanke <hanke@dkrz.de>
 *                                 Thomas Jahns <jahns@dkrz.de>
 *
 * @author Jörg Behrens <behrens@dkrz.de>
 *         Moritz Hanke <hanke@dkrz.de>
 *         Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Keywords:
 * Maintainer: Jörg Behrens <behrens@dkrz.de>
 *             Moritz Hanke <hanke@dkrz.de>
 *             Thomas Jahns <jahns@dkrz.de>
 * URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
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
#ifndef XT_ARITHMETIC_LONG_H
#define XT_ARITHMETIC_LONG_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdbool.h>

#include "xt/xt_core.h"
#include "xt_arithmetic_util.h"

enum {
  xt_int_bits = sizeof (Xt_int) * CHAR_BIT,
};
/* since the computation of the next overlap of two stripes might
 * overflow Xt_int (and indeed need 3*xt_int_bits to represent), what
 * follows is an implementation of double- and triple-length integer
 * arithmetic, depending on whether a type twice the size of Xt_int
 * is available.
 *
 * The implementation follows recipes from Henry S. Warren, Jr.;
 * Hacker's Delight, 2-16 and 8-2 (first edition)
 */
#ifdef XT_LONG
typedef XT_LONG Xt_long;
typedef XT_ULONG Xt_ulong;
/**
 * @return -1 if x < 0, 1 otherwise
 */
static inline int
xlsign(Xt_long x)
{
  return (x >= 0) - (x < 0);
}

static inline Xt_long
xlabs(Xt_long x)
{
  return x < 0 ? -x : x;
}

static inline Xt_uint
xlhi(Xt_long x)
{
  return (Xt_uint)(x >> xt_int_bits);
}

static inline Xt_uint
xllo(Xt_long x)
{
  return (Xt_uint)x;
}


/* is a represntable as Xt_int? */
static inline bool
xl_is_in_xt_int_range(Xt_long a)
{
  return (a <= XT_INT_MAX) & (a >= XT_INT_MIN);
}

static inline Xt_long
xiimul(Xt_int a, Xt_int b)
{
  return (Xt_long)a * b;
}

typedef struct { Xt_uint hi; Xt_ulong midlo; } Xt_tword;

static inline Xt_uint
xthi(Xt_tword x)
{
  return x.hi;
}

static inline Xt_uint
xtmid(Xt_tword x)
{
  return (Xt_uint)(x.midlo >> xt_int_bits);
}

static inline Xt_uint
xtlo(Xt_tword x)
{
  return (Xt_uint)x.midlo;
}

static inline Xt_tword
xlimulu(Xt_long a, Xt_uint b)
{
  Xt_tword r = { 0U, 0U };

  Xt_ulong t = (Xt_ulong)(Xt_uint)a*b;
  Xt_uint lo = (Xt_uint)t;
  t = (Xt_ulong)((Xt_ulong)a >> xt_int_bits)*b + (t >> xt_int_bits);
  Xt_uint mid = (Xt_uint)t;
  r.hi = (Xt_uint)(t >> xt_int_bits);

  // Now r has the unsigned product.  Correct by
  // subtracting b*2**(xt_int_bits*2)     if a < 0

  if (a < 0) {
    r.hi = (Xt_uint)(r.hi - (Xt_uint)b);
  }
  r.midlo = ((Xt_ulong)mid << xt_int_bits) | lo;
  return r;
}

static inline Xt_tword
xlimul(Xt_long a, Xt_int b)
{
  Xt_tword r = { 0U, 0U };

  Xt_ulong t = (Xt_ulong)(Xt_uint)a*(Xt_uint)b;
  Xt_uint lo = (Xt_uint)t;
  t = (Xt_ulong)((Xt_ulong)a >> xt_int_bits)*(Xt_uint)b + (t >> xt_int_bits);
  Xt_uint mid = (Xt_uint)t;
  r.hi = (Xt_uint)(t >> xt_int_bits);

  // Now r has the unsigned product.  Correct by
  // subtracting b*2**(xt_int_bits*2)     if a < 0, and
  // subtracting a*2**xt_int_bits if b < 0.

  if (a < 0) {
    r.hi = (Xt_uint)(r.hi - (Xt_uint)b);
  }
  if (b < 0) {
    t = (Xt_ulong)mid - (Xt_uint)a;
    mid = (Xt_uint)t;
    Xt_uint borrow = (Xt_uint)t >> (xt_int_bits-1);
    r.hi = (Xt_uint)(r.hi - (Xt_uint)((Xt_ulong)a >> xt_int_bits) - borrow);
  }
  r.midlo = ((Xt_ulong)mid << xt_int_bits) | lo;
  return r;
}

static inline int
xttcmp_eq(Xt_tword a, Xt_tword b)
{
  return (a.hi == b.hi) & (a.midlo == b.midlo);
}

static inline int
xllcmp_eq(Xt_long a, Xt_long b)
{
  return a == b;
}


#else
/* stores in one-complement form */
typedef struct { Xt_uint hi, lo; } Xt_long;
typedef Xt_long Xt_ulong;

/**
 * @return -1 if x < 0, 1 otherwise
 */
static inline int
xlsign(Xt_long x)
{
  int sign_bit = (int)(x.hi >> (xt_int_bits - 1));
  return (sign_bit^1) - sign_bit;
}

static inline Xt_long
xi2l(Xt_int a)
{
  Xt_long r = { .lo = (Xt_uint)a,
                .hi = (Xt_uint)(Xt_isign_mask(a)) };
  return r;
}

static inline Xt_long
xlabs(Xt_long a)
{
  Xt_uint sign_mask = (Xt_uint)(Xt_isign_mask((Xt_int)a.hi));
  Xt_uint borrow = sign_mask & (Xt_uint)(-(Xt_int)a.lo != 0);
  Xt_long r = { .hi = ((a.hi + sign_mask) ^ sign_mask) - borrow,
                .lo = (a.lo + sign_mask) ^ sign_mask };
  return r;
}

static inline Xt_uint
xlhi(Xt_long x)
{
  return x.hi;
}

static inline Xt_uint
xllo(Xt_long x)
{
  return x.lo;
}

static inline Xt_long
xlnegate(Xt_long a, bool negate)
{
  Xt_uint borrow = (Xt_uint)((Xt_uint)negate & (Xt_uint)(-(Xt_int)a.lo != 0));
  Xt_long r = { .hi = (a.hi ^ (Xt_uint)-(Xt_int)negate) + negate - borrow,
                .lo = (a.lo ^ (Xt_uint)-(Xt_int)negate) + negate };
  return r;
}

/* is a represntable as Xt_int? */
static inline bool
xl_is_in_xt_int_range(Xt_long a)
{
  return Xt_isign_mask((Xt_int)a.lo) == (Xt_int)a.hi;
}

static inline Xt_long
xiiadd(Xt_int a, Xt_int b)
{
  Xt_uint al = (Xt_uint)a, bl = (Xt_uint)b,
    ah = (Xt_uint)(Xt_isign_mask(a)), bh = (Xt_uint)(Xt_isign_mask(b));
  Xt_uint carry = ((al & bl) | ((al | bl) & ~(al + bl))) >> (xt_int_bits-1);
  Xt_long r = { .lo = al + bl, .hi = ah + bh + carry };
  return r;
}

static inline Xt_long
xliadd(Xt_long a, Xt_int b)
{
  Xt_uint al = a.lo, bl = (Xt_uint)b,
    ah = a.hi, bh = (Xt_uint)(Xt_isign_mask(b));
  Xt_uint carry = ((al & bl) | ((al | bl) & ~(al + bl))) >> (xt_int_bits-1);
  Xt_long r = { .lo = al + bl, .hi = ah + bh + carry };
  return r;
}

/* returns a incremented by value of b, i.e. zero or one */
static inline Xt_long
xlinc(Xt_long a, bool b)
{
  Xt_uint al = a.lo, bl = (Xt_uint)b, ah = a.hi;
  Xt_uint carry = bl & !~al;
  Xt_long r = { .lo = al + bl, .hi = ah + carry };
  return r;
}

static inline Xt_long
xlladd(Xt_long a, Xt_long b)
{
  Xt_uint al = a.lo, bl = b.lo, ah = a.hi, bh = b.hi;
  Xt_uint carry = ((al & bl) | ((al | bl) & ~(al + bl))) >> (xt_int_bits-1);
  Xt_long r = { .lo = al + bl, .hi = ah + bh + carry };
  return r;
}

static inline Xt_long
xiisub(Xt_int a, Xt_int b)
{
  Xt_uint al = (Xt_uint)a, bl = (Xt_uint)b,
    ah = (Xt_uint)(Xt_isign_mask(a)), bh = (Xt_uint)(Xt_isign_mask(b));
  Xt_uint carry = ((~al & bl) | ((~(al ^ bl)) & (al - bl))) >> (xt_int_bits-1);
  Xt_long r = { .lo = al - bl, .hi = ah - bh - carry };
  return r;
}

static inline Xt_long
xllsub(Xt_long a, Xt_long b)
{
  Xt_uint al = a.lo, bl = b.lo, ah = a.hi, bh = b.hi;
  Xt_uint carry = ((~al & bl) | ((~(al ^ bl)) & (al - bl))) >> (xt_int_bits-1);
  Xt_long r = { .lo = al - bl, .hi = ah - bh - carry };
  return r;
}

static inline Xt_long
xlisub(Xt_long a, Xt_int b)
{
  Xt_uint al = a.lo, ah = a.hi;
  Xt_uint carry
    = ((~al & (Xt_uint)b) | ((~(al ^ (Xt_uint)b)) & (al - (Xt_uint)b)))
    >> (xt_int_bits-1);
  Xt_long r = { .lo = al - (Xt_uint)b, .hi = ah - carry };
  return r;
}

static inline Xt_long
xilsub(Xt_int a, Xt_long b)
{
  Xt_uint al = (Xt_uint)a, bl = b.lo,
    ah = (Xt_uint)(Xt_isign_mask(a)), bh = b.hi;
  Xt_uint carry = ((~al & bl) | ((~(al ^ bl)) & (al - bl))) >> (xt_int_bits-1);
  Xt_long r = { .lo = al - bl, .hi = ah - bh - carry };
  return r;
}

enum {
  xt_hint_bits = xt_int_bits/2,
};

/* this is implemented following H.S.Warren Hacker's Delight mulhs (8-2) */
static inline Xt_long
xiimul(Xt_int a, Xt_int b)
{
  const Xt_uint lo_mask = ((Xt_uint)1 << xt_hint_bits) - 1U;
  Xt_uint a_lo = (Xt_uint)a & lo_mask,
    b_lo = (Xt_uint)b & lo_mask,
    a_hi = (Xt_uint)(a >> xt_hint_bits),
    b_hi = (Xt_uint)(b >> xt_hint_bits),
    lo_prod = a_lo*b_lo;
  Xt_int t = (Xt_int)(a_hi*b_lo + (lo_prod >> xt_hint_bits));
  Xt_int w1 = (Xt_int)((Xt_uint)t & lo_mask),
    w2 = t >> xt_hint_bits;
  w1 = (Xt_int)(a_lo*b_hi + (Xt_uint)w1);
  Xt_long r = { .hi = a_hi * b_hi + (Xt_uint)w2 + (Xt_uint)(w1 >> xt_hint_bits),
                .lo = (Xt_uint)(a * b) };
  return r;
}

typedef struct { Xt_uint hi, mid, lo; } Xt_tword;

static inline Xt_uint
xthi(Xt_tword x)
{
  return x.hi;
}

static inline Xt_uint
xtmid(Xt_tword x)
{
  return x.mid;
}

static inline Xt_uint
xtlo(Xt_tword x)
{
  return x.lo;
}

static inline Xt_tword
xlimul(Xt_long a, Xt_int b)
{
  const Xt_uint lo_mask = ((Xt_uint)1 << xt_hint_bits) - 1U;
  Xt_tword r;

  Xt_uint bl = (Xt_uint)b & lo_mask, bh = (Xt_uint)b >> xt_hint_bits,
    a0 = a.lo & lo_mask, a1 = a.lo >> xt_hint_bits,
    a2 = a.hi & lo_mask, a3 = a.hi >> xt_hint_bits;
  Xt_uint t = a0*bl;
  Xt_uint accum = t & lo_mask;
  Xt_uint k;
  t = a1*bl + (t >> xt_hint_bits);
  r.lo = accum | (t << xt_hint_bits);
  k = t >> xt_hint_bits;
  t = a2*bl + k;
  accum = t & lo_mask;
  k = t >> xt_hint_bits;
  t = a3*bl + k;
  r.mid = accum | (t << xt_hint_bits);
  r.hi = t >> xt_hint_bits;
  t = a0*bh + (r.lo >> xt_hint_bits);
  r.lo = (r.lo & lo_mask) | (t << xt_hint_bits);
  k = t >> xt_hint_bits;
  t = a1*bh + (r.mid & lo_mask) + k;
  accum = t & lo_mask;
  k = t >> xt_hint_bits;
  t = a2*bh + (r.mid >> xt_hint_bits) + k;
  r.mid = accum | (t << xt_hint_bits);
  k = t >> xt_hint_bits;
  r.hi += a3*bh + k;

  // Now r has the unsigned product.  Correct by
  // subtracting b*2**(xt_int_bits*2)     if a < 0, and
  // subtracting a*2**xt_int_bits if b < 0.

  if ((Xt_int)a.hi < 0) {
    r.hi -= (Xt_uint)b;
  }
  if (b < 0) {
    t = r.mid - a.lo;
    r.mid = (Xt_uint)t;
    Xt_uint borrow = t >> (xt_int_bits-1);
    r.hi -= a.hi + borrow;
  }
  return r;
}

typedef struct { Xt_int quot, rem; } Xt_idiv;
typedef struct { Xt_long quot, rem; } Xt_ldiv;

/* divide long unsigned value a by b, from
 * H.S.Warren, Hacker's Delight,
 * this code is only valid if the result fits an Xt_uint
 */
/* This version is divlu1 except with outer loop unrolled, and array un
changed into local variables.  Several of the variables below could be
"short," but having them fullwords gives better code on gcc/Intel.
Statistics:  Based on one million executions of this program, with
uniformly distributed random values for the dividend and divisor, the
number of times in each loop per execution of the program is:

again1: 0.4060
again2: 0.3469

This is the version that's in the hacker book. */
static inline Xt_idiv
xlidivu(Xt_ulong a, Xt_uint b)
{
  const Xt_uint lo_mask = ((Xt_uint)1 << xt_hint_bits) - 1U;
  const Xt_uint base = (Xt_uint)1 << xt_hint_bits; // Number base
  Xt_uint un1, un0,          // Norm. dividend LSD's.
    ah = a.hi, al = a.lo,
    vn1, vn0,                // Norm. divisor digits.
    q1, q0,                  // Quotient digits.
    un32, un21, un10,        // Dividend digit pairs.
    rhat;                    // A remainder.

  if (ah >= b)       // overflow, set remainder to impossible value
    return (Xt_idiv){ .quot = (Xt_int)~(Xt_uint)0, .rem = (Xt_int)~(Xt_uint)0 };

  // Shift amount for norm.
  int s = xinlz(b);            // 0 <= s <= xt_int_bits.
  enum { int_bits = sizeof (int) * CHAR_BIT };
  b = b << s;               // Normalize divisor.
  vn1 = b >> xt_hint_bits;  // Break divisor up into
  vn0 = b & lo_mask;        // two half-words

  un32 = (ah << s)
    | ((al >> (xt_int_bits - s)) & (Xt_uint)(-s >> (int_bits-1)));
  un10 = al << s;           // Shift dividend left.

  un1 = un10 >> xt_hint_bits; // Break right half of
  un0 = un10 & lo_mask;      // dividend into two digits.

  q1 = un32/vn1;            // Compute the first
  rhat = un32 - q1*vn1;     // quotient digit, q1.
  while (q1 >= base || q1*vn0 > base*rhat + un1) {
    q1 -= 1U;
    rhat += vn1;
    if (rhat >= base) break;
  }

  un21 = un32*base + un1 - q1*b;  // Multiply and subtract.

  q0 = un21/vn1;            // Compute the second
  rhat = un21 - q0*vn1;     // quotient digit, q0.
  while (q0 >= base || q0*vn0 > base*rhat + un0) {
    q0 -= 1U;
    rhat += vn1;
    if (rhat >= base) break;
  }

  return (Xt_idiv){ .quot = (Xt_int)(q1*base + q0),
      .rem = (Xt_int)((un21*base + un0 - q0*b) >> s)};
}

typedef XT_SHORT Xt_short;
typedef unsigned XT_SHORT Xt_ushort;


/* q[0], r[0], u[0], and v[0] contain the LEAST significant halfwords.
(The sequence is in little-endian order).

This first version is a fairly precise implementation of Knuth's
Algorithm D, for a binary computer with base b = 2**16.  The caller
supplies
   1. Space q for the quotient, m - n + 1 halfwords (at least one).
   2. Space r for the remainder (optional), n halfwords.
   3. The dividend u, m halfwords, m >= 1.
   4. The divisor v, n halfwords, n >= 2.
The most significant digit of the divisor, v[n-1], must be nonzero.  The
dividend u may have leading zeros; this just makes the algorithm take
longer and makes the quotient contain more leading zeros.  A value of
NULL may be given for the address of the remainder to signify that the
caller does not want the remainder.
   The program does not alter the input parameters u and v.
   The quotient and remainder returned may have leading zeros.  The
function itself returns a value of 0 for success and 1 for invalid
parameters (e.g., division by 0).
   For now, we must have m >= n.  Knuth's Algorithm D also requires
that the dividend be at least as long as the divisor.  (In his terms,
m >= 0 (unstated).  Therefore m+n >= n.) */
static inline void
xdivmnu(Xt_ushort *restrict q, Xt_ushort *restrict r,
        const Xt_ushort *restrict u, const Xt_ushort *restrict v,
        int m, int n) {

  const Xt_uint lo_mask = ((Xt_uint)1 << xt_hint_bits) - 1U;
  const Xt_uint base = (Xt_uint)1 << xt_hint_bits; // Number base
  Xt_int s, t;

  if (m < n || n <= 0 || v[n-1] == 0)
    return;              // Return if invalid param.

  if (n == 1) {                                      // Take care of
    Xt_uint k = 0;                                    // the case of a
    for (int j = m - 1; j >= 0; j--) {               // single-digit
      q[j] = (Xt_ushort)((k*base + u[j])/v[0]);      // divisor here.
      k = (k*base + (Xt_int)u[j]) - (Xt_int)(q[j]*v[0]);
    }
    if (r != NULL) r[0] = (Xt_ushort)k;
    return;
  }

   // Normalize by shifting v left just enough so that
   // its high-order bit is on, and shift u left the
   // same amount.  We may have to append a high-order
   // digit on the dividend; we do that unconditionally.

   s = xinlz(v[n-1]) - xt_hint_bits;        // 0 <= s <= 15.
   Xt_ushort vn[n];                       /* normalized v */
   for (int i = n - 1; i > 0; i--)
     vn[i] = (v[i] << s) | (v[i-1] >> (xt_hint_bits-s));
   vn[0] = v[0] << s;

   Xt_ushort un[m + 1];                   /* normalized u */
   un[m] = u[m-1] >> (xt_hint_bits-s);
   for (int i = m - 1; i > 0; i--)
     un[i] = (u[i] << s) | (u[i-1] >> (xt_hint_bits-s));
   un[0] = u[0] << s;

   for (int j = m - n; j >= 0; j--) {       // Main loop.
    // Compute estimated quotient digit qhat of q[j].
    Xt_uint qhat = (un[j+n]*base + un[j+n-1])/vn[n-1],
      // and remainder
      rhat = (un[j+n]*base + un[j+n-1]) - qhat*vn[n-1];
    while (qhat >= base || qhat*vn[n-2] > base*rhat + un[j+n-2]) {
        --qhat;
        rhat += vn[n-1];
        if (rhat >= base) break;
      }

      // Multiply and subtract.
      Xt_int k = 0;
      for (int i = 0; i < n; i++) {
        Xt_uint p = qhat*vn[i];  // Product of two digits.
        t = un[i+j] - k - (Xt_int)(p & lo_mask);
        un[i+j] = (Xt_ushort)t;
        k = (Xt_int)(p >> xt_hint_bits) - (t >> xt_hint_bits);
      }
      t = un[j+n] - k;
      un[j+n] = (Xt_ushort)t;

      q[j] = (Xt_ushort)qhat;   // Store quotient digit.
      if (t < 0) {              // If we subtracted too
         q[j] = q[j] - 1;       // much, add back.
         k = 0;
         for (int i = 0; i < n; i++) {
            t = un[i+j] + vn[i] + k;
            un[i+j] = (Xt_ushort)t;
            k = t >> xt_hint_bits;
         }
         un[j+n] = (Xt_ushort)(un[j+n] + k);
      }
   } // End j.
  // If the caller wants the remainder, unnormalize
  // it and pass it back.
  if (r != NULL) {
    for (int i = 0; i < n; i++)
      r[i] = (un[i] >> s) | (un[i+1] << (xt_hint_bits-s));
  }
}

static inline void
xl2a(Xt_ushort a[4], Xt_long u)
{
  a[0] = (Xt_ushort)u.lo;
  a[1] = (Xt_ushort)(u.lo >> xt_hint_bits);
  a[2] = (Xt_ushort)u.hi;
  a[3] = (Xt_ushort)(u.hi >> xt_hint_bits);
}

static inline Xt_long
xa2l(Xt_ushort a[4])
{
  Xt_long u;
  u.lo = (Xt_uint)a[0] | (Xt_uint)a[1] << xt_hint_bits;
  u.hi = (Xt_uint)a[2] | (Xt_uint)a[3] << xt_hint_bits;
  return u;
}

static inline Xt_ldiv
xlldivu(Xt_ulong u, Xt_ulong v)
{
  const Xt_uint lo_mask = ((Xt_uint)1 << xt_hint_bits) - 1U;
  const Xt_uint base = (Xt_uint)1 << xt_hint_bits; // Number base

  Xt_uint qhat;            // Estimated quotient digit.
  Xt_uint rhat;            // A remainder.
  Xt_ldiv d;

  if (v.hi == 0U) {             // Take care of the case of a single-word
    d.quot.hi = u.hi/v.lo;      // divisor here.
    Xt_uint k = u.hi - d.quot.hi*v.lo;
    Xt_idiv dl = xlidivu((Xt_ulong){ .hi = k, .lo = u.lo }, v.lo);
    d.quot.lo = (Xt_uint)dl.quot;    // divisor here.
    d.rem = xllsub((Xt_long){ .hi = k, .lo = u.lo },
                   xiimul((Xt_int)d.quot.lo, (Xt_int)v.lo));
  } else if (u.hi < v.hi || (u.hi == v.hi && u.lo < v.lo)) {
    d.quot = (Xt_long){ 0U, 0U };
    d.rem = u;
  } else {
    // Normalize by shifting v left just enough so that
    // its high-order bit is on, and shift u left the
    // same amount.  We may have to append a high-order
    // digit on the dividend; we do that unconditionally.

    int s = xinlz(v.hi);       // 0 <= s <= xt_int_bits.
    Xt_ulong                   // Normalized form of u, v.
      vn = { .hi = (v.hi << s) | (v.lo >> (xt_int_bits-s)),
             .lo = v.lo << s };
    Xt_tword
      un = { .hi = u.hi >> (xt_int_bits-s),
             .mid = (u.hi << s) | (u.lo >> (xt_int_bits-s)),
             .lo = u.lo << s };
    // the normalized divisor has 4, the dividend 6 half-words
    if ((u.hi >= base) == (v.hi >= base)) {
      // Compute estimate qhat of q[j].
      qhat = un.hi/(vn.hi >> xt_hint_bits);
      rhat = un.hi - qhat * (vn.hi >> xt_hint_bits);
      while (qhat >= base || qhat*(vn.hi & lo_mask) > base*rhat + (un.mid >> xt_hint_bits)) {
        --qhat;
        rhat += (vn.hi >> xt_hint_bits);
        if (rhat >= base) break;
      }

      // Multiply and subtract.
      Xt_uint p = qhat * (vn.lo & lo_mask);
      Xt_int t = (Xt_int)((un.lo & lo_mask) - (p & lo_mask));
      Xt_uint accum = (Xt_uint)t & lo_mask;
      Xt_int k = (Xt_int)(p >> xt_hint_bits) - (Xt_int)(t >> xt_hint_bits);
      p = qhat * (vn.lo >> xt_hint_bits);
      t = (Xt_int)(un.lo >> xt_hint_bits) - k - (Xt_int)(p & lo_mask);
      un.lo = accum | ((Xt_uint)t << xt_hint_bits);
      k = (Xt_int)(p >> xt_hint_bits) - (t >> xt_hint_bits);
      p = qhat * (vn.hi & lo_mask);
      t = (Xt_int)(un.mid & lo_mask) - k - (Xt_int)(p & lo_mask);
      accum = (Xt_uint)t & lo_mask;
      k = (Xt_int)(p >> xt_hint_bits) - (t >> xt_hint_bits);
      p = qhat * (vn.hi >> xt_hint_bits);
      t = (Xt_int)(un.mid >> xt_hint_bits) - k - (Xt_int)(p & lo_mask);
      un.mid = accum | ((Xt_uint)t << xt_hint_bits);
      k = (Xt_int)(p >> xt_hint_bits) - (t >> xt_hint_bits);
      t = (Xt_int)((un.hi & lo_mask) - (Xt_uint)k);
      un.hi = (un.hi & ~lo_mask) + (Xt_uint)t;

      d.quot.hi = 0U;
      d.quot.lo = qhat;         // Store quotient digit.
      if (t < 0) {              // If we subtracted too
        --d.quot.lo;            // much, add back.
        t = (Xt_int)(un.lo & lo_mask) + (Xt_int)(vn.lo & lo_mask);
        accum = (Xt_uint)t & lo_mask;
        k = t >> xt_hint_bits;
        t = (Xt_int)((un.lo >> xt_hint_bits) + (vn.lo >> xt_hint_bits)
                     + (Xt_uint)k);
        un.lo = accum | ((Xt_uint)t << xt_hint_bits);
        k = t >> xt_hint_bits;
        t = (Xt_int)((un.mid & lo_mask) + (vn.hi & lo_mask) + (Xt_uint)k);
        accum = (Xt_uint)t & lo_mask;
        k = t >> xt_hint_bits;
        t = (Xt_int)((un.mid >> xt_hint_bits) + (un.hi >> xt_hint_bits)
                     + (Xt_uint)k);
        un.mid = accum | ((Xt_uint)t << xt_hint_bits);
        k = t >> xt_hint_bits;
        un.hi += (Xt_uint)k;
      }
      // unnormalize the remainder and pass it back.
      d.rem = (Xt_long){ .lo = (un.lo >> s) | (un.mid << (xt_int_bits - s)),
                         .hi = (un.mid >> s) | (un.hi << (xt_int_bits - s)) };
    } else {
      Xt_ushort ua[4], va[4], qa[4] = { 0U }, ra[4] = { 0U };
      xl2a(ua, u);
      xl2a(va, v);
      int m = 4, n = 4;
      while (ua[m-1] == 0U) --m;
      while (va[n-1] == 0U) --n;
      xdivmnu(qa, ra, ua, va, m, n);
      d.quot = xa2l(qa);
      d.rem = xa2l(ra);
    }
  }
  return d;
}


/* compute a divided by b */
static inline Xt_ldiv
xtidiv(Xt_tword a, Xt_int b)
{
  /* todo: implement */
  (void)a;
  (void)b;
  Xt_ldiv r = { { 0, 0 }, { 0, 0 } };
  return r;
}

static inline Xt_ldiv
xlldiv(Xt_long a, Xt_long b)
{
  Xt_long au = xlabs(a), bu = xlabs(b);
  Xt_ldiv r = xlldivu(au, bu);
  bool negate = (Xt_int)(a.hi ^ b.hi) < 0;
  r.quot = xlnegate(r.quot, negate);
  r.rem = xlnegate(r.rem, negate);
  return r;
}

/**
 * @return 1 if a > b or 0 if a <= b
 */
static inline int
xlicmp_gt(Xt_long a, Xt_int b)
{
  Xt_uint al = a.lo, bl = (Xt_uint)b;
  Xt_int ah = (Xt_int)a.hi, bh = Xt_isign_mask(b);
  int r = ah > bh || (ah == bh && al > bl);
  return r;
}

/**
 * @return 1 if a >= b or 0 if a < b
 */
static inline int
xlicmp_ge(Xt_long a, Xt_int b)
{
  Xt_uint al = a.lo, bl = (Xt_uint)b;
  Xt_int ah = (Xt_int)a.hi, bh = Xt_isign_mask(b);
  int r = ah > bh || (ah == bh && al >= bl);
  return r;
}

/**
 * @return 1 if a <= b or 0 if a > b
 */
static inline int
xlicmp_le(Xt_long a, Xt_int b)
{
  Xt_uint al = a.lo, bl = (Xt_uint)b;
  Xt_int ah = (Xt_int)a.hi, bh = Xt_isign_mask(b);
  int r = ah < bh || (ah == bh && al <= bl);
  return r;
}

/**
 * @return 1 if a < b or 0 if a >= b
 */
static inline int
xlicmp_lt(Xt_long a, Xt_int b)
{
  Xt_uint al = a.lo, bl = (Xt_uint)b;
  Xt_int ah = (Xt_int)a.hi, bh = Xt_isign_mask(b);
  int r = ah < bh || (ah == bh && al < bl);
  return r;
}

static inline int
xllcmp_eq(Xt_long a, Xt_long b)
{
  return ((a.hi == b.hi) & (a.lo == b.lo));
}

static inline int
xttcmp_eq(Xt_tword a, Xt_tword b)
{
  return (a.hi == b.hi) & (a.mid == b.mid) & (a.lo == b.lo);
}


#endif

#endif
/*
 * Local Variables:
 * coding: utf-8
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
 * license-default: "bsd"
 * End:
 */
