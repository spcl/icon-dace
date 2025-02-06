/* Cross Platform Floating Point Arithmetics

   This header file defines several platform-dependent macros that ensure
   equal and deterministic floating point behaviour across several platforms,
   compilers and architectures.

   The current macros are currently only used on x86 and x86_64 architectures,
   on every other architecture, these macros expand to NOPs. This assumes that
   other architectures do not have an internal precision and the operhand types
   define the computational precision of floating point operations. This
   assumption may be false, in that case, the author is interested in further
   details on the other platform.

   For further details, please visit:
   http://www.christian-seiler.de/projekte/fpmath/

   Author: Christian Seiler <webmaster@christian-seiler.de>
   Version: 20081026

   Modifications for ScalES-PPM: Thomas Jahns <jahns@dkrz.de>, 20110913

   This file is released under public domain - or - in countries where this is
   not possible under the following license:

      Permission is hereby granted, free of charge, to any person obtaining a
      copy of this software, to deal in the software without restriction,
      including without limitation the rights to use, copy, modify, merge,
      publish, distribute, sublicense, and/or sell copies of the software,
      and to permit persons to whom the software is furnished to do so, subject
     to no condition whatsoever.

      This software is provided AS IS, without warranty of any kind, express or
      implied. */

#ifndef XPFPA_FUNC_H
#define XPFPA_FUNC_H

/*
 Implementation notes:

 x86_64:
  - Since all x86_64 compilers use SSE by default, it is probably unecessary
    to use these macros there. We define them anyway since we are too lazy
    to differentiate the architecture. Also, the compiler option -mfpmath=i387
    justifies this decision.

 General:
  - It would be nice if one could detect whether SSE if used for math via some
    funky compiler defines and if so, make the macros go to NOPs. Any ideas
    on how to do that?

 MS Visual C:
  - Since MSVC users tipically don't use autoconf or CMake, we will detect
    MSVC via compile time define.
*/

// MSVC detection (MSVC people usually don't use autoconf)
#ifdef _MSC_VER
# if _MSC_VER >= 1500
   // Visual C++ 2008 or higher, supports _controlfp_s
#  define HAVE__CONTROLFP_S
# else
   // Visual C++ (up to 2005), supports _controlfp
#  define HAVE__CONTROLFP
# endif // MSC_VER >= 1500
  // Tell MSVC optimizer that we access FP environment
# pragma fenv_access (on)
#endif // _MSC_VER

#ifdef HAVE__CONTROLFP_S

// float.h defines _controlfp_s
# include <float.h>

static inline void
xpfpa_save(uint32_t *_xpfpa_fpu_oldcw)
{
  _controlfp_s(_xpfpa_fpu_cw, 0, 0);
}


static inline void
xpfpa_switch_double(uint32_t *_xpfpa_fpu_oldcw)
{
  uint32_t _xpfpa_fpu_cw;
  _controlfp_s(&_xpfpa_fpu_cw, 0, 0);
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _controlfp_s(&_xpfpa_fpu_cw, _PC_53, _MCW_PC);
}

static inline void
xpfpa_switch_single(uint32_t *_xpfpa_fpu_oldcw)
{
  uint32_t _xpfpa_fpu_cw;
  _controlfp_s(&_xpfpa_fpu_cw, 0, 0);
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _controlfp_s(&_xpfpa_fpu_cw, _PC_24, _MCW_PC);
}

// NOTE: This only sets internal precision. MSVC does NOT support double-
// extended precision!
static inline void
xpfpa_switch_double_extended(uint32_t *_xpfpa_fpu_oldcw)
{
  _controlfp_s(&_xpfpa_fpu_cw, 0, 0);
  _xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _controlfp_s(&_xpfpa_fpu_cw, _PC_64, _MCW_PC);
}

static inline void
xpfpa_restore(uint32_t _xpfpa_fpu_oldcw)
{
  uint32_t _xpfpa_fpu_cw;
  _controlfp_s(&_xpfpa_fpu_cw, _xpfpa_fpu_oldcw, _MCW_PC);
}

// We do NOT use the volatile return trick since _controlfp_s is a function
// call and thus FP registers are saved in memory anyway. However, we do use
// a variable to ensure that the expression passed into val will be evaluated
// *before* switching back contexts.

#elif defined(HAVE__CONTROLFP)

// float.h defines _controlfp
# include <float.h>

static inline void
xpfpa_save(uint32_t *_xpfpa_fpu_oldcw)
{
  *_xpfpa_fpu_oldcw = _controlfp(0, 0);
}

static inline void
xpfpa_switch_double(uint32_t *_xpfpa_fpu_oldcw)
{
  uint32_t _xpfpa_fpu_cw = _controlfp(0, 0);
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _controlfp(_PC_53, _MCW_PC);
}

static inline void
xpfpa_switch_single(uint32_t *_xpfpa_fpu_oldcw)
{
  uint32_t _xpfpa_fpu_cw = _controlfp(0, 0);
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _controlfp(_PC_24, _MCW_PC);
}

// NOTE: This will only work as expected on MinGW.
static inline void
xpfpa_switch_double_extended(uint32_t *_xpfpa_fpu_oldcw)
{
  uint32_t _xpfpa_fpu_cw = _controlfp(0, 0);
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _controlfp(_PC_64, _MCW_PC);
}

static inline void
xpfpa_restore(uint32_t _xpfpa_fpu_oldcw)
{
  _controlfp(_xpfpa_fpu_oldcw, _MCW_PC);
}

#elif defined(HAVE__FPU_SETCW) // glibc systems

// fpu_control.h defines _FPU_[GS]ETCW
# include <fpu_control.h>

static inline void
xpfpa_save(uint32_t *_xpfpa_fpu_oldcw)
{
  fpu_control_t _xpfpa_fpu_cw;
  _FPU_GETCW(_xpfpa_fpu_cw);
  *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
}

static inline void
xpfpa_switch_double(uint32_t *_xpfpa_fpu_oldcw)
{
  fpu_control_t _xpfpa_fpu_cw;
  _FPU_GETCW(_xpfpa_fpu_cw);
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _xpfpa_fpu_cw =
    (fpu_control_t)((_xpfpa_fpu_cw & ~_FPU_EXTENDED & ~_FPU_SINGLE)
                    | _FPU_DOUBLE);
  _FPU_SETCW(_xpfpa_fpu_cw);
}

static inline void
xpfpa_switch_single(uint32_t *_xpfpa_fpu_oldcw)
{
  fpu_control_t _xpfpa_fpu_cw;
  _FPU_GETCW(_xpfpa_fpu_cw);
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _xpfpa_fpu_cw =
    (fpu_control_t)((_xpfpa_fpu_cw & ~_FPU_EXTENDED & ~_FPU_DOUBLE)
                    | _FPU_SINGLE);
  _FPU_SETCW(_xpfpa_fpu_cw);
}

static inline void
xpfpa_switch_double_extended(uint32_t *_xpfpa_fpu_oldcw)
{
  fpu_control_t _xpfpa_fpu_cw;
  _FPU_GETCW(_xpfpa_fpu_cw);
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _xpfpa_fpu_cw =
    (_xpfpa_fpu_cw & ~_FPU_SINGLE & ~_FPU_DOUBLE) | _FPU_EXTENDED;
  _FPU_SETCW(_xpfpa_fpu_cw);
}

static inline void
xpfpa_restore(uint32_t _xpfpa_fpu_cw)
{
  _FPU_SETCW(_xpfpa_fpu_cw);
}

#elif defined(HAVE_FPSETPREC) // FreeBSD

// fpu_control.h defines _FPU_[GS]ETCW
# include <machine/ieeefp.h>

static inline void
xpfpa_save(uint32_t *_xpfpa_fpu_oldcw)
{
  *_xpfpa_fpu_oldcw = (uint32_t)fpgetprec();
}

static inline void
xpfpa_switch_double(uint32_t *_xpfpa_fpu_oldprec)
{
  if (_xpfpa_fpu_oldprec) *_xpfpa_fpu_oldprec = (uint32_t)fpgetprec();
  fpsetprec(FP_PD);
}

static inline void
xpfpa_switch_single(uint32_t *_xpfpa_fpu_oldprec)
{
  if (_xpfpa_fpu_oldprec) *_xpfpa_fpu_oldprec = (uint32_t)fpgetprec();
  fpsetprec(FP_PS);
}

static inline void
xpfpa_switch_double_extended(uint32_t *_xpfpa_fpu_oldprec)
{
  if (_xpfpa_fpu_oldprec) *_xpfpa_fpu_oldprec = (uint32_t)fpgetprec();
  fpsetprec(FP_PE);
}

static inline void
xpfpa_restore(uint32_t _xpfpa_fpu_oldprec)
{
  fpsetprec((fp_prec_t)_xpfpa_fpu_oldprec);
}

#elif defined(HAVE_FPU_INLINE_ASM_X86)

/*
  Custom x86 inline assembler implementation.

  This implementation does not use predefined wrappers of the OS / compiler
  but rather uses x86/x87 inline assembler directly. Basic instructions:

  fnstcw - Store the FPU control word in a variable
  fldcw  - Load the FPU control word from a variable

  Bits (only bits 8 and 9 are relevant, bits 0 to 7 are for other things):
     0x0yy: Single precision
     0x1yy: Reserved
     0x2yy: Double precision
     0x3yy: Double-extended precision

  We use an unsigned int for the datatype. glibc sources add __mode__ (__HI__)
  attribute to it (HI stands for half-integer according to docs). It is unclear
  what the does exactly and how portable it is.

  The assembly syntax works with GNU CC, Intel CC and Sun CC.
*/

static inline void
xpfpa_save(uint32_t *_xpfpa_fpu_oldcw)
{
  __asm__ __volatile__ ("fnstcw %0" : "=m" (*_xpfpa_fpu_oldcw));
}

static inline void
xpfpa_switch_double(uint32_t *_xpfpa_fpu_oldcw)
{
  uint32_t _xpfpa_fpu_cw;
  __asm__ __volatile__ ("fnstcw %0" : "=m" (*&_xpfpa_fpu_cw));
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _xpfpa_fpu_cw = (_xpfpa_fpu_cw & ~0x100) | 0x200;
  __asm__ __volatile__ ("fldcw %0" : : "m" (*&_xpfpa_fpu_cw));
}

static inline void
xpfpa_switch_single(uint32_t *_xpfpa_fpu_oldcw)
{
  uint32_t _xpfpa_fpu_cw;
  __asm__ __volatile__ ("fnstcw %0" : "=m" (*&_xpfpa_fpu_cw));
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _xpfpa_fpu_cw = (_xpfpa_fpu_cw & ~0x300);
  __asm__ __volatile__ ("fldcw %0" : : "m" (*&_xpfpa_fpu_cw));
}

static inline void
xpfpa_switch_double_extended(uint32_t *_xpfpa_fpu_oldcw)
{
  uint32_t _xpfpa_fpu_cw;
  __asm__ __volatile__ ("fnstcw %0" : "=m" (*&_xpfpa_fpu_cw));
  if (_xpfpa_fpu_oldcw) *_xpfpa_fpu_oldcw = _xpfpa_fpu_cw;
  _xpfpa_fpu_cw = _xpfpa_fpu_cw | 0x300;
  __asm__ __volatile__ ("fldcw %0" : : "m" (*&_xpfpa_fpu_cw));
}

static inline void
xpfpa_restore(uint32_t _xpfpa_fpu_oldcw)
{
  __asm__ __volatile__ ("fldcw %0" : : "m" (*&_xpfpa_fpu_oldcw));
}


#else // FPU CONTROL

static inline void
xpfpa_save(uint32_t *_xpfpa_fpu_oldcw) { (void)_xpfpa_fpu_oldcw; }

static inline void
xpfpa_switch_double(uint32_t *_xpfpa_fpu_oldcw) { (void)_xpfpa_fpu_oldcw; }

static inline void
xpfpa_switch_single(uint32_t *_xpfpa_fpu_oldcw) { (void)_xpfpa_fpu_oldcw; }

static inline void
xpfpa_switch_double_extended(uint32_t *_xpfpa_fpu_oldcw)
{
  (void)_xpfpa_fpu_oldcw;
}

static inline void
xpfpa_restore(uint32_t _xpfpa_fpu_oldcw) { (void)_xpfpa_fpu_oldcw; }

#endif // FPU CONTROL

#endif // XPFPA_H
