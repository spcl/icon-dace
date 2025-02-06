// SPDX-FileCopyrightText: David N. Williams
//
// SPDX-License-Identifier: CC0-1.0

/*   Title:  Floating-point exception handling example
    Author:  David N. Williams
      File:  fe-handlng-example.c
   License:  Public Domain
   Version:  0.5.0
   Started:  21-Sep-09
   Revised:  22-Sep-09
   Revised:  30-Sep-09 (comment typo)

This code is an example of alternate, nondefault handling of
IEEE 754 floating-point exceptions in OS X and Linux, based on
the GNU functions feenableexcept(), fedisableeexcept(), and
fegetexcept() [in libm], plus POSIX sigaction().

The GNU functions above are not implemented in OS X Leopard,
gcc 4.x, but are present in Linux.  We implement them here for
OS X, at least until the underlying mechanism is no longer
supported by Apple.

The mechanism is to use the POSIX functions fegetenv() and
fesetenv(), which *are* present in OS X, to manipulate the ppc
and intel floating-point control registers, after changing bits
in fields corresponding to those registers in the fenv_t data
type.

Assembly language code to directly access the floating-point
status and control registers for ppc and intel is also included.

This example grew out of an update to legacy code for Apple
ppc's.  The original legacy code is in Listing 7-1 in "PowerPC
Numerics", 2004:

http://lists.apple.com/archives/unix-porting/2003/May/msg00026.html

Another version of the ppc legacy code is here:

http://developer.apple.com/documentation/Performance/Conceptual/Mac_OSX_Numerics/Mac_OSX_Numerics.pdf

Terry Lambert pointed out that our naive update of the legacy
example to Mac OS X Leopard made egregious unsupported use of
system context structures in the handler.  See his reply to

http://lists.apple.com/archives/Darwin-dev/2009/Sep/msg00091.html

*/

#ifndef __APPLE__
/* BEGIN quote
http://graphviz.sourcearchive.com/documentation/2.16/gvrender__pango_8c-source.html
*/
/* _GNU_SOURCE is needed (supposedly) for the feenableexcept
 * prototype to be defined in fenv.h on GNU systems.
 * Presumably it will do no harm on other systems.
 */
#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif

/* We are not supposed to need __USE_GNU, but I can't see
 * how to get the prototype for fedisableexcept from
 * /usr/include/fenv.h without it.
 */
#ifndef __USE_GNU
#define __USE_GNU
#endif
/* END quote */
#endif // !__APPLE__

#include <fenv.h>

#ifdef __APPLE__

// PPC
#if (defined(__PPC64__) || defined(__ppc64__) || defined(_ARCH_PPC64))

#define FE_EXCEPT_SHIFT 22  // shift flags right to get masks
#define FM_ALL_EXCEPT    FE_ALL_EXCEPT >> FE_EXCEPT_SHIFT

static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = (excepts & FE_ALL_EXCEPT) >> FE_EXCEPT_SHIFT,
               old_excepts;  // all previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = (fenv & FM_ALL_EXCEPT) << FE_EXCEPT_SHIFT;

  fenv = (fenv & ~new_excepts) | new_excepts;
  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

// INTEL
#elif (defined(__x86_64__) || defined(_M_X64) || defined(i386) || defined(__i386__) || defined(__i386) || defined(_M_IX86))

static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  unsigned int new_excepts = excepts & FE_ALL_EXCEPT,
               old_excepts;  // previous masks

  if ( fegetenv (&fenv) ) return -1;
  old_excepts = fenv.__control & FE_ALL_EXCEPT;

  // unmask
  fenv.__control &= ~new_excepts;
  fenv.__mxcsr   &= ~(new_excepts << 7);

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

// ARM
#elif (defined(__arm) || defined(__arm64) || defined(__aarch64__))

#define FE_EXCEPT_SHIFT 8

static int
feenableexcept (unsigned int excepts)
{
  static fenv_t fenv;
  if ( fegetenv (&fenv) ) return -1;
  unsigned int old_excepts = ((fenv.__fpcr) >> FE_EXCEPT_SHIFT) & FE_ALL_EXCEPT;

  unsigned int new_excepts = excepts & FE_ALL_EXCEPT;
  fenv.__fpcr |= new_excepts << FE_EXCEPT_SHIFT;

  return ( fesetenv (&fenv) ? -1 : old_excepts );
}

#else

#error "feenableexcept.h: unknown hardware architecture"

#endif  // PPC, INTEL, or ARM enabling
#endif  // __APPLE__
