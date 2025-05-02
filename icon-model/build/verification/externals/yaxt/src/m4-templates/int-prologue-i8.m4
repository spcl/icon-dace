dnl int-prologue-i8.m4 --- macros for 8 byte integer to combine with
dnl                        f03-redist-gen.m4 to provide specific
dnl                        xt_redist exchange calls
dnl
dnl Copyright  (C)  2022   Jörg Behrens <behrens@dkrz.de>
dnl                        Moritz Hanke <hanke@dkrz.de>
dnl                        Thomas Jahns <jahns@dkrz.de>
dnl
dnl Version: 1.0
dnl Author: Jörg Behrens <behrens@dkrz.de>
dnl         Moritz Hanke <hanke@dkrz.de>
dnl         Thomas Jahns <jahns@dkrz.de>
dnl Maintainer: Jörg Behrens <behrens@dkrz.de>
dnl             Moritz Hanke <hanke@dkrz.de>
dnl             Thomas Jahns <jahns@dkrz.de>
dnl URL: https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/
dnl
dnl Redistribution and use in source and binary forms, with or without
dnl modification, are  permitted provided that the following conditions are
dnl met:
dnl
dnl Redistributions of source code must retain the above copyright notice,
dnl this list of conditions and the following disclaimer.
dnl
dnl Redistributions in binary form must reproduce the above copyright
dnl notice, this list of conditions and the following disclaimer in the
dnl documentation and/or other materials provided with the distribution.
dnl
dnl Neither the name of the DKRZ GmbH nor the names of its contributors
dnl may be used to endorse or promote products derived from this software
dnl without specific prior written permission.
dnl
dnl THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
dnl IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
dnl TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
dnl PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
dnl OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
dnl EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
dnl PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
dnl PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
dnl LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
dnl NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
dnl SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
dnl
  INTEGER, PARAMETER :: pi8 =  14
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(pi8)
  PUBLIC :: i8
