dnl f03-redist-gen.m4 --- type-invariant generator code for a module
dnl                       that provides routines to exchange data of
dnl                       a specific type
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
include(forloop2.m4)dnl
define({interface_gen_nd},
{forloop({dim},{1},{7},{{    MODULE PROCEDURE $2_$1_}dim{d}
})})dnl
dnl xt_redist_s_exchange_impl_gen(typenametag,type)
define({xt_redist_s_exchange_impl_gen},{
  ! see @ref xt_redist_s_exchange
  SUBROUTINE xt_redist_s_exchange_}{$1}{_1d_as(redist, src_size, src_data, &
       dst_size, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: src_size, dst_size
    $2, TARGET, INTENT(in) :: src_data(src_size)
    $2, TARGET, INTENT(inout) :: dst_data(dst_size)
    TYPE(c_ptr) :: src_data_cptr, dst_data_cptr
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
    src_data_cptr = C_LOC(src_data)
    dst_data_cptr = C_LOC(dst_data)
ifelse($2,{LOGICAL},
{#else
    CALL xt_slice_c_loc(src_data, src_data_cptr)
    CALL xt_slice_c_loc(dst_data, dst_data_cptr)
#endif
})dnl
    CALL xt_redist_s_exchange1(redist, src_data_cptr, dst_data_cptr)
  END SUBROUTINE xt_redist_s_exchange_}{$1}{_1d_as
forloop({dim},{1},{7},{
  ! see @ref xt_redist_s_exchange
  SUBROUTINE {xt_redist_s_exchange_}$1{_}dim{d(redist, src_data, dst_data)
    TYPE(xt_redist), INTENT(in) :: redist
    $2, TARGET, INTENT(in) :: src_data(:}forloop({pdim},{2},dim,{,:}){)
    $2, TARGET, INTENT(inout) :: dst_data(:}forloop({pdim},{2},dim,{,:}){)

    $2, POINTER :: src_p(:}forloop({pdim},{2},dim,{,:}){), dst_p(:}forloop({pdim},{2},dim,{,:}){)
    $2, TARGET :: dummy(1}forloop({pdim},{2},dim,{,1}){)
    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
    IF (src_size > 0) THEN
      src_p => src_data
    ELSE
      src_p => dummy
      src_size = 1
    END IF
    IF (dst_size > 0) THEN
      dst_p => dst_data
    ELSE
      dst_p => dummy
      dst_size = 1
    END IF
    CALL xt_redist_s_exchange_}{$1}{_1d_as(redist, src_size, src_p, dst_size, dst_p)
  END SUBROUTINE xt_redist_s_exchange_}{$1}{_}dim{d
}})})dnl
dnl
define({xt_redist_a_exchange_impl_gen},{
  ! see @ref xt_redist_a_exchange
  SUBROUTINE xt_redist_a_exchange_}{$1}{_1d_as(redist, src_size, src_data, &
       dst_size, dst_data, request)
    TYPE(xt_redist), INTENT(in) :: redist
    INTEGER, INTENT(in) :: src_size, dst_size
    $2, TARGET, INTENT(in) :: src_data(src_size)
    $2, TARGET, INTENT(inout) :: dst_data(dst_size)
    TYPE(xt_request), INTENT(out) :: request

    $2, TARGET :: dummy(1)
    TYPE(c_ptr) :: src_data_cptr, dst_data_cptr
    IF (src_size > 0) THEN
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
      src_data_cptr = C_LOC(src_data)
ifelse($2,{LOGICAL},
{#else
      CALL xt_slice_c_loc(src_data, src_data_cptr)
#endif
})dnl
    ELSE
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
      src_data_cptr = C_LOC(dummy)
ifelse($2,{LOGICAL},
{#else
      CALL xt_slice_c_loc(dummy, src_data_cptr)
#endif
})dnl
    END IF
    IF (dst_size > 0) THEN
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
      dst_data_cptr = C_LOC(dst_data)
ifelse($2,{LOGICAL},
{#else
      CALL xt_slice_c_loc(dst_data, dst_data_cptr)
#endif
})dnl
    ELSE
ifelse($2,LOGICAL,
{#ifdef HAVE_FC_LOGICAL_INTEROP
})dnl
      dst_data_cptr = C_LOC(dummy)
ifelse($2,{LOGICAL},
{#else
      CALL xt_slice_c_loc(dummy, dst_data_cptr)
#endif
})dnl
    END IF
    CALL xt_redist_a_exchange1(redist, src_data_cptr, dst_data_cptr, request)
  END SUBROUTINE xt_redist_a_exchange_}{$1}{_1d_as
forloop({dim},{1},{7},{
  ! see @ref xt_redist_a_exchange
  SUBROUTINE {xt_redist_a_exchange_}$1{_}dim{d(redist, src_data, dst_data, &
       request)
    TYPE(xt_redist), INTENT(in) :: redist
    $2, TARGET, INTENT(in) :: src_data(:}forloop({pdim},{2},dim,{,:}){)
    $2, TARGET, INTENT(inout) :: dst_data(:}forloop({pdim},{2},dim,{,:}){)
    TYPE(xt_request), INTENT(out) :: request

    INTEGER :: src_size, dst_size
    src_size = SIZE(src_data)
    dst_size = SIZE(dst_data)
#ifdef HAVE_FC_IS_CONTIGUOUS
    IF (.NOT. (IS_CONTIGUOUS(src_data) .AND. IS_CONTIGUOUS(dst_data))) &
      CALL xt_abort('arguments to xt_redist_a_exchange must be contiguous!',&
      filename, __LINE__)
#endif
    CALL xt_redist_a_exchange_}{$1}{_1d_as(redist, src_size, src_data, dst_size, &
         dst_data, request)
  END SUBROUTINE xt_redist_a_exchange_}{$1}{_}dim{d
}})})dnl
dnl
dnl
dnl
dnl Local Variables:
dnl f90-continuation-indent: 5
dnl coding: utf-8
dnl mode: f90
dnl indent-tabs-mode: nil
dnl show-trailing-whitespace: t
dnl require-trailing-newline: t
dnl license-project-url: "https://dkrz-sw.gitlab-pages.dkrz.de/yaxt/"
dnl license-default: "bsd"
dnl End:
dnl
