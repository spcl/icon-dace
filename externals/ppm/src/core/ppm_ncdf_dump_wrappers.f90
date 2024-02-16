!
! Copyright  (C)  2015  Thomas Jahns <jahns@dkrz.de>
!
! Version: 1.0
! Keywords: NetCDF
! Author: Thomas Jahns <jahns@dkrz.de>
! Maintainer: Thomas Jahns <jahns@dkrz.de>
! URL: https://www.dkrz.de/redmine/projects/scales-ppm
!
! Redistribution and use in source and binary forms, with or without
! modification, are  permitted provided that the following conditions are
! met:
!
! Redistributions of source code must retain the above copyright notice,
! this list of conditions and the following disclaimer.
!
! Redistributions in binary form must reproduce the above copyright
! notice, this list of conditions and the following disclaimer in the
! documentation and/or other materials provided with the distribution.
!
! Neither the name of the DKRZ GmbH nor the names of its contributors
! may be used to endorse or promote products derived from this software
! without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!
!> @file ppm_ncdf_dump_wrappers.f90
!! helper functions for ppm_ncdf_dump.f90, needed to work around
!! different handling of scalar and array arguments
#define filename 'ppm_ncdf_dump_wrappers.f90'
FUNCTION ppm_ncdf_dump_wrap_pais(ncid, varid, name, xtype, slen, s) RESULT(r)
  USE netcdf, ONLY: nf90_put_att, nf90_char
  USE ppm_base, ONLY: assertion
  INTEGER, INTENT(in) :: ncid, varid, xtype, slen
  CHARACTER(len=*), INTENT(in) :: name
  INTEGER, INTENT(in) :: s(*)
  INTEGER :: r
  CALL assertion(xtype == nf90_char, filename, __LINE__, &
       'incorrect storage type')
  r = nf90_put_att(ncid, varid, name, s(1:slen))
END FUNCTION  ppm_ncdf_dump_wrap_pais

FUNCTION ppm_ncdf_dump_wrap_pads(ncid, varid, name, xtype, slen, s) RESULT(r)
  USE netcdf, ONLY: nf90_put_att, nf90_double
  USE ppm_base, ONLY: assertion
  INTEGER, INTENT(in) :: ncid, varid, xtype, slen
  CHARACTER(len=*), INTENT(in) :: name
  DOUBLE PRECISION, INTENT(in) :: s(*)
  INTEGER :: r
  CALL assertion(xtype == nf90_double, filename, __LINE__, &
       'incorrect storage type')
  r = nf90_put_att(ncid, varid, name, s(1:slen))
END FUNCTION  ppm_ncdf_dump_wrap_pads

FUNCTION ppm_ncdf_dump_wrap_pars(ncid, varid, name, xtype, slen, s) RESULT(r)
  USE netcdf, ONLY: nf90_put_att, nf90_real
  USE ppm_base, ONLY: assertion
  INTEGER, INTENT(in) :: ncid, varid, xtype, slen
  CHARACTER(len=*), INTENT(in) :: name
  REAL, INTENT(in) :: s(*)
  INTEGER :: r
  CALL assertion(xtype == nf90_real, filename, __LINE__, &
       'incorrect storage type')
  r = nf90_put_att(ncid, varid, name, s(1:slen))
END FUNCTION  ppm_ncdf_dump_wrap_pars
!
! Local Variables:
! license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
! license-markup: "doxygen"
! license-default: "bsd"
! End:
