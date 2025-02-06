! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

extern ttest_thresh(p1:number,p2:number) "fortran90" inline
!
! Adrian Tompkins 12/6/2009
! 
! simple program to calculate T-test limit 
! Adapted from Thomas Jung
!
      Program Main
      
      real(8) dsig, ddof,tcrit     ! communication with METVIEW
      
!* - NAG stuff
      external G01FBF
      real(8) G01FBF
      
      call MGETN2(dsig)
      call MGETN2(ddof)
      tcrit=G01FBF('L',dsig,ddof,ifail)
      call MSETN2(tcrit)
      
      end
      
end inline
