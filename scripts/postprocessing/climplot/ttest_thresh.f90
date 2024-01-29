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
