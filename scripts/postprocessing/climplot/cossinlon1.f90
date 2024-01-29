extern cossinlon1(f1:fieldset,f2:number) "fortran90" inline
      Program Main
!---------------------------------------------------------------------------------
!  function cossinlon(Field,T)
!
!           input: Field contains field, just for GRIB info
!                  T is time in hours (UTC, 0 to 24)
!           output: 2 fields containing cos(2*pi*((T-12)/24+LON/360))
!                                   and sin(2*pi*((T-12)/24+LON/360))
!
!          pgf90 -o cossinlon1  cossinlon1.f90 -lmars -lemos   # on linux
!
!--------------------------------------------------------------------------------
! reads GRIB field 
!
!

      use grib_api      
      parameter (NLONM=1024,NLATM=512)
      integer NLON(NLATM)      
      integer NLON_, length      
      real    ZINC_I, ZLON_FIRST      
      integer, allocatable :: IFIELD(:)      
      real    ZDLON(NLATM)      
      integer :: igrib_in, igrib_out                         ! grib_api 'handles'      
      integer :: iret                                        ! grib_api return code      
      integer (kind=kindOfSize) :: byte_size                 ! size of grib message for grib_api      
      character(len=1), dimension(:), allocatable :: message ! grib message container      
      integer :: numberOfPlValues, numberOfValues            ! size of 'pl' and 'vaues' arrays      
            
      real   zt     ! communication with METVIEW      
            
      real, allocatable :: valuesA(:), valuesB(:), pl(:), AVG(:)      
            
      character*50 filename      
            
! --- Get first argument for an input fieldset. ICNT is the number of fields.      
      call MGETG (IGRIB1,ICNT)      
            
! --- Get time ZT as further argument      
      call MGETN2 (ZT)      
!     write(*,*) "ZT: ", ZT      
! --- Create a new output fieldset      
      call MNEWG (IGRIB2)      
            
! --- field size      
      ISIZE = NLONM*NLATM      
      IPNTS = ISIZE      
               
! --- allocate some memory      
      allocate(IFIELD(ISIZE + 10000),stat=ierr)     ! +10k for grib header      
       if (ierr.ne.0) stop '#E to allocate memory'      
!      
            
! --- init      
      IFIELD = 0      
      ZSEC4 = 0.0      
      ZSEC4A = -9999.      
! --- load field      
      call MLOADG (IGRIB1,IFIELD,ISIZE)      
            
            
! ---       transfer the integer grib message into a character array for grib_api      
      length = size(transfer(IFIELD, message))      
      allocate(message(length))      
      message = transfer(IFIELD, message)      
            
      call grib_new_from_message(igrib_in, message)  ! get a grib_api handle      
            
      call grib_get (igrib_in, 'Nj', NLAT)      
            
      IF (NLAT .GT. NLATM) THEN      
         WRITE(*,*) '#E ERROR in SINCOSLON; NLAT > NALTM'      
         WRITE(*,*) ' NLAT, NLATM: ', NLAT,NLATM      
         STOP      
      ENDIF      
            
            
! --- check whether this is a reduced grid or not - if so, then 'numberOfPointsAlongAParallel'      
! --- will be missing.      
            
            
      is_missing=0;      
      call grib_is_missing(igrib_in,'Ni', is_missing);      
            
      IF ( is_missing /= 1 ) THEN  ! regular grid      
        call grib_get (igrib_in, 'Ni', NLON_)      
        call grib_get (igrib_in, 'iDirectionIncrementInDegrees', ZINC_I)      
            
        DO JL=1,NLAT      
          NLON(JL)=NLON_      
          ZDLON(JL)=ZINC_I      
        ENDDO      
            
      ELSE  ! irregular grid      
        call grib_get_size(igrib_in, 'pl', numberOfPlValues) ! get the size of the 'pl' array
        allocate(pl(numberOfPlValues), stat=iret)            ! allocate memory for it      
        call grib_get_real8_array(igrib_in, 'pl', pl, iret)  ! get the 'pl' values themselves
            
        DO JL=1,NLAT      
          NLON(JL)=pl(JL)      
          ZDLON(JL)=360./NLON(JL)      
        ENDDO      
      ENDIF           
            
            
            
! --- allocate some memory for the output 'values' arrays      
            
      call grib_get_size(igrib_in, 'values', numberOfValues) ! get the size of the 'values' array  
      allocate(valuesA(numberOfValues),stat=ierr)
      if (ierr.ne.0) stop '#E to allocate memory'
      allocate(valuesB(numberOfValues),stat=ierr)      
      if (ierr.ne.0) stop '#E to allocate memory'
            
            
            
! -- calculate sin and cos of longitude function      
      call grib_get (igrib_in, 'longitudeOfFirstGridPointInDegrees', ZLON_FIRST)
      ZPI=4.*atan(1.)      
      ZLONW=ZLON_FIRST      
!     write(*,*) "ZLONW: ", ZLONW, "ZDLON:", ZDLON(1)      
      inum = 0      
      do ii = 1, nlat       
        do jj = 1, nlon(ii)      
          inum = inum + 1      
          valuesA(inum)=cos( 2*zpi*( (ZT-12.)/24. +ZLONW/360.&
           &               +float(jj-1)*zdlon(ii)/360. ) )  
          valuesB(inum)=sin( 2*zpi*( (ZT-12.)/24. +ZLONW/360.&
           &               +float(jj-1)*zdlon(ii)/360. ) )      
        enddo      
      enddo      
            
            
! ---     create a new output grib with the 'A' data      
      call grib_get_message_size(igrib_in, byte_size)      
      allocate(message(byte_size), stat=ierr)      
            
      call grib_copy_message(igrib_in, message)      
      call grib_new_from_message(igrib_out, message)      
            
! ---     set the values array and other parameters      
      call grib_set_real8_array (igrib_out, 'values', valuesA, iret)      
            
            
      call grib_get_message_size(igrib_out, byte_size)      
      deallocate(message)      
      allocate(message(byte_size), stat=ierr)      
      call grib_copy_message(igrib_out, message)      
            
! --- add to output filedset      
      CALL MSAVEG (IGRIB2,message,byte_size/4)      
            
            
! ---     create a new output grib with the 'A' data      
      call grib_get_message_size(igrib_in, byte_size)      
      deallocate(message)      
      allocate(message(byte_size), stat=ierr)      
      call grib_copy_message(igrib_in, message)      
      call grib_new_from_message(igrib_out, message)      
            
! ---     set the values array and other parameters      
      call grib_set_real8_array (igrib_out, 'values', valuesB, iret)      
            
      call grib_get_message_size(igrib_out, byte_size)      
      deallocate(message)      
      allocate(message(byte_size), stat=ierr)      
      call grib_copy_message(igrib_out, message)      
            
! --- add to output filedset      
      CALL MSAVEG (IGRIB2,message,byte_size/4)      
            
            
! --- send sesult      
      CALL MSETG(IGRIB2)      
            
!      
      end      
end inline      
