extern avgtime(p1:fieldset,p2:number,p3:number,p4:number,p5:number,p6:number) "fortran90" inline

      Program Main
      use grib_api
!
! reads GRIB data and determines temporal average
! converted to GRIB_API by Iain Russels
!
! Compile: f90 -o avgtime avgtime.f90 -lmars -lemos
!
      integer :: IFIELDSET_IN, IFIELDSET_OUT                 ! Metview fieldset ids
      integer :: igrib_in, igrib_out                         ! grib_api 'handles'
      integer :: iret                                        ! grib_api return code
      integer(kind=kindOfSize) :: byte_size                  ! size of grib message for grib_api

      real*8 dlon, dlat, dpara, dlev, dtime     ! communication with METVIEW

      real, allocatable :: values(:), AVG(:,:,:)
      real, allocatable :: LEVELS(:)
      integer, allocatable :: PARAM(:)

      character*50 filename

! --- infos
!      write(*,*)'I message from avgtime()'

! --- Get first argument for an input fieldset. ICNT is the number of fields.
      call mfi_get_fieldset (IFIELDSET_IN,ICNT)

! --- Get field dimensions as further arguments
      call mfi_get_number (DLON)
      nlon = int(dlon)

      call mfi_get_number (DLAT)
      nlat = int(dlat)

      call mfi_get_number (DPARA)
      npara = int(dpara)

      call mfi_get_number (DLEV)
      nlev = int(dlev)

      call mfi_get_number (DTIME)
      ntime = int(dtime)

! --- field size
      ISIZE = nlon*nlat
      IPNTS = ISIZE
      
! --- allocate some memory
      allocate(values(IPNTS),stat=ierr)
       if (ierr.ne.0) stop 'E to allocate memory'
      allocate(AVG(IPNTS,nlev,npara),stat=ierr)
       if (ierr.ne.0) stop 'E to allocate memory'
      allocate(levels(nlev),param(npara),stat=ierr)
       if (ierr.ne.0) stop 'E to allocate memory'

! --- init
      values = 0.0
      AVG = 0.0

      do i = 1, ntime
        do j = 1, nlev
          do k = 1, npara

            call mfi_load_one_grib (IFIELDSET_IN, igrib_in) ! get a grib_api handle
!            //printf ("x: %f, DBL: %f , ", x, DBL_MAX);


! ---       get the array of values, plus some grib parameteters
            call grib_get_real8_array (igrib_in, 'values', values, iret)
            call grib_get (igrib_in, 'paramId', param(k))
            call grib_get (igrib_in, 'level',                levels(j))


! ---       perform the averaging
            do l = 1, nlon*nlat
              avg(l,j,k) = avg(l,j,k) + values(l)/real(ntime)
            enddo

! ---       free the memory for this field, but not if it's the very last one
            if (i.ne.ntime.or.j.ne.nlev.or.k.ne.npara) then
                call grib_release(igrib_in)
            endif

          enddo ! end of loop over parameters
        enddo   ! end of loop over levels
      enddo     ! end of loop over time


! --- Create a new output fieldset
      call mfi_new_fieldset (IFIELDSET_OUT)

      do j = 1, nlev
        do k = 1, npara

          do l = 1, nlon*nlat
            values(l) = avg(l,j,k)
          enddo

! ---     create a new output grib with this data
          call grib_clone (igrib_in, igrib_out)

! ---     set the values array and other parameters
          call grib_set(igrib_out, 'paramId', param(k))
          call grib_set(igrib_out, 'level', levels(j))
          call grib_set_real8_array (igrib_out, 'values', values, iret)


! ---     pass the grib message back to Metview
          call mfi_save_grib (IFIELDSET_OUT,igrib_out)
          call grib_release(igrib_out)

        enddo ! end of loop over parameters
      enddo   ! end of loop over levels

      CALL mfi_return_fieldset (IFIELDSET_OUT)

      end

end inline
