! RTTOV interface for the case, that RTTOV is not available
!
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

MODULE mo_rtifc_nort

!------------------------------------------------------------------------------
!
! Description:
!   This module contains version specific stuff for the mo_rtifc module for
!   the case that RTTOV is not available. All subroutines just set
!   the exit status to "ERR_NO_RTTOV_LIB".
!
!------------------------------------------------------------------------------

!---------------
! MACRO SETTINGS
!---------------

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
!
!+ Preprocessor macros for mo_rtifc_*.f90
!
!
! Description:
!   This file contains the Macro definitions
!
! Current Code Owner: DWD, Robin Faulwetter
!    phone: +49 69 8062 2746
!    fax:   +49 69 8062 3721
!    email:  robin.faulwetter@dwd.de
!
! History:
! Version      Date       Name
! ------------ ---------- ----
! @VERSION@    @DATE@     Robin Faulwetter
!  Initial release
!
! Code Description:
! Language: Fortran 2003.
! Software Standards:
!
! Heritage: mo_rtifc_10.f90
!=======================================================================


!---------------
! MACRO SETTINGS
!---------------

! Take into account _DACE_ and __DACE__

! Icon without RTTOV

! Set macros for RTTOV coefficient distribution

! Select mpi routines to be used for coefficient distribution





!-------------------------
!--------------------------------

!========================
!========================

!-------------
! Modules used
!-------------
  use mo_rtifc_base



  implicit none

  private


  ! subroutines
  public :: rtifc_version          ! Version string
  public :: rtifc_init             ! Initialise RTTOV modules, read coeffs
  public :: rtifc_coef_index       ! Returns index of coeffs for given satid/instr
  public :: rtifc_cleanup          ! frees memory allocated by rtifc_init
  public :: rtifc_fill_input       ! fills the profile-dependent part for RTTOV
  public :: rtifc_direct           ! calls RTTOV direct routine
  public :: rtifc_k                ! calls RTTOV K routine
  public :: rtifc_l2c_god          ! god-corrected l2c
  public :: rtifc_print_profiles
  public :: rtifc_coef_prop        ! get properties of coefs

  ! Fake RTTOV options
  public :: rttov_options

  ! parameters/variables
  public :: rtifc_vers
  public :: gas_id_ozone
  public :: gas_id_co2

  interface rtifc_fill_input
    module procedure rtifc_fill_input_var
  end interface

  type rttov_options
    ! Fake RTTOV options
  end type rttov_options

  ! ifc version
  integer, parameter :: rtifc_vers = 0

  integer, parameter :: jprb = wp

  ! Fake gas_id
  integer, parameter :: gas_id_ozone = -1
  integer, parameter :: gas_id_co2   = -1


contains


  function rtifc_version() result(vers)
    character(len=17) :: vers

    vers = rttov_version()
    write(vers(13:),'("IFC",I2.2)') rtifc_vers
  end function rtifc_version


  subroutine rtifc_init(instruments,channels,nchans_inst,ch_id,iopts,my_proc_id,n_proc,  &
                        io_proc_id,mpi_comm_type,status,path_coefs)
    integer,             intent(in)          :: instruments(:,:)
    integer,             intent(in)          :: channels(:,:)
    integer,             intent(in)          :: nchans_inst(:)
    integer,             intent(out)         :: ch_id(:,:)        ! channel indices (for RTTOV calls)
    integer,             intent(in)          :: iopts(:)          ! rttov_options for each instrument
    integer,             intent(in)          :: my_proc_id
    integer,             intent(in)          :: n_proc
    integer,             intent(in)          :: io_proc_id
    integer,             intent(in)          :: mpi_comm_type
    integer,             intent(out)         :: status                ! exit status
    character(*),        intent(in), optional:: path_coefs

    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_init

  subroutine rtifc_l2c_god(iopt, chans, valid, l2c, l2c_corr, transm, opdep, p_l2c, ideb)
    integer,             intent(in)           :: iopt
    integer,             intent(in)           :: chans(:)
    logical,             intent(in)           :: valid(:,:)
    real,                intent(in)           :: l2c(:,:)
    real,                intent(out)          :: l2c_corr(:,:)
    real(kind=jprb),     intent(in), optional :: transm(:,:,:)
    real(kind=jprb),     intent(in), optional :: opdep(:,:,:)
    real(kind=jprb),     intent(in), optional :: p_l2c(:,:)
    integer,             intent(in), optional :: ideb(:)

    call finish('rtifc_l2c_god', err_msg(ERR_NO_RTTOV_LIB))

  end subroutine rtifc_l2c_god

  function rtifc_coef_index(satid, platf, instr) result(idx)
    integer             :: idx
    integer, intent(in) :: satid
    integer, intent(in) :: platf
    integer, intent(in) :: instr

    idx = -ERR_NO_RTTOV_LIB

  end function rtifc_coef_index


  subroutine rtifc_cleanup(lprof, lcoef, latlas)
    logical, intent(in), optional :: lprof
    logical, intent(in), optional :: lcoef
    logical, intent(in), optional :: latlas

!   call finish('rtifc_cleanup', err_msg(ERR_NO_RTTOV_LIB))

  end subroutine rtifc_cleanup






  subroutine rtifc_fill_input_var (status,iopts,press,temp,humi,t2m,q2m,psurf,hsurf,&
                                   u10m,v10m,stemp,stype,lat,lon,sat_zen,sun_zen,   &
                                   sat_azi,sun_azi,cloud,cfrac,snw_frc,id,ctp,cfraction,&
                                   ivect,istart,iend,pe,wr_profs, wr_profs_fmt)
    integer,             intent(out)          :: status
    integer,             intent(in), target   :: iopts(:)        ! options index
    real(wp),            intent(in)           :: press     (:,:)
    real(wp),            intent(in)           :: temp      (:,:)
    real(wp),            intent(in)           :: humi      (:,:)
    real(wp),            intent(in)           :: t2m         (:)
    real(wp),            intent(in)           :: q2m         (:)
    real(wp),            intent(in)           :: psurf       (:)
    real(wp),            intent(in)           :: hsurf       (:)
    real(wp),            intent(in)           :: u10m        (:)
    real(wp),            intent(in)           :: v10m        (:)
    real(wp),            intent(in)           :: stemp       (:)
    integer,             intent(in)           :: stype       (:)
    real(wp),            intent(in)           :: lat         (:)
    real(wp),            intent(in)           :: lon         (:)
    real(wp),            intent(in)           :: sat_zen     (:)
    real(wp),            intent(in), optional :: sun_zen     (:)
    real(wp),            intent(in), optional :: ctp         (:)
    real(wp),            intent(in), optional :: cfraction   (:)
    real(wp),            intent(in), optional :: sat_azi     (:)
    real(wp),            intent(in), optional :: sun_azi     (:)
    real(wp),            intent(in), optional :: cloud   (:,:,:)
    real(wp),            intent(in), optional :: cfrac     (:,:)
    real(wp),            intent(in), optional :: snw_frc     (:) ! snow fraction
    integer,             intent(in), optional :: id          (:)
    integer,             intent(in), optional :: ivect
    integer,             intent(in), optional :: istart
    integer,             intent(in), optional :: iend
    integer,             intent(in), optional :: pe
    integer,             intent(in), optional :: wr_profs(:)
    character(len=*),    intent(in), optional :: wr_profs_fmt

    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_fill_input_var



  subroutine rtifc_direct (iopt,lprofs,chans,emissiv,t_b,status,t_b_clear,rad, &
                             radclear,radupclear,raddnclear,refdnclear,radovercast,radtotal,    &
                             transm,transmtotal,opdep,height,istore,errorstatus,reg_lim,rflag,  &
                             dealloc,iprint,rad_out_flg,pe,l_pio)
    integer,             intent(in)          :: iopt
    integer,             intent(in)          :: lprofs         (:)
    integer,             intent(in)          :: chans          (:)
    real(wp),            intent(inout)       :: emissiv      (:,:)
    real(wp),            intent(out)         :: t_b          (:,:)
    integer,             intent(out)         :: status
    real(wp),            intent(out),optional:: t_b_clear    (:,:)
    real(wp),            intent(out),optional:: rad          (:,:)
    real(wp),            intent(out),optional:: radclear     (:,:)
    real(wp),            intent(out),optional:: radupclear   (:,:)
    real(wp),            intent(out),optional:: raddnclear   (:,:)
    real(wp),            intent(out),optional:: refdnclear   (:,:)
    real(wp),            intent(out),optional:: radovercast(:,:,:)
    real(wp),            intent(out),optional:: radtotal     (:,:)
    real(wp),            intent(out),optional:: transm     (:,:,:)
    real(wp),            intent(out),optional:: transmtotal  (:,:)
    real(wp),            intent(out),optional:: opdep      (:,:,:)
    real(wp),            intent(out),optional:: height       (:,:)
    integer,             intent(in) ,optional:: istore       (:,:)
    integer,             intent(out),optional:: errorstatus    (:)
    integer,             intent(out),optional:: reg_lim    (:,:,:)
    integer,             intent(out),optional:: rflag        (:,:)
    integer,             intent(in) ,optional:: iprint         (:)
    logical,             intent(in) ,optional:: dealloc
    integer,             intent(in) ,optional:: rad_out_flg
    integer,             intent(in) ,optional:: pe
    logical,             intent(in) ,optional:: l_pio

    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_direct


  subroutine rtifc_k (iopt,lprofs,chans,emissiv,emissiv_k,temp_k,     &
                      humi_k,t2m_k,q2m_k,stemp_k,t_b,status,t_b_clear,rad,       &
                      radclear,radtotal,radovercast,transm,opdep,psurf_k,u10m_k,v10m_k,&
                      o3_surf_k,wfetc_k,ctp_k,cfraction_k,clw_k,o3_k,     &
                      co2_k,n2o_k,co_k,ch4_k,istore,reg_lim,rflag,dealloc, &
                      iprint,rad_out_flg,pe,l_pio)
    integer ,            intent(in)          :: iopt
    integer ,            intent(in)          :: lprofs         (:)
    integer ,            intent(in)          :: chans          (:)
    real(wp),            intent(inout)       :: emissiv      (:,:)
    real(wp),            intent(inout)       :: emissiv_k    (:,:)
    real(wp),            intent(out)         :: temp_k     (:,:,:)
    real(wp),            intent(out)         :: humi_k     (:,:,:)
    real(wp),            intent(out)         :: t2m_k        (:,:)
    real(wp),            intent(out)         :: q2m_k        (:,:)
    real(wp),            intent(out)         :: stemp_k      (:,:)
    real(wp),            intent(out)         :: t_b          (:,:)
    integer,             intent(out)         :: status
    real(wp),            intent(out),optional:: t_b_clear    (:,:)
    real(wp),            intent(out),optional:: rad          (:,:)
    real(wp),            intent(out),optional:: radclear     (:,:)
    real(wp),            intent(out),optional:: radtotal     (:,:)
    real(wp),            intent(out),optional:: radovercast(:,:,:)
    real(wp),            intent(out),optional:: transm     (:,:,:)
    real(wp),            intent(out),optional:: opdep      (:,:,:)
    real(wp),            intent(out),optional:: psurf_k      (:,:)
    real(wp),            intent(out),optional:: u10m_k       (:,:)
    real(wp),            intent(out),optional:: v10m_k       (:,:)
    real(wp),            intent(out),optional:: o3_surf_k    (:,:)
    real(wp),            intent(out),optional:: wfetc_k      (:,:)
    real(wp),            intent(out),optional:: ctp_k        (:,:)
    real(wp),            intent(out),optional:: cfraction_k  (:,:)
    real(wp),            intent(out),optional:: clw_k      (:,:,:)
    real(wp),            intent(out),optional:: o3_k       (:,:,:)
    real(wp),            intent(out),optional:: co2_k      (:,:,:)
    real(wp),            intent(out),optional:: n2o_k      (:,:,:)
    real(wp),            intent(out),optional:: co_k       (:,:,:)
    real(wp),            intent(out),optional:: ch4_k      (:,:,:)
    integer,             intent(in) ,optional:: istore       (:,:)
    integer,             intent(out),optional:: reg_lim    (:,:,:)
    integer,             intent(out),optional:: rflag        (:,:)
    logical,             intent(in) ,optional:: dealloc
    integer,             intent(in) ,optional:: rad_out_flg
    integer,             intent(in) ,optional:: iprint         (:)
    integer,             intent(in) ,optional:: pe
    logical,             intent(in) ,optional:: l_pio

    status = ERR_NO_RTTOV_LIB

  end subroutine rtifc_k

  subroutine rtifc_coef_prop(iopt, version, sign_opdep, preslev, nlevs, igas, gas_bkg, &
                             satid, grid, instr)
    integer, intent(in)             :: iopt     ! options index
    integer, intent(out),  optional :: version
    integer, intent(out),  optional :: sign_opdep
    real(wp), intent(out), optional :: preslev(:)
    integer,  intent(out), optional :: nlevs
    integer,  intent(in),  optional :: igas
    real(wp), intent(out), optional :: gas_bkg(:)
    integer,  intent(out), optional :: satid
    integer,  intent(out), optional :: grid
    integer,  intent(out), optional :: instr

    call finish('rtifc_l2c_god', err_msg(ERR_NO_RTTOV_LIB))

  end subroutine rtifc_coef_prop


  subroutine rtifc_print_profiles(unit)
    integer, intent(in) :: unit
  end subroutine rtifc_print_profiles


!===============================
!===============================

end module mo_rtifc_nort
