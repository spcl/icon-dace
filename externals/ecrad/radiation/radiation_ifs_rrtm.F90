! # 1 "radiation/radiation_ifs_rrtm.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_ifs_rrtm.f90"
! this file has been modified for the use in icon

! radiation_ifs_rrtm.f90 - interface to ifs implementation of rrtm-g
!
! (c) copyright 2015- ecmwf.
!
! this software is licensed under the terms of the apache licence version 2.0
! which can be obtained at http://www.apache.org/licenses/license-2.0.
!
! in applying this licence, ecmwf does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.
!
! author:  robin hogan
! email:   r.j.hogan@ecmwf.int
!
! modifications
!   2017-04-11  r. hogan  receive "surface" dummy argument
!   2017-09-08  r. hogan  reverted some changes
!   2017-10-18  r. hogan  added planck_function public function
!   2018-01-11  r. hogan  added optional spectral scaling of incoming solar radiation
!   2018-02-22  r. hogan  optimized reverse indexing of heights
!   2018-05-05  r. hogan  gas_optics can be called for reduced number of levels
!   2019-01-02  r. hogan  initialize shortwave props to zero in case sun below horizon

module radiation_ifs_rrtm

  implicit none

  public  :: setup_gas_optics, gas_optics, planck_function, set_gas_units

contains

  !---------------------------------------------------------------------
  ! setup the ifs implementation of rrtm-g gas absorption model
  subroutine setup_gas_optics(config, directory)

    use yoerrtm,   only : jpglw
    use yoesrtm,   only : jpgsw
    use yoerrtftr, only : ngb_lw => ngb
    use yoesrtm,   only : ngb_sw => ngbsw
    use ecradhook,   only : lhook, dr_hook, jphook

    use radiation_config
    use radiation_spectral_definition, only &
         &  : solarreferencetemperature, terrestrialreferencetemperature

    type(config_type), intent(inout), target :: config
    character(len=*), intent(in)     :: directory

    integer :: irep ! for implied do

    integer, parameter :: rrtm_gpoint_reordering_lw(140) = (/ &
          &   89, 90, 139, 77, 137, 69, 131, 97, 91, 70, 78, 71, 53, 72, 123, 54, 79, 98,  &
          &   92, 55, 80, 132, 124, 81, 73, 56, 99, 82, 57, 23, 125, 100, 24, 74, 93, 58, 25,  &
          &   83, 126, 75, 26, 11, 101, 133, 59, 27, 76, 140, 12, 84, 102, 94, 28, 127, 85,  &
          &   13, 39, 60, 86, 103, 87, 109, 14, 29, 115, 40, 95, 15, 61, 88, 41, 110, 104, 1,  &
          &   116, 42, 30, 134, 128, 138, 96, 62, 16, 43, 117, 63, 111, 44, 2, 64, 31, 65,  &
          &   105, 17, 45, 66, 118, 32, 3, 33, 67, 18, 129, 135, 46, 112, 34, 106, 68, 35, 4,  &
          &   119, 36, 47, 107, 19, 37, 38, 113, 48, 130, 5, 120, 49, 108, 20, 50, 51, 114,  &
          &   21, 121, 52, 136, 122, 6, 22, 7, 8, 9, 10 &
          & /)
    integer, parameter :: rrtm_gpoint_reordering_sw(112) = (/ &
          &   35, 45, 19, 27, 36, 57, 20, 46, 58, 21, 28, 67, 55, 68, 37, 1, 69, 22, 29, 59,  &
          &   78, 101, 79, 77, 70, 76, 47, 75, 30, 81, 60, 102, 80, 82, 23, 2, 83, 84, 85,  &
          &   86, 103, 61, 31, 87, 56, 38, 71, 48, 88, 3, 62, 89, 24, 7, 49, 32, 104, 72, 90,  &
          &   63, 39, 4, 8, 50, 91, 64, 40, 33, 25, 51, 95, 96, 73, 65, 9, 41, 97, 92, 105,  &
          &   52, 5, 98, 10, 42, 99, 100, 66, 11, 74, 34, 53, 26, 6, 106, 12, 43, 13, 54, 93,  &
          &   44, 107, 94, 14, 108, 15, 16, 109, 17, 18, 110, 111, 112 &
          & /)

    logical :: do_sw, do_lw
    
    real(jphook) :: hook_handle

!#include "surdi.intfb.h"

! # 1 "./include/surrtab.intfb.h" 1
interface
subroutine surrtab
end subroutine surrtab
end interface
! # 78 "radiation/radiation_ifs_rrtm.f90" 2

! # 1 "./include/surrtpk.intfb.h" 1
interface
subroutine surrtpk
end subroutine surrtpk
end interface
! # 79 "radiation/radiation_ifs_rrtm.f90" 2

! # 1 "./include/surrtrf.intfb.h" 1
interface
subroutine surrtrf
end subroutine surrtrf
end interface
! # 80 "radiation/radiation_ifs_rrtm.f90" 2

! # 1 "./include/rrtm_init_140gp.intfb.h" 1
interface
subroutine rrtm_init_140gp(directory)
character(len=*), intent(in) :: directory
end subroutine rrtm_init_140gp
end interface
! # 81 "radiation/radiation_ifs_rrtm.f90" 2

! # 1 "./include/srtm_init.intfb.h" 1
interface
subroutine srtm_init(directory, nwvcontinuum)
use parkind1, only : jpim
character(len=*), intent(in) :: directory
integer(kind=jpim), intent(in), optional :: nwvcontinuum
end subroutine srtm_init
end interface
! # 82 "radiation/radiation_ifs_rrtm.f90" 2

    if (lhook) call dr_hook('radiation_ifs_rrtm:setup_gas_optics',0,hook_handle)

    do_sw = (config%do_sw .and. config%i_gas_model_sw == igasmodelifsrrtmg)
    do_lw = (config%do_lw .and. config%i_gas_model_lw == igasmodelifsrrtmg)
    
    ! the ifs implementation of rrtmg uses many global variables.  in
    ! the ifs these will have been set up already; otherwise set them
    ! up now.
    if (config%do_setup_ifsrrtm) then
      !call surdi
      call surrtab
      call surrtpk
      call surrtrf
      if (do_lw) then
        call rrtm_init_140gp(directory)
      end if
      if (do_sw) then
        call srtm_init(directory)
      end if
    end if

    if (do_sw) then
      
      ! cloud and aerosol properties can only be defined per band
      config%do_cloud_aerosol_per_sw_g_point = .false.
      config%n_g_sw = jpgsw
      config%n_bands_sw = 14
      ! wavenumber ranges of each band may be needed so that the user
      ! can compute uv and photosynthetically active radiation for a
      ! particular wavelength range
      call config%gas_optics_sw%spectral_def%allocate_bands_only(solarreferencetemperature, &
           &  [2600.0_jprb, 3250.0_jprb, 4000.0_jprb, 4650.0_jprb, 5150.0_jprb, 6150.0_jprb, 7700.0_jprb, &
           &   8050.0_jprb, 12850.0_jprb, 16000.0_jprb, 22650.0_jprb, 29000.0_jprb, 38000.0_jprb, 820.0_jprb], &
           &  [3250.0_jprb, 4000.0_jprb, 4650.0_jprb, 5150.0_jprb, 6150.0_jprb, 7700.0_jprb, 8050.0_jprb, &
           &   12850.0_jprb, 16000.0_jprb, 22650.0_jprb, 29000.0_jprb, 38000.0_jprb, 50000.0_jprb, 2600.0_jprb])
      allocate(config%i_band_from_g_sw          (config%n_g_sw))
      allocate(config%i_band_from_reordered_g_sw(config%n_g_sw))
      allocate(config%i_g_from_reordered_g_sw   (config%n_g_sw))
      ! shortwave starts at 16: need to start at 1
      config%i_band_from_g_sw = ngb_sw - ngb_sw(1)+1

      if (config%i_solver_sw == isolverspartacus) then
        ! spartacus requires g points ordered in approximately
        ! increasing order of optical depth
        config%i_g_from_reordered_g_sw = rrtm_gpoint_reordering_sw
      else
        ! implied-do for no reordering
        !      config%i_g_from_reordered_g_sw = rrtm_gpoint_reordering_sw
        config%i_g_from_reordered_g_sw = (/ (irep, irep=1,config%n_g_sw) /)
      end if

      config%i_band_from_reordered_g_sw &
           = config%i_band_from_g_sw(config%i_g_from_reordered_g_sw)

      ! the i_spec_* variables are used solely for storing spectral
      ! data, and this can either be by band or by g-point
      if (config%do_save_spectral_flux .or. config%do_toa_spectral_flux) then
        if (config%do_save_gpoint_flux) then
          config%n_spec_sw = config%n_g_sw
          config%i_spec_from_reordered_g_sw => config%i_g_from_reordered_g_sw
        else
          config%n_spec_sw = config%n_bands_sw
          config%i_spec_from_reordered_g_sw => config%i_band_from_reordered_g_sw
        end if
      else
        config%n_spec_sw = 0
        nullify(config%i_spec_from_reordered_g_sw)
      end if
      
    end if

    if (do_lw) then
      ! cloud and aerosol properties can only be defined per band
      config%do_cloud_aerosol_per_lw_g_point = .false.
      config%n_g_lw = jpglw
      config%n_bands_lw = 16
      call config%gas_optics_lw%spectral_def%allocate_bands_only(terrestrialreferencetemperature, &
           &  [10.0_jprb, 350.0_jprb, 500.0_jprb, 630.0_jprb, 700.0_jprb, 820.0_jprb, 980.0_jprb, 1080.0_jprb, &
           &   1180.0_jprb, 1390.0_jprb, 1480.0_jprb, 1800.0_jprb, 2080.0_jprb, 2250.0_jprb, 2380.0_jprb, 2600.0_jprb], &
           &  [350.0_jprb, 500.0_jprb, 630.0_jprb, 700.0_jprb, 820.0_jprb, 980.0_jprb, 1080.0_jprb, 1180.0_jprb, &
           &   1390.0_jprb, 1480.0_jprb, 1800.0_jprb, 2080.0_jprb, 2250.0_jprb, 2380.0_jprb, 2600.0_jprb, 3250.0_jprb])
      allocate(config%i_band_from_g_lw          (config%n_g_lw))
      allocate(config%i_band_from_reordered_g_lw(config%n_g_lw))
      allocate(config%i_g_from_reordered_g_lw   (config%n_g_lw))
      config%i_band_from_g_lw = ngb_lw

      if (config%i_solver_lw == isolverspartacus) then
        ! spartacus requires g points ordered in approximately
        ! increasing order of optical depth
        config%i_g_from_reordered_g_lw = rrtm_gpoint_reordering_lw
      else
        ! implied-do for no reordering
        config%i_g_from_reordered_g_lw = (/ (irep, irep=1,config%n_g_lw) /)
      end if

      config%i_band_from_reordered_g_lw &
           = config%i_band_from_g_lw(config%i_g_from_reordered_g_lw)

      ! the i_spec_* variables are used solely for storing spectral
      ! data, and this can either be by band or by g-point
      if (config%do_save_spectral_flux .or. config%do_toa_spectral_flux) then
        if (config%do_save_gpoint_flux) then
          config%n_spec_lw = config%n_g_lw
          config%i_spec_from_reordered_g_lw => config%i_g_from_reordered_g_lw
        else
          config%n_spec_lw = config%n_bands_lw
          config%i_spec_from_reordered_g_lw => config%i_band_from_reordered_g_lw
        end if
      else
        config%n_spec_lw = 0
        nullify(config%i_spec_from_reordered_g_lw)
      end if

    end if
    
    if (lhook) call dr_hook('radiation_ifs_rrtm:setup_gas_optics',1,hook_handle)

  end subroutine setup_gas_optics


  !---------------------------------------------------------------------
  ! scale gas mixing ratios according to required units
  subroutine set_gas_units(gas)

    use radiation_gas,           only : gas_type, imassmixingratio
    type(gas_type),    intent(inout) :: gas

    call gas%set_units(imassmixingratio)

  end subroutine set_gas_units


  !---------------------------------------------------------------------
  ! compute gas optical depths, shortwave scattering, planck function
  ! and incoming shortwave radiation at top-of-atmosphere
  subroutine gas_optics(ncol,nlev,istartcol,iendcol, &
       &  config, single_level, thermodynamics, gas, & 
       &  od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
       &  incoming_sw)

    use parkind1,                 only : jprb, jpim




    use parrrtm  , only : jpband, jpxsec, jpinpx 
    use yoerrtm  , only : jpgpt_lw => jpgpt
    use yoesrtm  , only : jpgpt_sw => jpgpt  
    use ecradhook  , only : lhook, dr_hook, jphook

    use radiation_config,         only : config_type, isolverspartacus, igasmodelifsrrtmg
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_single_level,   only : single_level_type
    use radiation_gas

    integer, intent(in) :: ncol               ! number of columns
    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(gas_type),           intent(in) :: gas

    ! longwave albedo of the surface
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &  intent(in), optional :: lw_albedo

    ! gaseous layer optical depth in longwave and shortwave, and
    ! shortwave single scattering albedo (i.e. fraction of extinction
    ! due to rayleigh scattering) at each g-point
    real(jprb), dimension(config%n_g_lw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_lw
    real(jprb), dimension(config%n_g_sw,nlev,istartcol:iendcol), intent(out) :: &
         &   od_sw, ssa_sw

    ! the planck function (emitted flux from a black body) at half
    ! levels at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), &
         &   intent(out), optional :: planck_hl
    ! planck function for the surface (w m-2)
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &   intent(out), optional :: lw_emission

    ! the incoming shortwave flux into a plane perpendicular to the
    ! incoming radiation at top-of-atmosphere in each of the shortwave
    ! g-points
    real(jprb), dimension(config%n_g_sw,istartcol:iendcol), &
         &   intent(out), optional :: incoming_sw

    real(jprb) :: incoming_sw_scale(istartcol:iendcol)

    ! the variables in capitals are used in the same way as the
    ! equivalent routine in the ifs

    real(jprb) :: zod_lw(jpgpt_lw,nlev,istartcol:iendcol) ! note ordering of dimensions
    real(jprb) :: zod_sw(istartcol:iendcol,nlev,jpgpt_sw)
    real(jprb) :: zssa_sw(istartcol:iendcol,nlev,jpgpt_sw)
    real(jprb) :: zincsol(istartcol:iendcol,jpgpt_sw)

    real(jprb) :: zcolmol(istartcol:iendcol,nlev)
    real(jprb) :: zcoldry(istartcol:iendcol,nlev)
    real(jprb) :: zwbrodl(istartcol:iendcol,nlev) !broadening gases,column density (mol/cm2)
    real(jprb) :: zcolbrd(istartcol:iendcol,nlev) !broadening gases, column amount
    real(jprb) :: zwkl(istartcol:iendcol,jpinpx,nlev)

    real(jprb) :: zwx(istartcol:iendcol,jpxsec,nlev) ! amount of trace gases
    
    real(jprb) :: zfluxfac, zpi

    ! - from aer
    real(jprb) :: ztauaerl(istartcol:iendcol,nlev,jpband)

    !- from intfac      
    real(jprb) :: zfac00(istartcol:iendcol,nlev)
    real(jprb) :: zfac01(istartcol:iendcol,nlev)
    real(jprb) :: zfac10(istartcol:iendcol,nlev)
    real(jprb) :: zfac11(istartcol:iendcol,nlev)
    
    !- from for
    real(jprb) :: zforfac(istartcol:iendcol,nlev)
    real(jprb) :: zforfrac(istartcol:iendcol,nlev)
    integer    :: indfor(istartcol:iendcol,nlev) 

    !- from minor
    integer    :: indminor(istartcol:iendcol,nlev) 
    real(jprb) :: zscaleminor(istartcol:iendcol,nlev) 
    real(jprb) :: zscaleminorn2(istartcol:iendcol,nlev) 
    real(jprb) :: zminorfrac(istartcol:iendcol,nlev) 
    
    real(jprb)     :: &                 
         &  zrat_h2oco2(istartcol:iendcol,nlev),zrat_h2oco2_1(istartcol:iendcol,nlev), &
         &  zrat_h2oo3(istartcol:iendcol,nlev) ,zrat_h2oo3_1(istartcol:iendcol,nlev), & 
         &  zrat_h2on2o(istartcol:iendcol,nlev),zrat_h2on2o_1(istartcol:iendcol,nlev), &
         &  zrat_h2och4(istartcol:iendcol,nlev),zrat_h2och4_1(istartcol:iendcol,nlev), &
         &  zrat_n2oco2(istartcol:iendcol,nlev),zrat_n2oco2_1(istartcol:iendcol,nlev), &
         &  zrat_o3co2(istartcol:iendcol,nlev) ,zrat_o3co2_1(istartcol:iendcol,nlev)
    
    !- from intind
    integer :: jp(istartcol:iendcol,nlev)
    integer :: jt(istartcol:iendcol,nlev)
    integer :: jt1(istartcol:iendcol,nlev)

    !- from precise             
    real(jprb) :: zoneminus, zoneminus_array(istartcol:iendcol)

    !- from profdata             
    real(jprb) :: zcolh2o(istartcol:iendcol,nlev)
    real(jprb) :: zcolco2(istartcol:iendcol,nlev)
    real(jprb) :: zcolo3(istartcol:iendcol,nlev)
    real(jprb) :: zcoln2o(istartcol:iendcol,nlev)
    real(jprb) :: zcolch4(istartcol:iendcol,nlev)
    real(jprb) :: zcolo2(istartcol:iendcol,nlev)
    real(jprb) :: zco2mult(istartcol:iendcol,nlev)
    integer    :: ilaytrop(istartcol:iendcol)
    integer    :: ilayswtch(istartcol:iendcol)
    integer    :: ilaylow(istartcol:iendcol)

    !- from profile             
    real(jprb) :: zpavel(istartcol:iendcol,nlev)
    real(jprb) :: ztavel(istartcol:iendcol,nlev)
    real(jprb) :: zpz(istartcol:iendcol,0:nlev)
    real(jprb) :: ztz(istartcol:iendcol,0:nlev)
    
    !- from self             
    real(jprb) :: zselffac(istartcol:iendcol,nlev)
    real(jprb) :: zselffrac(istartcol:iendcol,nlev)
    integer :: indself(istartcol:iendcol,nlev)

    !- from sp             
    real(jprb) :: zpfrac(istartcol:iendcol,jpgpt_lw,nlev)
    
    !- from surface             
    integer :: ireflect(istartcol:iendcol)

    real(jprb) :: pressure_fl(ncol, nlev), temperature_fl(ncol, nlev)

    ! if nlev is less than the number of heights at which gas mixing
    ! ratios are stored, then we assume that the lower part of the
    ! atmosphere is required. this enables nlev=1 to be passed in to
    ! the routine, in which case the gas properties of the lowest
    ! layer are provided, useful for canopy radiative transfer.
    integer :: istartlev, iendlev

    logical :: do_sw, do_lw
    
    integer :: jlev, jgreorder, jg, ig, iband, jcol

    real(jphook) :: hook_handle


! # 1 "./include/rrtm_prepare_gases.intfb.h" 1
interface
subroutine rrtm_prepare_gases&
 & ( kidia, kfdia, klon, klev,&
 & paph , pap ,&
 & pth , pt ,&
 & pq , pco2 , pch4, pn2o , pno2, pc11, pc12, pc22, pcl4, pozn,&
 & pcoldry, pwbrodl, pwkl, pwx ,&
 & pavel , ptavel , pz , ptz , kreflect) 
use parkind1 , only : jpim, jprb
use parrrtm , only : jpband, jpxsec, jpinpx
integer(kind=jpim),intent(in) :: klon
integer(kind=jpim),intent(in) :: klev
integer(kind=jpim),intent(in) :: kidia, kfdia
real(kind=jprb) ,intent(in) :: paph(klon,klev+1)
real(kind=jprb) ,intent(in) :: pap(klon,klev)
real(kind=jprb) ,intent(in) :: pth(klon,klev+1)
real(kind=jprb) ,intent(in) :: pt(klon,klev)
real(kind=jprb) ,intent(in) :: pq(klon,klev)
real(kind=jprb) ,intent(in) :: pco2(klon,klev)
real(kind=jprb) ,intent(in) :: pch4(klon,klev)
real(kind=jprb) ,intent(in) :: pn2o(klon,klev)
real(kind=jprb) ,intent(in) :: pno2(klon,klev)
real(kind=jprb) ,intent(in) :: pc11(klon,klev)
real(kind=jprb) ,intent(in) :: pc12(klon,klev)
real(kind=jprb) ,intent(in) :: pc22(klon,klev)
real(kind=jprb) ,intent(in) :: pcl4(klon,klev)
real(kind=jprb) ,intent(in) :: pozn(klon,klev)
real(kind=jprb) ,intent(out) :: pcoldry(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: pwbrodl(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: pwkl(kidia:kfdia,jpinpx,klev)
real(kind=jprb) ,intent(out) :: pwx(kidia:kfdia,jpxsec,klev)
real(kind=jprb) ,intent(out) :: pavel(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: ptavel(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: pz(kidia:kfdia,0:klev)
real(kind=jprb) ,intent(out) :: ptz(kidia:kfdia,0:klev)
integer(kind=jpim),intent(out) :: kreflect(kidia:kfdia)
end subroutine rrtm_prepare_gases
end interface
! # 373 "radiation/radiation_ifs_rrtm.f90" 2

! # 1 "./include/rrtm_setcoef_140gp.intfb.h" 1
interface
subroutine rrtm_setcoef_140gp (kidia,kfdia,klev,p_coldry,p_wbroad,p_wkl,&
 & p_fac00,p_fac01,p_fac10,p_fac11,p_forfac,p_forfrac,k_indfor,k_jp,k_jt,k_jt1,&
 & p_colh2o,p_colco2,p_colo3,p_coln2o,p_colch4, p_colo2,p_co2mult, p_colbrd,&
 & k_laytrop,k_layswtch,k_laylow,pavel,p_tavel,p_selffac,p_selffrac,k_indself,&
 & k_indminor,p_scaleminor,p_scaleminorn2,p_minorfrac,&
 & prat_h2oco2, prat_h2oco2_1, prat_h2oo3, prat_h2oo3_1,&
 & prat_h2on2o, prat_h2on2o_1, prat_h2och4, prat_h2och4_1,&
 & prat_n2oco2, prat_n2oco2_1, prat_o3co2, prat_o3co2_1) 
use parkind1 , only : jpim, jprb
use parrrtm , only : jpinpx
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: p_coldry(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_wbroad(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_colbrd(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_wkl(kidia:kfdia,jpinpx,klev)
real(kind=jprb) ,intent(out) :: p_fac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_fac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_fac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_fac11(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_forfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_forfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(out) :: k_jp(kidia:kfdia,klev)
integer(kind=jpim),intent(out) :: k_jt(kidia:kfdia,klev)
integer(kind=jpim),intent(out) :: k_jt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_colh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_colco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_colo3(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_coln2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_colch4(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_colo2(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_co2mult(kidia:kfdia,klev)
integer(kind=jpim),intent(out) :: k_laytrop(kidia:kfdia)
integer(kind=jpim),intent(out) :: k_layswtch(kidia:kfdia)
integer(kind=jpim),intent(out) :: k_laylow(kidia:kfdia)
real(kind=jprb) ,intent(in) :: pavel(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_tavel(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_selffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_selffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(out) :: k_indself(kidia:kfdia,klev)
integer(kind=jpim),intent(out) :: k_indfor(kidia:kfdia,klev)
integer(kind=jpim),intent(out) :: k_indminor(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_scaleminor(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_scaleminorn2(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: p_minorfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) ::&
 & prat_h2oco2(kidia:kfdia,klev),prat_h2oco2_1(kidia:kfdia,klev),&
 & prat_h2oo3(kidia:kfdia,klev) ,prat_h2oo3_1(kidia:kfdia,klev),&
 & prat_h2on2o(kidia:kfdia,klev),prat_h2on2o_1(kidia:kfdia,klev),&
 & prat_h2och4(kidia:kfdia,klev),prat_h2och4_1(kidia:kfdia,klev),&
 & prat_n2oco2(kidia:kfdia,klev),prat_n2oco2_1(kidia:kfdia,klev),&
 & prat_o3co2(kidia:kfdia,klev) ,prat_o3co2_1(kidia:kfdia,klev) 
end subroutine rrtm_setcoef_140gp
end interface
! # 374 "radiation/radiation_ifs_rrtm.f90" 2

! # 1 "./include/rrtm_gas_optical_depth.intfb.h" 1
interface
subroutine rrtm_gas_optical_depth(kidia,kfdia,klev,pod,pavel, pcoldry,pcolbrd,pwx,&
 & ptauaerl,pfac00,pfac01,pfac10,pfac11,pforfac,pforfrac,kindfor,kjp,kjt,kjt1,poneminus,&
 & pcolh2o,pcolco2,pcolo3,pcoln2o,pcolch4,pcolo2,p_co2mult,&
 & klaytrop,klayswtch,klaylow,pselffac,pselffrac,kindself,pfrac,&
 & kindminor,pscaleminor,pscaleminorn2,pminorfrac,&
 & prat_h2oco2, prat_h2oco2_1, prat_h2oo3, prat_h2oo3_1,&
 & prat_h2on2o, prat_h2on2o_1, prat_h2och4, prat_h2och4_1,&
 & prat_n2oco2, prat_n2oco2_1, prat_o3co2, prat_o3co2_1) 
use parkind1 ,only : jpim ,jprb
use parrrtm , only : jpband ,jpxsec
use yoerrtm , only : jpgpt
integer(kind=jpim),intent(in) :: kidia
integer(kind=jpim),intent(in) :: kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(out) :: pod(jpgpt,klev,kidia:kfdia)
real(kind=jprb) ,intent(in) :: pavel(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcoldry(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pwx(kidia:kfdia,jpxsec,klev)
real(kind=jprb) ,intent(in) :: ptauaerl(kidia:kfdia,klev,jpband)
real(kind=jprb) ,intent(in) :: pfac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pfac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pfac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pfac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kjp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kjt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kjt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: poneminus
real(kind=jprb) ,intent(in) :: pcolh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolo3(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcoln2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolch4(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolo2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: p_co2mult(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: klaytrop(kidia:kfdia)
integer(kind=jpim),intent(in) :: klayswtch(kidia:kfdia)
integer(kind=jpim),intent(in) :: klaylow(kidia:kfdia)
real(kind=jprb) ,intent(in) :: pselffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pselffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindself(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: pfrac(kidia:kfdia,jpgpt,klev)
real(kind=jprb) ,intent(in) :: pforfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pforfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pminorfrac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pscaleminor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pscaleminorn2(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindminor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolbrd(kidia:kfdia,klev)
real(kind=jprb) , intent(in) ::&
 & prat_h2oco2(kidia:kfdia,klev),prat_h2oco2_1(kidia:kfdia,klev),&
 & prat_h2oo3(kidia:kfdia,klev),prat_h2oo3_1(kidia:kfdia,klev),&
 & prat_h2on2o(kidia:kfdia,klev),prat_h2on2o_1(kidia:kfdia,klev),&
 & prat_h2och4(kidia:kfdia,klev),prat_h2och4_1(kidia:kfdia,klev),&
 & prat_n2oco2(kidia:kfdia,klev),prat_n2oco2_1(kidia:kfdia,klev),&
 & prat_o3co2(kidia:kfdia,klev),prat_o3co2_1(kidia:kfdia,klev) 
end subroutine rrtm_gas_optical_depth
end interface
! # 375 "radiation/radiation_ifs_rrtm.f90" 2

! # 1 "./include/srtm_setcoef.intfb.h" 1
! this file has been modified for the use in icon

interface
subroutine srtm_setcoef&
 & ( kidia , kfdia , klev ,&
 & pavel , ptavel ,&
 & pcoldry , pwkl ,&
 & klaytrop,&
 & pcolch4 , pcolco2 , pcolh2o , pcolmol , pcolo2 , pcolo3 ,&
 & pforfac , pforfrac , kindfor , pselffac, pselffrac, kindself ,&
 & pfac00 , pfac01 , pfac10 , pfac11 ,&
 & kjp , kjt , kjt1 , prmu0&
 & ) 
use parkind1 , only : jpim, jprb
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: pavel(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: ptavel(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcoldry(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pwkl(kidia:kfdia,35,klev)
integer(kind=jpim),intent(inout) :: klaytrop(kidia:kfdia)
real(kind=jprb) ,intent(inout) :: pcolch4(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pcolco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pcolh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pcolmol(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pcolo2(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pcolo3(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pforfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pforfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(inout) :: kindfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pselffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pselffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(inout) :: kindself(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(inout) :: pfac11(kidia:kfdia,klev)
integer(kind=jpim),intent(inout) :: kjp(kidia:kfdia,klev)
integer(kind=jpim),intent(inout) :: kjt(kidia:kfdia,klev)
integer(kind=jpim),intent(inout) :: kjt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
end subroutine srtm_setcoef
end interface
! # 376 "radiation/radiation_ifs_rrtm.f90" 2

! # 1 "./include/srtm_gas_optical_depth.intfb.h" 1
interface
subroutine srtm_gas_optical_depth&
 & ( kidia , kfdia , klev , poneminus,&
 & prmu0,&
 & klaytrop,&
 & pcolch4 , pcolco2 , pcolh2o , pcolmol , pcolo2 , pcolo3 ,&
 & pforfac , pforfrac , kindfor , pselffac, pselffrac, kindself ,&
 & pfac00 , pfac01 , pfac10 , pfac11 ,&
 & kjp , kjt , kjt1 ,&
 & pod, pssa, pincsol) 
use parkind1 , only : jpim, jprb
use yoesrtm , only : jpgpt
integer(kind=jpim),intent(in) :: kidia, kfdia
integer(kind=jpim),intent(in) :: klev
real(kind=jprb) ,intent(in) :: poneminus(kidia:kfdia)
real(kind=jprb) ,intent(in) :: prmu0(kidia:kfdia)
integer(kind=jpim),intent(in) :: klaytrop(kidia:kfdia)
real(kind=jprb) ,intent(in) :: pcolch4(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolco2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolh2o(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolmol(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolo2(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pcolo3(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pforfac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pforfrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindfor(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pselffac(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pselffrac(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kindself(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pfac00(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pfac01(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pfac10(kidia:kfdia,klev)
real(kind=jprb) ,intent(in) :: pfac11(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kjp(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kjt(kidia:kfdia,klev)
integer(kind=jpim),intent(in) :: kjt1(kidia:kfdia,klev)
real(kind=jprb) ,intent(out) :: pod(kidia:kfdia,klev,jpgpt)
real(kind=jprb) ,intent(out) :: pssa(kidia:kfdia,klev,jpgpt)
real(kind=jprb) ,intent(out) :: pincsol(kidia:kfdia,jpgpt)
end subroutine srtm_gas_optical_depth
end interface
! # 377 "radiation/radiation_ifs_rrtm.f90" 2

    if (lhook) call dr_hook('radiation_ifs_rrtm:gas_optics',0,hook_handle)

    do_sw = (config%do_sw .and. config%i_gas_model_sw == igasmodelifsrrtmg)
    do_lw = (config%do_lw .and. config%i_gas_model_lw == igasmodelifsrrtmg)

! # 393 "radiation/radiation_ifs_rrtm.f90"

    !$acc data create(incoming_sw_scale, &
    !$acc             zod_lw, zod_sw, zssa_sw, zincsol, &
    !$acc             zcolmol, zcoldry, zwbrodl, zcolbrd, zwkl, &
    !$acc             zwx, &
    !$acc             ztauaerl, &
    !$acc             zfac00, zfac01, zfac10, zfac11, &
    !$acc             zforfac, zforfrac, indfor, &
    !$acc             indminor, zscaleminor, zscaleminorn2, zminorfrac, &
    !$acc             zrat_h2oco2,zrat_h2oco2_1, &
    !$acc             zrat_h2oo3 ,zrat_h2oo3_1, & 
    !$acc             zrat_h2on2o,zrat_h2on2o_1, &
    !$acc             zrat_h2och4,zrat_h2och4_1, &
    !$acc             zrat_n2oco2,zrat_n2oco2_1, &
    !$acc             zrat_o3co2 ,zrat_o3co2_1, &
    !$acc             jp, jt, jt1, &
    !$acc             zoneminus_array, &
    !$acc             zcolh2o, zcolco2, zcolo3, zcoln2o, zcolch4, zcolo2, &
    !$acc             zco2mult, &
    !$acc             ilaytrop, ilayswtch, ilaylow, &
    !$acc             zpavel, ztavel, zpz, ztz, &
    !$acc             zselffac, zselffrac, &
    !$acc             indself, &
    !$acc             zpfrac, &
    !$acc             ireflect, &
    !$acc             pressure_fl, temperature_fl) &
    !$acc     present(config, single_level, thermodynamics, gas, & 
    !$acc             od_lw, od_sw, ssa_sw, lw_albedo, planck_hl, lw_emission, &
    !$acc             incoming_sw)

    ! compute start and end levels for indexing the gas mixing ratio
    ! and thermodynamics arrays
    iendlev   = ubound(gas%mixing_ratio,2)
    istartlev = iendlev - nlev + 1

    zpi = 2.0_jprb*asin(1.0_jprb)
    zfluxfac = zpi * 1.e+4
    zoneminus = 1.0_jprb - 1.0e-6_jprb
    !$acc parallel default(none) async(1)
    !$acc loop gang vector
    do jcol= istartcol,iendcol
      zoneminus_array(jcol) = zoneminus
    end do
    !$acc end parallel

    ! are full level temperature and pressure available in thermodynmics? if not, interpolate.
    if (thermodynamics%rrtm_pass_temppres_fl) then
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2)
      do jlev=1,nlev
        do jcol= istartcol,iendcol
          pressure_fl   (jcol,jlev) = thermodynamics%pressure_fl   (jcol,jlev)
          temperature_fl(jcol,jlev) = thermodynamics%temperature_fl(jcol,jlev)
        end do
      end do
      !$acc end parallel
    else
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2)
      do jlev=1,nlev
        do jcol= istartcol,iendcol
          pressure_fl(jcol,jlev) &
              &  = 0.5_jprb * (thermodynamics%pressure_hl(jcol,jlev+istartlev-1) &
              &               +thermodynamics%pressure_hl(jcol,jlev+istartlev))
          temperature_fl(jcol,jlev) &
              &  = 0.5_jprb * (thermodynamics%temperature_hl(jcol,jlev+istartlev-1) &
              &               +thermodynamics%temperature_hl(jcol,jlev+istartlev))
        end do
      end do
      !$acc end parallel
    end if
    
    ! check we have gas mixing ratios in the right units
    call gas%assert_units(imassmixingratio)

    ! warning: o2 is hard-coded within the following function so the
    ! user-provided concentrations of this gas are ignored for both
    ! the longwave and shortwave
    call rrtm_prepare_gases &
         & ( istartcol, iendcol, ncol, nlev, &
         &   thermodynamics%pressure_hl(:,istartlev:iendlev+1), &
         &   pressure_fl, &
         &   thermodynamics%temperature_hl(:,istartlev:iendlev+1), &
         &   temperature_fl, &
         &   gas%mixing_ratio(:,istartlev:iendlev,ih2o), &
         &   gas%mixing_ratio(:,istartlev:iendlev,ico2), &
         &   gas%mixing_ratio(:,istartlev:iendlev,ich4), &
         &   gas%mixing_ratio(:,istartlev:iendlev,in2o), &
         &   gas%mixing_ratio(:,istartlev:iendlev,ino2), &
         &   gas%mixing_ratio(:,istartlev:iendlev,icfc11), &
         &   gas%mixing_ratio(:,istartlev:iendlev,icfc12), &
         &   gas%mixing_ratio(:,istartlev:iendlev,ihcfc22), &
         &   gas%mixing_ratio(:,istartlev:iendlev,iccl4), &
         &   gas%mixing_ratio(:,istartlev:iendlev,io3), &
         &  zcoldry, zwbrodl,zwkl, zwx, &
         &  zpavel , ztavel , zpz , ztz, ireflect)  

    if (do_lw) then
    
      call rrtm_setcoef_140gp &
           &( istartcol, iendcol, nlev , zcoldry  , zwbrodl , zwkl , &
           &  zfac00 , zfac01   , zfac10 , zfac11 , zforfac,zforfrac,indfor, jp, jt, jt1 , &
           &  zcolh2o, zcolco2  , zcolo3 , zcoln2o, zcolch4, zcolo2,zco2mult , zcolbrd, & 
           &  ilaytrop,ilayswtch, ilaylow, zpavel , ztavel , zselffac, zselffrac, indself, &
           &  indminor,zscaleminor,zscaleminorn2,zminorfrac,&
           &  zrat_h2oco2, zrat_h2oco2_1, zrat_h2oo3, zrat_h2oo3_1, &
           &  zrat_h2on2o, zrat_h2on2o_1, zrat_h2och4, zrat_h2och4_1, &
           &  zrat_n2oco2, zrat_n2oco2_1, zrat_o3co2, zrat_o3co2_1)   

      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(3)
      do jg = 1,jpband
        do jlev = 1,nlev
          do jcol = istartcol,iendcol
            ztauaerl(jcol,jlev,jg) = 0.0_jprb
          end do
        end do
      end do
      !$acc end parallel

      call rrtm_gas_optical_depth &
           &( istartcol, iendcol, nlev, zod_lw, zpavel, zcoldry, zcolbrd, zwx ,&
           &  ztauaerl, zfac00 , zfac01, zfac10 , zfac11 , zforfac,zforfrac,indfor, &
           &  jp, jt, jt1, zoneminus ,&
           &  zcolh2o , zcolco2, zcolo3, zcoln2o, zcolch4, zcolo2,zco2mult ,&
           &  ilaytrop, ilayswtch,ilaylow, zselffac, zselffrac, indself, zpfrac, &
           &  indminor,zscaleminor,zscaleminorn2,zminorfrac,&
           &  zrat_h2oco2, zrat_h2oco2_1, zrat_h2oo3, zrat_h2oo3_1, &
           &  zrat_h2on2o, zrat_h2on2o_1, zrat_h2och4, zrat_h2och4_1, &
           &  zrat_n2oco2, zrat_n2oco2_1, zrat_o3co2, zrat_o3co2_1)      
    
      if (present(lw_albedo)) then
    
        call planck_function_atmos(nlev, istartcol, iendcol, config, &
             &                     thermodynamics, zpfrac, planck_hl)

        if (single_level%is_simple_surface) then
          call planck_function_surf(istartcol, iendcol, config, &
               &                    single_level%skin_temperature, zpfrac(:,:,1), &
               &                    lw_emission)
          
          ! the following can be used to extract the parameters defined at
          ! the top of the planck_function routine below:
          !write(*,'(a,140(e12.5,","),a)') 'zpfrac_surf=[', &
          !&  sum(zpfrac(istartcol:iendcol,:,1),1) / (iendcol+1-istartcol), ']'
        
          ! lw_emission at this point is actually the planck function of
          ! the surface
          !$acc parallel default(none) async(1)
          !$acc loop gang vector collapse(2)
          do jcol = istartcol,iendcol
            do jg= 1,config%n_g_lw
              lw_emission(jg,jcol) = lw_emission(jg,jcol) * (1.0_jprb - lw_albedo(jg,jcol))
            end do
          end do
          !$acc end parallel
        else
          ! longwave emission has already been computed
          if (config%use_canopy_full_spectrum_lw) then
            lw_emission = transpose(single_level%lw_emission(istartcol:iendcol,:))
          else
            lw_emission = transpose(single_level%lw_emission(istartcol:iendcol, &
                 & config%i_emiss_from_band_lw(config%i_band_from_reordered_g_lw)))
          end if
        end if

      end if

      if (config%i_solver_lw == isolverspartacus) then
        ! we need to rearrange the gas optics info in memory: reordering
        ! the g points in order of approximately increasing optical
        ! depth (for efficient 3d processing on only the regions of the
        ! spectrum that are optically thin for gases) and reorder in
        ! pressure since the the functions above treat pressure
        ! decreasing with increasing index.  note that the output gas
        ! arrays have dimensions in a different order to the inputs,
        ! so there is some inefficiency here.
        do jgreorder = 1,config%n_g_lw
          iband = config%i_band_from_reordered_g_lw(jgreorder)
          ig = config%i_g_from_reordered_g_lw(jgreorder)
          
          ! top-of-atmosphere half level
          do jlev = 1,nlev
            do jcol = istartcol,iendcol
              ! some g points can return negative optical depths;
              ! specifically original g points 54-56 which causes
              ! unphysical single-scattering albedo when combined with
              ! aerosol
              od_lw(jgreorder,jlev,jcol) &
                   &   = max(config%min_gas_od_lw, zod_lw(ig,nlev+1-jlev,jcol))
            end do
          end do
        end do
      else
        ! g points have not been reordered 
        !$acc parallel default(none) async(1)
        !$acc loop gang collapse(3)
        do jcol = istartcol,iendcol
          do jlev = 1,nlev
            do jg= 1,config%n_g_lw
              ! check for negative optical depth
              od_lw(jg,jlev,jcol) = max(config%min_gas_od_lw, zod_lw(jg,nlev+1-jlev,jcol))
            end do
          end do
        end do
        !$acc end parallel
      end if

    end if

    if (do_sw) then
    
      call srtm_setcoef &
           & ( istartcol, iendcol, nlev,&
           & zpavel  , ztavel,&
           & zcoldry , zwkl,&
           & ilaytrop,&
           & zcolch4  , zcolco2 , zcolh2o , zcolmol  , zcolo2 , zcolo3,&
           & zforfac , zforfrac , indfor  , zselffac, zselffrac, indself, &
           & zfac00  , zfac01   , zfac10  , zfac11,&
           & jp      , jt       , jt1     , single_level%cos_sza(istartcol:iendcol)  &
           & )  
    
      ! srtm_gas_optical_depth will not initialize profiles when the sun
      ! is below the horizon, so we do it here
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(3)
      do jg = 1, jpgpt_sw
        do jlev = 1,nlev
          do jcol = istartcol,iendcol
            zod_sw(jcol,jlev,jg)  = 0.0_jprb
            zssa_sw(jcol,jlev,jg) = 0.0_jprb
          end do
        end do
      end do
      !$acc end parallel
      !$acc parallel default(none) async(1)
      !$acc loop gang vector collapse(2)
      do jg = 1, jpgpt_sw
        do jcol = istartcol,iendcol
          zincsol(jcol,jg)   = 0.0_jprb
        end do
      end do
      !$acc end parallel

      call srtm_gas_optical_depth &
           &( istartcol, iendcol , nlev  , zoneminus_array,&
           & single_level%cos_sza(istartcol:iendcol), ilaytrop,&
           & zcolch4 , zcolco2  , zcolh2o, zcolmol , zcolo2   , zcolo3,&
           & zforfac , zforfrac , indfor , zselffac, zselffrac, indself,&
           & zfac00  , zfac01   , zfac10 , zfac11  ,&
           & jp      , jt       , jt1    ,&
           & zod_sw  , zssa_sw  , zincsol )
    
      ! scale the incoming solar per band, if requested
      if (config%use_spectral_solar_scaling) then
        !$acc parallel default(none) async(1)
        !$acc loop gang vector collapse(2)
        do jg = 1,jpgpt_sw
          do jcol = istartcol,iendcol 
            zincsol(jcol,jg) = zincsol(jcol,jg) * &
                 &   single_level%spectral_solar_scaling(config%i_band_from_reordered_g_sw(jg))
          end do
        end do
        !$acc end parallel
      end if

      ! scaling factor to ensure that the total solar irradiance is as
      ! requested.  note that if the sun is below the horizon then
      ! zincsol will be zero.
      if (present(incoming_sw)) then
        !$acc parallel default(none) async(1)
        !$acc loop gang vector
        do jcol = istartcol,iendcol
          if (single_level%cos_sza(jcol) > 0.0_jprb) then
! added for dwd (2020)
!nec$ nounroll
            incoming_sw_scale(jcol) = single_level%solar_irradiance / sum(zincsol(jcol,:))
          else
            incoming_sw_scale(jcol) = 1.0_jprb
          end if
        end do
        !$acc end parallel
      end if

      if (config%i_solver_sw == isolverspartacus) then
!      if (.true.) then
        ! account for reordered g points
        do jgreorder = 1,config%n_g_sw
          ig = config%i_g_from_reordered_g_sw(jgreorder)
          do jlev = 1,nlev
            do jcol = istartcol,iendcol
              ! check for negative optical depth
              od_sw (jgreorder,nlev+1-jlev,jcol) &
                   &  = max(config%min_gas_od_sw, zod_sw (jcol,jlev,ig))
              ssa_sw(jgreorder,nlev+1-jlev,jcol) = zssa_sw(jcol,jlev,ig)
            end do
          end do
          if (present(incoming_sw)) then
            incoming_sw(jgreorder,:) &
                 &  = incoming_sw_scale(:) * zincsol(:,ig)
          end if
        end do
      else
        ! g points have not been reordered
        !$acc parallel default(none) num_gangs(iendcol-istartcol+1) num_workers((config%n_g_sw-1)/32+1) &
        !$acc   vector_length(32) async(1)
        !$acc loop gang
        do jcol = istartcol,iendcol
          !$acc loop seq
          do jlev = 1,nlev
            !$acc loop worker vector
            do jg = 1,config%n_g_sw
              ! check for negative optical depth
              od_sw (jg,nlev+1-jlev,jcol) = max(config%min_gas_od_sw, zod_sw(jcol,jlev,jg))
              ssa_sw(jg,nlev+1-jlev,jcol) = zssa_sw(jcol,jlev,jg)
            end do
          end do
          if (present(incoming_sw)) then
            !$acc loop worker vector
            do jg = 1,config%n_g_sw
              incoming_sw(jg,jcol) = incoming_sw_scale(jcol) * zincsol(jcol,jg)
            end do
          end if
        end do
        !$acc end parallel
      end if

    end if

    !$acc wait
    !$acc end data
    
    if (lhook) call dr_hook('radiation_ifs_rrtm:gas_optics',1,hook_handle)
    
  end subroutine gas_optics
  

  !---------------------------------------------------------------------
  ! compute planck function of the atmosphere
  subroutine planck_function_atmos(nlev,istartcol,iendcol, &
       config, thermodynamics, pfrac, &
       planck_hl)

    use parkind1,                 only : jprb, jpim

    use yoerrtm  , only : jpgpt_lw => jpgpt
    use yoerrtwn, only : totplnk, delwave

    use ecradhook, only : lhook, dr_hook, jphook

    use radiation_config,         only : config_type, isolverspartacus
    use radiation_thermodynamics, only : thermodynamics_type
    !use radiation_gas

    integer, intent(in) :: nlev               ! number of model levels
    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    type(thermodynamics_type),intent(in) :: thermodynamics
    real(jprb), intent(in) :: pfrac(istartcol:iendcol,jpgpt_lw,nlev)

    ! the planck function (emitted flux from a black body) at half
    ! levels at each longwave g-point
    real(jprb), dimension(config%n_g_lw,nlev+1,istartcol:iendcol), intent(out) :: &
         &   planck_hl

    ! planck function values per band
    real(jprb), dimension(istartcol:iendcol, config%n_bands_lw) :: planck_store

    ! look-up table variables for planck function
    real(jprb), dimension(istartcol:iendcol) :: frac
    integer,    dimension(istartcol:iendcol) :: ind

    ! temperature (k) of a half-level
    real(jprb) :: temperature

    real(jprb) :: factor, planck_tmp(istartcol:iendcol,config%n_g_lw)
    real(jprb) :: zfluxfac

    integer :: jlev, jgreorder, jg, ig, iband, jband, jcol, ilevoffset

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_atmos',0,hook_handle)

    zfluxfac = 2.0_jprb*asin(1.0_jprb) * 1.0e4_jprb
    
    ! nlev may be less than the number of original levels, in which
    ! case we assume that the user wants the lower part of the
    ! atmosphere
    ilevoffset = ubound(thermodynamics%temperature_hl,2)-nlev-1

    ! work out interpolations: for each half level, the index of the
    ! lowest interpolation bound, and the fraction into interpolation
    ! interval
    !$acc parallel default(none) create(planck_store, frac, ind, planck_tmp) present(config, thermodynamics, pfrac, &
    !$acc   planck_hl) async(1)
    !$acc loop seq
    do jlev = 1,nlev+1
      !$acc loop gang(static:1) vector private(temperature)
      do jcol = istartcol,iendcol
        temperature = thermodynamics%temperature_hl(jcol,jlev+ilevoffset)
        if (temperature < 339.0_jprb .and. temperature >= 160.0_jprb) then
          ! linear interpolation between -113 and 66 degc
          ind(jcol)  = int(temperature - 159.0_jprb)
          frac(jcol) = temperature - int(temperature)
        else if(temperature >= 339.0_jprb) then
          ! extrapolation above 66 degc
          ind(jcol)  = 180
          frac(jcol) = temperature - 339.0_jprb
        else
          ! cap below -113 degc (to avoid possible negative planck
          ! function values)
          ind(jcol)  = 1
          frac(jcol) = 0.0_jprb
        end if
      end do

      ! calculate planck functions per band
      !$acc loop seq private(factor)
      do jband = 1,config%n_bands_lw
        factor = zfluxfac * delwave(jband)
        !$acc loop gang(static:1) vector
        do jcol = istartcol,iendcol
          planck_store(jcol,jband) = factor &
               &  * (totplnk(ind(jcol),jband) &
               &  + frac(jcol)*(totplnk(ind(jcol)+1,jband)-totplnk(ind(jcol),jband)))
        end do
      end do


      if (config%i_solver_lw == isolverspartacus) then
        ! we need to rearrange the gas optics info in memory:
        ! reordering the g points in order of approximately increasing
        ! optical depth (for efficient 3d processing on only the
        ! regions of the spectrum that are optically thin for gases)
        ! and reorder in pressure since the the functions above treat
        ! pressure decreasing with increasing index.
        if (jlev == 1) then
          ! top-of-atmosphere half level - note that pfrac is on model
          ! levels not half levels
          do jgreorder = 1,config%n_g_lw
            iband = config%i_band_from_reordered_g_lw(jgreorder)
            ig = config%i_g_from_reordered_g_lw(jgreorder)
            planck_hl(jgreorder,1,:) = planck_store(:,iband) &
                 &   * pfrac(:,ig,nlev)
          end do
        else
          do jgreorder = 1,config%n_g_lw
            iband = config%i_band_from_reordered_g_lw(jgreorder)
            ig = config%i_g_from_reordered_g_lw(jgreorder)
            planck_hl(jgreorder,jlev,:) &
                   &   = planck_store(:,iband) &
                   &   * pfrac(:,ig,nlev+2-jlev)
          end do
        end if
      else

        ! g points have not been reordered 
        if (jlev == 1) then
          ! top-of-atmosphere half level - note that pfrac is on model
          ! levels not half levels
          !$acc loop seq private(iband)
          do jg = 1,config%n_g_lw
            iband = config%i_band_from_g_lw(jg)
            !$acc loop gang(static:1) vector
            do jcol = istartcol,iendcol
              planck_hl(jg,1,jcol) = planck_store(jcol,iband) * pfrac(jcol,jg,nlev)
            end do
          end do
        else
          !$acc loop seq private(iband)
          do jg = 1,config%n_g_lw
            iband = config%i_band_from_g_lw(jg)
            !$acc loop gang(static:1) vector
            do jcol = istartcol,iendcol
              planck_tmp(jcol,jg) = planck_store(jcol,iband) * pfrac(jcol,jg,nlev+2-jlev)
            end do
          end do
          !$acc loop gang(static:1) vector
          do jcol = istartcol,iendcol
            do jg = 1,config%n_g_lw
              planck_hl(jg,jlev,jcol) = planck_tmp(jcol,jg)
            end do
          end do
        end if

      end if

    end do
    !$acc end parallel

    !$acc wait

    if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_atmos',1,hook_handle)

  end subroutine planck_function_atmos


  !---------------------------------------------------------------------
  ! compute planck function of the surface
  subroutine planck_function_surf(istartcol, iendcol, config, temperature, pfrac, &
       &  planck_surf)

    use parkind1,                 only : jprb, jpim

    use yoerrtm  , only : jpgpt_lw => jpgpt
    use yoerrtwn, only : totplnk, delwave

    use ecradhook, only : lhook, dr_hook, jphook

    use radiation_config,         only : config_type, isolverspartacus
    !    use radiation_gas

    integer, intent(in) :: istartcol, iendcol ! range of columns to process
    type(config_type), intent(in) :: config
    real(jprb), intent(in) :: temperature(:)

    real(jprb), intent(in) :: pfrac(istartcol:iendcol,jpgpt_lw)

    ! planck function of the surface (w m-2)
    real(jprb), dimension(config%n_g_lw,istartcol:iendcol), &
         &  intent(out) :: planck_surf

    ! planck function values per band
    real(jprb), dimension(istartcol:iendcol, config%n_bands_lw) :: planck_store

    ! look-up table variables for planck function
    real(jprb), dimension(istartcol:iendcol) :: frac
    integer,    dimension(istartcol:iendcol) :: ind

    ! temperature (k)
    real(jprb) :: tsurf

    real(jprb) :: factor
    real(jprb) :: zfluxfac

    integer :: jgreorder, jg, ig, iband, jband, jcol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_surf',0,hook_handle)

    zfluxfac = 2.0_jprb*asin(1.0_jprb) * 1.0e4_jprb

    ! work out surface interpolations
    !$acc parallel default(none) create(planck_store, frac, ind) present(config, temperature, pfrac, planck_surf) &
    !$acc   async(1)
    !$acc loop gang(static:1) vector private(tsurf)
    do jcol = istartcol,iendcol
      tsurf = temperature(jcol)
      if (tsurf < 339.0_jprb .and. tsurf >= 160.0_jprb) then
        ! linear interpolation between -113 and 66 degc
        ind(jcol)  = int(tsurf - 159.0_jprb)
        frac(jcol) = tsurf - int(tsurf)
      else if(tsurf >= 339.0_jprb) then
        ! extrapolation above 66 degc
        ind(jcol)  = 180
        frac(jcol) = tsurf - 339.0_jprb
      else
        ! cap below -113 degc (to avoid possible negative planck
        ! function values)
        ind(jcol)  = 1
        frac(jcol) = 0.0_jprb
      end if
    end do

    ! calculate planck functions per band
    !$acc loop seq private(factor)
    do jband = 1,config%n_bands_lw
      factor = zfluxfac * delwave(jband)
      !$acc loop gang(static:1) vector
      do jcol = istartcol,iendcol
        planck_store(jcol,jband) = factor &
             &  * (totplnk(ind(jcol),jband) &
             &  + frac(jcol)*(totplnk(ind(jcol)+1,jband)-totplnk(ind(jcol),jband)))
      end do
    end do

    if (config%i_solver_lw == isolverspartacus) then
      ! we need to rearrange the gas optics info in memory: reordering
      ! the g points in order of approximately increasing optical
      ! depth (for efficient 3d processing on only the regions of
      ! the spectrum that are optically thin for gases) and reorder
      ! in pressure since the the functions above treat pressure
      ! decreasing with increasing index.
      !$acc loop seq private(iband, ig)
      do jgreorder = 1,config%n_g_lw
        iband = config%i_band_from_reordered_g_lw(jgreorder)
        ig = config%i_g_from_reordered_g_lw(jgreorder)
        !$acc loop gang(static:1) vector
        do jcol = istartcol,iendcol
          planck_surf(jgreorder,jcol) = planck_store(jcol,iband) * pfrac(jcol,ig)
        end do
      end do
    else
      ! g points have not been reordered 
      !$acc loop seq private(iband)
      do jg = 1,config%n_g_lw
        iband = config%i_band_from_g_lw(jg)
        !$acc loop gang(static:1) vector
        do jcol = istartcol,iendcol
          planck_surf(jg,jcol) = planck_store(jcol,iband) * pfrac(jcol,jg)
        end do
      end do
    end if
    !$acc end parallel

    !$acc wait

    if (lhook) call dr_hook('radiation_ifs_rrtm:planck_function_surf',1,hook_handle)
    
  end subroutine planck_function_surf


  !---------------------------------------------------------------------
  ! externally facing function for computing the planck function
  ! without reference to any gas profile; typically this would be used
  ! for computing the emission by facets of a complex surface.  note
  ! that this uses fixed "pfrac" values, obtained by averaging over
  ! those derived from rrtm-g for near-surface conditions over a line
  ! of meridian from the ecmwf model.
  subroutine planck_function(config, temperature, planck_surf)

    use parkind1,                 only : jprb, jpim

    use radiation_config,         only : config_type

    type(config_type), intent(in) :: config
    real(jprb), intent(in) :: temperature

    ! planck function of the surface (w m-2)
    real(jprb), dimension(config%n_g_lw), &
         &  intent(out) :: planck_surf

    ! fraction of each band contributed by each g-point within
    ! it. since there are 16 bands, this array sums to 16
    real(jprb), parameter, dimension(1,140) :: frac &
         = reshape( (/ 0.21227e+00, 0.18897e+00, 0.25491e+00, 0.17864e+00, 0.11735e+00, 0.38298e-01, 0.57871e-02, &
         &    0.31753e-02, 0.53169e-03, 0.76476e-04, 0.16388e+00, 0.15241e+00, 0.14290e+00, 0.12864e+00, &
         &    0.11615e+00, 0.10047e+00, 0.80013e-01, 0.60445e-01, 0.44918e-01, 0.63395e-02, 0.32942e-02, &
         &    0.54541e-03, 0.15380e+00, 0.15194e+00, 0.14339e+00, 0.13138e+00, 0.11701e+00, 0.10081e+00, &
         &    0.82296e-01, 0.61735e-01, 0.41918e-01, 0.45918e-02, 0.37743e-02, 0.30121e-02, 0.22500e-02, &
         &    0.14490e-02, 0.55410e-03, 0.78364e-04, 0.15938e+00, 0.15146e+00, 0.14213e+00, 0.13079e+00, &
         &    0.11672e+00, 0.10053e+00, 0.81566e-01, 0.61126e-01, 0.41150e-01, 0.44488e-02, 0.36950e-02, &
         &    0.29101e-02, 0.21357e-02, 0.19609e-02, 0.14134e+00, 0.14390e+00, 0.13913e+00, 0.13246e+00, &
         &    0.12185e+00, 0.10596e+00, 0.87518e-01, 0.66164e-01, 0.44862e-01, 0.49402e-02, 0.40857e-02, &
         &    0.32288e-02, 0.23613e-02, 0.15406e-02, 0.58258e-03, 0.82171e-04, 0.29127e+00, 0.28252e+00, &
         &    0.22590e+00, 0.14314e+00, 0.45494e-01, 0.71792e-02, 0.38483e-02, 0.65712e-03, 0.29810e+00, &
         &    0.27559e+00, 0.11997e+00, 0.10351e+00, 0.84515e-01, 0.62253e-01, 0.41050e-01, 0.44217e-02, &
         &    0.36946e-02, 0.29113e-02, 0.34290e-02, 0.55993e-03, 0.31441e+00, 0.27586e+00, 0.21297e+00, &
         &    0.14064e+00, 0.45588e-01, 0.65665e-02, 0.34232e-02, 0.53199e-03, 0.19811e+00, 0.16833e+00, &
         &    0.13536e+00, 0.11549e+00, 0.10649e+00, 0.93264e-01, 0.75720e-01, 0.56405e-01, 0.41865e-01, &
         &    0.59331e-02, 0.26510e-02, 0.40040e-03, 0.32328e+00, 0.26636e+00, 0.21397e+00, 0.14038e+00, &
         &    0.52142e-01, 0.38852e-02, 0.14601e+00, 0.13824e+00, 0.27703e+00, 0.22388e+00, 0.15446e+00, &
         &    0.48687e-01, 0.98054e-02, 0.18870e-02, 0.11961e+00, 0.12106e+00, 0.13215e+00, 0.13516e+00, &
         &    0.25249e+00, 0.16542e+00, 0.68157e-01, 0.59725e-02, 0.49258e+00, 0.33651e+00, 0.16182e+00, &
         &    0.90984e-02, 0.95202e+00, 0.47978e-01, 0.91716e+00, 0.82857e-01, 0.77464e+00, 0.22536e+00 /), (/ 1,140 /) )

    call planck_function_surf(1, 1, config, spread(temperature,1,1), &
         &                    frac, planck_surf)

  end subroutine planck_function

end module radiation_ifs_rrtm
! #define __atomic_acquire 2
! #define __char_bit__ 8
! #define __float_word_order__ __order_little_endian__
! #define __order_little_endian__ 1234
! #define __order_pdp_endian__ 3412
! #define __gfc_real_10__ 1
! #define __finite_math_only__ 0
! #define __gnuc_patchlevel__ 0
! #define __gfc_int_2__ 1
! #define __sizeof_int__ 4
! #define __sizeof_pointer__ 8
! #define __gfortran__ 1
! #define __gfc_real_16__ 1
! #define __stdc_hosted__ 0
! #define __no_math_errno__ 1
! #define __sizeof_float__ 4
! #define __pic__ 2
! #define _language_fortran 1
! #define __sizeof_long__ 8
! #define __gfc_int_8__ 1
! #define __dynamic__ 1
! #define __sizeof_short__ 2
! #define __gnuc__ 13
! #define __sizeof_long_double__ 16
! #define __biggest_alignment__ 16
! #define __atomic_relaxed 0
! #define _lp64 1
! #define __ecrad_little_endian 1
! #define __gfc_int_1__ 1
! #define __order_big_endian__ 4321
! #define __byte_order__ __order_little_endian__
! #define __sizeof_size_t__ 8
! #define __pic__ 2
! #define __sizeof_double__ 8
! #define __atomic_consume 1
! #define __gnuc_minor__ 3
! #define __gfc_int_16__ 1
! #define __lp64__ 1
! #define __atomic_seq_cst 5
! #define __sizeof_long_long__ 8
! #define __atomic_acq_rel 4
! #define __atomic_release 3
! #define __version__ "13.3.0"

