! # 1 "radiation/radiation_save.f90"
! # 1 "<built-in>"
! # 1 "<command-line>"
! # 1 "/users/pmz/gitspace/icon-model/externals/ecrad//"
! # 1 "radiation/radiation_save.f90"
! radiation_save.f90 - save data to netcdf files
!
! (c) copyright 2014- ecmwf.
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
!   2017-04-22  r. hogan  adapt for new way of describing longwave properties
!   2019-01-02  r. hogan  only save cloud properties if do_clouds==.true.

module radiation_save

  use parkind1, only : jprb

  implicit none

  ! save final fluxes and save intermediate radiative properties
  public :: save_fluxes, save_radiative_properties, save_inputs
  ! save net fluxes ifs style, where upwelling fluxes are actually net down
  public :: save_net_fluxes

contains

  !---------------------------------------------------------------------
  ! save fluxes in "flux" to netcdf file_name, plus pressure from the
  ! thermodynamics object
  subroutine save_fluxes(file_name, config, thermodynamics, flux, &
       &                 iverbose, is_hdf5_file, experiment_name, &
       &                 is_double_precision)

    use ecradhook,                  only : lhook, dr_hook, jphook

    use easy_netcdf

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, igasmodelmonochromatic
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_flux,           only : flux_type

    character(len=*),           intent(in) :: file_name
    type(config_type),          intent(in) :: config
    type(thermodynamics_type),  intent(in) :: thermodynamics
    type(flux_type),            intent(in) :: flux
    integer,          optional, intent(in) :: iverbose
    logical,          optional, intent(in) :: is_hdf5_file
    logical,          optional, intent(in) :: is_double_precision
    character(len=*), optional, intent(in) :: experiment_name

    type(netcdf_file)                      :: out_file
    integer                                :: ncol, n_lev_plus1
    character(5), parameter                :: default_lw_units_str = 'w m-2'
    character(5)                           :: lw_units_str
    integer                                :: i_local_verbose

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_save:save_fluxes',0,hook_handle)
    
    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

    ! work out array dimensions
    if (config%do_sw) then
      ncol = size(flux%sw_up,1)
      n_lev_plus1 = size(flux%sw_up,2)
    elseif (config%do_lw) then
      ncol = size(flux%lw_up,1)
      n_lev_plus1 = size(flux%lw_up,2)
    else
      if (i_local_verbose >= 1) then
        write(nulout,'(a,a,a)') 'warning: neither longwave nor shortwave computed so ', &
             &                  trim(file_name),' not written'
      end if
      return
    end if

    if (config%i_gas_model_lw == igasmodelmonochromatic &
         .and. config%mono_lw_wavelength > 0.0_jprb) then
      lw_units_str = 'w m-3'
    else
      lw_units_str = default_lw_units_str
    end if

    ! open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose, is_hdf5_file=is_hdf5_file)

    ! variables stored internally with column varying fastest, but in
    ! output file column varies most slowly so need to transpose
    call out_file%transpose_matrices(.true.)

    ! set default precision for file, if specified
    if (present(is_double_precision)) then
      call out_file%double_precision(is_double_precision)
    end if

    ! spectral fluxes in memory are dimensioned (nband,ncol,nlev), but
    ! are reoriented in the output file to be (nband,nlev,ncol), where
    ! the convention here is first dimension varying fastest
    call out_file%permute_3d_arrays( (/ 1, 3, 2 /) )

    ! define dimensions
    call out_file%define_dimension("column", ncol)
    call out_file%define_dimension("half_level", n_lev_plus1)

    if (config%do_save_spectral_flux .or. config%do_toa_spectral_flux) then
      if (config%do_lw) then
        call out_file%define_dimension("band_lw", config%n_spec_lw)
      end if
      if (config%do_sw) then
        call out_file%define_dimension("band_sw", config%n_spec_sw)
      end if
    else if (config%do_surface_sw_spectral_flux) then
      if (config%do_sw) then
        call out_file%define_dimension("band_sw", config%n_bands_sw)
      end if
    end if

    if (config%do_lw .and. config%do_canopy_fluxes_lw) then
      call out_file%define_dimension("canopy_band_lw", &
           &  size(flux%lw_dn_surf_canopy, 1))
    end if
    if (config%do_sw .and. config%do_canopy_fluxes_sw) then
      call out_file%define_dimension("canopy_band_sw", &
           &  size(flux%sw_dn_diffuse_surf_canopy, 1))
    end if

    ! put global attributes
    call out_file%put_global_attributes( &
         &   title_str="radiative flux profiles from the ecrad offline radiation model", &
         &   references_str="hogan, r. j., and a. bozzo, 2018: a flexible and efficient radiation " &
         &   //"scheme for the ecmwf model. j. adv. modeling earth sys., 10, 1990–2008", &
         &   source_str="ecrad offline radiation model")

    ! save "experiment" global attribute if present and not empty
    if (present(experiment_name)) then
      if (experiment_name /= " ") then
        call out_file%put_global_attribute("experiment", experiment_name)
      end if
    end if

    ! define variables
    call out_file%define_variable("pressure_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="pa", long_name="pressure", &
         &   standard_name="air_pressure")

    if (config%do_lw) then
      call out_file%define_variable("flux_up_lw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str=lw_units_str, long_name="upwelling longwave flux", &
           &   standard_name="upwelling_longwave_flux_in_air")
      call out_file%define_variable("flux_dn_lw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str=lw_units_str, long_name="downwelling longwave flux", &
           &   standard_name="downwelling_longwave_flux_in_air")
      if (config%do_clear) then
        call out_file%define_variable("flux_up_lw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str=lw_units_str, &
             &   long_name="upwelling clear-sky longwave flux")
        call out_file%define_variable("flux_dn_lw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str=lw_units_str, &
             &   long_name="downwelling clear-sky longwave flux")
      end if

      if (config%do_lw_derivatives) then
        call out_file%define_variable("lw_derivative", &
             &  dim2_name="column", dim1_name="half_level", &
             &  units_str="1", &
             &  long_name="derivative of upwelling lw flux w.r.t. surface value")
      end if

      if (config%do_save_spectral_flux) then
        call out_file%define_variable("spectral_flux_up_lw", &
             &   dim3_name="column", dim2_name="half_level", &
             &   dim1_name="band_lw", units_str=lw_units_str, &
             &   long_name="spectral upwelling longwave flux")
        call out_file%define_variable("spectral_flux_dn_lw", &
             &   dim3_name="column", dim2_name="half_level", &
             &   dim1_name="band_lw", units_str=lw_units_str, &
             &   long_name="spectral downwelling longwave flux")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_lw_clear", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_lw", units_str=lw_units_str, &
               &   long_name="spectral upwelling clear-sky longwave flux")
          call out_file%define_variable("spectral_flux_dn_lw_clear", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_lw", units_str=lw_units_str, &
               &   long_name="spectral downwelling clear-sky longwave flux")
        end if
      end if
   
      if (config%do_toa_spectral_flux) then
        call out_file%define_variable("spectral_flux_up_lw_toa", &
             &   dim2_name="column", dim1_name="band_lw", units_str="w m-2", &
             &   long_name="spectral upwelling longwave flux at top-of-atmosphere")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_lw_toa_clear", &
               &   dim2_name="column", dim1_name="band_lw", units_str="w m-2", &
               &   long_name="spectral upwelling clear-sky longwave flux at top-of-atmosphere")
        end if
      end if
   
      if (config%do_canopy_fluxes_lw) then
        call out_file%define_variable("canopy_flux_dn_lw_surf", &
             &   dim2_name="column", dim1_name="canopy_band_lw", units_str=lw_units_str, &
             &   long_name="surface downwelling longwave flux in canopy bands")
      end if

    end if

    if (config%do_sw) then
      call out_file%define_variable("flux_up_sw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str="w m-2", long_name="upwelling shortwave flux", &
           &   standard_name="upwelling_shortwave_flux_in_air")
      call out_file%define_variable("flux_dn_sw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str="w m-2", long_name="downwelling shortwave flux", &
           &   standard_name="downwelling_shortwave_flux_in_air")
      if (config%do_sw_direct) then
        call out_file%define_variable("flux_dn_direct_sw", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str="w m-2", &
             &   long_name="downwelling direct shortwave flux")
      end if
      if (config%do_clear) then
        call out_file%define_variable("flux_up_sw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str="w m-2", &
             &   long_name="upwelling clear-sky shortwave flux")
        call out_file%define_variable("flux_dn_sw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str="w m-2", &
             &   long_name="downwelling clear-sky shortwave flux")
        if (config%do_sw_direct) then
          call out_file%define_variable("flux_dn_direct_sw_clear", &
               &   dim2_name="column", dim1_name="half_level", &
               &   units_str="w m-2", &
               &   long_name="downwelling clear-sky direct shortwave flux")
        end if
      end if

      if (config%do_save_spectral_flux) then
        call out_file%define_variable("spectral_flux_up_sw", &
             &   dim3_name="column", dim2_name="half_level", &
             &   dim1_name="band_sw", units_str="w m-2", &
             &   long_name="spectral upwelling shortwave flux")
        call out_file%define_variable("spectral_flux_dn_sw", &
             &   dim3_name="column", dim2_name="half_level", &
             &   dim1_name="band_sw", units_str="w m-2", &
             &   long_name="spectral downwelling shortwave flux")
        if (config%do_sw_direct) then
          call out_file%define_variable("spectral_flux_dn_direct_sw", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_sw", units_str="w m-2", &
               &   long_name="spectral downwelling direct shortwave flux")
        end if
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_sw_clear", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_sw", units_str="w m-2", &
               &   long_name="spectral upwelling clear-sky shortwave flux")
          call out_file%define_variable("spectral_flux_dn_sw_clear", &
               &   dim3_name="column", dim2_name="half_level", &
               &   dim1_name="band_sw", units_str="w m-2", &
               &   long_name="spectral downwelling clear-sky shortwave flux")
          if (config%do_sw_direct) then
            call out_file%define_variable("spectral_flux_dn_direct_sw_clear", &
                 &   dim3_name="column", dim2_name="half_level", &
                 &   dim1_name="band_sw", units_str="w m-2", &
                 &   long_name="spectral downwelling clear-sky direct shortwave flux")
          end if
        end if
      else if (config%do_surface_sw_spectral_flux) then
        call out_file%define_variable("spectral_flux_dn_sw_surf", &
             &   dim2_name="column", dim1_name="band_sw", units_str="w m-2", &
             &   long_name="spectral downwelling shortwave flux at surface")
        call out_file%define_variable("spectral_flux_dn_direct_sw_surf", &
             &   dim2_name="column", dim1_name="band_sw", units_str="w m-2", &
             &   long_name="spectral downwelling direct shortwave flux at surface")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_dn_sw_surf_clear", &
               &   dim2_name="column", dim1_name="band_sw", units_str="w m-2", &
               &   long_name="spectral downwelling clear-sky shortwave flux at surface")
          call out_file%define_variable("spectral_flux_dn_direct_sw_surf_clear", &
               &   dim2_name="column", dim1_name="band_sw", units_str="w m-2", &
               &   long_name="spectral downwelling clear-sky direct shortwave flux at surface")
        end if
      end if

      if (config%do_toa_spectral_flux) then
        call out_file%define_variable("spectral_flux_dn_sw_toa", &
             &   dim2_name="column", dim1_name="band_sw", units_str="w m-2", &
             &   long_name="spectral downwelling shortwave flux at top-of-atmosphere")
        call out_file%define_variable("spectral_flux_up_sw_toa", &
             &   dim2_name="column", dim1_name="band_sw", units_str="w m-2", &
             &   long_name="spectral upwelling shortwave flux at top-of-atmosphere")
        if (config%do_clear) then
          call out_file%define_variable("spectral_flux_up_sw_toa_clear", &
               &   dim2_name="column", dim1_name="band_sw", units_str="w m-2", &
               &   long_name="spectral upwelling clear-sky shortwave flux at top-of-atmosphere")
        end if
      end if
   
      if (config%do_canopy_fluxes_sw) then
        call out_file%define_variable("canopy_flux_dn_diffuse_sw_surf", &
             &   dim2_name="column", dim1_name="canopy_band_sw", units_str="w m-2", &
             &   long_name="surface downwelling diffuse shortwave flux in canopy bands")
        call out_file%define_variable("canopy_flux_dn_direct_sw_surf", &
             &   dim2_name="column", dim1_name="canopy_band_sw", units_str="w m-2", &
             &   long_name="surface downwelling direct shortwave flux in canopy bands")
      end if

    end if
   
    if (config%do_lw .and. config%do_clouds) then
      call out_file%define_variable("cloud_cover_lw", &
           &  dim1_name="column", units_str="1", &
           &  long_name="total cloud cover diagnosed by longwave solver", &
           &  standard_name="cloud_area_fraction")
    end if
    if (config%do_sw .and. config%do_clouds) then
      call out_file%define_variable("cloud_cover_sw", &
           &  dim1_name="column", units_str="1", &
           &  long_name="total cloud cover diagnosed by shortwave solver", &
           &  standard_name="cloud_area_fraction")
    end if

    ! write variables

    call out_file%put("pressure_hl", thermodynamics%pressure_hl)

    if (config%do_lw) then
      call out_file%put("flux_up_lw", flux%lw_up)
      call out_file%put("flux_dn_lw", flux%lw_dn)
      if (config%do_clear) then
        call out_file%put("flux_up_lw_clear", flux%lw_up_clear)
        call out_file%put("flux_dn_lw_clear", flux%lw_dn_clear)
      end if

      if (config%do_lw_derivatives) then
        call out_file%put("lw_derivative", flux%lw_derivatives)
      end if

      if (config%do_save_spectral_flux) then
        call out_file%put("spectral_flux_up_lw", flux%lw_up_band)
        call out_file%put("spectral_flux_dn_lw", flux%lw_dn_band)
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_lw_clear", flux%lw_up_clear_band)
          call out_file%put("spectral_flux_dn_lw_clear", flux%lw_dn_clear_band)
        end if
      end if

      if (config%do_toa_spectral_flux) then
        call out_file%put("spectral_flux_up_lw_toa", flux%lw_up_toa_band, &
               &   do_transp=.false.)
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_lw_toa_clear", flux%lw_up_toa_clear_band, &
               &   do_transp=.false.)
        end if
      end if
      
      if (config%do_canopy_fluxes_lw) then
        call out_file%put("canopy_flux_dn_lw_surf", flux%lw_dn_surf_canopy, &
             &            do_transp = .false.)
      end if

    end if

    if (config%do_sw) then
      call out_file%put("flux_up_sw", flux%sw_up)
      call out_file%put("flux_dn_sw", flux%sw_dn)
      if (config%do_sw_direct) then
        call out_file%put("flux_dn_direct_sw", flux%sw_dn_direct)
      end if
      if (config%do_clear) then
        call out_file%put("flux_up_sw_clear", flux%sw_up_clear)
        call out_file%put("flux_dn_sw_clear", flux%sw_dn_clear)
        if (config%do_sw_direct) then
          call out_file%put("flux_dn_direct_sw_clear", flux%sw_dn_direct_clear)
        end if
      end if

      if (config%do_save_spectral_flux) then
        call out_file%put("spectral_flux_up_sw", flux%sw_up_band)
        call out_file%put("spectral_flux_dn_sw", flux%sw_dn_band)
        if (config%do_sw_direct) then
          call out_file%put("spectral_flux_dn_direct_sw", &
               &   flux%sw_dn_direct_band)
        end if
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_sw_clear", flux%sw_up_clear_band)
          call out_file%put("spectral_flux_dn_sw_clear", flux%sw_dn_clear_band)
          if (config%do_sw_direct) then
            call out_file%put("spectral_flux_dn_direct_sw_clear", &
                 &   flux%sw_dn_direct_clear_band)
          end if
        end if
      else if (config%do_surface_sw_spectral_flux) then
        call out_file%put("spectral_flux_dn_sw_surf", flux%sw_dn_surf_band, &
               &   do_transp=.false.)
        call out_file%put("spectral_flux_dn_direct_sw_surf", flux%sw_dn_direct_surf_band, &
               &   do_transp=.false.)
        if (config%do_clear) then
          call out_file%put("spectral_flux_dn_sw_surf_clear", flux%sw_dn_surf_clear_band, &
               &   do_transp=.false.)
          call out_file%put("spectral_flux_dn_direct_sw_surf_clear", &
               &            flux%sw_dn_direct_surf_clear_band, do_transp=.false.)
        end if
      end if

      if (config%do_toa_spectral_flux) then
        call out_file%put("spectral_flux_dn_sw_toa", flux%sw_dn_toa_band, &
               &   do_transp=.false.)
        call out_file%put("spectral_flux_up_sw_toa", flux%sw_up_toa_band, &
               &   do_transp=.false.)
        if (config%do_clear) then
          call out_file%put("spectral_flux_up_sw_toa_clear", flux%sw_up_toa_clear_band, &
               &   do_transp=.false.)
        end if
      end if
      
      if (config%do_canopy_fluxes_sw) then
        call out_file%put("canopy_flux_dn_diffuse_sw_surf", flux%sw_dn_diffuse_surf_canopy, &
             &            do_transp = .false.)
        call out_file%put("canopy_flux_dn_direct_sw_surf",  flux%sw_dn_direct_surf_canopy, &
             &            do_transp = .false.)
      end if

    end if

    if (config%do_lw .and. config%do_clouds) then
      call out_file%put("cloud_cover_lw", flux%cloud_cover_lw)
    end if
    if (config%do_sw .and. config%do_clouds) then
      call out_file%put("cloud_cover_sw", flux%cloud_cover_sw)
    end if

    ! close file
    call out_file%close()

    if (lhook) call dr_hook('radiation_save:save_fluxes',1,hook_handle)

  end subroutine save_fluxes
  

  !---------------------------------------------------------------------
  ! save ifs-style net fluxes in "flux" to netcdf file_name, plus
  ! pressure from the thermodynamics object
  subroutine save_net_fluxes(file_name, config, thermodynamics, flux, &
       &                     iverbose, is_hdf5_file, experiment_name, &
       &                     is_double_precision)

    use ecradhook,                  only : lhook, dr_hook, jphook

    use easy_netcdf

    use radiation_io,             only : nulout
    use radiation_config,         only : config_type, igasmodelmonochromatic
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_flux,           only : flux_type

    character(len=*),           intent(in) :: file_name
    type(config_type),          intent(in) :: config
    type(thermodynamics_type),  intent(in) :: thermodynamics
    type(flux_type),            intent(in) :: flux
    integer,          optional, intent(in) :: iverbose
    logical,          optional, intent(in) :: is_hdf5_file
    logical,          optional, intent(in) :: is_double_precision
    character(len=*), optional, intent(in) :: experiment_name

    type(netcdf_file)                      :: out_file
    integer                                :: ncol, n_lev_plus1
    character(5), parameter                :: default_lw_units_str = 'w m-2'
    character(5)                           :: lw_units_str
    integer                                :: i_local_verbose

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_save:save_net_fluxes',0,hook_handle)
    
    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

    ! work out array dimensions
    if (config%do_sw) then
      ncol = size(flux%sw_up,1)
      n_lev_plus1 = size(flux%sw_up,2)
    elseif (config%do_lw) then
      ncol = size(flux%lw_up,1)
      n_lev_plus1 = size(flux%lw_up,2)
    else
      if (i_local_verbose >= 1) then
        write(nulout,'(a,a,a)') 'warning: neither longwave nor shortwave computed so ', &
             &                  file_name,' not written'
      end if
      return
    end if

    if (config%i_gas_model_lw == igasmodelmonochromatic &
         .and. config%mono_lw_wavelength > 0.0_jprb) then
      lw_units_str = 'w m-3'
    else
      lw_units_str = default_lw_units_str
    end if

    ! open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose, is_hdf5_file=is_hdf5_file)

    ! variables stored internally with column varying fastest, but in
    ! output file column varies most slowly so need to transpose
    call out_file%transpose_matrices(.true.)

    ! set default precision for file, if specified
    if (present(is_double_precision)) then
      call out_file%double_precision(is_double_precision)
    end if

    ! spectral fluxes in memory are dimensioned (nband,ncol,nlev), but
    ! are reoriented in the output file to be (nband,nlev,ncol), where
    ! the convention here is first dimension varying fastest
    call out_file%permute_3d_arrays( (/ 1, 3, 2 /) )

    ! define dimensions
    call out_file%define_dimension("column", ncol)
    call out_file%define_dimension("half_level", n_lev_plus1)

    if (config%do_lw .and. config%do_canopy_fluxes_lw) then
      call out_file%define_dimension("canopy_band_lw", &
           &  size(flux%lw_dn_surf_canopy, 1))
    end if
    if (config%do_sw .and. config%do_canopy_fluxes_sw) then
      call out_file%define_dimension("canopy_band_sw", &
           &  size(flux%sw_dn_diffuse_surf_canopy, 1))
    end if

    ! put global attributes
    call out_file%put_global_attributes( &
         &   title_str="radiative flux profiles from the ecrad offline radiation model", &
         &   references_str="hogan, r. j., and a. bozzo, 2018: a flexible and efficient radiation " &
         &   //"scheme for the ecmwf model. j. adv. modeling earth sys., 10, 1990–2008", &
         &   source_str="ecrad offline radiation model")

    ! save "experiment" global attribute if present and not empty
    if (present(experiment_name)) then
      if (experiment_name /= " ") then
        call out_file%put_global_attribute("experiment", experiment_name)
      end if
    end if

    ! define variables
    call out_file%define_variable("pressure_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="pa", long_name="pressure", &
         &   standard_name="air_pressure")

    if (config%do_lw) then
      call out_file%define_variable("flux_net_lw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str=lw_units_str, long_name="net downward longwave flux", &
           &   standard_name="net_downward_longwave_flux_in_air")
      call out_file%define_variable("flux_dn_lw_surf", &
           &   dim1_name="column", &
           &   units_str=lw_units_str, long_name="surface downwelling longwave flux", &
           &   standard_name="surface_downwelling_longwave_flux_in_air")
      if (config%do_clear) then
        call out_file%define_variable("flux_net_lw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str=lw_units_str, &
             &   long_name="net downward clear-sky longwave flux", &
             &   standard_name="net_downward_longwave_flux_in_air_assuming_clear_sky")
        call out_file%define_variable("flux_dn_lw_clear_surf", &
             &   dim1_name="column", &
             &   units_str=lw_units_str, &
             &   long_name="surface downwelling clear-sky longwave flux", &
             &   standard_name="surface_downwelling_longwave_flux_in_air_assuming_clear_sky")
      end if

      if (config%do_lw_derivatives) then
        call out_file%define_variable("lw_derivative", &
             &  dim2_name="column", dim1_name="half_level", &
             &  units_str="1", &
             &  long_name="derivative of upwelling lw flux w.r.t. surface value")
      end if
   
      if (config%do_canopy_fluxes_lw) then
        call out_file%define_variable("canopy_flux_dn_lw_surf", &
             &   dim2_name="column", dim1_name="canopy_band_lw", units_str=lw_units_str, &
             &   long_name="surface downwelling longwave flux in canopy bands")
      end if

    end if

    if (config%do_sw) then
      call out_file%define_variable("flux_net_sw", &
           &   dim2_name="column", dim1_name="half_level", &
           &   units_str="w m-2", long_name="net downward shortwave flux", &
           &   standard_name="net_downward_shortwave_flux_in_air")
      call out_file%define_variable("flux_dn_sw_surf", &
           &   dim1_name="column", &
           &   units_str="w m-2", long_name="surface downwelling shortwave flux", &
           &   standard_name="surface_downwelling_shortwave_flux_in_air")
      call out_file%define_variable("flux_dn_sw_toa", &
           &   dim1_name="column", &
           &   units_str="w m-2", long_name="top-of-atmosphere downwelling shortwave flux", &
           &   standard_name="toa_incoming_shortwave_flux")
      if (config%do_sw_direct) then
        call out_file%define_variable("flux_dn_direct_sw_surf", &
             &   dim1_name="column", &
             &   units_str="w m-2", &
             &   long_name="surface downwelling direct shortwave flux")
      end if
      if (config%do_clear) then
        call out_file%define_variable("flux_net_sw_clear", &
             &   dim2_name="column", dim1_name="half_level", &
             &   units_str="w m-2", &
             &   long_name="net downward clear-sky shortwave flux")
        call out_file%define_variable("flux_dn_sw_clear_surf", &
             &   dim1_name="column", &
             &   units_str="w m-2", &
             &   long_name="surface downwelling clear-sky shortwave flux")
        if (config%do_sw_direct) then
          call out_file%define_variable("flux_dn_direct_sw_clear_surf", &
               &   dim1_name="column", &
               &   units_str="w m-2", &
               &   long_name="surface downwelling clear-sky direct shortwave flux")
        end if
      end if
   
      if (config%do_canopy_fluxes_sw) then
        call out_file%define_variable("canopy_flux_dn_diffuse_sw_surf", &
             &   dim2_name="column", dim1_name="canopy_band_sw", units_str="w m-2", &
             &   long_name="surface downwelling diffuse shortwave flux in canopy bands")
        call out_file%define_variable("canopy_flux_dn_direct_sw_surf", &
             &   dim2_name="column", dim1_name="canopy_band_sw", units_str="w m-2", &
             &   long_name="surface downwelling direct shortwave flux in canopy bands")
      end if

    end if
   
    ! write variables

    call out_file%put("pressure_hl", thermodynamics%pressure_hl)

    if (config%do_lw) then
      call out_file%put("flux_net_lw", flux%lw_dn-flux%lw_up)
      call out_file%put("flux_dn_lw_surf", flux%lw_dn(:,n_lev_plus1))
      if (config%do_clear) then
        call out_file%put("flux_net_lw_clear", flux%lw_dn_clear-flux%lw_up_clear)
        call out_file%put("flux_dn_lw_clear_surf", flux%lw_dn_clear(:,n_lev_plus1))
      end if

      if (config%do_lw_derivatives) then
        call out_file%put("lw_derivative", flux%lw_derivatives)
      end if

      if (config%do_canopy_fluxes_lw) then
        call out_file%put("canopy_flux_dn_lw_surf", flux%lw_dn_surf_canopy, &
             &            do_transp = .false.)
      end if

    end if

    if (config%do_sw) then
      call out_file%put("flux_net_sw", flux%sw_dn-flux%sw_up)
      call out_file%put("flux_dn_sw_surf", flux%sw_dn(:,n_lev_plus1))
      call out_file%put("flux_dn_sw_toa", flux%sw_dn(:,1))
      if (config%do_sw_direct) then
        call out_file%put("flux_dn_direct_sw_surf", flux%sw_dn_direct(:,n_lev_plus1))
      end if
      if (config%do_clear) then
        call out_file%put("flux_net_sw_clear", flux%sw_dn_clear-flux%sw_up_clear)
        call out_file%put("flux_dn_sw_clear_surf", flux%sw_dn_clear(:,n_lev_plus1))
        if (config%do_sw_direct) then
          call out_file%put("flux_dn_direct_sw_clear_surf", flux%sw_dn_direct_clear(:,n_lev_plus1))
        end if
      end if

      if (config%do_canopy_fluxes_sw) then
        call out_file%put("canopy_flux_dn_diffuse_sw_surf", flux%sw_dn_diffuse_surf_canopy, &
             &            do_transp = .false.)
        call out_file%put("canopy_flux_dn_direct_sw_surf",  flux%sw_dn_direct_surf_canopy, &
             &            do_transp = .false.)
      end if

    end if

    ! close file
    call out_file%close()

    if (lhook) call dr_hook('radiation_save:save_net_fluxes',1,hook_handle)

  end subroutine save_net_fluxes
  

  !---------------------------------------------------------------------
  ! save intermediate radiative properties, specifically the
  ! scattering and absorption properties at each g-point/band
  subroutine save_radiative_properties(file_name, nlev, &
       &  istartcol, iendcol, &
       &  config, single_level, thermodynamics, cloud, &
       &  planck_hl, lw_emission, lw_albedo, &
       &  sw_albedo_direct, sw_albedo_diffuse, &
       &  incoming_sw, &
       &  od_lw, ssa_lw, g_lw, &
       &  od_sw, ssa_sw, g_sw, &
       &  od_lw_cloud, ssa_lw_cloud, g_lw_cloud, &
       &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud)

    use radiation_config,        only : config_type
    use radiation_single_level,  only : single_level_type
    use radiation_thermodynamics,only : thermodynamics_type
    use radiation_cloud,         only : cloud_type
    use easy_netcdf

    character(len=*),         intent(in) :: file_name
    type(config_type),        intent(in) :: config
    type(single_level_type),  intent(in) :: single_level
    type(thermodynamics_type),intent(in) :: thermodynamics
    type(cloud_type),         intent(in) :: cloud

    integer, intent(in) :: nlev, istartcol, iendcol

    ! input variables, as defined in radiation_interface.f90

    ! layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each shortwave g-point
    real(jprb), intent(in), dimension(config%n_g_sw,nlev,istartcol:iendcol) :: od_sw, ssa_sw, g_sw

   ! layer optical depth, single scattering albedo and asymmetry factor of
    ! hydrometeors in each shortwave band
    real(jprb), intent(in), dimension(config%n_bands_sw,nlev,istartcol:iendcol)   :: &
         &  od_sw_cloud, ssa_sw_cloud, g_sw_cloud

    ! direct and diffuse surface albedo, and the incoming shortwave
    ! flux into a plane perpendicular to the incoming radiation at
    ! top-of-atmosphere in each of the shortwave g-points
    real(jprb), intent(in), dimension(config%n_g_sw,istartcol:iendcol) &
         &  :: sw_albedo_direct, sw_albedo_diffuse, incoming_sw

    ! layer optical depth, single scattering albedo and asymmetry factor of
    ! gases and aerosols at each longwave g-point, where the latter
    ! two variables are only defined if aerosol longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), intent(in), dimension(config%n_g_lw,nlev,istartcol:iendcol) :: od_lw
    real(jprb), intent(in), dimension(config%n_g_lw_if_scattering,nlev,istartcol:iendcol) :: &
         &  ssa_lw, g_lw

    ! layer optical depth, single scattering albedo and asymmetry factor of
    ! hydrometeors in each longwave band, where the latter two
    ! variables are only defined if hydrometeor longwave scattering is
    ! enabled (otherwise both are treated as zero).
    real(jprb), intent(in), dimension(config%n_bands_lw,nlev,istartcol:iendcol) :: od_lw_cloud
    real(jprb), intent(in), dimension(config%n_bands_lw_if_scattering,nlev,istartcol:iendcol) :: &
         &  ssa_lw_cloud, g_lw_cloud

    ! the planck function (emitted flux from a black body) at half
    ! levels and at the surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw,nlev+1,istartcol:iendcol) :: planck_hl

    ! emission (planck*emissivity) and albedo (1-emissivity) at the
    ! surface at each longwave g-point
    real(jprb), intent(in), dimension(config%n_g_lw, istartcol:iendcol) :: lw_emission, lw_albedo

    ! local variables

    integer :: n_col_local ! number of columns from istartcol to iendcol

    ! object for output netcdf file
    type(netcdf_file) :: out_file

    n_col_local = iendcol + 1 - istartcol

    ! alas the netcdf library is not thread-safe for writing, so we
    ! must write radiative-property files serially

    !$omp critical

    ! open the file
    call out_file%create(trim(file_name), iverbose=config%iverbose)

    ! configure matrix and 3d-array orientation
    call out_file%transpose_matrices(.true.)

    ! sometimes the planck function values are very large or small
    ! even if the fluxes are within a manageable range
    call out_file%double_precision(.true.)

    ! define dimensions
    !    call out_file%define_dimension("column", n_col_local)
    call out_file%define_dimension("column", 0) ! "unlimited" dimension
    call out_file%define_dimension("level", nlev)
    call out_file%define_dimension("half_level", nlev+1)
    if (config%do_clouds) then
      call out_file%define_dimension("level_interface", nlev-1)
    end if

    if (config%do_lw) then
      call out_file%define_dimension("gpoint_lw", config%n_g_lw) 
      if (config%do_clouds) then
        call out_file%define_dimension("band_lw", config%n_bands_lw)
      end if
    end if
    if (config%do_sw) then
      call out_file%define_dimension("gpoint_sw", config%n_g_sw) 
      if (config%do_clouds) then
        call out_file%define_dimension("band_sw", config%n_bands_sw)
      end if
    end if

    ! put global attributes
    call out_file%put_global_attributes( &
         &   title_str="spectral radiative properties from the ecrad offline radiation model", &
         &   references_str="hogan, r. j., and a. bozzo, 2018: a flexible and efficient radiation " &
         &   //"scheme for the ecmwf model. j. adv. modeling earth sys., 10, 1990–2008", &
         &   source_str="ecrad offline radiation model")

    ! define variables
    call out_file%define_variable("pressure_hl", &
         &  dim2_name="column", dim1_name="half_level", &
         &  units_str="pa", long_name="pressure on half-levels")

    if (allocated(thermodynamics%h2o_sat_liq) .and. config%use_aerosols) then
      call out_file%define_variable("q_sat_liquid", &
           &  dim2_name="column", dim1_name="level", &
           &  units_str="kg kg-1", long_name="specific humidity at liquid saturation")
    end if

    if (config%do_sw) then
      call out_file%define_variable("cos_solar_zenith_angle", &
           &  dim1_name="column", units_str="1", &
           &  long_name="cosine of the solar zenith angle")
    end if

    if (config%do_clouds) then
      call out_file%define_variable("cloud_fraction", &
           &  dim2_name="column", dim1_name="level", &
           &  units_str="1", long_name="cloud fraction")
      call out_file%define_variable("overlap_param", &
           &  dim2_name="column", dim1_name="level_interface", &
           &  units_str="1", long_name="cloud overlap parameter")
    end if

    if (config%do_lw) then
      call out_file%define_variable("planck_hl", &
           &  dim3_name="column", dim2_name="half_level", dim1_name="gpoint_lw", &
           &  units_str="w m-2", long_name="planck function on half-levels")
      call out_file%define_variable("lw_emission", &
           &  dim2_name="column", dim1_name="gpoint_lw", &
           &  units_str="w m-2", long_name="longwave surface emission")
      call out_file%define_variable("lw_emissivity", &
           &  dim2_name="column", dim1_name="gpoint_lw", &
           &  units_str="1", long_name="surface longwave emissivity")

      call out_file%define_variable("od_lw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_lw", &
           &  units_str="1", long_name="clear-sky longwave optical depth")
      if (config%do_lw_aerosol_scattering) then
        call out_file%define_variable("ssa_lw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_lw", &
           &  units_str="1", long_name="clear-sky longwave single scattering albedo")
        call out_file%define_variable("asymmetry_lw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_lw", &
           &  units_str="1", long_name="clear-sky longwave asymmetry factor")
      end if

      if (config%do_clouds) then
        call out_file%define_variable("od_lw_cloud", &
             &  dim3_name="column", dim2_name="level", dim1_name="band_lw", &
             &  units_str="1", long_name="in-cloud longwave optical depth")
        if (config%do_lw_cloud_scattering) then
          call out_file%define_variable("ssa_lw_cloud", &
               &  dim3_name="column", dim2_name="level", dim1_name="band_lw", &
               &  units_str="1", long_name="cloud longwave single scattering albedo")
          call out_file%define_variable("asymmetry_lw_cloud", &
               &  dim3_name="column", dim2_name="level", dim1_name="band_lw", &
               &  units_str="1", long_name="cloud longwave asymmetry factor")
        end if
      end if ! do_clouds
    end if ! do_lw
    
    if (config%do_sw) then
      call out_file%define_variable("incoming_sw", &
           &  dim2_name="column", dim1_name="gpoint_sw", &
           &  units_str="w m-2", long_name="incoming shortwave flux at top-of-atmosphere in direction of sun")

      call out_file%define_variable("sw_albedo", &
           &  dim2_name="column", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="surface shortwave albedo to diffuse radiation")
      call out_file%define_variable("sw_albedo_direct", &
           &  dim2_name="column", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="surface shortwave albedo to direct radiation")

      call out_file%define_variable("od_sw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="clear-sky shortwave optical depth")
      call out_file%define_variable("ssa_sw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="clear-sky shortwave single scattering albedo")
      call out_file%define_variable("asymmetry_sw", &
           &  dim3_name="column", dim2_name="level", dim1_name="gpoint_sw", &
           &  units_str="1", long_name="clear-sky shortwave asymmetry factor")

      if (config%do_clouds) then
        call out_file%define_variable("od_sw_cloud", &
             &  dim3_name="column", dim2_name="level", dim1_name="band_sw", &
             &  units_str="1", long_name="in-cloud shortwave optical depth")
        call out_file%define_variable("ssa_sw_cloud", &
             &  dim3_name="column", dim2_name="level", dim1_name="band_sw", &
             &  units_str="1", long_name="cloud shortwave single scattering albedo")
        call out_file%define_variable("asymmetry_sw_cloud", &
             &  dim3_name="column", dim2_name="level", dim1_name="band_sw", &
             &  units_str="1", long_name="cloud shortwave asymmetry factor")
      end if
    end if
   
    if (config%do_clouds) then
      if (allocated(cloud%fractional_std)) then
        call out_file%define_variable("fractional_std", &
             &  dim2_name="column", dim1_name="level", units_str="1", &
             &  long_name="fractional standard deviation of cloud optical depth")
      end if
      if (allocated(cloud%inv_cloud_effective_size)) then
        call out_file%define_variable("inv_cloud_effective_size", &
             &  dim2_name="column", dim1_name="level", units_str="m-1", &
             &  long_name="inverse of cloud effective horizontal size")
      end if
      if (allocated(cloud%inv_inhom_effective_size)) then
        call out_file%define_variable("inv_inhom_effective_size", &
             &  dim2_name="column", dim1_name="level", units_str="m-1", &
             &  long_name="inverse of cloud inhomogeneity effective horizontal size")
      end if
   end if

    ! write variables
    call out_file%put("pressure_hl", thermodynamics%pressure_hl(istartcol:iendcol,:))

    if (allocated(thermodynamics%h2o_sat_liq) .and. config%use_aerosols) then
      call out_file%put("q_sat_liquid", thermodynamics%h2o_sat_liq(istartcol:iendcol,:))
    end if

    if (config%do_clouds) then
      call out_file%put("cloud_fraction", cloud%fraction(istartcol:iendcol,:))
      call out_file%put("overlap_param", cloud%overlap_param(istartcol:iendcol,:))
    end if

    if (config%do_sw) then
      call out_file%put("cos_solar_zenith_angle", single_level%cos_sza(istartcol:iendcol))
      call out_file%put("sw_albedo", sw_albedo_diffuse, do_transp=.false.)
      call out_file%put("sw_albedo_direct", sw_albedo_direct, do_transp=.false.)
    end if

    if (config%do_lw) then
      call out_file%put("lw_emissivity", 1.0_jprb - lw_albedo, do_transp=.false.)
      call out_file%put("planck_hl", planck_hl)
      call out_file%put("lw_emission", lw_emission, do_transp=.false.)
      
      call out_file%put("od_lw", od_lw)
      if (config%do_lw_aerosol_scattering) then
        call out_file%put("ssa_lw", ssa_lw)
        call out_file%put("asymmetry_lw", g_lw)
      end if
      
      if (config%do_clouds) then
        call out_file%put("od_lw_cloud", od_lw_cloud)
        if (config%do_lw_cloud_scattering) then
          call out_file%put("ssa_lw_cloud", ssa_lw_cloud)
          call out_file%put("asymmetry_lw_cloud", g_lw_cloud)
        end if
      end if
    end if
    
    if (config%do_sw) then
      call out_file%put("incoming_sw", incoming_sw, do_transp=.false.)
      
      call out_file%put("od_sw", od_sw)
      call out_file%put("ssa_sw", ssa_sw)
      call out_file%put("asymmetry_sw", g_sw)
      
      if (config%do_clouds) then
        call out_file%put("od_sw_cloud", od_sw_cloud)
        call out_file%put("ssa_sw_cloud", ssa_sw_cloud)
        call out_file%put("asymmetry_sw_cloud", g_sw_cloud)
      end if
    end if
    
    if (config%do_clouds) then
      if (allocated(cloud%fractional_std)) then
        call out_file%put("fractional_std", cloud%fractional_std(istartcol:iendcol,:))
      end if
      if (allocated(cloud%inv_cloud_effective_size)) then
        call out_file%put("inv_cloud_effective_size", cloud%inv_cloud_effective_size(istartcol:iendcol,:))
      end if
      if (allocated(cloud%inv_inhom_effective_size)) then
        call out_file%put("inv_inhom_effective_size", cloud%inv_inhom_effective_size(istartcol:iendcol,:))
      end if
    end if

    ! close the file
    call out_file%close()

    !$omp end critical
    
  end subroutine save_radiative_properties
  

  !---------------------------------------------------------------------
  ! save inputs to the radiation scheme
  subroutine save_inputs(file_name, config, single_level, thermodynamics, &
       &                 gas, cloud, aerosol, lat, lon, iverbose)
    use ecradhook,                  only : lhook, dr_hook, jphook

    use radiation_config,         only : config_type
    use radiation_single_level,   only : single_level_type
    use radiation_thermodynamics, only : thermodynamics_type
    use radiation_gas
    use radiation_cloud,          only : cloud_type
    use radiation_aerosol,        only : aerosol_type
    use easy_netcdf

    character(len=*),             intent(in)   :: file_name
    type(config_type),            intent(in)   :: config
    type(single_level_type),      intent(in)   :: single_level
    type(thermodynamics_type),    intent(in)   :: thermodynamics
    type(gas_type),               intent(inout):: gas
    type(cloud_type),             intent(in)   :: cloud
    type(aerosol_type), optional, intent(in)   :: aerosol
    real(jprb),         optional, intent(in)   :: lat(:), lon(:)
    integer,            optional, intent(in)   :: iverbose

    real(jprb), allocatable :: mixing_ratio(:,:)
    real(jprb), allocatable :: aerosol_mmr(:,:,:)
    real(jprb), allocatable :: seed(:)
    integer       :: i_local_verbose
    integer       :: ncol, nlev
    integer       :: jgas
    character(32) :: var_name, long_name

    ! object for output netcdf file
    type(netcdf_file) :: out_file

    logical :: do_aerosol

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_save:save_inputs',0,hook_handle)

    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

    ! work out array dimensions
    ncol = size(thermodynamics%pressure_hl,1)
    nlev = size(thermodynamics%pressure_hl,2)
    nlev = nlev - 1
    
    do_aerosol = config%use_aerosols .and. present(aerosol)

    ! open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose)

    ! variables stored internally with column varying fastest, but in
    ! output file column varies most slowly so need to transpose
    call out_file%transpose_matrices(.true.)

    ! in the case of aerosols we convert dimensions (ncol,nlev,ntype)
    ! in memory to (nlev,ntype,ncol) in file (in both cases the first
    ! dimension varying fastest).
    call out_file%permute_3d_arrays( (/ 2, 3, 1 /) ) ! for aerosols

    ! define dimensions
    call out_file%define_dimension("column",     ncol)
    call out_file%define_dimension("level",      nlev)
    call out_file%define_dimension("half_level", nlev+1)
    if (allocated(cloud%overlap_param)) then
      call out_file%define_dimension("level_interface", nlev-1)
    end if
    call out_file%define_dimension("sw_albedo_band", &
         &                         size(single_level%sw_albedo,2))
    call out_file%define_dimension("lw_emissivity_band", &
         &                         size(single_level%lw_emissivity,2))

    if (do_aerosol) then
      call out_file%define_dimension("aerosol_type", size(aerosol%mixing_ratio,3))
    end if

    ! put global attributes
    call out_file%put_global_attributes( &
         &   title_str="input profiles to the ecrad offline radiation model", &
         &   references_str="hogan, r. j., and a. bozzo, 2018: a flexible and efficient radiation " &
         &   //"scheme for the ecmwf model. j. adv. modeling earth sys., 10, 1990–2008", &
         &   source_str="ecrad offline radiation model")

    ! define single-level variables
    call out_file%define_variable("solar_irradiance", &
         &   units_str="w m-2", long_name="solar irradiance at earth's orbit")
    if (present(lat)) then
      call out_file%define_variable("lat", &
           &   dim1_name="column", units_str="degrees_north", long_name="latitude")
    end if
    if (present(lon)) then
      call out_file%define_variable("lon", &
           &   dim1_name="column", units_str="degrees_east", long_name="longitude")
    end if
    call out_file%define_variable("skin_temperature", &
         &   dim1_name="column", units_str="k", long_name="skin_temperature")
    if (config%do_sw) then
      call out_file%define_variable("cos_solar_zenith_angle", &
           &   dim1_name="column", units_str="1", &
           &   long_name="cosine of the solar zenith angle")
    end if

    if (allocated(single_level%sw_albedo_direct)) then
      call out_file%define_variable("sw_albedo", &
           &   dim2_name="column", dim1_name="sw_albedo_band", &
           &   units_str="1", long_name="shortwave surface albedo to diffuse radiation")
            call out_file%define_variable("sw_albedo_direct", &
           &   dim2_name="column", dim1_name="sw_albedo_band", &
           &   units_str="1", long_name="shortwave surface albedo to direct radiation")
    else
      call out_file%define_variable("sw_albedo", &
           &   dim2_name="column", dim1_name="sw_albedo_band", &
           &   units_str="1", long_name="shortwave surface albedo")
      
    end if
    call out_file%define_variable("lw_emissivity", &
         &   dim2_name="column", dim1_name="lw_emissivity_band", &
         &   units_str="1", long_name="longwave surface emissivity")

    if (allocated(single_level%iseed)) then
      call out_file%define_variable("iseed", &
           &   dim1_name="column", units_str="1", is_double=.true., &
           &   long_name="seed for random-number generator")
    end if

    ! define thermodynamic variables on half levels
    call out_file%define_variable("pressure_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="pa", long_name="pressure")
    call out_file%define_variable("temperature_hl", &
         &   dim2_name="column", dim1_name="half_level", &
         &   units_str="k", long_name="temperature")

    ! define gas mixing ratios on full levels
    call out_file%define_variable("q", &
         &   dim2_name="column", dim1_name="level", &
         &   units_str="1", long_name="specific humidity")
    call out_file%define_variable("o3_mmr", &
         &   dim2_name="column", dim1_name="level", &
         &   units_str="1", long_name="ozone mass mixing ratio")
    do jgas = 1,nmaxgases
      if (gas%is_present(jgas) .and. jgas /= ih2o .and. jgas /= io3) then
        write(var_name,'(a,a)') trim(gaslowercasename(jgas)), '_vmr'
        write(long_name,'(a,a)') trim(gasname(jgas)), ' volume mixing ratio'
        call out_file%define_variable(trim(var_name), &
             &   dim2_name="column", dim1_name="level", &
             &   units_str="1", long_name=trim(long_name))
      end if
    end do

    if (config%do_clouds) then
      ! define cloud variables on full levels
      call out_file%define_variable("cloud_fraction", &
           &   dim2_name="column", dim1_name="level", &
           &   units_str="1", long_name="cloud fraction")
      call out_file%define_variable("q_liquid", &
           &   dim2_name="column", dim1_name="level", &
           &   units_str="1", long_name="gridbox-mean liquid water mixing ratio")
      call out_file%define_variable("q_ice", &
           &   dim2_name="column", dim1_name="level", &
           &   units_str="1", long_name="gridbox-mean ice water mixing ratio")
      call out_file%define_variable("re_liquid", &
           &   dim2_name="column", dim1_name="level", &
           &   units_str="m", long_name="ice effective radius")
      if (associated(cloud%re_ice)) then
        call out_file%define_variable("re_ice", &
             &   dim2_name="column", dim1_name="level", &
             &   units_str="m", long_name="ice effective radius")
      end if
      if (allocated(cloud%overlap_param)) then
        call out_file%define_variable("overlap_param", &
             &  dim2_name="column", dim1_name="level_interface", &
             &  units_str="1", long_name="cloud overlap parameter")
      end if
      if (allocated(cloud%fractional_std)) then
        call out_file%define_variable("fractional_std", &
             &  dim2_name="column", dim1_name="level", units_str="1", &
             &  long_name="fractional standard deviation of cloud optical depth")
      end if
      if (allocated(cloud%inv_cloud_effective_size)) then
        call out_file%define_variable("inv_cloud_effective_size", &
             &   dim2_name="column", dim1_name="level", units_str="m-1", &
             &   long_name="inverse of cloud effective horizontal size")
      end if
      if (allocated(cloud%inv_inhom_effective_size)) then
        call out_file%define_variable("inv_inhom_effective_size", &
             &  dim2_name="column", dim1_name="level", units_str="m-1", &
             &  long_name="inverse of cloud inhomogeneity effective horizontal size")
      end if
    end if ! do_clouds

    ! define aerosol mass mixing ratio
    if (do_aerosol) then
      call out_file%define_variable("aerosol_mmr", &
           &   dim3_name="column", dim2_name="aerosol_type", &
           &   dim1_name="level", units_str="kg kg-1", &
           &   long_name="aerosol mass mixing ratio")
    end if

    ! write variables
    call out_file%put("solar_irradiance", single_level%solar_irradiance)
    if (present(lat)) then
      call out_file%put("lat", lat)
    end if
    if (present(lon)) then
      call out_file%put("lon", lon)
    end if
    call out_file%put("skin_temperature", single_level%skin_temperature)
    if (config%do_sw) then
      call out_file%put("cos_solar_zenith_angle", single_level%cos_sza)
    end if
    call out_file%put("sw_albedo", single_level%sw_albedo)
    if (allocated(single_level%sw_albedo_direct)) then
      call out_file%put("sw_albedo_direct", single_level%sw_albedo_direct)
    end if
    call out_file%put("lw_emissivity", single_level%lw_emissivity)
    if (config%do_clouds .and. allocated(single_level%iseed)) then
      allocate(seed(ncol))
      seed = single_level%iseed
      call out_file%put("iseed", seed)
      deallocate(seed)
    end if

    call out_file%put("pressure_hl", thermodynamics%pressure_hl)
    call out_file%put("temperature_hl", thermodynamics%temperature_hl)

    allocate(mixing_ratio(ncol,nlev))
    call gas%get(ih2o, imassmixingratio, mixing_ratio)
    call out_file%put("q", mixing_ratio)
    call gas%get(io3, imassmixingratio, mixing_ratio)
    call out_file%put("o3_mmr", mixing_ratio)
    do jgas = 1,nmaxgases
      if (gas%is_present(jgas) .and. jgas /= ih2o .and. jgas /= io3) then
        write(var_name,'(a,a)') trim(gaslowercasename(jgas)), '_vmr'
        call gas%get(jgas, ivolumemixingratio, mixing_ratio)
        call out_file%put(trim(var_name), mixing_ratio)
      end if
    end do
    deallocate(mixing_ratio)

    if (config%do_clouds) then
      call out_file%put("cloud_fraction", cloud%fraction)
      call out_file%put("q_liquid", cloud%q_liq)
      call out_file%put("q_ice", cloud%q_ice)
      call out_file%put("re_liquid", cloud%re_liq)
      if (associated(cloud%re_ice)) then
        call out_file%put("re_ice", cloud%re_ice)
      end if
      if (allocated(cloud%overlap_param)) then
        call out_file%put("overlap_param", cloud%overlap_param)
      end if
      if (allocated(cloud%fractional_std)) then
        call out_file%put("fractional_std", cloud%fractional_std)
      end if
      if (allocated(cloud%inv_cloud_effective_size)) then
        call out_file%put("inv_cloud_effective_size", cloud%inv_cloud_effective_size)
      end if
      if (allocated(cloud%inv_inhom_effective_size)) then
        call out_file%put("inv_inhom_effective_size", cloud%inv_inhom_effective_size)
      end if
    end if

    if (do_aerosol) then
      allocate(aerosol_mmr(ncol, nlev, size(aerosol%mixing_ratio,3)))
      aerosol_mmr = 0.0_jprb
      aerosol_mmr(:,aerosol%istartlev:aerosol%iendlev,:) = aerosol%mixing_ratio
      call out_file%put("aerosol_mmr", aerosol_mmr)
      deallocate(aerosol_mmr)
    end if

    ! close the file
    call out_file%close()

    if (lhook) call dr_hook('radiation_save:save_inputs',1,hook_handle)
    
  end subroutine save_inputs


  !---------------------------------------------------------------------
  ! save shortwave diagnostics computed from "flux" to netcdf
  ! file_name.  the "mapping" matrix maps from fluxes in bands or
  ! g-points to user specified spectral intervals, and should have
  ! been produced by config%get_sw_mapping. see the example in
  ! ecrad_driver.f90.
  subroutine save_sw_diagnostics(file_name, config, wavelength_bound, mapping, flux, &
       &                         iverbose, is_hdf5_file, experiment_name, &
       &                         is_double_precision)

    use ecradhook,                  only : lhook, dr_hook, jphook

    use easy_netcdf

    use radiation_io,             only : nulout
    use radiation_flux,           only : flux_type
    use radiation_config,         only : config_type
    use radiation_matrix,         only : sparse_x_dense

    character(len=*),           intent(in) :: file_name
    type(config_type),          intent(in) :: config
    real(jprb),                 intent(in) :: wavelength_bound(:) ! m
    real(jprb),                 intent(in) :: mapping(:,:)
    type(flux_type),            intent(in) :: flux
    integer,          optional, intent(in) :: iverbose
    logical,          optional, intent(in) :: is_hdf5_file
    logical,          optional, intent(in) :: is_double_precision
    character(len=*), optional, intent(in) :: experiment_name

    integer                                :: nwav ! number of wavelength intervals
    real(jprb), allocatable                :: flux_out(:,:)
    type(netcdf_file)                      :: out_file
    integer                                :: ncol, n_lev_plus1
    integer                                :: i_local_verbose

    real(jphook) :: hook_handle

    if (lhook) call dr_hook('radiation_save:save_sw_diagnostics',0,hook_handle)
    
    if (present(iverbose)) then
      i_local_verbose = iverbose
    else
      i_local_verbose = config%iverbose
    end if

        ! open the file
    call out_file%create(trim(file_name), iverbose=i_local_verbose, is_hdf5_file=is_hdf5_file)

    ! set default precision for file, if specified
    if (present(is_double_precision)) then
      call out_file%double_precision(is_double_precision)
    end if

    ! define dimensions
    ncol = size(flux%sw_up,1)
    call out_file%define_dimension("column", ncol)
    nwav = size(wavelength_bound)-1
    call out_file%define_dimension("wavelength", nwav)
    
    ! put global attributes
    call out_file%put_global_attributes( &
         &   title_str="shortwave spectral diagnostics from the ecrad offline radiation model", &
         &   references_str="hogan, r. j., and a. bozzo, 2018: a flexible and efficient radiation " &
         &   //"scheme for the ecmwf model. j. adv. modeling earth sys., 10, 1990–2008", &
         &   source_str="ecrad offline radiation model")

    ! save "experiment" global attribute if present and not empty
    if (present(experiment_name)) then
      if (experiment_name /= " ") then
        call out_file%put_global_attribute("experiment", experiment_name)
      end if
    end if

    ! define variables
    call out_file%define_variable("wavelength1", dim1_name="wavelength", &
         units_str="m", long_name="wavelength lower bound")
    call out_file%define_variable("wavelength2", dim1_name="wavelength", &
         units_str="m", long_name="wavelength upper bound")

    ! these are always available if
    ! config%do_surface_sw_spectral_flux=true, which is required for
    ! this routine
    call out_file%define_variable("flux_dn_sw_surf", &
         &  dim2_name="column", dim1_name="wavelength", &
         &  units_str="w m-2", long_name="surface downwelling shortwave flux")
    call out_file%define_variable("flux_dn_direct_sw_surf", &
         &  dim2_name="column", dim1_name="wavelength", &
         &  units_str="w m-2", long_name="surface downwelling direct shortwave flux")
    if (config%do_clear) then
      call out_file%define_variable("flux_dn_sw_surf_clear", &
           &  dim2_name="column", dim1_name="wavelength", &
           &  units_str="w m-2", long_name="surface downwelling clear-sky shortwave flux")
      call out_file%define_variable("flux_dn_direct_sw_surf_clear", &
           &  dim2_name="column", dim1_name="wavelength", &
           &  units_str="w m-2", long_name="surface downwelling clear-sky direct shortwave flux")
    end if

    ! the following are only availble if
    ! config%do_save_spectral_flux=true
    if (allocated(flux%sw_up_band)) then
      call out_file%define_variable("flux_up_sw_toa", &
           &  dim2_name="column", dim1_name="wavelength", &
           &  units_str="w m-2", long_name="top-of-atmosphere upwelling shortwave flux")
      call out_file%define_variable("flux_dn_sw_toa", &
           &  dim2_name="column", dim1_name="wavelength", &
           &  units_str="w m-2", long_name="top-of-atmosphere downwelling shortwave flux")
      call out_file%define_variable("flux_up_sw_surf", &
           &  dim2_name="column", dim1_name="wavelength", &
           &  units_str="w m-2", long_name="surface upwelling shortwave flux")
      if (allocated(flux%sw_up_clear_band)) then
        call out_file%define_variable("flux_up_sw_toa_clear", &
             &  dim2_name="column", dim1_name="wavelength", &
             &  units_str="w m-2", long_name="top-of-atmosphere upwelling clear-sky shortwave flux")
        call out_file%define_variable("flux_up_sw_surf_clear", &
             &  dim2_name="column", dim1_name="wavelength", &
             &  units_str="w m-2", long_name="surface upwelling clear-sky shortwave flux")
      end if
    end if

    ! write variables
    call out_file%put("wavelength1", wavelength_bound(1:nwav))
    call out_file%put("wavelength2", wavelength_bound(2:nwav+1))

    n_lev_plus1 = size(flux%sw_up,2)

    ! the mapping matrix is usually sparse, in which case we can check
    ! its elements before multiplying a column of bandwise fluxes by
    ! it.






    
    flux_out = sparse_x_dense(mapping, flux%sw_dn_surf_band)
    call out_file%put("flux_dn_sw_surf", flux_out)
    flux_out = sparse_x_dense(mapping, flux%sw_dn_direct_surf_band)
    call out_file%put("flux_dn_direct_sw_surf", flux_out)
    if (config%do_clear) then
      flux_out = sparse_x_dense(mapping, flux%sw_dn_surf_clear_band)
      call out_file%put("flux_dn_sw_surf_clear", flux_out)
      flux_out = sparse_x_dense(mapping, flux%sw_dn_direct_surf_clear_band)
      call out_file%put("flux_dn_direct_sw_surf_clear", flux_out)
    end if

    if (allocated(flux%sw_up_band)) then
      flux_out = sparse_x_dense(mapping, flux%sw_up_band(:,:,n_lev_plus1))
      call out_file%put("flux_up_sw_surf", flux_out)
      flux_out = sparse_x_dense(mapping, flux%sw_up_band(:,:,1))
      call out_file%put("flux_up_sw_toa", flux_out)
      flux_out = sparse_x_dense(mapping, flux%sw_dn_band(:,:,1))
      call out_file%put("flux_dn_sw_toa", flux_out)
      if (allocated(flux%sw_up_clear_band)) then
        flux_out = sparse_x_dense(mapping, flux%sw_up_clear_band(:,:,1))
        call out_file%put("flux_up_sw_toa_clear", flux_out)
        flux_out = sparse_x_dense(mapping, flux%sw_up_clear_band(:,:,n_lev_plus1))
        call out_file%put("flux_up_sw_surf_clear", flux_out)
      end if
    end if

    call out_file%close()
    
    if (lhook) call dr_hook('radiation_save:save_sw_diagnostics',1,hook_handle)

  end subroutine save_sw_diagnostics
  
end module radiation_save
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
! #define use_sparse_matmul 1
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
! #define my_matmul sparse_x_dense
! #define __sizeof_long_long__ 8
! #define __atomic_acq_rel 4
! #define __atomic_release 3
! #define __version__ "13.3.0"

