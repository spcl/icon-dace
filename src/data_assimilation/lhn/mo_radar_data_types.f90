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

MODULE mo_radar_data_types
  USE mo_kind,               ONLY: wp
  USE mo_var_list,           ONLY: t_var_list_ptr
  USE mtime,                 ONLY: datetime

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: t_radar_fields
  PUBLIC :: t_lhn_diag
  PUBLIC :: t_radar_td_fields
  PUBLIC :: t_radar_ct_fields


TYPE t_radar_td_fields

  REAL(wp), POINTER      &



    & ::               &
    & obs(:,:,:)      ,& ! observations on model grid at six observation time levels
    & spqual(:,:,:)   ,& ! spatial quality function on model grid at two observation time levels
    & radheight(:,:,:)   ! DX radar heights

  TYPE(datetime),POINTER ::  obs_date(:) ! reference date of observations

END TYPE t_radar_td_fields


TYPE t_radar_ct_fields

  REAL(wp), POINTER      &



    & :: blacklist(:,:) ! blacklist for DX radar data

END TYPE t_radar_ct_fields


TYPE t_radar_fields

  TYPE (t_radar_td_fields) :: radar_td
  TYPE (t_var_list_ptr) :: radar_td_list

  TYPE (t_radar_ct_fields) :: radar_ct
  TYPE (t_var_list_ptr) :: radar_ct_list

END TYPE t_radar_fields


TYPE t_lhn_diag
!
  REAL(wp), POINTER      &



    & ::                    &
    & ttend_lhn(:,:,:)     ,& ! temperature increment due to LHN
    & qvtend_lhn(:,:,:)    ,& ! moisture increment due to LHN
    & pr_obs_sum(:,:)      ,& ! cumulated precipitation (hourly)
    & pr_mod_sum(:,:)      ,& ! cumulated precipitation (hourly)
    & pr_ref_sum(:,:)         ! cumulated precipitation (hourly)

  REAL(wp)              &
    & ::                    &
    & ref_bias                ! value of bias correction for reference precipitation

  LOGICAL, POINTER          &



    & ::                    &
    & brightband(:,:)         ! bright band mask field

END TYPE t_lhn_diag

!===============================================================================

END MODULE mo_radar_data_types


