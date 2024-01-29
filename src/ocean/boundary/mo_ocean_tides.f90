! Implementation of tides by computation of the Sun's and Moon's full tidal potential
! This will be used in the pressure gradient calculation
!
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

!----------------------------
#include "icon_definitions.inc"
#include "omp_definitions.inc"
#include "iconfor_dsl_definitions.inc"
!----------------------------
MODULE mo_ocean_tides
  !-------------------------------------------------------------------------
  USE mo_kind,                   ONLY: wp, dp
   USE mtime,                    ONLY: datetime, datetimeToString, deallocateDatetime,              &
       &                               timedelta, newTimedelta, deallocateTimedelta,                &
       &                               MAX_DATETIME_STR_LEN, newDatetime,                           &
       &                               MAX_MTIME_ERROR_STR_LEN, no_error, mtime_strerror,           &
       &                               OPERATOR(-), OPERATOR(+), OPERATOR(>), OPERATOR(*),          &
       &                               ASSIGNMENT(=), OPERATOR(==), OPERATOR(>=), OPERATOR(/=),     &
       &                               event, eventGroup, newEvent,                                 &
       &                               addEventToEventGroup, isCurrentEventActive
  USE mo_exception,              ONLY: message, finish, warning, message_text
  USE mo_model_domain,           ONLY: t_patch, t_patch_3d
  USE mo_ocean_nml,              ONLY: n_zlev, tides_esl_damping_coeff, use_tides,   &
    & tides_mod, use_tides_SAL, tides_SAL_coeff, OceanReferenceDensity_inv, ab_beta,               &
    & tides_smooth_iterations, OceanReferenceDensity
  USE mo_physical_constants,     ONLY: grav
  USE mo_grid_subset,            ONLY: t_subset_range, get_index_range
  USE mo_parallel_config,        ONLY: nproma
  USE mo_impl_constants,         ONLY: sea_boundary
  USE mo_ocean_math_operators,   ONLY: grad_fd_norm_oce_2d_3d
  USE mo_util_dbg_prnt,          ONLY: dbg_print, debug_printValue

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: calculate_tides_potential
  
  CHARACTER(LEN=12)  :: str_module = 'tides'  ! Output of module for 1 line debug
  INTEGER            :: idt_src    = 1               ! Level of detail for 1 line debug
  INTEGER            :: load_iteration = 0
  INTEGER            :: tide_iteration = 0
  !-------------------------------------------------------------------------
  
CONTAINS

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_tides_potential(patch_3d,current_time,rho,h,tides_potential, SelfAtrractionLoad)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    TYPE(datetime), POINTER              :: current_time
    onCells                              :: rho
    onCells_2D                           :: h
    onCells_2D, INTENT(inout)            :: tides_potential
    onCells_2D, INTENT(inout)            :: SelfAtrractionLoad

    IF (use_tides) THEN
      SELECT CASE(tides_mod)
      CASE(1) 
        CALL tide(patch_3d,current_time,tides_potential)
    !  CASE(2) 
    !    CALL tide_mpi(patch_3d,current_time,tides_potential)
      CASE default
        CALL finish(str_module, "Unknown tides_mod")
      END SELECT
   
    ENDIF

    IF (use_tides_SAL) THEN
      CALL calculate_tides_SAL_simple(patch_3d,rho,h,SelfAtrractionLoad)
    ENDIF
  
    
!     CALL dbg_print('tides_potential', tides_potential, str_module, 1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
!     CALL dbg_print('SelfAtrractionLoad',      SelfAtrractionLoad,      str_module, 1, in_subset=patch_3d%p_patch_2d(1)%cells%owned)
        
  END SUBROUTINE calculate_tides_potential
  !-------------------------------------------------------------------------

  !-------------------------------------------------------------------------
  SUBROUTINE calculate_tides_SAL_simple(patch_3d,rho,h,SelfAtrractionLoad)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    onCells                              :: rho
    onCells_2D                           :: h
    onCells_2D, INTENT(inout)            :: SelfAtrractionLoad
    
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: jb,jc,jk,start_index,end_index, levels
    REAL(wp) :: grav_eps_explicit, smooth, barotropic_part, z_grav_rho_inv_eps
    onCells :: cell_thickness
   
    patch_2D        => patch_3d%p_patch_2d(1)
    all_cells       => patch_2D%cells%ALL
   cell_thickness  => patch_3D%p_patch_1d(1)%prism_thick_c
!     cell_thickness  => patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c
  
    IF (load_iteration >= tides_smooth_iterations) THEN 
      smooth = 1.0_wp
    ELSE
      smooth = REAL(load_iteration,wp) / REAL(tides_smooth_iterations,wp)
      load_iteration = load_iteration + 1
    ENDIF
 
    z_grav_rho_inv_eps  = grav * tides_SAL_coeff * smooth
!     grav_eps_explicit = grav * tides_SAL_coeff * (1.0_wp-ab_beta) * smooth

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, jk, start_index,end_index) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      DO jc = start_index, end_index
         SelfAtrractionLoad(jc,jb) = h(jc,jb) * z_grav_rho_inv_eps
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO     
 
  END SUBROUTINE calculate_tides_SAL_simple
  !-------------------------------------------------------------------------
  
  !-------------------------------------------------------------------------
  SUBROUTINE calculate_tides_SAL_Thomas(patch_3d,rho,h,SelfAtrractionLoad)
    TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
    onCells                              :: rho
    onCells_2D                           :: h
    onCells_2D, INTENT(inout)            :: SelfAtrractionLoad
    
    TYPE(t_patch), POINTER :: patch_2D
    TYPE(t_subset_range), POINTER :: all_cells
    INTEGER :: jb,jc,jk,start_index,end_index, levels
    REAL(wp) :: grav_eps_explicit, smooth, barotropic_part, z_grav_rho_inv_eps
    onCells :: cell_thickness
   
    patch_2D        => patch_3d%p_patch_2d(1)
    all_cells       => patch_2D%cells%ALL
   cell_thickness  => patch_3D%p_patch_1d(1)%prism_thick_c
!     cell_thickness  => patch_3D%p_patch_1d(1)%prism_thick_flat_sfc_c
  
    IF (load_iteration >= tides_smooth_iterations) THEN 
      smooth = 1.0_wp
    ELSE
      smooth = REAL(load_iteration,wp) / REAL(tides_smooth_iterations,wp)
      load_iteration = load_iteration + 1
    ENDIF
 
    z_grav_rho_inv_eps  = OceanReferenceDensity_inv * grav * tides_SAL_coeff * smooth
!     grav_eps_explicit = grav * tides_SAL_coeff * (1.0_wp-ab_beta) * smooth

!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, jk, start_index,end_index, levels, &
!ICON_OMP barotropic_part) ICON_OMP_DEFAULT_SCHEDULE
    DO jb = all_cells%start_block, all_cells%end_block
      CALL get_index_range(all_cells, jb, start_index, end_index)
      DO jc = start_index, end_index
         SelfAtrractionLoad(jc,jb) = 0.0_wp
         DO jk = 1, patch_3d%p_patch_1d(1)%dolic_c(jc,jb)
           SelfAtrractionLoad(jc,jb) = SelfAtrractionLoad(jc,jb) + cell_thickness(jc,jk,jb) * rho(jc,jk,jb)
        ENDDO   
        SelfAtrractionLoad(jc,jb) = SelfAtrractionLoad(jc,jb) * z_grav_rho_inv_eps
      ENDDO
    ENDDO
!ICON_OMP_END_PARALLEL_DO     
 
  END SUBROUTINE calculate_tides_SAL_Thomas
  !-------------------------------------------------------------------------
 
  !-------------------------------------------------------------------------
  SUBROUTINE tide(patch_3d,mtime_current,tides_potential)
  
  IMPLICIT NONE
  
  TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
  TYPE(datetime), POINTER              :: mtime_current
  REAL(wp), INTENT(inout)              :: tides_potential(:,:)
  
  TYPE(t_patch), POINTER :: patch_2D
  TYPE(t_subset_range), POINTER :: all_cells
  INTEGER :: jb,jc,jk,start_index,end_index
  CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: datestri
  INTEGER  :: jahr,monat,tag,stunde,minute
  REAL(dp) :: sekunde,gezhochfahr,rs,rm
  REAL(dp) :: dtim,gst,dcl_m,alp_m,h_m,dcl_s,alp_s,h_s,mpot,spot,longitude,latitude
  REAL(wp) :: smooth
  
  
  
  patch_2D        => patch_3d%p_patch_2d(1)
  all_cells       => patch_2D%cells%ALL
  
  
  CALL datetimeToString(mtime_current, datestri)
  
  READ(datestri(1:4),'(i4)') jahr
  READ(datestri(6:7),'(i2.2)') monat
  READ(datestri(9:10),'(i2.2)') tag
  READ(datestri(12:13),'(i2.2)') stunde
  READ(datestri(15:16),'(i2.2)') minute
  READ(datestri(18:23),'(f6.3)') sekunde
  
  IF (tide_iteration >= tides_smooth_iterations) THEN 
    smooth = 1.0_wp
  ELSE
    smooth = REAL(tide_iteration,wp) / REAL(tides_smooth_iterations,wp)
    tide_iteration = tide_iteration + 1
  ENDIF
   
  call timing(jahr,monat,tag,stunde,minute,sekunde,dtim)  ! compute the number of days since 2010 January 0.0
  call get_gst(jahr,monat,tag,stunde,minute,sekunde,gst) ! compute the Greenwich siderial time gst
  call moon_declination(dtim,gst,dcl_m,alp_m,h_m) ! declination (deg) and right ascension (hours) and Greenwich hour angle h (hours) of the Moon after Duffet and Zwart (1979)
  call sun_declination(dtim,gst,dcl_s,alp_s,h_s) ! declination (deg) dcl_s, right ascension (hours) alp_s and Greenwich hour angle h (hours) of the Sun after Duffet and Zwart (1979)
  
  call dist_earth_sun_moon(jahr,monat,tag,stunde+minute*1.0_dp/60.0_dp,rs,rm)  ! distance between rs = Earth - Sun and rm = Earth - Moon  (m)
  
  
   gezhochfahr = 1.0_dp  ! no spin-up 
!  gezhochfahr = min(1.0_dp,(dtim-tide_start)/30.0_dp)  ! spin-up over the first 30 days, when starting at "tide_start"
!  gezhochfahr = min(1.0_dp,(dtim+3286.0_dp)/30.0_dp)  ! spin-up over the first 30 days, when starting at 2001-01-01 00:00
!  gezhochfahr = min(1.0_dp,(dtim+2160.0_dp)/30.0_dp)  ! spin-up over the first 30 days, when starting at 2004-02-01 00:00
  
! Compute tidal potential ==================================================
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, start_index,end_index,longitude,latitude,mpot,spot) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, start_index, end_index)

    DO jc = start_index, end_index

      longitude = patch_2d%cells%center(jc,jb)%lon
      latitude = patch_2d%cells%center(jc,jb)%lat           
      call Moon_pot(dcl_m,h_m,latitude,longitude,rm,mpot) ! computation of the Moon's tidal potential
      call Sun_pot(dcl_s,h_s,latitude,longitude,rs,spot)  ! computation of the Sun's tidal potential

      tides_potential(jc,jb) = tides_esl_damping_coeff * gezhochfahr * (mpot + spot) * smooth ! combined tidal potential 
  
    END DO
  END DO
!ICON_OMP_END_PARALLEL_DO
  
  
  RETURN
END SUBROUTINE tide
  
!-------------------------------------------------------------------------
! find the number of days since 2010 January 0.0
subroutine timing(year,mon,day,hour,minute,sekunde,d)
implicit none

REAL(dp) :: d,d0,sekunde
INTEGER  :: year,mon,day,hour,minute

INTEGER ::  molen(12),ylen(1904:2500),y,m

data molen /31,28,31,30,31,30,31,31,30,31,30,31/



do y = 1904,2500,1
  ylen(y) = 365
enddo
do y = 1904,2500,4
  ylen(y) = 366
  if(year.eq.y) molen(2) = 29
enddo
do y = 2100,2500,100
  ylen(y) = 365
  if(year.eq.y) molen(2) = 28
enddo


d0 = 0.
do y = 1904,2009,1
  d0 = d0 + ylen(y)
enddo

d = 0.
do y = 1904,year-1,1
  d = d + ylen(y)
enddo


do m = 1,mon-1,1
  d = d + molen(m)
enddo

d = d + day
d = d + hour*1.0_dp/24.0_dp
d = d + minute*1.0_dp/(24.0_dp*60)
d = d + sekunde*1.0_dp/(24.0_dp*60*60)

d = d - d0

return
end subroutine timing

! ==========================================================================
! declination (deg) and right ascension (hours) of the Moon after Duffet and Zwart (1979)

subroutine moon_declination(d,gst,dcl,alp,h)
implicit none

REAL(dp) :: d,dcl,alp,gst,h
REAL(dp) :: ma,lama
REAL(dp) :: epsg,omg,e,pi,rad,l,l0,mm,n,p0,n0,ev,ae,a3,c
REAL(dp) :: mms,ec,a4,ls,v,lss,ns,i,lamm,betm,y,x,at,t,eps

pi = 3.141592653589793_dp
rad = pi/180.0_dp


epsg = 279.557208_dp  ! ecliptic longitude at epoch 2010.0
omg = 283.112438_dp   ! ecliptic longitude of perigee at epoch 2010.0
e = 0.016705_dp       ! eccentricity of orbit at epoch 2010.0

ma = 360.0_dp/365.242191_dp*d + epsg - omg

ma = mod(ma,360.0_dp)
if(ma.lt.0.0_dp) ma = ma + 360.0_dp   ! Sun's mean anomaly


lama = 360.0_dp/365.242191_dp*d + 360.0_dp/pi*e*dsin((360.0_dp/365.242191_dp*d + epsg - omg)*rad) + epsg

lama = mod(lama,360.0_dp)
if(lama.lt.0.0_dp) lama = lama + 360.0_dp   ! Sun's ecliptic longitude


! ----------------------------------------------------------------------------------

l0 = 91.929336_dp  ! Moon's mean longitude at the epoch
p0 = 130.143076_dp ! mean longitude of the perigee at the epoch
n0 = 291.682547_dp ! mean longitude of the node at the epoch

l = 13.1763966_dp*d + l0  ! Moon's mean longitude
l = mod(l,360.0_dp)
if(l.lt.0.0_dp) l = l + 360.0_dp

mm = l - 0.1114041_dp*d - p0  ! Moon's mean anomaly
mm = mod(mm,360.0_dp)
if(mm.lt.0.0_dp) mm = mm + 360.0_dp

n = n0 - 0.0529539_dp*d  ! ascending node's mean longitude
n = mod(n,360.0_dp)
if(n.lt.0.0_dp) n = n + 360.0_dp


! --------------------------------------------------------------

C = l-lama

Ev = 1.2739_dp*dsin(rad*(2*C-Mm))  ! corrections for evection

Ae = 0.1858_dp*dsin(rad*Ma) ! annual equation

A3 = 0.37_dp*dsin(rad*Ma)  ! third correction


! -----------------------------------------------------------

Mms = Mm + Ev - Ae -A3   ! Moon's corrected anomaly

Ec = 6.2886_dp*dsin(rad*Mms)  ! equation of the centre

A4 = 0.214_dp*dsin(rad*2*Mms) ! another correction term

ls = l + Ev + Ec - Ae + A4 ! Moon's corrected longitude

V = 0.6583_dp*dsin(rad*2*(ls-lama)) ! variation

lss = ls + v  ! Moon's true orbital longitude


! ---------------------------------------------------------

Ns = N - 0.16_dp*dsin(rad*Ma) ! corrected longitude of the node
i = 5.145396_dp    ! inclination of Moon's orbit

y = dsin(rad*(lss-Ns))*dcos(rad*i)
x = dcos(rad*(lss-Ns)) + 1.e-8_dp


at = datan(y/x)*1.0_dp/rad

call conv(at,x,y)


lamm = at + Ns  ! ecliptic longitude

if(lamm.gt.360.0_dp) lamm = lamm-360.0_dp


betm = dasin( dsin(rad*(lss-Ns))*dsin(rad*i))*1.0_dp/rad  ! ecliptic latitude


! -------------------------------------------------------------

T = (d + 3651.5_dp)/36525.0_dp ! number of Julian centuries since epoch 2000 January 1 12:00

eps = 23.43929166663333336_dp - 0.013004166666666666_dp*T - 1.6666666666666665e-7_dp*T**2 + 5.027777777777778e-7_dp*T**3 ! mean obliquity of the ecliptic

y = dsin(rad*lamm)*dcos(rad*eps) - dtan(rad*betm)*dsin(rad*eps)
x = dcos(rad*lamm) + 1.e-8_dp


alp = datan(y/x)*1.0_dp/rad

call conv(alp,x,y)

alp = alp*0.06666666666666667_dp ! right ascension (hours)

call adj2(alp)

dcl = dasin(dsin(rad*betm)*dcos(rad*eps) + dcos(rad*betm)*dsin(rad*eps)*dsin(rad*lamm))*1.0_dp/rad ! declination

h = gst - alp  ! hour angle (hours)

call adj2(h)


return
end subroutine moon_declination

! ============================================================
subroutine conv(at,x,y)
implicit none

REAL(dp) :: at,x,y

if(at.lt.0.0_dp) at = at + 360.0_dp

if((x.ge.0.0_dp).and.(y.ge.0.0_dp)) then
  if(at.gt.90.0_dp) at = at - 180.0_dp
endif

if((x.lt.0.0_dp).and.(y.ge.0.0_dp)) then
  if(at.gt.270.0_dp) at = at - 180.0_dp
endif

if((x.lt.0.0_dp).and.(y.lt.0.0_dp)) then
  if(at.lt.90.0_dp) at = at + 180.0_dp
endif

if((x.ge.0.0_dp).and.(y.lt.0.0_dp)) then
  if(at.lt.180.0_dp) at = at + 180.0_dp
endif

return
end subroutine conv

! ==========================================================================
! declination (deg) and right ascension (hours) of the Sun after Duffet and Zwart (1979)

subroutine sun_declination(d,gst,dcl,alp,h)
implicit none

REAL(dp) :: d,dcl,alp,gst,h
REAL(dp) :: eps,t,ma,pi,rad,lam,n,ec,x,y,bet,azi

pi = 3.141592653589793_dp
rad = pi/180.0_dp


n = 360.0_dp/365.242191_dp*d  

call adj(n)

Ma = n + 279.557208_dp-283.112438_dp  !mean anomaly

call adj(ma)

ec = 360.0_dp/pi*0.016705_dp*dsin(rad*Ma)

lam = n + ec + 279.557208_dp  ! longitude of the Sun
call adj(lam)


! ===================================================================

bet = 0.0_dp

T = (d + 3651.5_dp)/36525.0_dp ! number of Julian centuries since epoch 2000 January 1 12:00

eps = 23.43929166663333336_dp - 0.013004166666666666_dp*T - 1.6666666666666665e-7_dp*T**2 +5.027777777777778e-7_dp*T**3 ! mean obliquity of the ecliptic

y = dsin(rad*lam)*dcos(rad*eps) - dtan(rad*bet)*dsin(rad*eps)
x = dcos(rad*lam) + 1.e-8_dp


alp = datan(y/x)*1.0_dp/rad

call conv(alp,x,y)

alp = alp*0.06666666666666667_dp ! right ascension (hours)

dcl = dasin(dsin(rad*bet)*dcos(rad*eps) + dcos(rad*bet)*dsin(rad*eps)*dsin(rad*lam))*1.0_dp/rad ! declination

h = gst - alp

call adj2(h)



return
end subroutine sun_declination

! =======================================================================
! Umformung zu postivem Wert zwischen 0 und 360.

subroutine adj(x)
implicit none

REAL(dp) :: x

1 continue
if(x.lt.0.0_dp) then
  x = x + 360.0_dp
  goto 1
endif

2 continue
if(x.gt.360.0_dp) then
  x = x - 360.0_dp
  goto 2
endif

return
end subroutine  adj

! =======================================================================
! Umformung zu positivem Wert zwischen 0. und 24.

subroutine adj2(x)
implicit none

REAL(dp) :: x

1 continue
if(x.lt.0.0_dp) then
  x = x + 24.0_dp
  goto 1
endif

2 continue
if(x.gt.24.0_dp) then
  x = x - 24.0_dp
  goto 2
endif

return
end subroutine adj2

! =========================================================================
subroutine get_gst(year,mon,day,hour,minute,sekunde,gst) ! computation of Greenwich siderial time
implicit none

INTEGER  :: year,mon,day,hour,minute
REAL(dp) :: sekunde,gst,d,s,t0,t,ut

call timing(year,mon,day,0,0,0.0_dp,d)

s = d + 3651.5_dp 
t = s/36525.0_dp
t0 = 6.697374558_dp + 2400.051336_dp*t + 0.000025862_dp*t**2

call adj2(t0)

ut = hour*1.0_dp + minute*1.0_dp/60.0_dp + sekunde*1.0_dp/3600.0_dp 

ut = ut*1.002737909_dp

gst = t0 + ut

call adj2(gst)

return
end subroutine get_gst

! ==========================================================================

subroutine Moon_pot(dcl,h,lat,lon,r,pot) ! computation of the Moon's tidal potential
implicit none
     
REAL(dp) :: pot
REAL(dp) :: dcl,h

REAL(dp) :: pi,rad,gam,m,r,pot0,lat,lon,er,tet,phi,codec,latd,lond
REAL(dp) :: cosg,wrz


pi = 3.141592653589793_dp
rad = pi/180.0_dp


gam = 6.67408e-11_dp ! gravitational constant [m^3/kg/s]
m = 7.349e22_dp ! Moon's mass [kg]
er = 6371.e3_dp ! earth radius [m]

pot0 = gam*m/r

latd = lat/rad
lond = lon/rad
  
tet = (90.0_dp - latd)*rad
codec = (90.0_dp - dcl)*rad
phi = (h*360.0_dp/24.0_dp + lond)*rad

  
cosg = dcos(tet)*dcos(codec) + dsin(tet)*dsin(codec)*dcos(phi)  
wrz = dsqrt(1.0_dp - 2.0_dp*er/r*cosg + (er/r)**2)
  
pot = pot0*(1.0_dp + er/r*cosg - 1.0_dp/wrz)
  

return
end subroutine Moon_pot


! ==========================================================================

subroutine Sun_pot(dcl,h,lat,lon,r,pot) ! computation of the Sun's tidal potential
implicit none
     
REAL(dp) :: pot
REAL(dp) :: dcl,h

REAL(dp) :: pi,rad,gam,m,r,pot0,lat,lon,er,tet,phi,codec,latd,lond
REAL(dp) :: cosg,wrz


pi = 3.141592653589793_dp
rad = pi/180.0_dp


gam = 6.67408e-11_dp ! gravitaional constant [m^3/kg/s]
m = 1.9884e30_dp ! Sun's mass [kg]
er = 6371.e3_dp ! earth radius [m]

pot0 = gam*m/r

latd = lat/rad
lond = lon/rad
  
tet = (90.0_dp - latd)*rad
codec = (90.0_dp - dcl)*rad
phi = (h*360.0_dp/24.0_dp + lond)*rad

  
cosg = dcos(tet)*dcos(codec) + dsin(tet)*dsin(codec)*dcos(phi)  
wrz = dsqrt(1.0_dp - 2.0_dp*er/r*cosg + (er/r)**2)
  
pot = pot0*(1.0_dp + er/r*cosg - 1.0_dp/wrz)
  

return
end subroutine Sun_pot

! ===============================================================================================================
subroutine dist_earth_sun_moon(y,mon,day,h,rs,rm) ! distance between rs = Earth - Sun and rm = Earth - Moon  (m)
implicit none

INTEGER  :: y,mon,day,ys,mons,a,b,c,d
REAL(dp) :: h
REAL(dp) :: jd,t,epsg,omg,e,ms,ee,ny,pi,rad,rs,r0
REAL(dp) :: lams,p0,l0,l,dd,cc,A3,Ae,Ev,Mm,Mms,ec,rm,n
INTEGER  :: ih,im


pi = 4.e0_dp*datan(1.0_dp)
rad = pi/180.0_dp

r0 = 1.495985e11_dp   ! semi-major axis (m)

if(mon.lt.3) then
  ys = y-1
  mons = mon+12
else
  ys = y
  mons = mon
endif

a = dint(0.01_dp*ys)
b = 2 - a + dint(0.25_dp*a)
c = dint(365.25_dp*ys)
d = dint(30.6001_dp*(mons +1))

jd = b + c + d + day + h/24.0_dp + 1720994.5_dp
t = (jd - 2415020.0_dp)/36525.0_dp

epsg = 279.6966778_dp + 36000.76892_dp*t + 0.0003025_dp*t**2 ! ecliptic longitute (degrees)
call adj(epsg)

omg = 281.2208444_dp + 1.719175_dp*t + 0.000452778_dp*t**2 ! ecliptic longitude of perigee (degrees)
call adj(omg)

e = 0.01675104_dp - 0.0000418_dp*t - 0.000000126_dp*t**2 ! eccentricity of orbit

ms = epsg - omg
call adj(ms)

ms = pi/180.0_dp*ms

call findEE(e,ms,ee)

ny = 2.0_dp*datan( dsqrt((1.0_dp+e)/(1.0_dp-e))*dtan(0.5_dp*ee))
ny = 180.0_dp/pi*ny
call adj(ny)

rs = r0 *(1.0_dp - e**2)/(1.0_dp+e*dcos(ny*rad))

! =========================================================================
! distance Earth - Moon (m)

ih = dint(h)
im = dint((h-ih)*1./60 + 0.5)

call timing(y,mon,day,ih,im,0.0_dp,DD)  ! DD = number of days since 2010 January 0.0

a = 384401.e3_dp ! semi-major axis (m)
e = 0.054900_dp ! eccentricity

n = 360.0_dp/365.242191_dp*DD
call adj(n)

ms = ms*180.0_dp/pi  ! Sun's mean anomaly (deg)
call adj(ms)

lams = ny + omg   ! Sun's ecliptic longitude (deg)
call adj(lams)

l0 = 91.929336_dp   ! Moon's mean longitude at the epoch (deg)
p0 = 130.143076_dp  ! mean longitude of the perigee at the epoch (deg)

l = 13.1763966_dp*DD + l0  ! Moon's mean longitude (deg)
call adj(l)

cc = l - lams
call adj(cc)

A3 = 0.37_dp*dsin(ms*rad)
Ae = 0.1858_dp*dsin(ms*rad)
Mm = l - 0.111404_dp*DD-p0  ! Moon's mean anomaly (deg)
call adj(Mm)

Ev = 1.2739_dp*dsin((2.0_dp*cc - mm)*rad) ! evection

Mms = Mm + Ev - Ae - A3  ! Moon's corrected anomaly (deg)

Ec = 6.2886_dp*dsin(Mms*rad) 

rm = a*(1.0_dp - e**2)/(1.0_dp + e*dcos((Mms + Ec)*rad))


return
end subroutine dist_earth_sun_moon

! ================================================

subroutine findee(e,ms,ee)
implicit none

INTEGER  :: it
REAL(dp) :: e,ms,ee,ee1,ee2,incr,pi,v1,v2,fkt

INTEGER :: ie, ie1,ie2,iincr

pi = 4.0_dp*datan(1.0_dp)
it = 0


ie1 = 1
ie2 = 10
iincr = 1

ee1  = 0.0_dp
ee2  = 2.0_dp*pi
incr = 0.2_dp*pi

2 continue

it = it + 1
ee=ee1
! do ee = ee1,ee2,incr
do ie = ie1,ie2,iincr
  ee = ee + incr
  fkt = ee - e*dsin(ee)
  if(fkt.gt.ms) then
    v1 = ee-incr
    v2 = ee
    goto 1
  endif
enddo

1 continue
ee1 = v1
ee2 = v2
incr = 0.1_dp*incr

if(it.lt.8) goto 2

ee = 0.5_dp*(v1+v2)

return
end subroutine findee  
! =============================================================================================
! =============================================================================================

  SUBROUTINE tide_mpi(patch_3d,mtime_current,tides_potential)
  
  IMPLICIT NONE
  
  TYPE(t_patch_3d ),TARGET, INTENT(in) :: patch_3d
  TYPE(datetime), POINTER              :: mtime_current
  REAL(wp), INTENT(inout)              :: tides_potential(:,:)
  
  TYPE(t_patch), POINTER :: patch_2D
  TYPE(t_subset_range), POINTER :: all_cells
  INTEGER :: jb,jc,jk,start_index,end_index
  CHARACTER(LEN=MAX_DATETIME_STR_LEN)    :: datestri
  INTEGER  :: jahr,monat,tag,stunde,minute
  REAL(dp) :: sekunde,gezhochfahr,rs,rm
  REAL(dp) :: dtim,gst,dcl_m,alp_m,h_m,dcl_s,alp_s,h_s,mpot,spot,longitude,latitude
  
  REAL(dp) :: pic,pi2,dt,T,sidt,ecl,nutob,nutl,res(3,2),res2(3,2)
  REAL(dp) :: rkomp,rkosp,erdrad,rekts,dekls,cris3,rektm,deklm,crim3,deklm2,dekls2,sidm,sidmq
  REAL(dp) :: codm,codmq,sids,sidsq,cods,codsq,sidm2,sids2,lon,lat,argp,alatr,hamp,hasp,tipoto
  REAL(dp) :: silato,colato
  INTEGER  :: mmccdt,fnut
  REAL(wp) :: smooth
  
  
  
  patch_2D        => patch_3d%p_patch_2d(1)
  all_cells       => patch_2D%cells%ALL

  IF (tide_iteration >= tides_smooth_iterations) THEN 
    smooth = 1.0_wp
  ELSE
    smooth = REAL(tide_iteration,wp) / REAL(tides_smooth_iterations,wp)
    tide_iteration = tide_iteration + 1
  ENDIF
 
  
  CALL datetimeToString(mtime_current, datestri)
  
  READ(datestri(1:4),'(i4)') jahr
  READ(datestri(6:7),'(i2.2)') monat
  READ(datestri(9:10),'(i2.2)') tag
  READ(datestri(12:13),'(i2.2)') stunde
  READ(datestri(15:16),'(i2.2)') minute
  READ(datestri(18:23),'(f6.3)') sekunde
  
  call timing(jahr,monat,tag,stunde,minute,sekunde,dtim)  ! compute the number of days since 2010 January 0.0
  
  gezhochfahr = 1.0_dp  ! no spin-up 
! gezhochfahr = min(1.0_dp,(dtim-tide_start)/30.0_dp)  ! spin-up over the first 30 days, when starting at "tide_start"
!  gezhochfahr = min(1.0_dp,(dtim+3286.0_dp)/30.0_dp)  ! spin-up over the first 30 days, when starting at 2001-01-01 00:00
!  gezhochfahr = min(1.0_dp,(dtim+2160.0_dp)/30.0_dp)  ! spin-up over the first 30 days, when starting at 2004-02-01 00:00

  t = dtim+3651.5_dp ! Julian days since 2000-01-01 12:00
  t = t/36525.0_dp  !fractional julian centuries t since 2000-01-01 12:00
   
   
  pi2 = dacos(-1.0_dp) * 2.0_dp
  pic = dacos(-1.0_dp)/DBLE(180.0_dp)
   
  call sidt2(pic,pi2,t,sidt) ! corresponding sidereal time Greenwich sidt
  
  
  fnut=0  ! set fnut (perform nutation -> 1; don't -> 0)

  CALL obliq(fnut,pic,T,ecl,nutob,nutl)  ! obliquity of the ecliptic
  
  CALL sun_n(fnut,pic,pi2,T,ecl,nutl,res)!  calculation of position of the Sun according to Duffett, 1990 
  CALL moon(fnut,pic,pi2,T,ecl,nutl,res) !  calculation of position of the Moon according to Duffett, 1990 
  
  CALL aufb2(sidt,res,res2) ! modifications in preparaion of calculation of potentials
  
  rkomp = -4.113e-07_dp ! factor of the tidal potential due to the moon
                     ! attention the factor is defined negative (contrary to the standard).
  rkosp = 0.46051_dp * rkomp ! FIXME: replace with radius from mo_planetary constants
  erdrad = 6371000.0_dp

  rekts=res2(1,1)
  dekls=res2(2,1)
  cris3=res2(3,1)

  rektm=res2(1,2)
  deklm=res2(2,2)
  crim3=res2(3,2)


  deklm2 = deklm * 2.0_dp
  dekls2 = dekls * 2.0_dp
  sidm   = dsin(deklm)
  sidmq  = sidm*sidm
  codm   = dcos(deklm)
  codmq  = codm*codm
  sids   = dsin(dekls)
  sidsq  = sids*sids
  cods   = dcos(dekls)
  codsq  = cods*cods
  sidm2 = dsin(deklm2)
  sids2 = dsin(dekls2)
   
  
! Compute tidal potential ==================================================
!ICON_OMP_PARALLEL_DO PRIVATE(jb, jc, start_index,end_index,longitude,latitude,argp,alatr, silato,colato, &
!ICON_OMP hamp, hasp) ICON_OMP_DEFAULT_SCHEDULE
  DO jb = all_cells%start_block, all_cells%end_block
    CALL get_index_range(all_cells, jb, start_index, end_index)

    DO jc = start_index, end_index

      longitude = patch_2d%cells%center(jc,jb)%lon
      latitude = patch_2d%cells%center(jc,jb)%lat 

      argp = longitude
      alatr = latitude

      silato = sin(latitude)
      colato = cos(latitude)

      hamp = rektm + argp
      hasp = rekts + argp

      ! attention : mpiom uses a negative tidal potential due to negative factor rkomp
      ! FIXME: move 1./3. to parameter
      tides_potential(jc,jb) = smooth * tides_esl_damping_coeff * gezhochfahr*(erdrad * rkomp * crim3 &
         * (3.0_dp * (silato**2 - 1.0_dp/3.0_dp) * (sidmq - 1.0_dp/3.0_dp)&
      &  + DSIN(2.0_dp * alatr) * sidm2 * DCOS(hamp) &
      &  + colato**2 * codmq * DCOS(2.0_dp * hamp))        &
      &  + erdrad * rkosp * cris3 &
      &  * (3.0_dp * (silato**2 - 1.0_dp/3.0_dp) &
      &  * (sidsq - 1.0_dp/3.0_dp) &
      &  + DSIN(2.0_dp * alatr) * sids2 * DCOS(hasp) &
      &  + colato**2 * codsq * DCOS(2.0_dp * hasp)))      

    END DO
  END DO
!ICON_OMP_END_PARALLEL_DO
  
  
  RETURN
END SUBROUTINE tide_mpi

! ========================================================================
!  convertion of Julian Centuries since J2000 T to siderical time sidt
!  according to Duffett, 1990

  SUBROUTINE sidt2(pic,pi2,T,sidt)
! convertion of Julian Centuries since J2000 T to siderical time sidt
! according to Duffett, 1990
!
  REAL(dp) :: pic,pi2,T,sidt,JD,T2,T3
!
  T2=T*T
  T3=T2*T
!
! Julian days
  jd = t * 36525.0_dp + 2451545.0_dp
!
!  mean siderial time of Greenwich in rad
      sidt = (280.46061837_dp + 360.98564736629_dp * (jd - 2451545.0_dp) &
           + 0.000387933_dp * T2 - T3/38710000.0_dp)*pic
!  between 0 and 2pi :
  IF (sidt .LT. 0.0_dp) THEN
    Call negangle2(pi2,sidt)
  ENDIF
  IF(sidt.Ge.pi2) THEN
    CALL langle2(pi2,sidt)
  ENDIF
!
  END SUBROUTINE sidt2
  
! =============================================================
  
  SUBROUTINE negangle2(pi2,x)
! transformation of negative angles
! to angles in interval [0 degree; 360 degree)
!
  Logical endvar
  REAL(dp) :: pi2,x

  endvar=.False.
1 If(.Not.endvar) Then
    x=x+pi2
    endvar = (x .GE. 0.0_dp)
  Goto 1
  Endif
!
  Return
  END  SUBROUTINE  negangle2
      
      
! ===========================================================

  SUBROUTINE langle2(pi2,x)
! transformation of large angles
! to angles in interval [0 degree; 360 degree)
!
  Logical endvar
  REAL(dp) :: pi2,x
!
  endvar=.False.
1 If(.Not.endvar) Then
    x=x-pi2
    endvar=(x.Lt.pi2)
    Goto 1
  Endif
!
  Return
  END SUBROUTINE    langle2
  
! ===============================================================

  Subroutine obliq(fnut,pic,T,ecl,nutob,nutl)
  !  calculation of obliquity of ecliptic
  !  according to Duffett, 1990

  INTEGER :: fnut
  REAL(dp) :: pic,T,ecl,nutob,nutl,A,B,C,T1,T2,T3
  REAL(dp) :: L1,L2,D1,D2,M1,M2,N1,N2
!
!  see page 57
! correction terms if requested
  If(fnut.Eq.1) Then
    t1 = t + 1.0_dp     ! adding one century to "shift" reference time to J1900
    T2=T1*T1
!
    a = 100.0021358_dp * t1
    ! FIXME: consider using intrinsic AINT here
    b = 360.0_dp * (a - DBLE(DINT(a)))
    l1 = 279.6967_dp + 0.000303_dp * t2 + b
    l2 = 2.0_dp * l1 * pic
    a = 1336.855231_dp * t1
    ! FIXME: consider using intrinsic AINT here
    b = 360.0_dp * (a - DBLE(DINT(a)))
    d1 = 270.4342_dp - 0.001133_dp * t2 + b
    d2 = 2.0_dp * d1 * pic
    a = 99.99736056_dp * t1
    ! FIXME: consider using intrinsic AINT here
    b = 360.0_dp * (a - DBLE(DINT(a)))
    m1 = (358.4758_dp - 0.00015_dp * t2 + b) * pic
    a = 1325.552359_dp * t1
    b = 360.0_dp * (a - DBLE(DINT(a)))
    m2 = (296.1046_dp + 0.009192_dp * t2 + b) * pic
    a = 5.372616667_dp * t1
    b = 360.0_dp * (a - DBLE(DINT(a)))
    n1 = (259.1833_dp + 0.002078_dp * t2 - b) * pic
    n2 = 2.0_dp * n1
!   correction term for nutation in longitude
    nutl = ((-17.2327_dp - 0.01737_dp * t1) * DSIN(n1)               &
         + (-1.2729_dp - 0.00013_dp * t1) * DSIN(l2) + 0.2088_dp * DSIN(n2) &
         - 0.2037_dp * DSIN(d2) + (0.1261_dp - 0.00031_dp * t1) * DSIN(m1)  &
         + 0.0675_dp * DSIN(m2) - (0.0497_dp - 0.00012_dp * t1) * DSIN(l2 + m1) &
         - 0.0342_dp * DSIN(d2 - n1) - 0.0261_dp * DSIN(d2 + m2)      &
         + 0.0214_dp * DSIN(l2 - m1) - 0.0149_dp * DSIN(l2 - d2 + m2) &
         + 0.0124_dp * DSIN(l2 - n1) + 0.0114_dp * DSIN(d2 - m2)) &
         / 3600.0_dp * pic
! correction term for nutation in obliquity of the ecliptic
    nutob = ((9.21_dp + 0.00091_dp * t1) * DCOS(n1)                  &
             + (0.5522_dp - 0.00029_dp * t1) * DCOS(l2) - 0.0904_dp * COS(n2) &
             + 0.0884_dp * DCOS(d2) + 0.0216_dp * DCOS(l2 + m1)           &
             + 0.0183_dp * DCOS(d2 - n1) + 0.0113_dp * DCOS(d2 + m2)      &
             - 0.0093_dp * DCOS(l2 - m1) - 0.0066_dp * DCOS(l2 - n1)) &
             / 3600.0_dp * pic
!
  Else
    nutob = 0.0_dp
    nutl = 0.0_dp
  Endif
   ! obliquity of the ecliptic
   ! FIXME : adding one century to "shift" reference time to J1900 was missing
   !t1 = t + 1._wp     ! adding one century to "shift" reference time to J1900
   !T2=T1*T1
   !T3=T2*T1
   !c = 46.815_wp * t1 + 0.0006_wp * t2 - 0.00181_wp * t3
   !ecl = (23.43929167_wp - c/3600.0_dp) * pic + nutob

    T2=T*T
    T3=T2*T
    c = 46.815_dp * t + 0.0006_dp * t2 - 0.00181_dp * t3
    ecl = (23.43929167_dp - c/3600.0_dp) * pic + nutob
!
  END SUBROUTINE obliq
  
! ==================================================================

      SUBROUTINE Sun_n(fnut,pic,pi2,T,ecl,nutl,res)
!  calculation of position of the Sun
!  according to Duffett, 1990
!
      Integer fnut
      REAL(dp) :: pic,pi2,T,T1,T2,T3,ecl,A,B,nutl
      REAL(dp) :: L,M1,EC,AT,AE
      REAL(dp) :: A1,B1,C1,D1,E1,H1,D2,D3,S1,S2,S3,SW,X1,X2,res(3,2)
!
!  see page 116
      t1 = t + 1.0_dp ! adding one century to "shift" reference time to J1900
      T2=T1*T1
      T3=T2*T1
!
      a = 100.0021359_dp * t1
      ! FIXME: consider using intrinsic AINT here
      b = 360.0_dp * (a - DBLE(DINT(a)))
      l = (279.69668_dp + 0.0003025_dp * t2 + b) * pic
      a = 99.99736042_dp * t1
      ! FIXME: consider using intrinsic AINT here
      b = 360.0_dp * (a - DBLE(DINT(a)))
      m1 = (358.47583_dp - 0.00015_dp * t2 + 0.0000033_dp * t3 + b) * pic
      ec = 0.01675104_dp - 0.0000418_dp * t1 - 0.000000126_dp * t2
!
!  true and eccentric anomaly in rad
      Call anomaly(pi2,M1,EC,AT,AE)
!
!  various arguments in rad
      a = 62.55209472_dp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360.0_dp * (a - DBLE(DINT(a)))
      a1 = (153.23_dp + b) * pic
      a = 125.1041894_dp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360.0_dp * (a - DBLE(DINT(a)))
      b1 = (216.57_dp + b) * pic
      a = 91.56766028_dp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360.0_dp * (a - DBLE(DINT(a)))
      c1 = (312.69_dp + b) * pic
      a = 1236.853095_dp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360.0_dp * (a - DBLE(DINT(a)))
      d1 = (350.74_dp - 0.00144_dp * t2 + b) * pic
      e1 = (231.19_dp + 20.2_dp * t1) * pic
      a = 183.1353208_dp * t1
      !FIXME: consider using intrinsic AINT here
      b = 360.0_dp * (a - DBLE(DINT(a)))
      h1 = (353.4_dp + b) * pic
!
      d2 = (0.00134_dp * DCOS(a1) + 0.00154_dp * DCOS(b1) + 0.002_dp * DCOS(c1) &
           + 0.00179_dp * DSIN(d1) + 0.00178_dp * DSIN(e1)) * pic
      d3 = 0.00000543_dp * DSIN(a1) + 0.00001575_dp * DSIN(b1) &
           + 0.00001627_dp * DSIN(c1) + 0.00003076_dp * DCOS(d1) &
           + 0.00000927_dp * DSIN(h1)
!
!  geocentric ecliptic coordinates of the Sun
      S1=AT+L-M1+D2
      If(fnut.Eq.1) Then
        S1=S1+nutl
      Endif
      IF (s1 .LT. 0.0_dp) CALL negangle2(pi2,S1)
      If(S1.Ge.pi2) Call langle2(pi2,S1)
      s2 = 0.0_dp
      s3 = 1.0000002_dp * (1.0_dp - ec * COS(ae)) + d3
!
!  geocentric equatorial coordinates of the Sun
      SW  =  -1.0_dp
      Call eqecl(pi2,S1,S2,X1,X2,ecl,SW)
      res(1,1)=X1
      res(2,1)=X2
      res(3,1)=S3
!
      End Subroutine sun_n
      
! ===============================================================

  SUBROUTINE anomaly(pi2,AM,EC,AT,AE)
!  calculation of true anomaly AT and eccentric anomaly AE
!  given mean anomaly AM and eccentricity EC
!  according to Duffett, 1990
!
   REAL(dp) :: pi2,AM,EC,AT,AE,M,D,A
!
!  see page 113
      ! FIXME: consider using intrinsic AINT here
      m = am - pi2 * DBLE(DINT(am / pi2))
      AE=M
  1   D=AE-(EC*DSin(AE))-M
      IF (ABS(d) .GE. 0.000006_dp) THEN
        d = d/(1.0_dp - ec * DCOS(ae))
        AE=AE-D
        Goto 1
      Else
        a = DSQRT((1.0_dp + ec)/(1.0_dp - ec))*DTAN(ae/2.0_dp)
        at = 2.0_dp * DATAN(a)
      Endif
!
      END SUBROUTINE anomaly
      
!  =============================================================

  SUBROUTINE eqecl(pi2,X,Y,P,Q,ecl,SW)
!  conversion of ecliptic into equatorial coordinates
!  according to Duffett, 1990
!
! if SW=+1: equatorial (X,Y..alpha,delta) to ecliptic (P,Q..lambda,beta)
! if SW=-1: equatorial (X,Y..lambda,beta) to ecliptic (P,Q..alpha,delta)
!
      REAL(dp) :: pi2,ecl,P,Q,X,Y,SW
!
!  see page 62
      p = DATAN2((DSIN(x) * DCOS(ecl) + DTAN(y) * DSIN(ecl) * sw), DCOS(X))
      IF (p .LT. 0.) CALL negangle2(pi2,p)
      IF (P .GE. pi2) CALL langle2(pi2,P)
      Q=DAsin(DSin(Y)*DCos(ecl)-DCos(Y)*DSin(ecl)*DSin(X)*SW)
!
      END SUBROUTINE eqecl
      
! ================================================================
  Subroutine moon(fnut,pic,pi2,T,ecl,nutl,res)
  !  calculation of position of the Moon
  !  according to Duffett, 1990 
  !
  Integer fnut
  REAL(dp) :: pic,pi2,T,T1,T2,T3,ecl,nutl,A,B,C,SW
  REAL(dp) :: Q,M1,M2,M3,M4,M5,M6,ML,MS,MD,ME,MF,NA,S1,S2,S3,S4,E,E2
  REAL(dp) :: L,G,W1,W2,PM,MO1,MO2,MO3,X1,X2,res(3,2)
!
! see page 157
  t1 = t + 1.0_dp ! adding one century to "shift" reference time to J1900
  T2=T1*T1
  T3=T2*T1
!
  q = t1 * 36525.0_dp
  m1 = q/27.32158213_dp
  !FIXME: consider using intrinsic AINT here
  m1 = 360.0_dp * (m1 - DBLE(DINT(m1)))
  m2 = q / 365.2596407_dp
  !FIXME: consider using intrinsic AINT here
  m2 = 360.0_dp * (m2 - DBLE(DINT(m2)))
  m3 = q / 27.55455094_dp
  !FIXME: consider using intrinsic AINT here
  m3 = 360.0_dp * (m3 - DBLE(DINT(m3)))
  m4 = q / 29.53058868_dp
  !FIXME: consider using intrinsic AINT here
  m4 = 360.0_dp * (m4 - DBLE(DINT(m4)))
  m5 = q / 27.21222039_dp
  !FIXME: consider using intrinsic AINT here
  m5 = 360.0_dp * (m5 - DBLE(DINT(m5)))
  m6 = q / 6798.363307_dp
  !FIXME: consider using intrinsic AINT here
  m6 = 360.0_dp * (m6 - DBLE(DINT(m6)))
  ml = 270.434164_dp + m1 - 0.001133_dp * t2 + 0.0000019_dp * t3
  ms = 358.475833_dp + m2 - 0.00015_dp * t2 + 0.0000033_dp * t3
  md = 296.104608_dp + m3 + 0.009192_dp * t2 + 0.0000144_dp * t3
  me = 350.737486_dp + m4 - 0.001436_dp * t2 + 0.0000019_dp * t3
  mf = 11.250889_dp + m5 - 0.003211_dp * t2 - 0.0000003_dp * t3
  na = (259.183275_dp - m6 + 0.002078_dp * t2 + 0.0000022_dp * t3) * pic
  s2 = DSIN(na)
  a = (51.2_dp + 20.2_dp * t1) * pic
  s1 = DSIN(a)
  b = (346.56_dp + 132.87_dp * t1 - 0.0091731_dp * t2) * pic
  s3 = 0.003964_dp * DSIN(b)
  c = na + (275.05_dp - 2.3_dp * t1) * pic
  s4 = DSIN(c)
  ml = (ml + 0.000233_dp * s1 + s3 + 0.001964_dp * s2) * pic
  ms = (ms - 0.001778_dp * s1) * pic
  md = (md + 0.000817_dp * s1 + s3 + 0.002541_dp * s2) * pic
  mf = (mf + s3 - 0.024691_dp * s2 - 0.004328_dp * s4) * pic
  me = (me + 0.002011_dp * s1 + s3 + 0.001964_dp * s2) * pic
  e = 1.0_dp - 0.002495_dp * t1 + 0.00000752_dp * t2
  e2 = e * e

!  ecliptic longitude MO1
  l = 6.28875_dp * DSIN(md) + 1.274018_dp * DSIN(2.0_dp * me - md)        &
           + 0.658309_dp * DSIN(2.0_dp * me) + 0.213616_dp * DSIN(2.0_dp * md) &
           - e * 0.185596_dp * DSIN(ms) - 0.114336_dp * DSIN(2.0_dp * mf)     &
           + 0.058793_dp * DSIN(2.0_dp * (me - md))                          &
           + 0.057212_dp * e * DSIN(2.0_dp * me - ms - md) + 0.05332_dp      &
           & * DSIN(2.0_dp * me + md)                                        &
           + 0.045874_dp * e * DSIN(2.0_dp * me - ms)                        &
           + 0.041024_dp * e * DSIN(md - ms)                                &
           - 0.034718_dp * DSIN(me) - e * 0.030465_dp * DSIN(md + ms)        &
           + 0.015326_dp * DSIN(2.0_dp * (me - mf))                          &
           - 0.012528_dp * DSIN(2.0_dp * mf + md)                            &
           - 0.01098_dp * DSIN(2.0_dp * mf - md)                             &
           + 0.010674_dp * DSIN(4.0_dp * me - md)                            &
           + 0.010034_dp * DSIN(3.0_dp * md)                                 &
           + 0.008548_dp * DSIN(4.0_dp * me - 2.0_dp * md)                    &
           - e * 0.00791_dp * DSIN(ms - md + 2.0_dp * me)                    &
           - e * 0.006783_dp * DSIN(2.0_dp * me + ms)                        &
           + 0.005162_dp * DSIN(md - me) + e * 0.005_dp * DSIN(me + ms)      &
           + 0.003862_dp * DSIN(4.0_dp * me)                                 &
           + e * 0.004049_dp * DSIN(md - ms + 2.0_dp * me)                   &
           + 0.003996_dp * DSIN(2.0_dp * (md + me))                          &
           + 0.003665_dp * DSIN(2.0_dp * me - 3.0_dp * md)                    &
           + e * 0.002695_dp * DSIN(2.0_dp * md - ms)                        &
           + 0.002602_dp * DSIN(md - 2.0_dp * (mf + me))                     &
           + e * 0.002396_dp * DSIN(2.0_dp * (me - md) - ms)                 &
           - 0.002349_dp * DSIN(me + md)                                    &
           + e2 * 0.002249_dp * DSIN(2.0_dp * (me - ms))                     &
           - e * 0.002125_dp * DSIN(ms + 2.0_dp * md)                        &
           - e2 * 0.002079_dp * DSIN(2.0_dp * ms)                            &
           + e2 * 0.002059_dp * DSIN(2.0_dp * (me - ms) - md)                &
           - 0.001773_dp * DSIN(2.0_dp * (me - mf) + md)                     &
           - 0.001595_dp * DSIN(2.0_dp * (me + mf))                          &
           + e * 0.00122_dp * DSIN(4.0_dp * me - ms - md)                    &
           - 0.00111_dp * DSIN(2.0_dp * (md + mf))                           &
           + 0.000892_dp * DSIN(md - 3.0_dp * me)                            &
           - e * 0.000811_dp * DSIN(ms + md + 2.0_dp * me)                   &
           + e * 0.000761_dp * DSIN(4.0_dp * me - ms - 2.0_dp * md)           &
           + e2 * 0.000704_dp * DSIN(md - 2.0_dp * (ms + me))                &
           + e * 0.000693_dp * DSIN(ms - 2.0_dp * (md - me))                 &
           + e * 0.000598_dp * DSIN(2.0_dp * (me - mf) - ms)                 &
           + 0.00055_dp * DSIN(md + 4.0_dp * me)                             &
           + 0.000538_dp * DSIN(4.0_dp * md)                                 &
           + e * 0.000521_dp * DSIN(4.0_dp * me - ms)                        &
           + 0.000486_dp * DSIN(2.0_dp * md - me)                            &
           + e2 * 0.000717_dp * DSIN(md - 2.0_dp * ms)
      MO1=ML+L*pic
      If(fnut.Eq.1) Then
        MO1=MO1+nutl
      Endif
      IF (MO1 .LT. 0.0_dp) CALL negangle2(pi2,MO1)
      IF(MO1 .GE. pi2) CALL langle2(pi2,MO1)

!  ecliptic latitude MO2
      g = 5.128189_dp * DSIN(mf) + 0.280606_dp * DSIN(md + mf)                 &
           + 0.277693_dp * DSIN(md - mf) + 0.173238_dp * DSIN(2.0_dp * me - mf) &
           + 0.055413_dp * DSIN(2.0_dp * me + mf - md)                         &
           + 0.046272_dp * DSIN(2.0_dp * me - mf - md)                         &
           + 0.032573_dp * DSIN(2.0_dp * me + mf)                              &
           + 0.017198_dp * DSIN(2.0_dp * md + mf)                              &
           + 0.009267_dp * DSIN(2.0_dp * me - mf + md)                         &
           + 0.008823_dp * DSIN(2.0_dp * md - mf)                              &
           + e * 0.008247_dp * DSIN(2.0_dp * me - ms - mf)                     &
           + 0.004323_dp * DSIN(2.0_dp * (me + md) - mf)                       &
           + 0.0042_dp * DSIN(2.0_dp * me + md + mf)                           &
           + e * 0.003372_dp * DSIN(mf - ms - 2.0_dp * me)                     &
           + e * 0.002472_dp * DSIN(2.0_dp * me - md + mf - ms)                &
           + e * 0.002222_dp * DSIN(2.0_dp * me + mf - ms)                     &
           + e * 0.002072_dp * DSIN(2.0_dp * me - md - mf - ms)                &
           + e * 0.001877_dp * DSIN(mf - ms + md)                             &
           + 0.001828_dp * DSIN(4.0_dp * me - md - mf)                         &
           - e * 0.001803_dp * DSIN(ms + mf) - 0.00175_dp * DSIN(3.0_dp * mf)   &
           + e * 0.00157_dp * DSIN(md - mf - ms) - 0.001487_dp * DSIN(me + mf) &
           - e * 0.001481_dp * DSIN(mf + ms + md)                             &
           + e * 0.001417_dp * DSIN(mf - ms - md)                             &
           + e * 0.00135_dp * DSIN(mf - ms) + 0.00133_dp * DSIN(mf - me)       &
           + 0.001106_dp * DSIN(mf + 3.0_dp * md)                              &
           + 0.00102_dp * DSIN(4.0_dp * me - mf)                               &
           + 0.000833_dp * DSIN(mf + 4.0_dp * me - md)                         &
           + 0.000781_dp * DSIN(md - 3.0_dp * mf)                              &
           + 0.00067_dp * DSIN(mf + 3.0_dp * me - 2.0_dp * md)                  &
           + 0.000606_dp * DSIN(2.0_dp * me - 3.0_dp * mf)                      &
           + 0.000597_dp * DSIN(2.0_dp * (me + md) - mf)                       &
           + e * 0.000492_dp * DSIN(2.0_dp * me + md - ms - mf)                &
           + 0.00045_dp * DSIN(2.0_dp * (md - me) - mf)                        &
           + 0.000439_dp * DSIN(3.0_dp * me - mf)                              &
           + 0.000423_dp * DSIN(mf + 2.0_dp * (me + md))                       &
           + 0.000422_dp * DSIN(2.0_dp * me - 3.0_dp * md - mf)                 &
           - e * 0.000367_dp * DSIN(mf + ms + 2.0_dp * me - md)                &
           - e * 0.000353_dp * DSIN(mf + ms + 2.0_dp * me)                     &
           + 0.000331_dp * DSIN(mf + 4.0_dp * me)                              &
           + e * 0.000317_dp * DSIN(2.0_dp * me + md - ms + mf)                &
           + e2 * 0.000306_dp * DSIN(2.0_dp * (me - ms) - mf)                  &
           - 0.000283_dp *DSIN(md + 3.0_dp * mf)
      w1 = 0.0004664_dp * DCOS(na)
      w2 = 0.0000754_dp * DCOS(c)
      mo2 = g * pic * (1.0_dp - w1 - w2)

!  horizontal parallax PM
      pm = 0.950724_dp + 0.051818_dp * DCOS(md)                             &
           + 0.009531_dp * DCOS(2.0_dp * me - md)                            &
           + 0.007843_dp * DCOS(2.0_dp * me) + 0.002824_dp * DCOS(2.0_dp * md) &
           + 0.000857_dp * DCOS(2.0_dp * me + md)                            &
           + e * 0.000533_dp * DCOS(2.0_dp * me - ms)                        &
           + e * 0.000401_dp * DCOS(2.0_dp * me - md - ms)                   &
           + e * 0.00032_dp * DCOS(md - ms) - 0.000271_dp * DCOS(me)         &
           - e * 0.000264_dp * DCOS(md + ms)                                &
           - 0.000198_dp * DCOS(2.0_dp * mf - md)                            &
           + 0.000173_dp * DCOS(3.0_dp * md)                                 &
           + 0.000167_dp * DCOS(4.0_dp * me - md)- e * 0.000111_dp * DCOS(ms) &
           + 0.000103_dp * DCOS(4.0_dp * me - 2.0_dp * md)                    &
           - 0.000084_dp * DCOS(2.0_dp * md - 2.0_dp * me)                    &
           - e * 0.000083_dp * DCOS(2.0_dp * me + ms)                        &
           + 0.000079_dp * DCOS(2.0_dp * me + 2.0_dp * md)                    &
           + 0.000072_dp * DCOS(4.0_dp * me)                                 &
           + e * 0.000064_dp * DCOS(2.0_dp * me - ms + md)                   &
           - e * 0.000063_dp * DCOS(2.0_dp * me + ms - md)                   &
           + e * 0.000041_dp * DCOS(ms + me)                                &
           + e * 0.000035_dp * DCOS(2.0_dp * md - ms)                        &
           - 0.000033_dp * DCOS(3.0_dp * md - 2.0_dp * me)                    &
           - 0.00003_dp * DCOS(md + me)                                     &
           - 0.000029_dp * DCOS(2.0_dp * (mf - me))                          &
           - e * 0.000029_dp * DCOS(2.0_dp * md + ms)                        &
           + e2 * 0.000026_dp * DCOS(2.0_dp * (me - ms))                     &
           - 0.000023_dp * DCOS(2.0_dp * (mf - me) + md)                     &
           + e * 0.000019_dp * DCOS(4.0_dp * me - md - ms)
      PM=PM*pic

!  geocentric distance MO3 in km
      ! FIXME: how is this related to radius in mo_planetary_constants?
      mo3 = 6378.14_dp/DSIN(pm)

!  geocentric equatorial coordinates of the Moon
      sw = -1.0_dp
      Call eqecl(pi2,MO1,MO2,X1,X2,ecl,SW)
      res(1,2)=X1
      res(2,2)=X2
      res(3,2)=MO3
!
      END SUBROUTINE moon
      
! ================================================================

  SUBROUTINE aufb2(sidt,res,res2)

  !  modifications according to "ephaufb.f" by Maik Thomas
  ! (rekt(rad)->sid.time.green.-r.asc.; dekl(rad); cri3->(a/r)^3
  ! for Sun and Moon) 
  !
  REAL(dp) :: sidt,h(3)
  REAL(dp) :: res(3,2),res2(3,2)
!
  res2(1,1)=sidt-res(1,1)
  res2(1,2)=sidt-res(1,2)
  res2(2,1)=res(2,1)
  res2(2,2)=res(2,2)
  h(1) = 1.0_dp / res(3,1)
  h(2) = 384400.0_dp / res(3,2)
  res2(3,1)=h(1)*h(1)*h(1)
  res2(3,2)=h(2)*h(2)*h(2)
!
  Return
  END SUBROUTINE aufb2
  

END MODULE mo_ocean_tides
!=============================================================================
