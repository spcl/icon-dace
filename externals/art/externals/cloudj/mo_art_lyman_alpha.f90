!
! mo_art_lyman_alpha
! This module provides the extension for CloudJ to include
! the Lyman-alpha Line for photolysis rate calculation
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

MODULE mo_fjx_lyman

     USE FJX_CMN_MOD
     USE mo_kind,                 ONLY: wp
     
     
 IMPLICIT NONE
 PRIVATE

 PUBLIC ::   photlyman

 CONTAINS

     SUBROUTINE photlyman (tj, dj, sza, valj, lu)

      IMPLICIT NONE


      integer, intent(in) :: lu                              !< given levels

      real(wp),  intent(in), dimension(lu)  :: tj            !< temperature K
      real(wp),  intent(in), dimension(lu)  :: dj            !< air mass
      !real(wp), intent(in) :: solar_lyman
      real(wp), intent(inout), dimension(lu-1,JVN_)::  valj  !< jvalues 1/s
      real(wp), intent(in) :: sza                            !< solar zenith angle degree
      real(wp), parameter ::  solar_lyman = 4.84e15_wp       !< incoming lyman alpha radiation 
                                                             !  in 1 / m2 * s
                                                             !< as preliminary value, 
                                                             !  later is should be read in with time


      real(wp), parameter :: degree =acos(-1d0)/180.0_wp
      real(wp) ::  chi
      real(wp), dimension(L1_+1) :: column, flux, o2_col
      real(wp) ::  t1, t2
      integer :: jk

      real(wp) :: sigma_o2all, phi_o1d, sigma, sigma_h2o, sigma_co2, sigma_ch4

       DO jk=L1_,1,-1
            o2_col(jk) =  dj(jk)*0.20948_wp*1e4_wp
       ENDDO

       ! limit chi < 90 degree
       chi = min(sza, 90d0-1e-5_wp)

       ! slant o2 column

       DO jk=1,L1_-1
           column(jk) = o2_col(jk)/cos(chi*degree)
       ENDDO

       t1 = tj( minloc(abs(column(1:L1_)-1e24_wp),1))
       t2 = tj( minloc(abs(column(1:L1_)-6e24_wp),1))



       DO jk=1,lu-1

        flux(jk) = & 
    &        ( 1 + 0.07_wp * (cos(chi*degree) - 0.5_wp)* exp(-4e-24_wp*column(jk)) ) &
    &        * exp (        &   !* exp(-tau_bar) &
    &        - (2.48e-22_wp + 6e-26_wp*t1)*column(jk)**0.9d0 &
    &        + (2.6e-25_wp - 1e-27_wp*t2)*column(jk) &
    &        - 2.5e-51_wp*column(jk)**2 ) &
    &        + 0.005_wp * (cos(chi*degree) + 0.07_wp) & ! + Q_s
    &        * exp( -9.8e-25_wp * column(jk)  &
    &        * (cos(chi*degree) + 0.02_wp) ) 

          ! multiply with solar flux of Lyman-alpha

           flux(jk)= flux(jk) * solar_lyman

       ENDDO



        do jk = 1, lu-1
! Absorption cross sections from total O2(m2) and quantum yield from
! O2_1d, both from Reddmann and Uhl (2003)

    sigma_o2all = 0.765e-24_wp * &
& ( 1 + 0.35d0 * exp(-5e-25_wp*column(jk)) )     &
&        * ( 1.1_wp + 0.1_wp * tanh(3e-25_wp*column(jk)     &
&        *max(0.8_wp-cos(chi*degree),0d0) - 2.4_wp) )    &
&        * ( 1.16_wp- 0.0021_wp * tj(jk)               &
&        + 6e-6_wp * tj(jk)**2 )                          
    phi_o1d = 0.48_wp                    &! Abschnitt 3.4
&        * ( 1 + 0.2_wp * exp(-3e-25_wp*column(jk)) )      &
&        * ( 1.06_wp + 0.06_wp * tanh(3.5e-25_wp*column(jk) &
&        *max(0.8_wp-cos(chi*degree),0.0_wp) - 2.4_wp) )    

!Absorption cross sections (m2) for all Lyman alpha species

          sigma = sigma_o2all *  phi_o1d
!- not yet implemented
   !CASE
         sigma_h2o = 1.53e-21_wp ! Chabrillat & Kockarts (1997)
!- not yet implemented
         sigma_co2 = 6.2e-24_wp ! vgl. Yoshino et al. (1996)
!- not yet implemented
         sigma_ch4  = 1.9e-21_wp ! Brasseur & Solomon S.159

! Enhancement of photolysiys rates via Lyman-alpha influence
! rewrite jvalues

          valj(jk, 2) = valj(jk,2) + flux(jk)*sigma
          valj(jk, 71) = valj(jk, 71) + flux(jk) *sigma_h2o
          valj(jk, 72) = valj(jk, 72) + flux(jk) *sigma_co2
          valj(jk, 73) = valj(jk, 73) + flux(jk) *sigma_ch4
      


! this is just for temporary output
! write(116,'(i4,6es14.5)') jk, column(jk), & 
!      & flux(jk),valj(jk,1), valj(jk,63), valj(jk,64), valj(jk,65)
          

 end do



      END SUBROUTINE

         
 END MODULE


