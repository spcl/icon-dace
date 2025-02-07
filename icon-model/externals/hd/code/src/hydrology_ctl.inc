! hydrology_ctl.inc - Namelist that sets controls for the HD model
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann and Ha Ho-Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

!  ---------------
!   HYDROLOGY_CTL controls the hydrlology module
!  ---------------
!      
!    NAME           TYPE       PURPOSE                                       
!    ----           ----       -------                                       
!                                                                       
!    ldebughd      LOGICAL     debug HD model  
!    diag_water_budget LOGICAL additional prints to diagnose the global water budget 
!    lbase         LOGICAL     switch for baseflow calculations
!    locean        LOGICAL     closure of water budget for ocean coupling
!    nhd_diag      INTEGER     region number for outflow diagnostic (formerly isolog)
!                                 0   none
!                                 1   Bothnian Bay/Sea
!                                 2   Torneaelven
!                                 4   St.Lawrence
!                                 5   Paraguay
!                                 6   Oder
!                                 7   Elbe
!                                 8   Oranje
!                                 9   Amudarya
!                                10   Lena
!                                99   user defined (FBLOG1, FLLOG1, FBLOG2, FLLOG2)
!    lhd_highres   LOGICAL     swich for outflow diagnostic on HD model grid
!    fllog1        REAL(dp)    user defined grid cells for diagnostics (with nhd_diag=99)
!    fblog1        REAL(dp)      fllog1, fblog1: longitude, latitude of grid cell 1
!    fllog2        REAL(dp)      fllog2, fblog2: longitude, latitude of grid cell 2
!    fblog2        REAL(dp)
!    nremap        INTEGER     Type of Interpolation from input (atmospheric) grid to HD grid
!                                 0   Input = Output
!                                 1   using HDMAP routine by Veronika (default)
!                                 2   0.5 degree to 5 Min.
!    lhd_rout      LOGICAL     switch for original routing or via index arrays
!    fk_rfk        REAL(dp)    Modification factor for k value of riverflow
!    fk_lfk        REAL(dp)    Modification factor for k value of overland flow
!    fk_gfk        REAL(dp)    Modification factor for k value of baseflow
!    irf_vel       INTEGER     Index for type of discharge dependence of riverflow velocity
!                                 0    None
!                                 1    velocity ~ 4th squareroot of Q (riverflow)
!    qrf_ref       REAL(dp)    Reference discharge for discharge dependent riverflow velocity
!    

NAMELIST /HYDROLOGY_CTL/ &
     LDEBUGHD, DIAG_WATER_BUDGET, LBASE, LOCEAN, NHD_DIAG, LHD_HIGHRES, &
     FLLOG1, FBLOG1, FLLOG2, FBLOG2, NREMAP, LHD_ROUT, FK_RFK, FK_LFK, FK_GFK, &
     irf_vel, qrf_ref

!     ------------------------------------------------------------------
