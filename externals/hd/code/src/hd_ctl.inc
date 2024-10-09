! hd_ctl.inc - Namelist for controls of HD in stand alone runs.
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: Apache-2.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann and Ha Ho-Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

!  -----------
!  HD_CTL: Controls of HD if used as a separate model, e.g. in stand alone runs.
!  -----------
!
!  In echam runs, some of these parameters are set in RUNCTL. Original name: HDALONE_CTL
!                                                                 
!      
!    NAME           TYPE       PURPOSE                                       
!    ----           ----       -------                                       
!                                                                       
!    out_expname   CHARACTER   experiment name
!    out_datapath  CHARACTER   path to where the output data shall be written
!    year1         INTEGER     initial year of the run
!    month1        INTEGER     initial month of the run
!    date_start    CHARACTER   start date of the run, format YYYYMMDD or YYYY-MM-DD
!    date_end      CHARACTER   end date of the run, format YYYYMMDD or YYYY-MM-DD
!    nstep         INTEGER     number of time steps within the run if date_start & date_end are not provided
!    delta_time    REAL(dp)    model time step lenght in seconds
!    ufakru        REAL(dp)    unit factor for runoff and drainage input data
!    runoff_file   CHARACTER   file with input runoff data
!    drainage_file CHARACTER   file with input drainage data
!    forcing_freq  INTEGER     frequency of the forcing data (0: stepwise, 1: daily)
!    iout          INTEGER     averaging period of the HD output
!                                 1   30-day averages
!                                 2   decadal averages
!                                 3   weekly averages
!                                 4   monthly averages without leap years
!                                 5   monthly averages with leap years
!                                 6   daily averages
!                                 7   stepwise output
!    coupling_type  INTEGER     Coupling switch that also indicates type of coupling to ocean model
!                                 0   no coupling (default)
!                                 1   coupling with no interpolation of discharge within HD
!                                 2   coupling with allocation of target grid boxes in HD
!    lcoupling_atm  LOGICAL     Switch for coupling to atmosphere (default: .False.) 
!    lcoupling_oce  LOGICAL     Switch for coupling to ocean (default: .False.)
!    icpl_sinks     INTEGER     Redistribution of water in sinks for ocean coupling 
!                                 0=None (Def.), 1=Relatively to mouth boxes, 2= equally to all ocean boxes
!    icpl_mask_tohd INTEGER     Switch for the mask on which HD may receive input for atmosphere coupling 
!                                 0=HD land boxes (Def.), 1=all HD boxes
!    coupling_file  CHARACTER   Input file with coupling information for coupling_type 2
!    lcoupling_out  LOGICAL     Write discharge on ocean grid (no/yes) (coupling_type 2 only)
!    iform_input    INTEGER     Format Input files: 0 = SRV, 1 = NetCDCF
!    ltransport     LOGICAL     Switch for Transport on/off (default .false.)   
!    ibc_type       INTEGER     Bias correction type: 0=None, 1 = Mean Bias, 2 = Low, Mid and High Biases        
!    dn_bcpara      CHARACTER   File with Bias correction parameters
!    lbc_write      LOGICAL     Write bias corrected discharges to file

NAMELIST /HD_CTL/ &
     OUT_EXPNAME, OUT_DATAPATH, YEAR1, MONTH1, DATE_START, DATE_END, NSTEP, DELTA_TIME, &
     UFAKRU, RUNOFF_FILE, DRAINAGE_FILE, FORCING_FREQ, IOUT, &
     COUPLING_TYPE, LCOUPLING_ATM, LCOUPLING_OCE, ICPL_SINKS, ICPL_MASK_TOHD, &
     COUPLING_FILE, LCOUPLING_OUT, IFORM_INPUT, LTRANSPORT, IBC_TYPE, DN_BCPARA, LBC_WRITE

!     ------------------------------------------------------------------
