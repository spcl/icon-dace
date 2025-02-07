# Namelist settings of HD model

```
namelist_settings.md - HD model namelist settings 
 
Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
SPDX-License-Identifier: CC-BY-4.0
See ./LICENSES/ for license information

Authors: Stefan Hagemann
Contact: <stefan.hagemann@hereon.de>
_________________________________________
```

The namelists have to be written by the HD run script into the HD run directory.

### Abbreviations
- Op.S = Operational Setting, e.g. used in HD Vs. 5. (Hagemann and Stacke 2022, doi:10.1016/j.oceano.2022.07.003)

## Namelist HD_CTL - Controls of HD if used as a separate model, e.g. in stand alone runs
- File: namelist.hd

```
     NAME           TYPE        PURPOSE

     out_expname    CHARACTER   experiment name, Default: 'hd'
     out_datapath   CHARACTER   path to where the output data shall be written, Default: './'
     year1          INTEGER     initial year of the run, Default: 1900
     month1         INTEGER     initial month of the run, Default: 1
     date_start     CHARACTER   start date of the run, format YYYYMMDD or YYYY-MM-DD
     date_end       CHARACTER   end date of the run, format YYYYMMDD or YYYY-MM-DD
     nstep          INTEGER     number of time steps within the run (Def. 365) if date_start & date_end are not provided
     delta_time     REAL(dp)    model time step length in seconds, Default: 86400.
     ufakru         REAL(dp)    unit factor for runoff and drainage input data so that their unit becomes [m/s], Default: 1.
     runoff_file    CHARACTER   file with input runoff data if stand alone run, Default: "runoff.nc"
     drainage_file  CHARACTER   file with input drainage data if stand alone run, Default: "drainage.nc"
     forcing_freq   INTEGER     frequency of the forcing data (0: stepwise (default), 1: daily)
     iout           INTEGER     averaging period of the HD output
                                  1   30-day averages
                                  2   decadal averages
                                  3   weekly averages
                                  4   monthly averages without leap years
                                  5   monthly averages with leap years (default)
                                  6   daily averages (Op.S)
                                  7   stepwise output
     coupling_type  INTEGER     Type of coupling to oceanmodel
                                  0   no coupling (default)
                                  1   coupling with no interpolation of discharge within HD
                                  2   coupling with allocation of target grid boxes in HD
     lcoupling_atm  LOGICAL     Switch for coupling to atmosphere, Default: .FALSE.
     lcoupling_oce  LOGICAL     Switch for coupling to ocean, Default: .FALSE.
     icpl_sinks     INTEGER     Redistribution of water in sinks for ocean coupling 
                                  0=None (Def.), 1=Relatively to mouth boxes, 
                                  2= equally to all ocean boxes
     icpl_mask_tohd INTEGER     Switch for the mask on which HD may receive input from atmosphere coupling via YAC
                                  0=HD land boxes (Def.), 1=all HD boxes, 2=use hd_receive_mask 
                                  from file *hd_receive.nc*
                                  If option 2 is used with ocean coupling via YAC, 
                                  receiving boxes with FDIR=-1 are also send to the ocean model.  
     coupling_file  CHARACTER   Input file with coupling information for coupling_type 2, Default: "hdcouple.nc"
     lcoupling_out  LOGICAL     Write discharge on ocean grid (no/yes) (coupling_type 2 only), Default: .FALSE.
     iform_input    INTEGER     Format of input files if stand alone run: 0 = SRV, 1 = NetCDCF (default)
     ltransport     LOGICAL     Switch for tracer transport on/off (Default: .FALSE.), 
                                This transport is on-going development not available in operational HD version.
     ibc_type       INTEGER     Bias correction type: 0=None (default), 1=Mean Bias, 2=Low, Mid and High Biases        
     dn_bcpara      CHARACTER   File (path from directory <HDFILE>) with bias correction parameters
     lbc_write      LOGICAL     Switch for writing bias corrected discharges to file (Default: .FALSE.)

```
Note that when HD was formerly a part of ECHAM, some of these parameters were set in RUNCTL. Original name of namelist was HDALONE_CTL.


## Namelist HYDROLOGY_CTL:  General controls of the HD model
- File: namelist.hdset

```
     NAME           TYPE        PURPOSE
                                                                        
     ldebughd       LOGICAL     debug HD model, Default: .FALSE.
     diag_water_budget LOGICAL  Provide diagnostic global sums of surface and subsurfacerunoff from the 
                                forcing and of discharge entering the ocean, Default: .FALSE.
     lbase          LOGICAL     switch for baseflow calculations, Default: .TRUE.
     locean         LOGICAL     closure of water budget for ocean coupling that was used when HD was called within ECHAM, Default: .FALSE.
     nhd_diag       INTEGER     region number for outflow diagnostic. At 5 Min. resolution, only 7 and 99 are implemented.
                                  0   none (default)
                                  1   Bothnian Bay/Sea
                                  2   Torneaelven
                                  4   St.Lawrence
                                  5   Paraguay
                                  6   Oder
                                  7   Elbe
                                  8   Oranje
                                  9   Amudarya
                                 10   Lena
                                 99   user defined (FBLOG1, FLLOG1, FBLOG2, FLLOG2)
     lhd_highres    LOGICAL     switch for outflow diagnostic on HD model grid, Default: .FALSE.
     fllog1         REAL(dp)    user defined grid cells for diagnostics (with nhd_diag=99)
     fblog1         REAL(dp)      fllog1, fblog1: longitude, latitude of grid cell 1, Default: 0.
     fllog2         REAL(dp)      fllog2, fblog2: longitude, latitude of grid cell 2, Default: 0.
     fblog2         REAL(dp)
     nremap         INTEGER     Type of Interpolation from input (atmospheric) grid to HD grid
                                  0   Input = Output
                                  1   using HDMAP routine by Veronika (default)
                                  2   0.5 degree to 5 Min.
        !                         3   Input = Output + longitudinal shift by 180 degree.
     lhd_rout       LOGICAL     switch for original routing (F) or via index arrays (T), Default: .FALSE.
     fk_rfk         REAL(dp)    Modification factor for k value of riverflow, Default: 1.
     fk_lfk         REAL(dp)    Modification factor for k value of overland flow, Default: 1.
     fk_gfk         REAL(dp)    Modification factor for k value of baseflow, Default: 1.
     irf_vel        INTEGER     Index for type of discharge dependence of riverflow velocity, Op.S: 1
                                  0    None (default)
                                  1    velocity ~ 4th squareroot of Q (riverflow)
     qrf_ref        REAL(dp)    Reference discharge for discharge dependent riverflow velocity, Default/Op.S: 1000.
```

## Namelist HDUSER_CTL:  Provide user-related metainfo to the HD output files
- File: namelist.hduser

```
     NAME           TYPE        PURPOSE
                                                                        
     hd_user        CHARACTER   User name
     hd_cont        CHARACTER   Contact information, e.g. email and/or web address
     hd_inst        CHARACTER   The user's affiliation/institute/company
     hd_instid      CHARACTER   Institutional ID
```