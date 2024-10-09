! readme_hdmodel_vs5_0.txt - Readme to HD model version 5.0  
! 
! Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
! SPDX-License-Identifier: CC-BY-4.0
! See ./LICENSES/ for license information
!
! Authors: Stefan Hagemann
! Contact: <stefan.hagemann@hereon.de>
!_________________________________________

======================
Readme to run HD model
======================

0. Change to Directory ./scr

1. Compile HD model using compile script compile_HD_euro.sh
   You may need to adapt this script and related settings to your own computing system.
   
   compile_HD_euro.sh -r <resolution> -c <compiler> -o <oasis_coupling>

   Help-Listing: Run with Option -h or --help or without any option
 
   Options for resolution/domain settings with -r/--resolution. Please use Option -r <No.>
     -r <No. or domain>        Set Resolution/domain (must be set)
       -r 0  or -r 05deg          Global 0.5 degree domain
       -r 1  or -r 5min           Global 5 Min. domain
       -r 2  or -r euro5min       European 5 Min. domain:   -11°W-69°E, 27-72°N
       -r 3  or -r aus5min        Australian 5 Min. domain:  112-154°E, 10-45°S
 
   Options for compiler settings with -c/--compiler: Please use Option -c <No.> or <compiler name>
     -c <No.or name>           Set compiler (Default=1 if omitted)
       -c 1  or -c intel          Intel on Mistral at DKRZ (Default)
       -c 2  or -c bull           Bull compiler on Mistral at DKRZ  
       -c 3  or -c intel-strand   Intel on Strand at HZG"
 
   Options for using OASIS coupling with -o/--oasis. Please use Option -o <No.> or ON/OFF
     -o <No. or ON/OFF>        Coupling to OASIS 
       -o 0  or -o OFF            No (default)
       -o 1  or -o ON             Yes    


2. Set HD user directory names/pathes to your system in the file ./scr/setuserdir.com
   This is a ksh script that defines the four main HD directories via export statements.
   An example is provided in setuserdir_example.com. setuserdir.com should comprise:
     export HDMAIN=<HD main Directory that includes script, grid and log directories>
     export  HDDIR=<HD Run Directory>
     export HDFILE=<HD input data Directory, e.g. HD parameter and start files>
     export HDFORCING=<HD directory with forcing data in subdirectories>.
 
   Note that the folders in setuserdir.com have to exist. If they do not exist you must create them.

3. Unpack the data tar file hdvs5_input_data.tar within the directory HDFILE
 
4. Adapt run script run_hdmodel.ksh to your requirements

  a) Either copy script run_hdmodel.ksh into your <HDMAIN>-Directory (see under 2.) 
     and name it run_hdmodel_<EXP>.ksh or modify run_hdmodel.ksh directly. 
     (It will do the respective copying during the first run.)
     EXP is the chosen experiment ID you set in the script under d)

  b) Change Email and Project account in job header

  c) Adjust time limit and memory in job header, i.e. time=00:25:00 and mem-per-cpu=500. 
       Values for Mistral (DKRZ) using daily forcing data:
         Domain          Memory   Run time   NodeHours    MaxRSS 
         Global 0.5°     <500     < 1 Min.   0.0002        123 MB
         Global 5 Min.   <2000    < 55 Min.  0.019        3673 MB
         Euro5min        <500     < 5 Min.   0.0016        243 MB

       Memory = mem-per-cpu
       Run time = Typical run time per model year (without remap file generation)
       MaxRSS = Maximum Resident Set Size of all tasks in job. It shows how much memory 
                is allocated to that process and is in RAM. It includes all stack and heap memory.

     Note for running the HD model at the global 5 Min. resolution:
       To avoid a stacksize related segmentation fault, either set the stacksize to 
       unlimited (i.e. do not limit with ulimit -s 102400) or compile with option -heap-arrays. 
       On Mistral, method 1 is faster (0.019 node-h and 55 Min run time) than method 2 
       (0.023 node-h and 66 Min. run time).

  d) set an experiment ID "EXP=" that fulfils your needs

  e) Select forcing "IFORCE=" and HD Model resolution "HDRES="

     Note that in the settings for the chosen HD model resolution (see further below), 
     you may exchange the standard HD start file by a restart file from a previous simulations.
     An example is provided for HDRES = 2 dependent on the experiment number. 

  f) Specify the format of the forcing files in CFORM: 
       'nc' = NetCDCF. Variable names of surface runoff and drainage (Subsurface runoff)
                       must be runoff and drainage, respectively.
      'srv' = Service Format (Default)
              For easyness, the srv format has been often used in the past as input file format, 
              which is a simple binary format that cdos are able to handle. 

  g) Set the Coupling type (Def.: 0), i.e. whether the HD output is also provided on a specific
     ocean grid using a couping file (ICOUPLE=2) (Note that the coupling file is not part of this
     HD model distribution. Please contact Stefan Hagemann (stefan.hagemann@hereon.de) for 
     coupling purposes.

  h) Set first (YYYY) and last (YEND) year of your simulation.
     The first year "YYYY=" of the whole simulation should be kept on its initial value
     even if you add more years later by increasing the last year YEND. 
     Set final Year "YEND=" of your simulation

  i) The standard method of running the HD model is annual (IWORK=1), i.e. starting in 1 January. 
     However, it is possible to run the model for monthly time slices (IWORK=2). Note that
     the latter method is not fully explored. 

  j) Set forcing data resolution that is read by the HD model. If you interpolate the original
     forcing data before the data are read, this is the resolution after interpolation.
     This may be adapted to your own purposes. Using cdos to convert the 
     forcing data to the chosen HD model grid is usually the most efficent method. This
     can be done in the forcing data related section below (see under l). 

  k) Set user specific settings that will be put as NetCDF attributes in the HD output file
     via namelist HDUSER_CTL in the file namelist.hduser

  l) Preparation of the HD forcing data in script prepare_hdforcing.ksh.
     This script must also provide two parameters back to the run script:
        The number of forcing time steps per day (ndt_day) and the unit factor (UFAK) 
        that needs to be applied to the forcing data so that the unit becomes [m/s].

     Check settings (directory path, unit) and implement data handling (variable selection,
     interpolation) for your forcing and the selected IFORCE. Define new forcing nos. IFORCE
     if appropriate. Note that the number of forcing time steps per day (1 or more) determines the
     time step of the HD model.

     The current script prepare_hdforcing.ksh provides some examples for the preparation of HD forcing data.
     Note that the forcing data for these example are not part of this distribution.
     Currently, the HD scripts are written in a way that it expects the following file names 
     for the forcing data (Note that this can be modified in the run script, too.):

     NetCDF format:  hdforcing.nc with the variables runoff and drainage representing 
	             surface runoff and drainage (subsurface runoff), respectively.
     Service format: runoff.srv for surface runoff 
                     drainage.srv for drainage (subsurface runoff)  

  m) Set your email to which the finished job email should be send in the end of the script.

  n) Submit the runscript. On a slurm system, (e.g. mistral at DKRZ) use:
     sbatch run_hdmodel.ksh

======
Remark
======
Note that this HD model distribution should include everything that is necessary to run the 
standalone/offline version of the HD model. The HD model is also equipped with coupling capabilities
using the OASIS coupler that can be enabled with the compiling Option -o ON. Then, the HD model expects 
to receive input fields of surface runoff and drainage/subsurface runoff on the respective 
HD model resolution/domain. This means that an interpolation of these fields from the 
atmosphere/land resolution to the HD model domain has to be done by OASIS before. 

The interpolation of the discharge at the ocean inflow points (i.e. the river mouths) onto
the ocean model grid is conducted within the HD model as OASIS is not able to thoroughly conduct
such an interpolation of point-related fluxes (state of March 2021). In this respect, a 
coupling file is required for the HD model that is not part of this distribution. 
The program to generate such a coupling file is quite complicated and currently not in a very user
friendly shape (this is work in progress). Hence, this is not part of this distribution. 
However, please feel free to contact Stefan Hagemann at HZG (stefan.hagemann@hereon.de) who
is open to cooperate and generate such a file for you.

========================================================================
Preparation of 5 Min. subdomain from the global 5 min. HD parameter file 
========================================================================
The script hdpara_cut.ksh cuts out a regional subdomain from the global 5 Min. HD parameter file.
Gridboxes that may laterally flow out of the domain at the regional boundaries are 
converted into sinks to avoid numerical exceedances of array bounds.
Examples are provided for Europe and Australia. Further regional domains should be 
implemented analogously. 

The script needs to be run as
    hdpara_cut.ksh -f <file name> -r <No. or name>

You need to use the option -f or --file <file name> 
    <file name> specifies the path+name of the global 5 Min HD model parameter file
 
Help-Listing: Run with option -h or --help or without any option
 
Options for resolution/domain settings with -r/--region. Please use Option -r <No. or name>
 No.: 1  = euro5min  -->   European 5 Min. domain: -11°W-69°E, 27-72°N
      2  = aus5min   --> Australian 5 Min. domain:  112-154°E, 10-45°S


=======================
Copyright and Licenses
=======================

Copyright 2021: Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
Licensed under the Apache License, Version 2.0
   (http://www.apache.org/licenses/LICENSE-2.0)
Authors of HD Code: Stefan Hagemann, Ha Ho-Hagemann  
Authors of HD parameter data: Stefan Hagemann, Tobias Stacke
Main Contact: Stefan Hagemann, Helmholtz-Zentrum Hereon, Email: stefan.hagemann@hereon.de

=======================
Acknowledgements
=======================

We acknowledge contributions of colleagues (see code headers) from the 
Max Planck Institute of Meteorology (MPI-M) and collaborating institutions
who developed various infrastructure routines and who contributed to the implementation
of the HD model into the coupled system ECHAM5/MPIOM and its successor MPI-ESM. 
We especially appreciated the help of Otto Böhringer (MPI-M) who separated the HD code
from the fixed association with the MPI-ESM infrastructure. We are also grateful
to Veronika Gayler (MPI-M) for her adapatations to Fortran90 and the implementation of various diagnostics. 

=======================
Copyright and Licenses
=======================

Copyright 2021: Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon

HD model code is licensed under the Apache License, Version 2.0
   (http://www.apache.org/licenses/LICENSE-2.0)
Authors of HD Code: Stefan Hagemann, Ha Ho-Hagemann  

HD parameter data are licensed under the Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0; https://creativecommons.org/licenses/)
Authors of HD parameter data: Stefan Hagemann, Tobias Stacke

Main Contact: Stefan Hagemann, Helmholtz-Zentrum Hereon, Email: stefan.hagemann@hereon.de

=======================
Acknowledgements
=======================

We acknowledge contributions of colleagues (see code headers) from the 
Max Planck Institute of Meteorology (MPI-M) and collaborating institutions
who developed various infrastructure routines and who contributed to the implementation
of the HD model into the coupled system ECHAM5/MPIOM and its successor MPI-ESM. 
We especially appreciated the help of Otto Böhringer (MPI-M) who separated the HD code
from the fixed association with the MPI-ESM infrastructure. We are also grateful
to Veronika Gayler (MPI-M) for her adapatations to Fortran90 and the implementation of various diagnostics. 

=======================
References for HD model
=======================
Main publication for 5 Min. version of HD model (Vs. 4+):
Hagemann, S., T. Stacke and H. Ho-Hagemann (2020) High resolution discharge simulations over Europe and the Baltic Sea catchment. Front. Earth Sci., 8:12. doi: 10.3389/feart.2020.00012.

Further basic HD publications (Vs. 1, 0.5 degree)
Hagemann, S., L. Dümenil (1998) A parameterization of the lateral waterflow for the global scale. Clim. Dyn. 14 (1), 17-31
Hagemann, S., L. Dümenil Gates (2001) Validation of the hydrological cycle of ECMWF and NCEP reanalyses using the MPI hydrological discharge model, J. Geophys. Res. 106, 1503-1510

First publication using coupled 5 Min. version of HD model over Europe 
Ho-Hagemann, H.T.M., Hagemann, S., Grayek, S., Petrik, R., Rockel, B., Staneva, J., Feser, F., and Schrum, C. (2020) Internal variability in the regional coupled system model GCOAST-AHOI, Atmos., 11, 227, doi: 10.3390/atmos11030227


