
# Readme for the Hydrological Discharge (HD) model

```
readme_hd_model.md - Readme for setting up and running the HD model Vs 5.1 and higher
 
Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
SPDX-License-Identifier: CC-BY-4.0
See ./LICENSES/ for license information

Authors: Stefan Hagemann
Contact: <stefan.hagemann@hereon.de>
_________________________________________
```


Note that the HD model is available from different sources. Released versions with a version tag are available at Zenodo under https://zenodo.org/doi/10.5281/zenodo.4893098. Alternatively, you may access the master branch of the public gitlab repository at git@gitlab.dkrz.de:hd-model/hd-couple.git .


0. Content overview 

This readme file comprises a checklist to perform standalone model runs with the HD model. In order to perform simulations where the HD model is coupled with other Earth system model components, you may check the markdown file *./docu/readme_coupling.md* and the scripts stored in the directory ./scr_coupling. In addition, the markdown file *./docu/namelist_settings.md* comprises a description of the HD model namelist settings that are utilized in both standalone and coupled applications. Here, is an overview on the following six sections:

1. Steps to conduct to run the HD model
1.1. Compile the HD model
1.2. Set HD use direcory namesand pathes
1.3. Unpack HD data tar file
1.4. Adapt run script settings in *run_hdmodel.ksh*
1.5. Choose from where you start a new model run
1.6. Adapt run script settings in subscript *hd_run_settings.ksh*
2. Coupling the HD model
3. Necessary preparations for a new regional HD 5 min. domain
3.1. Preparation of 5 Min. subdomain from the global 5 min. HD parameter
file
3.2. Preparation of grid descriptor file and changes in hd_run_settings.ksh
3.3. Preparation of forcing data masks file
3.4. Preparation of remap file
4. Copyright and Licenses
5. Acknowledgements
6. References for HD model


## 1. Steps to conduct to run the HD model 

0. Change to Directory *./scr*

1. Compile HD model using compile script *compile_HD_euro.sh*
   You may need to adapt this script and related settings to your own computing system.
   
   `compile_HD_euro.sh -r <resolution> -c <compiler> -o <oasis_coupling>`

   To get a listing of script options, you can run the script *compile_HD_euro.sh* with the options *-h* or *--help* or *without any option*.
 
   Options for resolution/domain settings with -r/--resolution. Please use Option -r <No.>
   ``` sh
       -r <No. or domain>        Set Resolution/domain (must be set)
       -r 0  or -r 05deg          Global 0.5 degree domain
       -r 1  or -r 5min           Global or regional 5 Min. domains
   ```
       
   Options for compiler settings with -c/--compiler: Please use Option -c <No.> or <compiler name>
   ``` sh
       -c <No.or name>           Set compiler (Default=7 if omitted)
       -c 1  or -c intel           Intel on Mistral (Bull/Atos) at DKRZ 
       -c 2  or -c bull            Bull compiler on Mistral at DKRZ  
       -c 3  or -c intel-strand    Intel on Strand at Hereon"
       -c 4  or -c openmpi         Intel OpenMPI on Mistral at DKRZ
       -c 5  or -c bullxmpi        Updated Bull compiler setting on Mistral at DKRZ
       -c 6  or -c intel-levante   Intel on Levante (new Bull/Atos) at DKRZ
       -c 7  or -c openmpi-levante Intel OpenMPI on Levante at DKRZ (Default)
   ```

   Options for using OASIS coupling with -o/--oasis. Please use Option -o <No.> or ON/OFF
   ``` sh
       -o <No. or ON/OFF>        Coupling to OASIS 
       -o 0  or -o OFF            No (default)
       -o 1  or -o ON             Yes
   ```
   Options for using YAC coupling with -y/--yac. Please use Option -y <No.> or ON/OFF
   ``` ksh
       -y <No. or ON/OFF>        Coupling to YAC 
       -y 0  or -y OFF            No (default)
       -y 1  or -y ON             Yes 
   ```

2. Set HD user directory names/pathes to your system in the file *./scr/setuserdir.com*
   This is a ksh script that defines the four main HD directories via export statements.
   An example is provided in *setuserdir_example.com*. The script *setuserdir.com* should comprise:
   ``` ksh
       export HDMAIN=<HD main Directory that includes script, grid and log directories>
       export HDDIR=<HD Run Directory>
       export HDFILE=<HD input data Directory, e.g. HD parameter and start files>
       export HDFORCING=<HD directory with forcing data in subdirectories>.
    ```
 
   Note that the folders in *setuserdir.com* have to exist. If they do not exist you must create them.

3. Unpack the HD data tar file with the newest data version (e.g. *hdvs5_1_input_data.tar*) 
    within the directory HDFILE
 
4. Adapt run script settings in *./scr/run_hdmodel.ksh* to your requirements

    Note that the HD model (stand alone model run) can be run with the script *run_hdmodel.ksh*. 
    However, this scripts calls several sub-scripts:
    | Script                 | Purpose                                   | 
    | ------                 | -------                                   |
    | *hd_run_settings.ksh*    | Individual simulations settings           |
    | *prepare_hdforcing.ksh*  | Forcing data preparation                  |
    | *hd_post.ksh*            | Postprocessing after end of HD simulation |
   
    In order to run these scripts, the cdo (Climate Data Operators; https://code.mpimet.mpg.de/projects/cdo/) package needs to be installed and loaded. The postprocessing script *hd_post.ksh* also needs access to the nco (NetCDF Operators; https://nco.sourceforge.net/) package.   
    
    * Change Email and Project account in job header

    * Adjust time limit and memory in job header, i.e. time=00:25:00 and mem-per-cpu=500.  

        **Values for Mistral (DKRZ) using daily forcing data:**
        | Domain        | Memory | Run time  | NodeHours |  MaxRSS | 
        | ------        | ------ | --------  | --------- |  ------ |
        | Global 0.5°   |  <500  | < 1 Min.  | 0.0002    |  123 MB |
        | Global 5 Min. |  <2000 | < 55 Min. | 0.019     | 3673 MB |
        | Euro5min      |  <500  | < 5 Min.  | 0.0016    |  243 MB |

        Memory = mem-per-cpu  
        Run time = Typical run time per model year (without remap file generation)  
        MaxRSS = Maximum Resident Set Size of all tasks in job. It shows how much memory 
                is allocated to that process and is in RAM. It includes all stack and heap memory.

    Note for running the HD model at the global 5 Min. resolution:
        To avoid a stacksize related segmentation fault, either set the stacksize to 
        unlimited (i.e. do not limit with ulimit -s 102400) or compile with option -heap-arrays. 
        On Mistral, method 1 was faster (0.019 node-h and 55 Min run time) than method 2 
        (0.023 node-h and 66 Min. run time).

5. Choose from where you start a **new** model run. 
   * Select an experiment ID <EXP> that fulfils your needs.

   * You may either start a model run from the directory <HDMAIN>/scr or from 
     <HDMAIN>/<EXP> (see under 2.). In the first case, you need to 
     modify *hd_run_settings.ksh* directly. In the second case, 
     you have to copy and rename a few scripts (see below). In the first case, the 
     respective copying is done during the first run. 
     Continuation runs always need to be started from <HDMAIN>/<EXP>, and, hence, 
     the modifications of the run settings should be done there.
 
   * Necessary actions in case 2
     ``` ksh
         mkdir ${HDMAIN}/${EXP}  
         cp -p ${HDMAIN}/scr/run_hdmodel.ksh ${HDMAIN}/${EXP}/run_hdmodel_${EXP}.ksh  
         sed -i "s/EXPERIMENT/${EXP}/"  ${HDMAIN}/${EXP}/run_hdmodel_${EXP}.ksh
         cp -p ${HDMAIN}/scr/hd_run_settings.ksh ${HDMAIN}/${EXP}/hd_run_settings_${EXP}.ksh  
         cp -p ${HDMAIN}/scr/hd_post.ksh ${HDMAIN}/${EXP}/hd_post.ksh  
     ```

6. Adapt run script settings in *./scr/hd_run_settings.ksh* (or *<HDMAIN>/<EXP>/hd_run_settings_<EXP>.ksh*)

   a. Set an ID for the forcing: EXPINP=<Forcing ID>       

   b. Set the experiment ID: "EXP=<EXP>" 
         
   c. Set first (YYYY) and last (YEND) year of your simulation.
      - The first year "YYYY=" of the whole simulation should be kept on its initial value
        even if you add more years later by increasing the last year YEND. 

   d. Select forcing "IFORCE=" and HD Model resolution "HDRES="

      - Note that in the settings for the chosen HD model resolution (see further below), 
        you may exchange the standard HD start file by a restart file from a previous simulations.
        An example is provided for HDRES = 2 dependent on the experiment number. 

   e. Set original forcing data resolution FORCE_RES analog to HDRES      
         
      - This is the forcing data resolution without interpolation (e.g. with cdo). 
        It might be used in the forcing data preparation script *prepare_hdforcing.ksh*.

   f. Specify the format of the forcing files in CFORM:  
      - 'nc' = NetCDF (Default) 
         Variable names of surface runoff and drainage (Subsurface runoff)  
             must be runoff and drainage, respectively.
      - 'srv' = Service Format  
         For easyness, the srv format has been often used in the past as input file format,  
             which is a simple binary format that cdos are able to handle.

   g. Set the Coupling type (Default: 0)
      - This setting determines whether the HD output is also provided on a specific
        ocean grid using a coupling file (ICOUPLE=2).  
        Note that the coupling file is not part of this HD model distribution.
        Please contact Stefan Hagemann (stefan.hagemann@hereon.de) for
        coupling purposes.

   h. Set the run time mode IWORK
      - The standard method of running the HD model is annual (IWORK=1), i.e. starting in 1 January. 
        
        | IWORK | Run time mode                            |
        | ----- | -------------                            |
        |   1   | annual, standard method                  |
        |   2   | monthly, not fully explored              |
        |   3   | annual with 30-day months                |
        |   4   | as 1 but final year with nday_final days |

      - Set Variable nday_final if IWORK=4, e.g. nday_final=212 (Jan-July)

   i. Set forcing data resolution that is read by the HD model. 
      - If you interpolate the original forcing data before the data are read, 
        i.e. this is the resolution after interpolation.
        This may be adapted to your own purposes. Using cdos to convert the 
        forcing data to the chosen HD model grid is usually the most efficent method. This
        can be done in the forcing data related section below (see under k). 

   j. Set user specific settings
      - This information will be put as NetCDF attributes in the HD output file
        via the namelist **HDUSER_CTL** in the file *namelist.hduser*.

   k. Preparation of the HD forcing data in script prepare_hdforcing.ksh.
      - This script must also provide two parameters back to the run script:  
        The number of forcing time steps per day (ndt_day) and the unit factor (UFAK) 
        that needs to be applied to the forcing data so that the unit becomes [m/s].

      - Check settings (directory path, unit) and implement data handling (variable selection,
        interpolation) for your forcing and the selected IFORCE. Define new forcing nos. IFORCE
        if appropriate. Note that the number of forcing time steps per day (1 or more)
        determines the time step of the HD model.

      - The current script *prepare_hdforcing.ksh* provides some examples 
        for the preparation of HD forcing data.

      - Note that the forcing data for these example are not part of this distribution.
        Currently, the HD scripts are written in a way that they expect the following file names 
        for the forcing data (Note that this can be modified in the run scripts, too.):

        - NetCDF format: *hdforcing.nc* with the variables runoff and drainage representing   
	                    surface runoff and drainage (subsurface runoff), respectively.
        - Service format: *runoff.srv* for surface runoff  
                         *drainage.srv* for drainage (subsurface runoff)  

   l. Submit the runscript *run_hdmodel.ksh*. 
      - On a slurm system, (e.g. Levante at DKRZ) use:
        `sbatch run_hdmodel.ksh`

        
## 2. Coupling the HD model

Note that this HD model distribution should include everything that is necessary to run the 
standalone/offline version of the HD model. The HD model is also equipped with coupling capabilities using the couplers YAC or OASIS. Information on coupling issues is provided in the file *./docu/readme_coupling.md*.


## 3. Necessary preparations for a new regional HD 5 min. domain 
         
If you want to apply the HD model on a new regional domain, you need to prepare a number of files for this domain. These files are linked by the run script in HDDIR/<EXP> as:
 
| Link name   | File description    | Necessary preparation |       
| ---------   | ----------------    | --------------------- |
| hdpara.nc   | HD parameter file   | Cut out and adapt to the regional domain --> see Sect. 3.1 | 
| hdstart.nc  | HD start file       | Cut out --> see Sect. 3.2 |
| grid_hd.txt | HD grid description | Read from regional HD parameter file --> see Sect. 3.2 |
| masks.nc    | Masks of forcing    | Provide (IMAP = 1) or cut out --> see Sect. 3.3        |   
| rmp_hd.nc   | Remapping file      | Provide remapping weights if IMAP = 1 --> see Sect. 3.4 |
         

### 3.1 Preparation of 5 Min. subdomain from the global 5 min. HD parameter file 

The script *./scr/hdpara_cut.ksh* cuts out a regional subdomain from the global 5 Min. HD parameter file. Gridboxes that may laterally flow out of the domain at the regional boundaries are converted into sinks to avoid numerical exceedances of array bounds. Examples are provided for Europe, Australia and Southeastern Asia (incl. Australia). Further regional domains should be implemented analogously. 

The script needs to be run as
    `hdpara_cut.ksh -f <HD parameter file name> -r <No. or name>`

You need to use the option -f or --file <file name> 
    Here, <HD parameter file name> specifies the path+name of the global 5 Min HD model parameter file
 
Help-Listing: Run with option -h or --help or without any option
 
Options for resolution/domain settings with -r/--region. Please use Option -r <No. or name>
| No. | Tag      | Domain           | Domain bounds        |
| --- | ---      | ------           | -------------        |
|  2  | euro5min | Europe-5 Min.    | -11°W-69°E, 27-72°N  |
|  3  | aus5min  | Australia-5 Min. | 112-154°E, 10-45°S   |
|  4  | sea5min  | SE Asia-5 Min.   | 89-154°E, 37°N-45°S  |

The script will generate the regional HD parameter file as *HDFILE/<HD parameter file name without .nc>_<Tag>.nc*. Note that you have to move the file to the directory *HDFILE/5min/<Tag>* (see also Sect. 3.2) if you use the runscript ./scr/*run_hdmodel.ksh* to run the HD model.
         
### 3.2 Preparation of grid descriptor file and changes in *hd_run_settings.ksh*
         
If you select a new regional HD domain and use the runscript ./scr/*run_hdmodel.ksh*, then you need to prepare a grid descriptor file for this HD model domain. The grid descriptor file is used to interpolate forcing data in the subscript ./scr/*prepare_hdforcing.ksh*. The new grid descriptor file should be stored in ./grid/ 

After you generated the new HD parameter file for your regional domain, you may generate the grid descriptor file with:

1. cdo griddes hdpara_<new domain>.nc
2. From the screen output of the cdo command, copy the upper part related to gridID 1 into an ASCII file *grid_<new_domain>.nc* 
3. Store *grid_<new_domain>.nc* in ./util/
4. Define your new domain in *hd_run_settings.ksh* with an unused number *HDRES* in the inner case $HDRES block within the outer case $HDRES block under '*'.  In the respective part for the new number in the inner case $HDRES block, set:
   REG_TAG="<new domain>" 
   REG_TAG is equal to <Tag> in Sect. 3.1.
5. In addition, you need to set the domain bounds for the cdo statement *cdosel* in the same case block part:
         
```
cdosel='-sellonlatbox,<W longitude>,<E longitude>,<S latitude>,<N latitude>'
```
         
   This cdo statement is used to cut out the HD start file for the new regional domain from the global 5-Min start file.

### 3.3 Preparation of forcing data masks file
      
The masks file for the forcing data is used to provide the grid information on the forcing data that are read by the HD model. If you interpolate the forcing data to the HD grid before they are read by the HD model, the grid should match the used HD model grid generated above (Sect. 3.1/3.2).
In addition it reads four variables. These are relics from the ECHAM5/JSBACH coupling where HD was coupled as a subroutine. The separation from the main part of the HD model is somewhat more complicated so that this among the future plans. Currently, the file still needs to be provided with the following variable names.
  - SLM:   Land sea mask: (1=land, 0=sea/lakes[alake>=0.5]) 
  - SLF:   Land fraction (Read but not used)
  - ALAKE: Lake fraction of grid box
  - GLAC:  Fraction of land covered by glaciers

In the common setup, when the interpolation of the forcing data to the HD grid is done beforehand, you may use the run script *./scr/run_hdmodel.ksh*. In its subscript *hd_run_settings.ksh* with the correct setting of *cdosel* (see Sect. 3.2), the masks file is automatically provided on the new regional domain as *HDFILE/5min/REG_TAG/masks_${REG_TAG}.nc*
         
### 3.4 Preparation of remap file

If you are using HD internal interpolation of the forcing data to the HD grid via a mapping weight file (IMAP = 1), this weight file needs to be provided as *HDFILE*/rmp_hd_<RES_INP>.nc whereat the tag *RES_INP* for the forcing data resolution needs to be defined in *./scr./hd_run_settings.ksh*. The weight file may be generated by CDOs, e.g. using *'cdo genycon'* for conservative remapping.
         
## 4. Copyright and Licenses

Copyright 2021: Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon  
Licensed under the Apache License, Version 2.0    (http://www.apache.org/licenses/LICENSE-2.0)  
Authors of HD Model up to Vs 5.1: Stefan Hagemann, Ha Ho-Hagemann  
Authors of HD Model Vs. 5.2+: Stefan Hagemann, Ha Ho-Hagemann, Moritz Hanke
         
HD parameter data are licensed under the Creative Commons Attribution-ShareAlike 4.0 International License (CC BY-SA 4.0; https://creativecommons.org/licenses/)
Authors of HD parameter data: Stefan Hagemann, Tobias Stacke  

Main Contact: Stefan Hagemann, Helmholtz-Zentrum Hereon, Email: stefan.hagemann@hereon.de

## 5. Acknowledgements

We acknowledge contributions of colleagues (see code headers) from the Max Planck Institute of Meteorology (MPI-M) and collaborating institutions who developed various infrastructure routines and who contributed to the implementation of the HD model into the coupled system ECHAM5/MPIOM and its successor MPI-ESM. We especially appreciated the help of Otto Böhringer (MPI-M) who separated the HD code from the fixed association with the MPI-ESM infrastructure. We are also grateful to Veronika Gayler (MPI-M) for her adapatations to Fortran90 and the implementation of various diagnostics. 


## 6. References for HD model

Main publication for 5 Min. version of HD model (Vs. 4+):  
Hagemann, S., T. Stacke and H. Ho-Hagemann (2020) High resolution discharge simulations over Europe and the Baltic Sea catchment. Front. Earth Sci., 8:12. doi:10.3389/feart.2020.00012.

First publication using HD Vs. 5:  
Hagemann, S., Stacke, T. (2022) Complementing ERA5 and E-OBS with high-resolution river discharge over Europe. Oceanologia, doi:10.1016/j.oceano.2022.07.003.         
         
Further basic HD publications (Vs. 1, 0.5 degree):  
Hagemann, S., L. Dümenil (1998) A parameterization of the lateral waterflow for the global scale. Clim. Dyn. 14 (1), 17-31  
Hagemann, S., L. Dümenil Gates (2001) Validation of the hydrological cycle of ECMWF and NCEP reanalyses using the MPI hydrological discharge model, J. Geophys. Res. 106, 1503-1510  

First publication using coupled 5 Min. version of HD model over Europe:  
Ho-Hagemann, H.T.M., Hagemann, S., Grayek, S., Petrik, R., Rockel, B., Staneva, J., Feser, F., and Schrum, C. (2020) Internal variability in the regional coupled system model GCOAST-AHOI, Atmos., 11, 227, doi:10.3390/atmos11030227

