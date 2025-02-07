#  History of major HD model changes

```
history.md - History of major HD model changes 
 
Copyright (C) 2021, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
SPDX-License-Identifier: CC-BY-4.0
See ./LICENSES/ for license information

Authors: Stefan Hagemann
Contact: <stefan.hagemann@hereon.de>
_________________________________________
```

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## HD Model Version 5.2.3, xx September 2024

### Changed
+ All HD model routines, scripts and utility programs were equipped with a license boilerplate notice.
+ License files were moved to the new directory ./LICENSES

### Added

+ Several text files, especially in ./docu, were associated with the CC-BY-4.0 license.

## HD Model Version 5.2.2, 19 December 2023

### Changed
+ The license for the infrastructure routines inherited from MPI-ESM has changed from the MPI-M software license to BSD-C3.

### Added

+ The possibility to bias correct discharge at river mouths has been implemented. This is related to unpublished work so that the generation of bias correction parameters is currently not part of the HD model distribution.

### Fixed
+ The global attribut istep (designates the number of HD time steps) in the restart file was always one step too low. This has been corrected.

## HD Model Version 5.2.1, 14 November 2023

### Changed
+ The coupling interface to the YAC coupler was updated from YAC2 to YAC3.

### Added

+ It was overseen that the directory ./util was not included in the archived HD model tar file. The directory has been added now.
+ The HD model run can now also be steered by using date_start and date_end via the naemlist HD_CTL. The offline run script has been adapted accordingly.
+ The python routine *./util/pyutil/apply_weights.py* can be used to apply YAC mapping for **ICON to HD** or **HD to ICON**. It requires the YAML file *config_apply_weights.yml*. Examples for both mapping applications are in 
  + *./util/pyutil/example/config_apply_weights_hdtoicon.yml*
  + *./util/pyutil/example/config_apply_weights_icontohd.yml*
Surely, the paths and file names have to be adapted to the individual system settings. 
+ The HD restart file includes now current_date as a global attribute, i.e. the model day the restart file was written.


### Fixed
+ In the coupling of HD using YAC and icpl_mask_tohd = 2, the array hd_receive_mask was read as integer from the file *hd_receive.nc*, which will cause problems if the mask is fractional. Now it is read as double precision array.
+ Overlooked inconsistency that some compilers do not notice. However, the SX-nfort does.
  In *oas_hd.f90* - line 573 (routine decomp_def): remove the DIMENSION(id_size) so that the line becomes:
  INTEGER, INTENT(out) :: id_paral(OASIS_Box_Params)
+ if compiler directive __SX__ is set, the 
  CALL util_backtrace  in routine *mo_exception.f90* is not conducted. 
  On the NEC at DWD, this call was leading to a compilation error.
+ Minor bug corrected for the redistribution of sink inflows.

## HD Model Version 5.2.0, 3 May 2023

### Added

+ Moritz Hanke, DKRZ became a co-author. He did:
  + Reprogramming of interface for coupled application.
  + Implementation of interface to the YAC coupler.
+ Added utility *./util/pyutil/pl_weights.py* to plot/generate sending & receiving HD masks from YAC weight files. It also allows to plot the mapping arrows from HD to ICON.
+ Added several switches and related code to allow closure of water balance in a coupled system. Switches can be set in namelist *hd_ctl*:
   + If *coupling_type* !=0: separate switches were implemented for atmosphere and ocean coupling.
   + Switches were implemented to enable the conservation of all runoff/discharge fluxes when mapping from Atmosphere -> HD -> Ocean in coupled applications.
+ Added readme file *./docu/readme_coupling.md* for HD coupling exercises (e.g. YAC and OASIS couplers) including a description of the necessary steps to close the water balance in a coupled system.
+ Add a description of necessary preparations for using a new regional HD 5 min. domain into the readme file *./docu/readme_hd_model.md*. 
+ Settings for DKRZ HPC system Levante were implemented. OpenMP and MPI versions yield bit identical results.

### Fixed
+ Correct routine *mo_time_control.f90* for the rerun warning. The original ECHAM related programming was expecting that the first time step of the next day after the simulation period was run before the model stops. However, in HD the model stops in the end of the simulation period.
+ In *mo_exception*: Remove "#ifdef __SX__   CALL mesput('Traceback: ', 11, 1)"
+ Obsolete parameter *iswrit* was removed from namelist *hd_ctl*.

### Changed

#### Cleanups of code and scripts
+ *runscript run_hdmodel.ksh*: write info on forcing ID and date of exe file so that info appears in the job log output.
+ Write name of the start file as a global attribute into model output file *<EXP>_meanflow_<YYYY>.nc*>.nc* of the first year <YYYY>.
+ Separate experiment/run settings and postprocessing from HD run script *run_hdmodel.ksh* into the sub scripts *hd_run_settings.ksh* and *hd_post.ksh*, respectively. Note that if a job hangs on levante after the model has been run, one can call the postprocessing script also interactively.
+ Rename *namelist.echam* into *namelist.hdset* in run script and *mo_hydrology.f90*
+ Remove *hd_domain.inc* from HD code --> results are binary identical, but lon/lat info slightly differs as this is now taken from the parameter file instead of purely calculated. Therefore, the HD executables do not include the regional domain name anymore: *hd_05.exe* or *hd_5min.exe*
+ Move JSBACH related code parts to *mo_jsbach_to_hd.f90* (unused module)
+ Rename compile script for the offline version from *compile_HD_euro.sh* to *compile_HD_model.sh* 
+ Write discharge_on_ocean with float (F32), not double (F64) --> Variable conversion using cdo is not necessary   
+ Set defaults of *lhd_rout* to .TRUE. and *locean* to .FALSE.
+ Set default for *iform_input* = 1
+ All readme files are located now in ./docu. Hence, *history.md* is moved from ./scr to ./docu
+ Revise readme file, change it to markdown format and rename it to *./docu/readme_hd_model.md*.
  + Revise/expand desciption of forcing data preparation in readme file
  + Description of the preparation of a 5 Min. subdomain from the global 5 min. HD parameter file becomes now a subsection of the section 3 on necessary preparations for using a new regional HD 5 min. domain.



## HD Model Version 5.1.0, 17 Nov. 2021

### Added
- Implement Netcdf conform output
- Implemented dependency of flow velocity on discharge

  For our initial 5 Min. setup over Europe, reasonable results were yielded with the HD model for many European rivers (with KGE often larger than 0.4). However, a global simulation showed that for very large rivers, the simulated discharge was lagging behind. This is something that can also be seen for the Danube, the second largest river of Europe. 

  For small or medium-sized rivers, the present HD5 parameter settings seems to be appropriate, while for larger river it seems that a dependence of the flow velocity v on the discharge q must be regarded. Leopold and Maddock (1953) found an empirical dependency v = k Q^m. 
  We defined a reference discharge qref = 1000 m^3/s below which no correction is necessary. 
  For q=qref, vref = k * qref^m  
  v = vref * (Q/Qref)^m  

  Consequently, for q > qref, we correct the flow velocity vref determined by the HD model river flow parameters by a correction factor (Q/Qref)^m.

  Leopold and Maddock (1953) found an average value of m=0.34 for 20 semi-arid rivers at stations in the Great Plains, and an average value of m=0.1 was obtained from various US rivers downstream. For our setup we chose an intermediate m=0.25.

### Fixed
- Removal of two bugs in the parameter file for Vs. 5.0
  - Accidentally, there was no wetland impact on riverflow.
  - For the wetland fraction at 1 km derived from GLWD wetlands, GLWD type 10 (50-100%) was erroneously associated with a fraction of 0.875 instead of 0.75, and GLWD type 12 was erroneously omitted. The latter is now associated with 0.125. 

### Changed
- Corrected usage of GLWD wetlands
  - Removal of unspecified wetland areas (GLWD types 10-12) for the impact on riverflow that mainly only occur over North America. Now, GLWD wetland types 4-9 affect river flow.
  - Removal of GLWD rivers from ESA waterbodies. For the calculation of lake fraction from ESA water bodies, the fraction of GLWD rivers (type 3) was substracted from the ESA waterbody fraction.
- Scaling factors on 5 Min. model parameters for overland flow (2) and baseflow OLF (4) were implemented into the HD parameter file and taken out of the run script. 

## Reference HD model Version: 5.0, 10 June 2021
