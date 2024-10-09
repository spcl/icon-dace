# Coupling the HD model 

```
readme_coupling.md - Readme file for HD model coupling applications 
 
Copyright (C) 2023, Institute of Coastal Systems - Analysis and Modelling, Helmholtz-Zentrum Hereon
SPDX-License-Identifier: CC-BY-4.0
See ./LICENSES/ for license information

Authors: Stefan Hagemann
Contact: <stefan.hagemann@hereon.de>
_________________________________________
```

## 1. Overview

1. Overview
2. Coupling the HD model to an ocean model - Generally or via OASIS
2.1. Coupling file generation
2.2.  Coupling via file - prepare river runoff forcing for ocean model
2.3. Implementation of coupling in model code
3. Coupling the HD model via YAC
3.1. Implementation of coupling
3.2. Closing the water cycle in the coupled system
3.3. Plot characteristics of the YAC mapping between HD and ICON

Please note that there are also various example scripts in the directory ./scr_coupling/, which can be used for coupling purposes. For further information on those scripts, please contact Ha Hagemann at Hereon (ha.hagemann@hereon.de).

## 2. Coupling the HD model to an ocean model - Generally or via OASIS

Using the OASIS coupler can be enabled with the compiling Option -o ON. Then, the HD model expects to receive input fields of surface runoff and drainage/subsurface runoff on the respective HD model resolution/domain. This means that an interpolation of these fields from the atmosphere/land resolution to the HD model domain has to be done by OASIS before.

The interpolation of the discharge at the ocean inflow points (i.e. the river mouths) onto the ocean model grid is conducted within the HD model as OASIS is not able to thoroughly conduct such an interpolation of point-related fluxes (state of March 2021). In this respect, a coupling file is required for the HD model that can be generated as described in Sect. 2.1. However, you may also feel free to contact Stefan Hagemann at Hereon (stefan.hagemann@hereon.de) who is open to cooperate and generate such a file for you.


### 2.1. Coupling file generation using the script *generate_hd_coupling_file.ksh*
 
#### Running the script

The script is located in ./util and can be called under linux with the following options whereat source and target file must be provided:
```
   -i ,               ID = Experiment name tag or ID. Default if not set: hd_to_ocean 
   -s <Source>,       = HD model parameter file
   -o <Ocean target>, Ocean model land sea mask file 
   -m <Mode No.>,     Ocean mask type:  0: Sea=0, 1: Sea=1 (Def.)
```
Before the coupling file generation is possible, a file had to be created that comprises only the ocean land-sea mask (sea = 1, land = 0) and the respective coordinate information. 
* Name of longitude and latitude coordinate variables: lon, lat
* Name of coordinate dimensions: lon, lat OR x, y
 
Example:
```
generate_hd_coupling_file.ksh -s hdpara.nc -o sea_mask_nemo.nc -i hd_to_nemo
```
#### Settings within the script
Note that within the script, a search radius SRAD needs to be set to regard that the ocean coast line is often smoother than the HD coast line. The search radius is then provided in the namelist *hdtoocean_ctl.nml* set within the script.  Typical search radius values are:
* 200 km to map from the HD 5 arcmin (1/12°) resolution to the NEMO domain covering Baltic Sea, North Sea and NE Atlantic.
* 700 km to map from the HD 0.5° resolution to the global ICON R2B6 grid.  

In order to test which search radius is suitable, one has to look at the variable FMOU_HD_TO_NEMO in the coupling file and check those boxes that are Zero. Here, Zero means that there is a HD river mouth without a coastal ocean point within the search radius.
```
FMOU_HD_TO_NEMO may comprise values of:
 0  Mouth in ocean domain but no coastal mouth point nearby." ;
 1  Mouth pointing to coastal ocean mouth point with no other inflow source." ;
 2  Mouth pointing to coastal ocean mouth point with other inflow sources.
```
Note that when all HD river mouth boxes are associated with an ocean target box for a specific search radius, a larger search radius has no additional effect.
		
In addition, the directory path *DSRC* to the source code files must be edited.
Currently the program runs in the directory from where it is called (*DOUT* is currently only a dummy).

### 2.2. Coupling via file - prepare river runoff forcing for ocean model
The HD coupling file can be utilised in a coupled model system (see Sect. 2.3) or used to generate files with the river runoff forcing on the ocean model grid from the output of an existing HD model simulation. In order to prepare river runoff forcing for the ocean model, the script *convert_hd_discharge_to_ocean.ksh* from the *./util* directory can be used. The script can be executed with

```
    convert_hd_discharge_to_ocean.ksh -i <Exp. ID> -y <First Year> -z <Last Year> 

     -i <7-digit ID>,   ID = Experiment no.
     -y  <First Year>   e.g. YYYY, Default is 1979 if omitted
     -z  <Last Year>    e.g. YYYY, Default: <First Year>
```
e.g. *convert_hd_discharge_to_ocean.ksh -i 7055157 -y 1979* conducts the conversion of discharge from the HD model simulation with the ID 7055157 for the year 1979.

Before you run the script, you have to adapt the section of the script that is denoted as **"This needs to be EDITED"**. This section comprises directory paths, two tags for the naming of the converted output files and the name of the HD coupling file that must be generated beforehand using the script *generate_hd_coupling_file.ksh*.

### 2.3. Implementation of coupling
In order to implement the coupling in a coupled system you may use the existing coupling directive **COUP_OAS** when using the OASIS coupler (see especially in *hd_driver.f90*). However, if you want to implement the coupling in a different way, I recommend using the module file *mo_couple_to_ocean.f90* where the routines *dis_to_ocean* and *read_coupling_info* should be used in a coupled ESM setup. This can be done analogously to how it is done in *hd_driver.f90*. 

The two settings  
```
! coupling_type    type of coupling to oceanmodel (0=0, 1= with no interpol. in HD, 2=with direct alloc.
! coupling_file    input file with coupling information for coupling_type 2
```
should be provided within the namelist **HD_CTL**. Here, coupling_type should be set to 2. Setting 1 would be the normal coupling where the interpolation is done by the coupler. The coupling file is the file generated by the script *generate_hd_coupling_file.ksh*

After the HD model routine *init_hydrology* is called, you may add:

```
  IF (coupling_type .EQ. 2) THEN
    CALL read_coupling_info(coupling_file)
  ENDIF
```

After the HD model routine *hydrologly model* is called, you may add
```
     IF (coupling_type.EQ.2) THEN
       hd_outflow(:,:) = 0._dp
       WHERE (fmou_hd_to_nemo.GT.0.5)
         hd_outflow(:,:) = water_to_ocean(:,:)
       END WHERE
!    
!      ******** Transfer HD model river discharge to ocean grid
       CALL dis_to_ocean
       IF (ABS(SUM(hd_outflow) - SUM(discharge_on_ocean)).GT. 0.01_dp) THEN
         WRITE (message_text,*) 'Discharge sums differ between HD and ocean grid: ', &
               SUM(hd_outflow), ' != ', SUM(discharge_on_ocean)
         CALL message('hd_driver', message_text)
         CALL finish ('hd_driver', 'run terminated.')
       ENDIF
     ENDIF
```
However, then you need to take care that the array (public array of *mo_couple_to_ocean*) *discharge_on_ocean* is sent to your ocean model directly without interpolation.

## 3. Coupling the HD model via YAC

Using the YAC coupler can be enabled using the respective compiling option (see Sect. 3.1). Then, the HD model expects to receive input fields of surface runoff and drainage/subsurface runoff on the respective HD model resolution/domain. This means that an interpolation of these fields from the atmosphere/land resolution to the HD model domain has to be done by YAC before. The YAC interface in HD was implemented by Moritz Hanke (DKRZ).

In order to steer the coupling via YAC, an xml file needs to be generated and put into the run directory of the coupled model (see YAC documentation), e.g. named coupling_<Experiment ID or Name>.xml. With regard to the HD model, the following statements should be included:
```
<components>
    <!--HD-->
    <component id="1">
        <name>LAND</name>
        <model>LAND</model>
        <simulated>land</simulated>
        <transient_grid_refs>
            <!--Surface Runoff-->
            <transient_grid_ref id="1" transient_ref="1" grid_ref="3" collection_size="1"/>
            <!--Subsurface runoff-->
            <transient_grid_ref id="2" transient_ref="2" grid_ref="3" collection_size="1"/>
            <!--River runoff/Discharge-->
            <transient_grid_ref id="3" transient_ref="3" grid_ref="3" collection_size="1"/>
        </transient_grid_refs>
    </component>
    <!--Add the components ocean and HD analogously with component IDs 2 & 3, and grid_ref 2 & 3, respectively.-->
</components>
:
:
<!--This block is for the Atmosphere/Land to HD coupling-->
<couple>
    <!--land-->
    <component1 component_id="1"/>
    <!--HD-->
    <component2 component_id="3"/>

    <!--This block is for surface runoff. Repeat this block analogously for subsurface runoff
        with a transientID=2 -->
    <transient_couple transient_id="1">
        <source component_ref="1" transient_grid_ref="1"/>
        <target transient_grid_ref="1"/>
        <timestep>
            <!-- Example: land provides average fluxes with a time step of 450s. 
                 HD time step and coupling period are 1 day. Land lag 1 is necessary for ICON-->
            <source>PT450S</source>
            <target>P01D</target>
            <!--Note that you may use "accumulate" to sum up for daily values-->
            <coupling_period operation="average">P01D</coupling_period>
            <source_timelag>1</source_timelag>
            <target_timelag>0</target_timelag>
        </timestep>
        <mapping_on_source>true</mapping_on_source>
        <interpolation_requirements>
            <interpolation method="conservative" order="1"/>
        </interpolation_requirements>
        <enforce_write_weight_file filename="weight_runoff_s.nc">true</enforce_write_weight_file>
    </transient_couple>
</couple>
<!--This block is for the HD to ocean coupling-->
<couple>
    <!--Ocean-->
    <component1 component_id="2"/>
    <!--HD-->
    <component2 component_id="3"/>
    <transient_couple transient_id="3">
        <source component_ref="2" transient_grid_ref="3"/>
        <target transient_grid_ref="3"/>
        <timestep>
            <!-- Example: HD provides average fluxes with a time step of 1 day. 
                 Coupling period is 1 day. Ocean time step is 1 hour, Ocean lag 1 is necessary for ICON-->
            <source>P01D</source>
            <target>PT60M</target>
            <coupling_period operation="average">P01D</coupling_period>
            <source_timelag>0</source_timelag>
            <target_timelag>1</target_timelag>
        </timestep>
        <mapping_on_source>true</mapping_on_source>
        <interpolation_requirements>
            <interpolation method="source_to_target_map"  spread_distance="0.0" 
                 max_search_distance="3.6" weighted="ARITHMETIC_AVERAGE"/>
        </interpolation_requirements>
        <enforce_write_weight_file filename="weight_hd_to_ocean.nc">true</enforce_write_weight_file>
    </transient_couple>
</couple>			
```

### 3.1. Implementation of coupling

In order to implement the coupling in a coupled system you may use the existing coupling directive **COUP_YAC**. If you use the compile script *compile_HD_model.sh*, this can be invoked by 
```
compile_HD_model.sh -r 0 -y 1         # -r 0 is compilation for the global 0.5° domain
```
The namelist **HD_CTL** comprises several switches where the coupling via YAC can be steered. You may couple HD to the atmosphere (*lcoupling_atm=.True.*), ocean (*lcoupling_oce=.True.*) or both compartments. In case the water cycle should be closed within a global system, *icpl_sinks* (c) and *icpl_mask_tohd* (b, e) provide some options for this (see below under **3.2.**). 

### 3.2. Closing the water cycle in the coupled system

In a coupled system, the grid resolutions and land sea masks often differ between the atmosphere/land, HD hydrology and the ocean compartments. If these mismatches are not accounted for, some water gets lost in the coupling process. For regional or short-term applications, these water losses can usually be neglected. However, in long-term global climate simulations, a closed water cycle is desired. In the following, we provide some advice and measures how these water losses can be avoided. 

#### a) Atmosphere/Land 

The Atmosphere/land compartment usually comprises lake and glacier areas. In a standard climate model setup, precipitation minus evaporation (P-E) is often not regarded with respect to generating runoff. On the one hand, over glaciers with P > E, snow is generally accumulating. In reality this is compensated by glacier calving, but in a climate model that does not regard this calving, this leads to a continous piling up of the snowpack over glaciers. This is analogous to a continous sink of water in the global water cycle. Hence, introducing a glacier calving routine or putting P-E into runoff may be potential solutions to close the water cycle. 
On the other hands, often no runoff over lake areas is generated within a cimate model. As P < E for most lakes, this is analogous to a continous source of water. In the optimum case, a dynamical lake module should be used where lake storage and lake area can react to the P-E fluxes, and potentially may also generate runoff. If such a module is not available, putting P-E over lakes into surface runoff may be a potential solution. Note that the HD model can deal with negative runoff values. In the past it did not crash when receiving negative runoff from lake or wetland areas. 

#### b) Atmosphere/Land --> HD

Due to the mismatch between the grid resolution and land sea mask of the Atmosphere/Land (A/L) model and the HD model, the interpolation (even if conservative) may generate surface and subsurface runoff over HD ocean grid boxes. In the default coupling setup, runoff is only send to HD land boxes. However, this can be changed by setting the switch *icpl_mask_tohd* in the namelist **HD_CTL**, which determines the mask on which HD may receive input from A/L coupling via YAC. 

```
    icpl_mask_tohd = 0  only HD land boxes (Def.)
                   = 1  all HD boxes may receive runoff water from the A/L model
                   = 2  use hd_receive_mask from file *hd_receive.nc*
```
To avoid sending zero values to HD ocean boxes that do not receive water (if set to 1), set the switch to 2. Then, HD will read the mask of receiving HD boxes with the variable name hd_receive_mask from the file *hd_receive.nc*. Note that runoff water received over HD ocean boxes is just dumped into the ocean as discharge (freshwater flux) without any transport or delay. These discharge fluxes do also appear in the standard output file of the HD model. 
If option 2 is used together with ocean coupling via YAC, discharge values at receiving boxes with FDIR=-1 are also sent to the ocean model (see e) below).  

#### c) Interior drainage basins = Sinks in the HD model

Many interior drainage basins exist around the globe. These are basins without an outflow to one of the oceans. Examples are the Caspian Sea, Lake Chad, the Okavango basin or the Gobi desert. The terminal points of such basins are represented by sink cells, which are HD land grid boxes without an outflow (FDIR = 5). Runoff and discharge that enter these sinks cells is usually lost within a coupled system. In reality, the water entering those cells is evaporating. This might be partially accounted for by a dynamical lake model that is coupled to the discharge model. However, in dry regions this water even evaporates before the deepest point of a basin is reached. To account for this in a coupled system is more difficult and is usually not accounted for. For global long-term climate applications where no water should get lost, this water should be redistributed without leading to noticeable effect in any regions. In HD, the water in sink cells can be redistributed before the transport into the ocean by setting the switch *icpl_sinks*.
```
    icpl_sinks =  0   No redistribution (Default)
                  1   Distribute relatively equal to all mouth boxes with discharge > Zero, 
                      i.e. the added percentage value is the same at all mouth boxes, so that boxes with 
                      small discharge only get a small extra amount, and larger inflows get larger add ons.
                  2   Distribute equally to all ocean boxes.
```

#### d) HD --> Ocean

Some HD river mouths may be too far away from ICON-coastal ocean boxes so that the respective discharge is not mapped. Then you should extend the search radius (max_search_distance="<Distance in degree"; Example default is 3.6) for the YAC coupling in the coupling xml file (see above) until all river mouths are mapped (maybe with ignoring Antarctica). 

#### e) Atmosphere/Land --> HD --> Ocean

In the default setup setup, only discharge at HD river mouths (FDIR=0) is mapped onto coastal ocean points of the ocean model. Hence, some water fluxes as HD ocean points, which are not river mouths but that contain runoff due to the interpolation under b) Atmosphere/Land -> HD, would get lost. This water loss can be avoided by extending the coupling mask from where discharge is sent to the ocean model. Here, the switch *icpl_mask_tohd* (see b) above) can be utilised. If *icpl_mask_tohd* is set to 2 together with ocean coupling via YAC, the receiving HD boxes with FDIR=-1 are also sent to the ocean model. 

The file hd_receive.nc should comprise the mask *hd_receive_mask* indicating all HD boxes that receive runoff from the Atmosphere/Land compartment. Below you find an example on how to create the file hd_receive.nc for ICON (R2B4 grid): 
```
   ICON_GRID=icon_grid_0043_R02B04_G.nc     # ICON grid file, e.g. R2B4, incl. cell_sea_mask 
   HDPARA=hdpara.nc                         # HD Parameter File 
   EXP=<Experiment ID or Name>
   # 1. Select mask with all HD land points
   cdo gtc,0. -selvar,FDIR $HDPARA mask_hd.nc
   # 2. Interpolate a mask with all ICON land points to the HD grid
   cdo -L -remapycon,mask_hd.nc -gtc,0 -selvar,cell_sea_land_mask $ICON_GRID icon_land_r2b4_to_hd_05.nc
   # 3. Rename mask variable to hd_receive_mask and store in HD run directory  
   cdo setvar,hd_receive_mask icon_land_r2b4_to_hd_05.nc ${HDDIR}/<EXP>/hd_receive.nc
```

### 3.3. Plot characteristics of the YAC mapping between HD and ICON

The directory ./util/pyutil comprises the python script *pl_weights.py* that can be used to plot characteristics of the YAC mapping between HD and ICON. The script is steered using the YAML file *config_pl_weights.yml*. It reads and converts YAC weight files from 1D coordinate fields as scatter plots. It has two working modes (IWORK):

1. Plot weights *<DNOUT>.pdf* and write NetCDF file *<DNOUT>.nc* of sending HD boxes for HD --> ICON-O. Here, the script reads the variable src_address from the respective YAC weight file (e.g. weight11.nc) and plots the hd_to_ICON boxes. In addition, it plots the mapping arrows from the HD river mouths to the respective ICON coastal ocean boxes in *mapping_<PLTAG>.pdf*. For this, you also need to provide the ICON grid configuration file (DNICON) including the variable *cell_sea_land_mask*. 
2. Plot weights *<DNOUT>.pdf* and write NetCDF file *<DNOUT>.nc* of receiving HD boxes for ICON-A --> HD. Here, the script reads the variable dst_address from the respective YAC weight file (e.g. weight15.nc) and plots the HD_from_ICON boxes.
     
In order to run the script you need to copy a YAML file *config_pl_weights.yml* to the directory where you run the script as this file is necessary to steer the program. It must include the necessary file names and paths as well as the working mode. Example files are provided in ./util/pyutil/example:
* *config_pl_weights_hdtoicon.yml* can be used to plot mapping arrows and plot/generate mask of HD boxes sending output to ICON-O
* *config_pl_weights_icontohd.yml* can be used to plot/generate mask of HD boxes receiving input from ICON-A

The script runs on Levante at DKRZ. If you want to run the script on other machines, it requires the following python libraries:
* cartopy, datetime, fnmatch, math, matplotlib, numpy, netCDF4, os, psyplot, psy_maps, psy_simple, subprocess, sys, yaml
    
Even though they are not used for the current plotting, one module (taken from a larger project) also needs:
* pyproj, shapefile    

