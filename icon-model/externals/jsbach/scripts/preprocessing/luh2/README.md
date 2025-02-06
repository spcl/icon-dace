# LUH2 preprocessing

> LUH2 preprocessing is conducted in two steps:
- in the first step a preprocessing is conducted as was done for **CMIP5**
- in the second step the states are scaled with a given `veg_ratio_max` to conserve LUH crop and pasture area as close as possible ("**CMIP6** way") -- the transitions are scaled accordingly 

    *[justification: we know that deserts are assumed to be part of the natural vegetation in LUH2 -- see TR email 02.11.16]*

  LUH2 preprocessing is a collection of several bash scripts (`tools`) and small `fortran` programs (`src`).

## Pre-requires

- Downloaded LUH2 data from, e.g., [https://luh.umd.edu/] in a directory `orgdata` 
  > INFO: For downloading LUH2 data, you might adapt and use `get_raw_data_LUH2_scenarios.bash` in the `attic`.
  > NOTE: in addition to the transitions and the states, also the `staticData_quarterdeg.nc` might be downloaded, in order to extract the `fstnf` value required to assign rangelands to pasture or to natural vegetation
- Alternatively to downloading of the `fstnf` map, it is also possible to create an own `fstnf` map, e.g. in order to assign all rangelands to pasture (`fstnf = 1`) or all to natural vegetation (`fstnf = 0`), the latter is what has been done for CMIP6
- Paths correctly set-in the `MAIN_preprocess_LUH2.bash`
- Membership to project `mj0060` and/or access to a maximum vegetation fraction file, e.g., 
  >`vegRatioMaxFile="/pool/data/JSBACH/prepare/T${ires}/${vrmFileName}.nc"`

  
## Quick start

1. Execute the compilation script from the `luh2_preprocessing` main directory
    
    ```bash
    ./compile.sh
    ```
    or
    ``` bash
    ./compile.sh -e <YOUR-EMAIL-ADDRESS> -p <YOUR-PROJECT-NUMBER> -c <USER-CONFIG-FILE>
    ```
    > INFO: 
    > - Upon the first call of `compile.sh` without options will copy the MAIN runscript `preprocess_LUH2.bash` and has to be edited (-> see 3.)
    > - Calling `compile.sh` with the given options will insert email, project number and config file into the header of `preprocess_LUH2.bash (-> you can edit it afterwards -> see 3.).
    > - `preprocess_LUH2.bash` will not be overwritten by consecutive calls of `compile.sh`
    > - In case you want to clean up the `bin` directory execute `./compile.sh clean`

2. Copy one of the example configurations (`examples`) and adapt to your needs
   
    ```bash
    cp ./examples/user_config_test ./examples/user_config_xxx_yyyy
    ```

    Usually these need adapting:
  - Your contact information
    >`user_email`
  - Output directory (use absolute paths) 
    >`mainWorkingPath`
  - Years to process (note: at least two years are required, and transitions are usually only available for one year less than states are available)
    >`selStartYear` and `selEndYearStates`
  - Downloaded original data
    >`pathToOrgData`
  - Name of the `states` files
    > `statesFile`
  - Name of the `transistions` files
    > `transitionFile`
  - Add a contact for the data in
    > `calledFrom="MAIN_preprocess_LUH2.bash [contact:] -- using ${vegRatioMaxFile}"`


3. Open `preprocess_LUH2.bash` in an editor of your choice and 

- YOUR-PROJECT: accounting information
- yourname@YOUR-INSTITUTE-DOMAIN: email address
- YOUR-CONFIG-FILE: own `user_config`

4. Execution

- For **test** purposes (max. 2 years!) you can run

    ```bash
    nohup ./preprocess_LUH2.bash 2>&1 > log &
    ```

  on the login node. This will put the process into the background and write all messages to a `log` file.

- For **production** send it to the queueing system (**obs: right user- and project ID**)

    ```bash
    sbatch ./MAIN_preprocess_LUH2.bash
    ```

### For levante (formerly mistral)

The machine is recognized and the compile environment is set automatically.

- On `levante`: `ifort` set up is available
- On `mistral`: `nag` compiler was used (will be used on `levante` when it becomes available)

### locally

- On any other machine, you'd adapt `compile.sh` and add the necessary library paths
- Example: jessie (OBS: machine recognition not tested!)

## Note

> - Area fractions for the different vegetation types in the grid cells of the LUH data do not necessarily add up to one.
> - For LUH2v2h we know that the given area fractions add up to 1 - (ice + water),  i.e. ocean, glaciers, and water bodies on land lead to area fractions 
smaller than one (or in pure ocean/inland lake cells to a missing value).
> - Since JSBACH requires a vegetated area of 1, all vegetation types in a grid cell with an area fraction sum less than one were increased proportionally in the preprocessing for JSBACH in **CMIP5**.
> - This was justified by the assumption that some of the cells would be counted as ocean and some as land which globally would cancel out (?). 
> - However, for rivers and smaller lakes this procedure might not be appropriate and could lead to an overestimation of crop and pasture in the respective grid cells, 
particularly in combination with the veg_ratio_max scaling of crop and pasture introduced in the preprocessing for JSBACH in **CMIP6**. 
> - For **CMIP6** a more sophisticated treatment is not possible, because it would require using additional information.
> - Beyond **CMIP6**, for ICON-Land, it would be desirable to treat coastlines differently than inland water bodies, with the latter being treated in separate tiles, such that the above described proportional increase of the area fractions would not be required anymore. 
> - **Thus, the LUH processing will need to be revisited, particularly also the application of relative transitions, which then should only refer to the dry fraction of a grid cell.**

## HISTORY

> Note: in `aggregate_LUH2_transitons.bash`, `aggregate_LUH2_harvest.bash` and `flip_latdata_LUH2_transitions.bash` an old `cdo` version was used on `mistral`, because else issues with the dimensions occurred

- `cdo`: /sw/rhel6-x64/cdo/cdo-1.9.6-magicsxx-gcc64/bin/cdo
- `nco`: /sw/rhel6-x64/nco/nco-4.7.5-gcc64/bin
