# ICON Community Interface :: User Guide

The Community Interface (ComIn) organizes the data exchange and simulation events between the ICON model and "3rd party modules". 
The information in this document provides a starting point for new users (plugin developers).
A more detailed technical specification of the interfaces is given [here](icon_comin_doc.md).

@tableofcontents{html:2}

## Use existing ComIn plugins

ComIn is bundled with several example plugins. Here is a quick start guide on using these plugins with ICON as the host model.

The build process is located in the sub-directory `$BASEDIR` which you should set as an environment variable:
```bash
cd $BASEDIR
```
### Repository checkout
  
Clone the ICON repository:  
```bash  
module load git
git clone git@gitlab.dkrz.de:icon/icon.git
  
cd icon
git submodule update --init
```

### Example plugins

After cloning ICON, the example plugins can be found in `externals/comin/plugins`. In this folder there are different examples written in Fortran, C and Python:

- @ref simple_fortran_plugin.F90  
  Simple ComIn plugin written in the **Fortran** programming language
- @ref simple_c_plugin.c  
  Simple ComIn plugin written in the **C** programming language
- @ref simple_python_plugin.py  
  Simple ComIn plugin written in **Python**, using the ComIn Python adapter
- @ref point_source.py  
  Test plugin requesting a tracer that participates in ICON's turbulence and convection scheme.

In order to use the plugins with ICON, the first step involves building the ICON and plugins. The instruction of building and using the plugins is explained here for the *Levante_gcc* and the *DWD_nec* platforms:

### Levante platform

#### Build instructions for Levante_gcc 

- **Load the necessary Python packages.**  
  To make use of the Python adapter a proper Python installation must be available, in particular if the `mpi4py` package should be used by the plugins.  
We recommend to use the Python installation contained in Levante's `spack` install tree. However, the package `py-mpi4py` is not contained there, therefore we created a custom `spack` install tree which can be used by placing the following file in your home directory:
```bash
mkdir ~/.spack
cat <<EOF >> ~/.spack/upstreams.yaml
upstreams:
  community_spack:
    install_tree: /work/k20200/k202160/community-spack/install
  system_installs:
    install_tree: /sw/spack-levante
EOF
```
Load the package with
```bash  
spack load py-mpi4py  
```

- **Build ICON with ComIn.** Out-of-source build is recommended. 
```bash
cd icon
mkdir build && cd build
../config/dkrz/levante.gcc --enable-comin
make -j6
```

- **Build the plugins.**  
  Through the building process, a shared library of the plugin is generated, allowing for dynamic loading while running  ICON. Here we use the Python adapter that is shipped together with the ComIn library.  Note that the plugin has to be enabled first:
```bash
(cd externals/comin/build && cmake -DCOMIN_ENABLE_EXAMPLES=ON  -DCOMIN_ENABLE_PYTHON_ADAPTER=ON .)
(cd externals/comin/build && make)
```

#### Modify your ICON run script.  

It is recommended  to adjust your run script template before configuring ICON. This modified version of your template will be copied to the `build/run` subdirectory then.

Modifying your experiment's template involves two parts.

##### Add the path of your plugin's shared library in to the `LD_LIBRARY_PATH`. 

To accomplish this, there are two options available.

- Add it directly in to your run script:
```bash
export LD_LIBRARY_PATH="${path_to_plugin}:$LD_LIBRARY_PATH"
```
			
- Alternatively use the auxiliary function `add_comin_setup` in the run script which does the same automatically. To use this function your `basedir` variable must be set to the `build` directory of your ICON installation.
```bash
add_comin_setup "/externals/comin/build/plugins/simple_fortran"
```


##### Add the &comin_nml to the namelist "atmo_namelist".

Here you add one or several plugins in the order you want to use them.  `plugin_list` has different items (different plugins ) which each item could have 5 components to be set.  

- `name `:  the name of the plugin.
- `plugin_library  `:  the shared library file associated with the plugin.
- `primary_constructor  `: name of the primary constructor. It must be specified  if it is **not** ```comin_main```.
- `options  `: offers the possibility to pass a character string (e.g. a python script filename) to the plugin.
- `comm  `: denotes the name of the MPI communicator that is created for this particular plugin. This is useful when exchanging data with other running processes. The parameter `comm` can be left as an empty string if the application does not require a communicator for this plugin.<br>

 Note: The first two components are mandatory to set. 
```bash
&comin_nml
plugin_list(1)%name          = "simple_fortran_plugin"
plugin_list(1)%plugin_library = "libsimple_fortran_plugin.so"
```

#### Run the experiment on Levante.

The following commands copy the modified version of your template to the `build/run` subdirectory and launch the batch job

```bash
cd build
./make_runscripts --all
cd run
sbatch name_of_your_runscript.run
```

An alternative option is to run your experiment on interactive nodes.  Allocate resources for a suitable (cheap!) cluster queue on `Levante` and wait for the interactive job to start:
	  
```bash  
salloc -p interactive -A <account> --time=04:00:00 --tasks-per-node=128 --nodes=1 --exclusive  
```  

Then run the test interactively (remember to make your `$BASEDIR` known to the new shell: `export BASEDIR= ...`):  

```bash  
cd $BASEDIR/icon/build/run  
./name_of_your_runscript.run  
```

### DWD NEC platform

The entire process closely mirrors the steps outlined for Levante_gcc. However, it is important to note that in the present release, the changes including the building of ICON using ComIn interactively or on Buildbot on DWD_nec are not integrated into the source code. To address this, interested users can reach out to a developer to acquire the version that incorporates the necessary modifications for building on DWD_nec. Once the updated version is obtained, the subsequent instructions can be  followed.

* Create a build directory.  
```bash
mkdir build && cd build
mkdir VH VE
```

* Build ICON on vector host.  
```bash
cd VH
../../config/dwd/rcl.VH.gcc --enable-comin
make -j6
```

* Build plugins on vector host.  
```bash
cd externals/comin/build
module purge
module load 'apps sx/default gcc/11.2.0 mpi/2.22.0 netcdf4/4.7.3-VH-gnu hdf5/1.10.5-VH-	gnu eccodes/2.25.0-VH-gnu aec/1.0.3-VH-gnu aocl/2.1-VH-gnu szip/2.1.1-VH-gnu unsupported cmake/3.26.4'
sed -i 's/-static//g' CMakeCache.txt
cmake -DCMAKE_C_COMPILER=mpincc -DCMAKE_Fortran_COMPILER=mpinfort -DCMAKE_CXX_COMPILER=g++ -DCOMIN_ENABLE_EXAMPLES=ON .
make
```

* Build ICON on vector engine.  
```bash
cd ../VE
../../config/dwd/rcl.VE.nfort --enable-comin
make -j6
```

* Build plugins on vector engine.  
```bash
cd externals/comin/build
module purge
module load 'sx/default nfort/4.0.0 nc++/4.0.0 mpi/2.22.0 netcdf4/4.7.3-sx hdf5/1.10.5-sx eccodes/2.28.0-sx aec/1.0.3-sx nlc/2.3.0 libxml2/2.9.10-sx szip/2.1.1-sx zlib/1.2.11-sx unsupported cmake/3.26.4'
cmake -DCOMIN_ENABLE_EXAMPLES=ON . 
make
```

* Modify the run script.  
This step is almost the same as is explained for Levante_gcc except that one also must add the path of shared library of plugin(s) to `VE_LD_LIBRARY_PATH`.  
```bash
export LD_LIBRARY_PATH="${path_to_plugin_on_VH}:$LD_LIBRARY_PATH"
export VE_LD_LIBRARY_PATH="${path_to_plugin_on_VE}:$VE_LD_LIBRARY_PATH"
```
 or use the auxiliary function `add_comin_setup` . This function does the same for both vector host and vector engine automatically.
```bash
path_to_plugin=`/externals/comin/build/plugins/simple_fortran`
add_comin_setup "$path_to_plugin"
```

* Run the experiment.  
```bash
cd /build/VE
./make_runscripts --all
cd run
qsub name_of_your_runscript.run
```

## Write ComIn plugins

The initial stage involves choosing the preferred programming language. While the adapter library is coded in Fortran 2003, it offers interfaces for incorporating plugins developed in C/C++ and Python. 
As an illustration, provided here is a guide on creating a plugin using Fortran. Each plugin must have three parts:

* Primary constructor
* Secondary constructor
* Callback function(s)

### Primary constructor

The primary constructor registers the plugin. It basically contains the following steps:

- It is explicitly checked that the major **versions of the ComIn library** that is used by the ICON and the one that is used by the plugin match. 
```fortran
version = comin_setup_get_version()
```
As many components of the development are still in the testing phase, the initial public release is set to version number 0.1.0. 
There is a **finish subroutine** which can be called in different occasions. For example here it can be used to check if the component ``version_no_major`` in the data structure ``version`` is ``0``.
```fortran
IF (version%version_no_major > 1) THEN
	CALL comin_plugin_finish("comin_main (simple_fortran_plugin)", "incompatible version!")
END IF
```

- There is a library function for getting the ComIn-internal ID that is used to identify a specific plugin during the subsequent operations.
```fortran
CALL comin_current_get_plugin_info(this_plugin, ierr)
```
Here, `COMIN_SUCCESS` is a defined constant. If `ierr == COMIN_SUCCESS`, it indicates a successful outcome.

#### Registering additional variables.  

A list of to-be-created variables made known to the ICON via the function `comin_var_request_add`.
```fortran
CALL comin_var_request_add(var_descriptor, lmodexclusive, ierr)
```
	
- `var_descriptor` is required to describe (and uniquely identify) a model variable in ICON.  Encapsulation of this information into a (constant) data structure of the data type is necessary for two reasons: a) iterating over the list of available variable is simplified, and b) future extensions, e.g. to lat-lon variables, are possible without changing 3rd party code.  
```fortran
var_descriptor=t_comin_var_descriptor( id = domain_id, name = "variable_name")
```

- Flag `lmodexclusive`: Whenever a plugin calls `comin_var_request_add`, there is a check to determine if the requested variable is already registered. If the variable exists and is either exclusively requested in the current call or was exclusively requested before, the model aborts based on the `lmodexclusive` setting. Alternatively, if the variable is not registered or the exclusive conditions are not met, a new variable with the specified properties is added to the list of requested variables.

- Variables may also be appended to ICON's container of tracer variables through the `tracer` flag (part of the metadata). Apart from that aspect it is not possible to create additional variable containers via the adapter library. It cannot be assumed (if only because of the "sharing" of variables between multiple ComIn plugins) that the tracers generated by a module are stored consecutively.
```fortran
CALL comin_metadata_set(var_descriptor, "tracer", .TRUE., ierr)
```
 While it is possible to create variables only for certain domains, ICON has the restriction that tracer variables have to be present on _every_ domain. For this reason, it is necessary to choose domain id `-1` (meaning all domains) as part of the `var_descriptor` for variables with `tracer = .true.`
	
- Newly created fields can be added to ICON's set of restart variables.
```fortran
CALL comin_metadata_set(var_descriptor, "restart", .TRUE., ierr)
```

#### Registering callbacks.  

The primary constructor appends subroutines of the 3rd party module to the callback register via the adapter library subroutine `comin_callback_register`.
```fortran
CALL comin_callback_register(entry_point_id, fct_ptr, ierr)
```
	
- `entry_point_id`: entry points denote events during the ICON model simulation, which can trigger a subroutine call of the plugin. Entry points are denoted by named integer constants.
  The table of available entry points is available in the [technical documentation](icon_comin_doc.md). 
- `fct_ptr`: this is the callback function.

#### Getting descriptive data structures.  

The descriptive data structures contain information on the ICON setup (e.g. Fortran `KIND` values), the computational grid(s), and the simulation status.
All descriptive data structures are treated as read-only (seen from the perspective of the 3rd party plugins). However, this read-only nature is (currently) not enforced. For efficiency reasons, the adapter library directly uses pointers to ICON data structures where possible. This holds mostly for components of `p_patch`, while non `p_patch` descriptive data are copied from the host model.
- Global data is available for the plugins primary constructor and all subsequent subroutine callbacks. Global data is never changed or updated and invariant w.r.t. the computational grid (logical domain ID). The detailed table of the components of global data can be found in the [technical documentation](icon_comin_doc.md).
```fortran
TYPE(t_comin_descrdata_global), POINTER :: p_global
p_global => comin_descrdata_get_global()
```
- Grid information is available for the 3rd party module's primary constructor and all subsequent subroutine callbacks. Grid information is never changed or updated. The data structures in this section are replicated for each computational domain (logical domain ID).
```fortran
TYPE(t_comin_descrdata_domain), POINTER :: p_patch
 p_patch => comin_descrdata_get_domain(jg)
```
- Timing information on the simulation.
```fortran
TYPE(t_comin_descrdata_simulation_interval), POINTER :: p_simulation_interval
p_simulation_interval => comin_descrdata_get_simulation_interval()
```
- Time step length per domain.  
```fortran
dtime=comin_descrdara_get_timesteplength()
```

### Secondary constructor

A secondary constructor is called _after_ the allocation of ICON variable lists and fields and _before_ the time loop.
* Access to ICON data fields happens via an accessor function `comin_var_get`.
Basically, `comin_var_get(context, var_descriptor, flag, var_pointer)` returns a 5-dimensional `REAL(wp)` pointer `var_pointer`. A return value `var_pointer /= NULL` means "success".
	* `context`: the name of the entry point.
	* `var_descriptor`: same as described in primary constructor part.
	* `flag`: the optional argument `flag` provides information w.r.t. the data flow. Flags may be combined like `flag = IOR(COMIN_FLAG_READ, COMIN_FLAG_WRITE)`.  
  It is important to highlight that when the `comin_var_request_add` procedure is executed, a variable is not immediately created. This step only involves the registration of a new variable. To use this variable later, it must be queried, similar to the other variables, using the `comin_var_get` function with `flag=COMIN_FLAG_WRITE`.
	* `comin_var_get` registers the access to a variable and returns the variable handle.

Code example:
```fortran
TYPE(t_comin_var_ptr), POINTER :: p	
CALL comin_var_get(context, var_descriptor, flag, p)
```

Convenience function `comin_var_to_3d` for accessing 2D/3D fields: In practice, access to fields can be simplified, under the condition that the sequence of dimensions is `(jc,jk,jb)`. This exact dimension sequence is (currently) fulfilled by the ICON model. In this case, a 3D pointer variable `REAL(wp) :: slice(:,:,:)` can be generated directly from a variable of type `TYPE(t_comin_var_ptr)` using the function.
```fortran
tracer_slice => comin_var_to_3d()
```

### Callback function

The plugin allows users to write subroutines that can be called at predefined events (entry points) throughout the model simulation.

## Python ComIn plugins

Python plugins can be attached to ComIn via the Python adapter which is located in the `plugins` directory in the ComIn source code. It is compiled with ComIn if `COMIN_ENABLE_PYTHON_ADAPTER` is enabled in the CMake configuration, see the build instructions above. The Python adapter embeds a Python interpreter, which also has the `comin` Python module available. This module contains all the functions, variables, constants and data structures of the [Python language API](comin_python_api.md). When including the Python adapter in the namelist, the Python plugin script must be specified as the `options` which can be modified while the actual ComIn plugin Python adapter (`libpython_adapter.so`) remains unchanged. This script is executed in the primary constructor of the Python adapter. Further callbacks can then be registered by the `comin.register_callback` function decorator.

```bash
 plugin_list(2)%name           = "simple_python_plugin"
 plugin_list(2)%plugin_library = "libpython_adapter.so"
 plugin_list(2)%options        = "${basedir}/externals/comin/plugins/python_adapter/examples/simple_python_plugin.py"
```

##  Build  ComIn plugin

For building a ComIn plugin we recommend to use [CMake](https://www.cmake.org). In the first step you should create a separate CMake project and place your plugin there. 

In the next step, one must build ComIn. We strongly recommend the out-of-source build (the instruction can be found in the next section). Following that, ComIn offers a CMake config (`ComInConfig.cmake`) such that it can be easily found in your CMake project.  Then, establish a connection between your CMake project and ComIn.
```bash
cd your_project
export ComIn_DIR=path_to_the_comin_build_directory
```
In the next step,  generate a `CMakeLists.txt` in to your CMake project with the following lines:
```
project(name_of_your_project LANGUAGES Fortran)
find_package(ComIn)
add_library(your_plugin MODULE your_plugin.F90)
target_link_libraries(your_plugin ComIn::ComIn)
```
*Note*: In the example above is assumed you want to build a Fortran plugin. In case of C plugin, the `LANGUAGES` does not need to be specified. 
Afterwards, you can create a build directory and build your plugin:
```bash
mkdir build && cd build
cmake ..
make
```

## Using ComIn's testing mechanism

ComIn offers functionality to test your plugin with the `minimal_example` host model emulator using CTest. In particular this can be used in a CI/CD setup for validating that the plugin builds and can be executed. The `minimal_example` does not provide any physical meaningful data. (A 'replay' functionality where that values of variables and metadata are read from a file is planned for the future.)

To add and configure tests in your projects ComIn provides [utility functions](cmake.md).

To add a test, you can use the `comin_add_test` CMake function in `CMakeLists.txt`.

```
comin_add_test(TEST your_test)
```

This generates a CTest test with the name `your_test` and sets up everything to run the `minimal_example`. To add a plugin to the test use the function `comin_test_add_plugin` in `CMakeLists.txt` .

```
comin_test_add_plugin(TEST your_test
  NAME "your_plugin"
  PLUGIN_LIBRARY $<TARGET_FILE:your_plugin>
  PRIMARY_CONSTRUCTOR "your_plugin_main"
  OPTIONS "some options"
  COMM "your_comm")
```

The parameters correspond to the parameters in the namelist (`t_comin_plugin_description`) for configuring a plugin.

The tests can be executed with `ctest` or `make test`. Note that the CMake variable `BUILD_TESTING` must be set to `ON` to build the tests.

## ComIn-standalone setup on Levante (DKRZ)  

It is also possible to build ComIn without ICON as a host model and test plugins using the standalone emulator (`minimal_example`) distributed with ComIn. 

The following command loads all `spack` packages required. 
```bash  
spack load /jlxcfzu
spack load py-mpi4py
```  
*Actually it loads `netcdf-fortran`, but it depends on all required packages, so they will be loaded as well. You may verify this by the command `spack find --loaded`.*  
  
The necessary modules are loaded with the following command:  
```bash  
module load gcc/11.2.0-gcc-11.2.0 netcdf-c/4.8.1-gcc-11.2.0 netcdf-fortran/4.5.3-gcc-11.2.0
export MPI_ROOT='/sw/spack-levante/openmpi-4.1.2-mnmady'
```  
  
Clone the ComIn git repository. 
```bash  
module load git  
git clone git@gitlab.dkrz.de:icon-comin/comin.git
```  
  
Then follow the standard CMake workflow: create a build directory, configure and build.  
```bash  
mkdir comin/build  
cd comin/build
cmake -DCOMIN_ENABLE_EXAMPLES=ON -DCOMIN_ENABLE_PYTHON_ADAPTER=ON -DBUILD_TESTING=ON -DCMAKE_C_COMPILER="${MPI_ROOT}/bin/mpicc" -DCMAKE_Fortran_COMPILER="${MPI_ROOT}/bin/mpif90" ..
  
make -j6
```  
The above command line enables in the `CMake` call  
- the Python adapter by setting `-DCOMIN_ENABLE_PYTHON_ADAPTER=ON`,  
- the building of the standalone NWP emulator `minimal_example` by setting `-DCOMIN_ENABLE_EXAMPLES=ON`,  
- the generation of CI/CD tests (`ctest` command) by setting `-DBUILD_TESTING=ON`.  
  
Besides, for debugging purposes the `CMake` build option `VERBOSE=1` might be useful.  
  
You can run the tests in the build directory with  
```bash  
ctest  
```  
  
If the parallel tests fail, you might need to add the environment variables described [here](https://docs.dkrz.de/doc/levante/running-jobs/runtime-settings.html#openmpi). 
  
```bash  
export OMPI_MCA_osc="ucx"  
export OMPI_MCA_pml="ucx"  
export OMPI_MCA_btl="self"  
export UCX_HANDLE_ERRORS="bt"  
export OMPI_MCA_pml_ucx_opal_mem_hooks=1
```
