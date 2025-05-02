# Testing in ComIn
In this subdirectory the tests are defined. To make the definition of
tests as easy as possible, we provide the cmake function
`comin_add_test` that can be found in
`cmake/Scripts/CominAddTest.cmake`.  This function can also be used by
downstream projects if they use cmake and have called `find_package(ComIn)`.

## Executing the tests
To use the test facilities of ComIn it needs to be configured with
`BUILD_TESTING=ON`. This can be achieved for example by executing
```bash
cmake -DBUILD_TESTING=ON .
```
in your build directory. The `make` command then should build all the
necessary things like the `minimal-example`, which is a dummy host code using ComIn.
To run the tests just call
```bash
ctest
```
in your build directory.

If one or more tests fail, you can add the `--output-on-failure` flag to `ctest`.

## Add a new test
To add a test you need to provide a template for the `master.nml`,
usually it has the ending `nml.in`. This is passed to `comin_add_test`
in the named argument `NAMELIST_TEMPLATE`. An example can be found in
[test/simple_c/simple_c.nml.in](test/simple_c/simple_c.nml.in).

## External processes
If your test needs external processes to run you can add these with
the named argument `EXTERNAL_PROCESSES`. You can pass multiple
executables. Every executable is added to the `mpirun` command with
one process (`-n 1`).

## Reference output (optional)
If you have defined your test and it passes, you can add the output as
reference to the test. This ensures that the output does not change
accidentally in the future. To do this add the named argument
`OUTPUT_REFERNCE` to the `comin_add_test` command and pass a path for
placing the reference output. Then call `make update_test_references`
to copy the output directory to that reference path. Make sure that
there are no absolute paths in the reference output. This would lead
to failures on other systems.

Examples can be found in [test/CMakeLists.txt](test/CMakeLists.txt).

# Specified tests
## `simple_fortran`
This test tests the example plugin `simple_fortran_plugin` that can be found in the `plugins/simple_fortran` subdirectory of this project. It tests
- the basic Fortran infrastructure for defining ComIn plugins in Fortran
- registering callbacks
- querying some basic descriptive data structures
- requesting to add variables
- accessing variables via `var_get`
## `simple_c`
This test tests the example plugin `simple_c_plugin` that can be found in the `plugins/simple_c` subdirectory of this project. It tests
- the basic C infrastructure for defining ComIn plugins in C
- registering callbacks
- querying some basic descriptive data structures
- requesting to add variables
- accessing variables via `var_get`
## `mpi_communicator`
The test tests the host and plugin MPI communicators with an external process.
## `multi_plugin`
This test simply tests if multiple plugins can be defined in the master.nml and if they run correctly.
## `parallel`
This test tests the parallel execution of the minimal_example with a simple plugin attached.
