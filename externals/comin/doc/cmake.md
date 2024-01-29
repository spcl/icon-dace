# ComIn CMake utils

ComIn provides cmake functions to setup tests with CTest. An example
how to use this functions can be found in the [User Guide](user_guide.md).

## comin_add_test
This is a wrapper around [`add_test`](https://cmake.org/cmake/help/latest/command/add_test.html)
that adds a test to the projects and takes care of the configuration of the `minimal_example`.
```cmake
comin_add_test(NAME <name>
               [NUM_PROCS <num_procs>]
               [REFERENCE_OUTPUT <dir>])
```
- `NAME`: name of the test (forwarded to `add_test`)
- `NUM_PROCS`: Number of processes for the minimal_example

## comin_test_add_plugin
Adds a plugin in the `comin_nml` of a given test. The arguments are forwrded to the `comin_nml`
```cmake
comin_test_add_plugin(TEST <test>
                      NAME <name>
                      [PLUGIN_LIBRARY <filename>]
                      [PRIMARY_CONSTRUCTOR <functionname>]
                      [OPTION <string>]
                      [COMM <string>]
                    )
```
- `TEST`: name of the test the plugin should be added to
- `NAME`: name of the plugin
- `PLUGIN_LIBRARY`: filename of the shared object of the plugin (if any)
- `PRIMARY_CONSTRUCTOR`: name of the primary constructor (default: `comin_main`)
- `OPTIONS`: a options string passed to the plugin (default: "")
- `COMM`: name of the plugin communicator (defualt: "", meaning no communicator is created for the plugin)

## comin_test_add_external_process
Add external process to the test. The processes are appended the mpirun command.
```cmake
comin_add_test(TEST <test>
               [NUM_PROCS <n>]
               COMMAND <command>
             )
```
- `TEST`: name of the test the plugin should be added to
- `NUM_PROCS`: Number of processes for the minimal_example
- `COMMAND`: command to be executed on the additional processes
