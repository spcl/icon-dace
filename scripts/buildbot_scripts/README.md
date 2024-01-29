# Overview of how Buildbot builds ICON and runs the testcases
Builder and experiment definitions all reside in the ICON source code. While there are a few variants of these files, only the defaults are mentioned here. The file `create_all_builders` defines the builders, including the names of the build wrapper scripts as well as meta data describing things like the builder name, configure flags or the machine on which the builder operates. The file `create_list_icon-dev` extends this definition to add experiments to each builder and/or machine. Experiments can be added to multiple builders and or machines and can carry meta data like run flags and dependencies which are potentially different for different builders/machines. This constitutes the following configuration matrix where each "x" represents a distinct experiment object (defined in `experiment.py`) carrying meta data:

|                            | daint             |                   |                   | mistral         | ... |
| -                          | :-:               | :-:               | :-:               | :-:             | --- |
|                            | **DAINT_CPU_cce** | **DAINT_CPU_pgi** | **DAINT_GPU_pgi** | **MISTRAL_gcc** | ... |
| **atm_amip_test**          | x                 | x                 | x                 | x               | ... |
| **mch_ch_lowres**          | x                 | x                 | x                 |                 | ... |
| **test_ocean_omip_10days** |                   |                   |                   | x               | ... |
| ...

Each builder creates the full configuration using the `create_all_builders` and `create_list_icon-dev`. Subsequently, it will look up the appropriate configure wrapper script to build ICON. Finally, it will take the builder-specific subset of experiments that are to be run. This workflow is implemented in the master configuration residing on the buildbot server and thus *not* apparent from within ICON.  

> **Nomenclature:** What is called "*configuration matrix*" here is sometimes also referred to as "*experiment list*" of simply "*list*" due to its creation by the `create_list_<list>` scripts. The implementation of the "*configuration matrix*" is done by an object called `BuildbotConfig` and is serialized to a file with the name of the list (e.g. `icon-dev`).

# Interface functions
## Utility functions used by server code
The server code is kept as general as possible so that most of the buildbot infrastructure code can be maintained by the ICON community. The server will first create the configuration matrix from specified by the list (`create_list_<list>`) and report the generated list as a string to the buildbot logfile (`lsexperiments`). Then the build wrapper is invoked, runscripts are generated from templates (`build`) and experiments are submitted to the batch system and exit codes are checked (`runexp`). The following tables provides more details:
| Operation | Description
| -         | -
| create_list_\<list\> | The server calls `create_list_<list>` to set up the configuration matrix. The parameter `<list>` can be set from the buildbot website's UI `name:value` pair. By default, the `icon-dev` list is used. This means that the `create_list_<list>` scripts need to specify the full configuration matrix, including definitions for all machines, builders and experiments.
| lsexperiments | Buildbot writes a summary of the configuration matrix to stdout. The `lsexperiments` script should provide exactly this.
| runexp | The runexp script submits the experiments for a particular builder and list to the batch system of the builders machine. It then waits for all the jobs to complete and checks the returncodes to generates the `STATUS` file which is shown on the Buildbot web UI to inform the user if the experiments were `OK` or `FAILED`. In the current implementation of the server-side code, `runexp` is not aware of 1) the buildbot configuration and 2) what builder it is executed by. This information is passed from the `build` script by writing a `buildbot_config.json` file.
| build | The build script configures, builds and creates the ICON experiments from the templates in a single step. It also writes the name of the builder as well as the list name to  the `buildbot_config.json` file which is read by `runexp`.

## Utility functions provided to list maintainers
Lists are usually built by a two-step process. First, all builders and associated machines need to be defined. Second, experiments (including meta data) can be added for individual builders, machines, or combinations and subsets of both, controlled by builder meta data. For example, it is allowed to add an experiment to all builders with the `--with-mpi` flags. One can even specify a list of builders, a set of flags and machines and the experiment will be added to the smallest set of builders satisfying all requirements. The implemented list operations are:
| Operation | Description
| -         | -
| addexp    | Add an experiment to the list by builder, machine, matching or excluding configure flags. Addidional flags for the runscript generator can be passed.
| rmexp     | remove an experiment from the list. The same selection criteria as for `addexp` can be used.
| set_builder_flags | Builders can be in 3 states: `Active`, `Inactive` und `build_only`. This can be set either when adding the builder or using this method
| adddep | Add dependencies between experiments. Different modes are supported: 1) a single experiment depends on another single experiment, 2) a single experiment depends on multiple experiments, 3) multiple experiments depend on a single experiment. For which builders and machines this dependency applies is specified the same way as for `addexp`.

# Implementation
The buildbot infrastructure code is written in Python and bash. Python takes care of handling the buildbot configuration object as well as building and submitting the jobs. The list creation is written in bash for backwards compatibility with the previous infrastructure.

## Experiment
The experiment class stores experiment meta data:
- The name of the experiment
- The flags passed to the runscript generator
- The name of the associated builder
- The name of the batch script that submits the experiment
- Dependencies by keeping track of parent and child experiments
- The appropriate batch job implementation

Furthermore, it handles the creation of experiment runscripts from templates, submitting jobs to the batch system and experiment dependencies.

## Batch Job
The batch job class stores the following data:
- The submit command
- The directory from where the job should be submitted
- The ID of the submitted job
- The subprocess object associated with the job
- The returncode of the job
- A list of job ID's it depends on

Since different machines can have different batch systems, interaction is abstracted by a separate class providing functions to `submit` jobs.

## Builder
The builder class stores the following data:
- The name of the builder
- The name of the associated machine
- The name of the configure wrapper script to be used for building
- The configure flags (they have no effect on the build and are purely used for matching builders by configure flags. This is a legacy from the old build system still present in the `create_list_<list>` files and thus needs to be supported)
- The builder flag (`Active`, `Inactive`, `build_only`)

The builder class builds ICON.

## Configuration Matrix
With it's tabular layout, the configuration matrix is a perfect candidate to be handled by `pandas.DataFrame`. Subset creation is straight forward with the `DataFrame.loc` function. The custom class `BuildbotConfig` stores all the information needed about machines, builders and experiments. There are three main storage objects: 
- a python `dict` storing `Builder` objects by builder name
- a python `dict` storing machine meta data as a dictionary (machines so far only have a single meta data type `queue` which is not used in current `create_list_<list>` scripts or member functions)
- and the configuration matrix `DataFrame` which stores `Experiment` objects in a machine-builder-experiment `MultiIndex` `DataFrame` that looks like the one at the beginning of this document.

# Requirements
The implementation makes use of some non-standard python libraries:
| Library | Use
| - | - 
| click | Click is a widely used command line interface library for python and the successor of the standard `argparse`. It is used here to provide the interface to the bash scripts.
| pandas | Pandas is a very powerful data analytics library. We use it here only for it's nice utility functions for multi-dimensional tables
| numpy | Only used for `np.nan`. (Numpy is also a requirement of pandas.)
| pathlib | Facilitates filesystem interactions