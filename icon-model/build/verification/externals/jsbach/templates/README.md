# Adding a new process

This short guide is mostly copied from the jsbach wiki and can be found [here](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/processes#how-to-register-a-process).

## Creating a new process
Scientific processes are usually capsuled within an own sub-folder within the ICON-L [code structure](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/code-structure). For example: `hydrology`, `phenology`, `srf_energy_bal`, ... . A process is usually realised by a couple of obligatory and optional [process modules](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/processes#process-modules) (depending on the type of process and process requirements). A [usecase](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/usecases) specifies which process is running on which tiles and thereby determines the [process actions](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/processes#process-actions), which further can be modified (processes switched on and off) via the namelist.

When implementing a new process in ICON-LAND, mandatory [process modules](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/processes#process-modules) need to be added, and the process needs to be made known in several modules of the ICON-L infrastructure (see [How to "register" a process](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/processes#how-to-register-a-process)). The templates for the mandatory modules are found in the `templates` directory.

## How to register a new process
In order to use a process within ICON-LAND, the process needs to be made known to several modules - depending on the type of the process and process requirements.

- All processes need a process ID which needs to be added to the according enumerator in [mo_jsb_process_class](https://gitlab.dkrz.de/jsbach/jsbach/-/blob/dev/src/base/mo_jsb_process_class.f90). In this module, also a mapping onto a process name and back needs to be specified.
- The process itself needs to be specified in [mo_jsb_model_usecases](https://gitlab.dkrz.de/jsbach/jsbach/-/blob/dev/src/base/mo_jsb_model_usecases.f90) for the [usecases](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/usecases) for which it should be calculated.
- The tasks of a process need to be added to the [task queue](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/task-queue) in [mo_jsb_model_init](https://gitlab.dkrz.de/jsbach/jsbach/-/blob/dev/src/base/mo_jsb_model_init.f90) [Note: mind task order requirements].
- The interface and the config class of a process need to be added to [mo_jsb_process_factory](https://gitlab.dkrz.de/jsbach/jsbach/-/blob/dev/src/base/mo_jsb_process_factory.f90). Don't forget to update the `max_no_of_processes` variable accordingly.
- If a process has an init module, the init routine needs to be called in `jsbach_init` of  [mo_jsb_model_init](https://gitlab.dkrz.de/jsbach/jsbach/-/blob/dev/src/base/mo_jsb_model_init.f90).
- If a process has its own memory the memory module also needs to be added to [mo_jsb_process_factory](https://gitlab.dkrz.de/jsbach/jsbach/-/blob/dev/src/base/mo_jsb_process_factory.f90) - if a process **does not have an own memory** the `has_memory` flag needs to be set to `.FALSE.` for this process in the `Create_process` function of the mo_jsb_process_factory.
- In the special case of an land cover change process:
    * the init lcc routine needs to be called in `jsbach_init` of  [mo_jsb_model_init](https://gitlab.dkrz.de/jsbach/jsbach/-/blob/dev/src/base/mo_jsb_model_init.f90).
    * the update routine of a task of an lcc process needs to call certain lcc infrastructure routines of [mo_jsb_lcc](https://gitlab.dkrz.de/jsbach/jsbach/-/blob/dev/src/base/mo_jsb_lcc.f90) - see [LCC infrastructure](https://gitlab.dkrz.de/jsbach/jsbach/-/wikis/Documentation/ICON-Land-framework/lcc-infrastructure).

Tip: usually it is helpful to identify an existing process with a similar type and purpose and to search how this is used in the different ICON-LAND infrastructure modules.

