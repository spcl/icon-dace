<!--
This file is written using Markdown language, which might make it difficult to
read it in a plain text editor. Please, visit ICON project page on DKRZ GitLab
(https://gitlab.dkrz.de/icon/icon/-/tree/master/config/cscs) to see this file
rendered or use a Markdown viewer of your choice
(https://www.google.com/search?q=markdown+viewer).
-->

This document aims to cover some peculiarities of
[dynamic linking](../../README.md#dynamic-libraries) on
[Piz Daint](https://www.cscs.ch/computers/piz-daint/).

[Piz Daint](https://www.cscs.ch/computers/piz-daint/) is a Cray machine. Cray
environments usually imply that the required compiler modules are set both at
the build and run times. However, it is still possible to build an executable in
a way that it does not require environment modifications before running it, i.e.
we have to neither set `LD_LIBRARY_PATH` variable nor load additional modules.
This might be not so straightforward though.

Let's consider two cases: in the first case we build with Intel compiler, in the
second one, we build with PGI compiler. In both cases, we build with the Cray's
Fortran wrapper `ftn` and link to Cray's `MPICH` libraries
(`module load cray-mpich`).

1. When we build with Intel compiler (`module load PrgEnv-intel`), the compiler
wrapper `ftn` links to `libmpich_intel.so.3` and `libmpichf90_intel.so.3` (in
the designated order), which reside in `/opt/cray/pe/lib64`, and adds a linker
flag that sets `-rpath` to the directory with the Intel runtime libraries,
e.g. `/opt/intel/compilers_and_libraries_2019.1.144/linux/compiler/lib/intel64`.
By default, `-rpath` is added as `RUNPATH` to our executable, which means that
it will be used only for location of the immediate dependencies of our
executable. Now, we have the following dependencies (simplified example):
    ```mermaid
    graph LR
    A(executable) --> B(libmpich_intel.so.3)
    A --> C(libmpichf90_intel.so.3)
    C --> B
    B --> D(libifport.so.5)
    C --> D
    ```
    where `libifport.so.5` is Intel runtime library residing in
`/opt/intel/compilers_and_libraries_2019.1.144/linux/compiler/lib/intel64`. Note
that our executable does not need this runtime library directly. When we run
the executable, the dynamic linker will try to resolve the dependencies. First,
it will check whether our executable contains `RUNPATH` entries. Since it does,
the dependencies will be resolved in the *new way*. This means that since
`libifport.so.5` is not an immediate dependency of the executable, the linker
will not try to find it in the `RUNPATH` entries of the executable.
Additionally, none of the `RUNPATH` entries of `libmpichf90_intel.so.3` is point
to the right directory and `libmpich_intel.so.3` does not have either `RUNPATH`
or `RPATH` entries at all. Therefore, the dynamic linker will not find
`libifport.so.5` and produce the following error:
    ```
    error while loading shared libraries: libifport.so.5: cannot open shared object file: No such file or directory
    ```
    There are two ways to handle this:
    *  Make `libifport.so.5` an immediate dependency of the executable (i.e.
*overlink*) by prepending `-lifport` to `LIBS`. This way, when the dynamic
linker figures out that it needs `libifport.so.5` for `libmpich_intel.so.3`, the
former will already be in the linker's registry and the dependency will be
resolved successfully. Although overlinking is generally a bad idea (e.g. see
[Better understanding Linux secondary dependencies solving with examples](http://www.kaizou.org/2015/01/linux-libraries.html))
it's not such a big problem in our case since we don't install our executable
and when the libraries are updated we can simply rebuild it. However, our
executable might indirectly depend on other runtime libraries. For example,
ICON might indirectly depend on other Intel runtime libraries (e.g.
`libifcoremt.so.5`, `libifcore.so.5`, `libintlc.so.5`, `libirng.so`, and
`libsvml.so`) and we would need to overlink to them too, i.e. additionally
prepend `-lifcore -lintlc -lirng -lsvml` to `LIBS`.
    * We can tell the linker to add `-rpath` as `RPATH` to our executable. This
will force the dynamic linker to resolve the dependencies in the *old way*. In
particular, this means that the `RPATH` entry will be used for location of the
secondary (a.k.a transitive, a.k.a non-immediate) dependencies. This can be done
by adding `-Wl,--disable-new-dtags` to `LIBS`. This is easier to maintain but
this works only because `libmpich_intel.so.3` and `libmpichf90_intel.so.3`
depend on the same set of runtime libraries, `libmpich_intel.so.3` does not have
`RUNPATH` entries and `ftn` links to `libmpich_intel.so.3` before
`libmpichf90_intel.so.3`. Basically, since `libmpich_intel.so.3` does not have
`RUNPATH` entries, the dynamic linker resolves its dependencies using `RPATH` of
the executable. When the linker starts resolving the dependencies of
`libmpichf90_intel.so.3` it finds a `RUNPATH` entry in it and switches to the
*new way*. Luckily, all the dependencies of `libmpichf90_intel.so.3` are already
resolved at this point (since they are the same as for `libmpich_intel.so.3`)
and we don't get the error.

2. It is not so relevant for the following example but it should be noted that
in contrast to Intel, in the case of PGI, it is not the Cray's wrapper `ftn`
that adds `-rpath` flags pointing to the directory with the compiler's runtime
libraries (e.g. `/opt/pgi/19.7.0/linux86-64/19.7/lib`) but the compiler itself.
Therefore each library built with PGI will likely have `RUNPATH` entries and the
approach with `-Wl,--disable-new-dtags` will not work. Unless, of course,
`RUNPATH` is set properly, which, unfortunately, is not the case for
`libmpich_pgi.so.3`: its `RUNPATH` entry points to a non-existing directory
`/notbackedup/users/sko/opt_pgi/19.1.0/linux86-64-llvm/19.1/lib`. Now, if we
compile the same executable with PGI (`module load PrgEnv-pgi`), the simplified
version of our dependency tree looks like this:
    ```mermaid
    graph LR
    A(executable) --> B(libmpich_pgi.so.3)
    A --> C(libmpichf90_pgi.so.3)
    C --> B
    C --> D(libomptarget.so)
    ```
    where `libomptarget.so` is PGI runtime library residing in
`/opt/pgi/19.7.0/linux86-64/19.7/lib`. Note that `libomptarget.so` is a
dependency only for `libmpichf90_pgi.so.3`, which has a `RUNPATH` entry and the
value of the entry is invalid. As a result, even if we link our executable with
`-Wl,--disable-new-dtags`, the dynamic linker will not use its `RPATH`:
    ```
    error while loading shared libraries: libomptarget.so: cannot open shared object file: No such file or directory
    ```
    The only solution here is to overlink to `libomptarget.so` by prepending
`-lomptarget` to `LIBS`.


> **_NOTE:_** In both cases, the dynamic linker manages to locate `MPICH`
libraries because the directory `/opt/cray/pe/lib64`, where they reside, is
listed in `/etc/ld.so.conf.d/cray-pe.conf`.