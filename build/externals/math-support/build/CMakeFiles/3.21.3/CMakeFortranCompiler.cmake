set(CMAKE_Fortran_COMPILER "/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/bin/gfortran")
set(CMAKE_Fortran_COMPILER_ARG1 "")
set(CMAKE_Fortran_COMPILER_ID "GNU")
set(CMAKE_Fortran_COMPILER_VERSION "10.2.0")
set(CMAKE_Fortran_COMPILER_WRAPPER "")
set(CMAKE_Fortran_PLATFORM_ID "")
set(CMAKE_Fortran_SIMULATE_ID "")
set(CMAKE_Fortran_SIMULATE_VERSION "")




set(CMAKE_AR "/usr/bin/ar")
set(CMAKE_Fortran_COMPILER_AR "/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/bin/gcc-ar")
set(CMAKE_RANLIB "/usr/bin/ranlib")
set(CMAKE_Fortran_COMPILER_RANLIB "/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/bin/gcc-ranlib")
set(CMAKE_COMPILER_IS_GNUG77 1)
set(CMAKE_Fortran_COMPILER_LOADED 1)
set(CMAKE_Fortran_COMPILER_WORKS TRUE)
set(CMAKE_Fortran_ABI_COMPILED TRUE)
set(CMAKE_COMPILER_IS_MINGW )
set(CMAKE_COMPILER_IS_CYGWIN )
if(CMAKE_COMPILER_IS_CYGWIN)
  set(CYGWIN 1)
  set(UNIX 1)
endif()

set(CMAKE_Fortran_COMPILER_ENV_VAR "FC")

set(CMAKE_Fortran_COMPILER_SUPPORTS_F90 1)

if(CMAKE_COMPILER_IS_MINGW)
  set(MINGW 1)
endif()
set(CMAKE_Fortran_COMPILER_ID_RUN 1)
set(CMAKE_Fortran_SOURCE_FILE_EXTENSIONS f;F;fpp;FPP;f77;F77;f90;F90;for;For;FOR;f95;F95)
set(CMAKE_Fortran_IGNORE_EXTENSIONS h;H;o;O;obj;OBJ;def;DEF;rc;RC)
set(CMAKE_Fortran_LINKER_PREFERENCE 20)
if(UNIX)
  set(CMAKE_Fortran_OUTPUT_EXTENSION .o)
else()
  set(CMAKE_Fortran_OUTPUT_EXTENSION .obj)
endif()

# Save compiler ABI information.
set(CMAKE_Fortran_SIZEOF_DATA_PTR "8")
set(CMAKE_Fortran_COMPILER_ABI "")
set(CMAKE_Fortran_LIBRARY_ARCHITECTURE "")

if(CMAKE_Fortran_SIZEOF_DATA_PTR AND NOT CMAKE_SIZEOF_VOID_P)
  set(CMAKE_SIZEOF_VOID_P "${CMAKE_Fortran_SIZEOF_DATA_PTR}")
endif()

if(CMAKE_Fortran_COMPILER_ABI)
  set(CMAKE_INTERNAL_PLATFORM_ABI "${CMAKE_Fortran_COMPILER_ABI}")
endif()

if(CMAKE_Fortran_LIBRARY_ARCHITECTURE)
  set(CMAKE_LIBRARY_ARCHITECTURE "")
endif()





set(CMAKE_Fortran_IMPLICIT_INCLUDE_DIRECTORIES "/scratch/wbenjami/dafy/sdfg_gpu/spack/opt/spack/linux-centos8-zen2/gcc-10.2.0/netcdf-fortran-4.6.0-p7jdeag2groulzasdfqnokx2qiris3zh/include;/scratch/wbenjami/dafy/sdfg_gpu/spack/opt/spack/linux-centos8-zen2/gcc-10.2.0/netcdf-c-4.9.2-zynxjqlltpkunevnlaltrm7kfqlwjgpu/include;/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/lib/gcc/x86_64-pc-linux-gnu/10.2.0/finclude;/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/lib/gcc/x86_64-pc-linux-gnu/10.2.0/include;/usr/local/include;/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/include;/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/lib/gcc/x86_64-pc-linux-gnu/10.2.0/include-fixed;/usr/include")
set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "gfortran;m;gcc_s;gcc;quadmath;m;gcc_s;gcc;c;gcc_s;gcc")
set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "/scratch/pmazumde/gitspace/icon-artifacts/velocity/integrate_meh;/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/lib/gcc/x86_64-pc-linux-gnu/10.2.0;/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/lib64;/lib64;/usr/lib64;/scratch/wbenjami/dafy/sdfg_gpu/spack/opt/spack/linux-centos8-zen2/gcc-10.2.0/netcdf-fortran-4.6.0-p7jdeag2groulzasdfqnokx2qiris3zh/lib;/scratch/wbenjami/dafy/sdfg_gpu/spack/opt/spack/linux-centos8-zen2/gcc-10.2.0/netcdf-c-4.9.2-zynxjqlltpkunevnlaltrm7kfqlwjgpu/lib;/scratch/wbenjami/dafy/sdfg_gpu/spack/opt/spack/linux-centos8-zen2/gcc-10.2.0/openblas-0.3.23-ara3inpfwh3myjbgqpaiqye5qc2b6wnc/lib;/apps/ault/spack/opt/spack/linux-centos8-zen/gcc-8.4.1/gcc-10.2.0-fqf4ze2mfclmx6e7ehrjckznpoqnevmi/lib")
set(CMAKE_Fortran_IMPLICIT_LINK_FRAMEWORK_DIRECTORIES "")
