# Copyright (c) 2018-2024, MPI-M
#
# Author: Sergey Kosukhin <sergey.kosukhin@mpimet.mpg.de>
#
# SPDX-License-Identifier: BSD-3-Clause
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
#    this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its
#    contributors may be used to endorse or promote products derived from
#    this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

# NAG Fortran Compiler emits the output of its preprocessor to the standard
# output when called with the following flags:
#   1) '-fpp' enables the preprocessing;
#   2) '-F' tells the compiler to stop after the preprocessing;
#   3) '-o -' tells the compiler to emit the output to stdout;
#   4) '-Wp,-w' prevents 'fpp', which is internally used by the
#      compiler, from emitting warning messages, including the
#      one that says "fpp: warning: bad option: - , ignored").
# 
# However, the compiler might not be able to compile the output of its own
# preprocessor. For example:
# 
#   $ cat ./conftest.f90
#   program main
#     print *, ("A very very very very very very very very very very very very long line of code. We need more than 132 characters on this line. &
#   &Finally!")
#   end program
#   $ nagfor -o - -Wp,-w -F -fpp ./conftest.f90 > ./conftest_preprocessed.f90
#   NAG Fortran Compiler Release 6.2(Chiyoda) Build 6223
#   $ nagfor -fpp ./conftest_preprocessed.f90
#   NAG Fortran Compiler Release 6.2(Chiyoda) Build 6223
#   Error: conftest.f90, line 5: Invalid continuation
#   Error: conftest.f90, line 5: Syntax error
#   Error: conftest.f90, line 3: Invalid continuation
#   [NAG Fortran Compiler pass 1 error termination, 3 errors]
#
# There are three known ways to solve this problem:
# 
#   1) compile the output of the preprocessing without the flag '-fpp', which
#      is apparently the way intended by the developers:
# 
#      $ nagfor ./conftest_preprocessed.f90
#      NAG Fortran Compiler Release 6.2(Chiyoda) Build 6223
#      [NAG Fortran Compiler normal termination]
# 
#      however, we might need to introduce additional preprocessor directives
#      to the preprocessed file before compiling it;
# 
#   2) disable the generation of line numbering directives by extending the
#      preprocessor command with the '-Wp,-P' flag:
#
#      $ nagfor -o - -Wp,-w -F -fpp -Wp,-P ./conftest.f90 > ./conftest_preprocessed.f90
#      NAG Fortran Compiler Release 6.2(Chiyoda) Build 6223
#      $ nagfor -fpp ./conftest_preprocessed.f90
#      NAG Fortran Compiler Release 6.2(Chiyoda) Build 6223
#      [NAG Fortran Compiler normal termination]
# 
#      however, we might want to keep the directives, e.g. in order to know the
#      files that were included by the preprocessor (there is flag '-Wp,-M'
#      that is meant for that but other compilers might have different flags
#      and/or different logic for this feature, which would require additional
#      customization of the dependency generation);
# 
#   3) clean the preprocessor output from the problematic directives, which
#      is implemented in this script.

# The first argument is supposed to be the compiler:
FC=$1; shift
# Set the flag enabling preprocessing:
pp_args='-fpp -F -o - -Wp,-w'

# Sed script replacing:
#
# banana&
# # 123
# &orange
# 
# with:
# 
# banana&
# &orange
# 
# The script should work with both GNU and BSD versions of sed:
sed_script='N;$q;s/&\n# [1-9][0-9]*/\&/;P;D'

# The following is a portable version of the following:
#
# #!/bin/bash
# set -e
# set -o pipefail
# $FC $pp_args "$@" | sed "$sed_script"
{ { { { $FC $pp_args "$@" || echo $? >&3; } | sed "$sed_script" >&4; echo $? >&3; } 3>&1; } | { read xs; exit $xs; } } 4>&1 || exit $?
