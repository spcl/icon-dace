#!/usr/bin/bash

set -e -o pipefail

mkdir build
cd build

export CFLAGS="-Wall -Wunused-parameter"
export CXXFLAGS="-Wall -Wunused-parameter -pedantic"
export FFLAGS="-Wall -Wno-maybe-uninitialized -Wl,--no-warn-execstack -std=f2008 -fcheck=bounds,do,mem,pointer,recursion -fmodule-private -fimplicit-none -fmax-identifier-length=63"

cmake -DCOMIN_ENABLE_PYTHON_ADAPTER=ON \
      -DBUILD_TESTING=ON \
      -DCOMIN_ENABLE_STATIC_LINKING_TEST=ON \
      -DCOMIN_ENABLE_REPLAY_TOOL=ON \
      ..

make -j8  2>&1 | tee build.log
