#!/bin/bash

# Copyright (c) 2024 The YAC Authors
#
# SPDX-License-Identifier: BSD-3-Clause

prev_log_name=

for template_name in \
  $(find . -maxdepth 1 -type f -a -name '*.sh.in' -exec grep -l '@MPI_LAUNCH@' {} \; | sort -r); do

  template_name="${template_name#'./'}"
  log_name="${template_name%'.sh.in'}.log"

  if test -n "${prev_log_name}"; then
    echo "${prev_log_name}: ${log_name}"
  fi

  prev_log_name=${log_name}

done
