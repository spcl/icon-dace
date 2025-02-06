// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <stdlib.h>
#include <stdio.h>

extern unsigned err_count__;

#define PUT_ERR(string) err_count__++, fputs((string "\n"), stderr);

#define INC_ERR (err_count__++)

#define TEST_EXIT_CODE ((!err_count__)?EXIT_SUCCESS:EXIT_FAILURE)
#define EXIT_SKIP_TEST (77)


