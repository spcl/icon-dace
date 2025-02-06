// Copyright (c) 2013-2024 MPI-M, Luis Kornblueh, Rahul Sinha and DWD, Florian Prill. All rights reserved.
//
// SPDX-License-Identifier: BSD-3-Clause
//
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef ENABLE_CHECK

#include <check.h>
#include <stdlib.h>

#include "mtime_calendar_test.h"
#include "mtime_date_test.h"
#include "mtime_time_test.h"
#include "mtime_datetime_test.h"
#include "mtime_julianDay_test.h"
#include "mtime_timedelta_test.h"

int
main(void)
{
  Suite *suite = suite_create("libmtime-check");
  SRunner *sr = srunner_create(suite);

  add_mtime_calendar_test_to_suite(suite);
  add_mtime_date_test_to_suite(suite);
  add_mtime_time_test_to_suite(suite);
  add_mtime_datetime_test_to_suite(suite);
  add_mtime_julianDay_test_to_suite(suite);
  add_mtime_timedelta_test_to_suite(suite);

  srunner_run_all(sr, CK_ENV);

  int number_failed = srunner_ntests_failed(sr);
  srunner_free(sr);

  return number_failed == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

#else

int
main(void)
{
  return 77;
}

#endif  // ENABLE_CHECK
