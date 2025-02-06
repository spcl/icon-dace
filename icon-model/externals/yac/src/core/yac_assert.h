// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YAC_ASSERT_H
#define YAC_ASSERT_H

// YAC PUBLIC HEADER START

void yac_abort_message(char const*, char const*, int);

#define die(msg) \
  yac_abort_message((msg), __FILE__, __LINE__)

#ifndef YAC_CODE_COVERAGE_TEST
#define YAC_ASSERT(exp, msg) \
  {if(!((exp))) yac_abort_message(((msg)), __FILE__, __LINE__);}

#define YAC_ASSERT_F(exp, format, ...) \
  { \
    if(!((exp))) { \
      char msg_buffer[1024]; \
      int ret = snprintf( \
        msg_buffer, sizeof(msg_buffer), ((format)), __VA_ARGS__); \
      if ((ret >= 0) && ((size_t)ret < sizeof(msg_buffer))) \
        yac_abort_message(((msg_buffer)), __FILE__, __LINE__); \
      else \
        yac_abort_message( \
          "an error occured, but error message could not be " \
          "generated", __FILE__, __LINE__); \
    } \
  }
#else
#define YAC_ASSERT(exp, msg) {(void)(exp);}
#define YAC_ASSERT_F(exp, format, ...) {(void)(exp);}
#endif

// YAC PUBLIC HEADER STOP

#endif // YAC_ASSERT_H

