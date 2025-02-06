// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <string.h>

#include "utils_core.h"
#include "interpolation_exchange.h"

enum exchange_state {
  EXCHANGE_IDLE,     //! exchange is currently not in use
  EXCHANGE_WAIT_PUT, //! contains valid put data, but is not yet communicated
  EXCHANGE_WAIT_GET, //! contains valid get data, but is not yet communicated
  EXCHANGE_ACTIVE,   //! contains valid data, which is being communicated
};

enum empty_exchange_state {
  EXCHANGE_INVALID, //! exchange is source and/or target
  EXCHANGE_UNSET,   //! exchange is neither source nor target
                    //! and not put/get has been called
  EXCHANGE_AT_PUT,  //! exchange is neither source nor target
                    //! and exchange is executed at put
  EXCHANGE_AT_GET,  //! exchange is neither source nor target
                    //! and exchange is executed at get
};

struct yac_interpolation_exchange {

  char * name;

  enum exchange_state state;
  enum empty_exchange_state empty_state;
  Xt_redist redist;
  Xt_request request;
  size_t count;

  int is_source, is_target;
  void const ** temp_send_buffer;
  void **       temp_recv_buffer;
  void ** dummy_buffer;
};

static Xt_redist combine_redists(
  Xt_redist * redists, size_t num_redists, size_t collection_size) {

  // check redists
  if (redists != NULL) {
    int redists_are_null = redists[0] == NULL;
    for (size_t i = 1; i < num_redists; ++i)
      YAC_ASSERT(
        redists_are_null == (redists[i] == NULL),
        "ERROR(combine_redists): invalid argument")
    if (redists_are_null) return NULL;
  } else {
    return NULL;
  }

  if (collection_size == 1) {
    if (num_redists == 1) return xt_redist_copy(redists[0]);
    else
      return
        xt_redist_collection_new(
          redists, (int)(num_redists * collection_size), -1,
          xt_redist_get_MPI_Comm(redists[0]));

  } else {
    Xt_redist * redists_buffer =
      xmalloc(collection_size * num_redists * sizeof(*redists_buffer));
    for (size_t i = 0; i < collection_size; ++i)
      for (size_t j = 0; j < num_redists; ++j)
        redists_buffer[i * num_redists + j] = redists[j];
    Xt_redist redist =
      xt_redist_collection_new(
        redists_buffer, (int)(num_redists * collection_size), -1,
        xt_redist_get_MPI_Comm(redists[0]));
    free(redists_buffer);
    return redist;
  }
}

static void do_empty_exchange(
  struct yac_interpolation_exchange * exchange, int is_put,
  char const * routine_name) {

  YAC_ASSERT_F(
    exchange->empty_state != EXCHANGE_INVALID,
    "ERROR(%s): state of exchange \"%s\" is inconsistent",
    routine_name, exchange->name);

  if (exchange->empty_state == EXCHANGE_UNSET)
    exchange->empty_state =
      is_put?EXCHANGE_AT_PUT:EXCHANGE_AT_GET;

  if ((is_put && (exchange->empty_state == EXCHANGE_AT_PUT)) ||
      (!is_put && (exchange->empty_state == EXCHANGE_AT_GET)))
    xt_redist_s_exchange(
      exchange->redist, exchange->count,
      (const void **)(exchange->dummy_buffer),
      (void **)(exchange->dummy_buffer));
}

static struct yac_interpolation_exchange * yac_interpolation_exchange_new_(
  Xt_redist redist, size_t count, char const * name) {

  struct yac_interpolation_exchange * exchange = xmalloc(1 * sizeof(*exchange));

  exchange->name = strdup(name);

  exchange->redist = redist;
  exchange->request = XT_REQUEST_NULL;
  exchange->state = EXCHANGE_IDLE;
  exchange->count = count;

  exchange->is_source =
    (redist != NULL) && (xt_redist_get_num_send_msg(redist) > 0);
  exchange->is_target =
    (redist != NULL) && (xt_redist_get_num_recv_msg(redist) > 0);
  exchange->temp_send_buffer =
    xmalloc(exchange->count * sizeof(*(exchange->temp_send_buffer)));
  exchange->temp_recv_buffer =
    xmalloc(exchange->count * sizeof(*(exchange->temp_recv_buffer)));
  exchange->dummy_buffer =
    xcalloc(exchange->count, sizeof(*(exchange->dummy_buffer)));

  exchange->empty_state =
    (exchange->is_source || exchange->is_target)?
      EXCHANGE_INVALID:EXCHANGE_UNSET;

  return exchange;
}

struct yac_interpolation_exchange * yac_interpolation_exchange_new(
  Xt_redist * redists, size_t num_fields, size_t collection_size,
  int with_frac_mask, char const * name) {

  if (with_frac_mask) collection_size *= 2;

  Xt_redist redist =
    combine_redists(redists, num_fields, collection_size);

  return
    yac_interpolation_exchange_new_(
      redist, num_fields * collection_size, name);
}

struct yac_interpolation_exchange * yac_interpolation_exchange_copy(
  struct yac_interpolation_exchange * exchange) {

  return yac_interpolation_exchange_new_(
    (exchange->redist == NULL)?NULL:xt_redist_copy(exchange->redist),
    exchange->count, exchange->name);
}

int yac_interpolation_exchange_is_source(
  struct yac_interpolation_exchange * exchange) {

  return exchange->is_source;
}

int yac_interpolation_exchange_is_target(
  struct yac_interpolation_exchange * exchange) {

  return exchange->is_target;
}

void yac_interpolation_exchange_wait(
  struct yac_interpolation_exchange * exchange, char const * routine_name) {

  // ensure that the exchange is not in waiting state
  YAC_ASSERT_F(
    exchange->state != EXCHANGE_WAIT_PUT, "ERROR(%s): "
    "the \"%s\"-exchange is currently in wait state (are you missing a get?)",
    routine_name, exchange->name);
  YAC_ASSERT_F(
    exchange->state != EXCHANGE_WAIT_GET, "ERROR(%s): "
    "the \"%s\"-exchange is currently in wait state (are you missing a put?)",
    routine_name, exchange->name);

  // if a previous put has not yet been completed
  if (exchange->state == EXCHANGE_ACTIVE) {
    xt_request_wait(&(exchange->request));
    exchange->state = EXCHANGE_IDLE;
  }
}

int yac_interpolation_exchange_put_test(
  struct yac_interpolation_exchange * exchange, char const * routine_name) {

  YAC_ASSERT_F(
    (exchange->state == EXCHANGE_IDLE) ||
    (exchange->state == EXCHANGE_WAIT_PUT) ||
    (exchange->state == EXCHANGE_WAIT_GET) ||
    (exchange->state == EXCHANGE_ACTIVE),
    "ERROR(%s): invalid exchange state", routine_name);

  switch (exchange->state) {
    default:
    case (EXCHANGE_IDLE):
    case (EXCHANGE_WAIT_GET):
      return 1;
    case (EXCHANGE_WAIT_PUT):
      return 0;
    case (EXCHANGE_ACTIVE): {
      int flag;
      xt_request_test(&(exchange->request), &flag);
      if (flag) exchange->state = EXCHANGE_IDLE;
      return flag;
    }
  }
}

int yac_interpolation_exchange_get_test(
  struct yac_interpolation_exchange * exchange, char const * routine_name) {

  YAC_ASSERT_F(
    (exchange->state == EXCHANGE_IDLE) ||
    (exchange->state == EXCHANGE_WAIT_PUT) ||
    (exchange->state == EXCHANGE_WAIT_GET) ||
    (exchange->state == EXCHANGE_ACTIVE),
    "ERROR(%s): invalid exchange state", routine_name);

  switch (exchange->state) {
    default:
    case (EXCHANGE_IDLE):
    case (EXCHANGE_WAIT_PUT):
      return 1;
    case (EXCHANGE_WAIT_GET):
      return 0;
    case (EXCHANGE_ACTIVE): {
      int flag;
      xt_request_test(&(exchange->request), &flag);
      if (flag) exchange->state = EXCHANGE_IDLE;
      return flag;
    }
  }
}

enum YAC_INTERP_EXCH_STATUS yac_interpolation_exchange_status(
  struct yac_interpolation_exchange * exchange, char const * routine_name) {

  YAC_ASSERT_F(
    (exchange->state == EXCHANGE_IDLE) ||
    (exchange->state == EXCHANGE_WAIT_PUT) ||
    (exchange->state == EXCHANGE_WAIT_GET) ||
    (exchange->state == EXCHANGE_ACTIVE),
    "ERROR(%s): invalid exchange state", routine_name);

  switch (exchange->state) {
    default:
    case (EXCHANGE_IDLE):
      return YAC_INTERP_EXCH_IDLE;
    case (EXCHANGE_WAIT_PUT):
      return YAC_INTERP_EXCH_WAIT_PUT;
    case (EXCHANGE_WAIT_GET):
      return YAC_INTERP_EXCH_WAIT_GET;
    case (EXCHANGE_ACTIVE):
      return YAC_INTERP_EXCH_ACTIVE;
  }
}

void yac_interpolation_exchange_execute_put(
  struct yac_interpolation_exchange * exchange, double const ** send_data,
  char const * routine_name) {

  // if we have to do an exchange
  if (exchange->redist != NULL) {

    switch ((exchange->is_source << 1) +
            (exchange->is_target)) {

      // neither source nor target
      default:
      case(0): {

        // we should be in idle state
        YAC_ASSERT_F(
          exchange->state == EXCHANGE_IDLE,
          "ERROR(%s): the \"%s\"-exchange is not in idle state, "
          "internal error?", routine_name, exchange->name);
        do_empty_exchange(exchange, 1, routine_name);
        break;
      }
      // only target
      case(1): {

        // we should be either in idle or active state
        YAC_ASSERT_F(
          (exchange->state == EXCHANGE_IDLE) ||
          (exchange->state == EXCHANGE_ACTIVE),
          "ERROR(%s): the \"%s\"-exchange is invalid state, "
          "internal error?", routine_name, exchange->name);
        yac_interpolation_exchange_get_test(exchange, routine_name);
        break;
      }
      // only source
      case(2): {

        // if a previous exchange operation has not yet been completed
        if (exchange->state == EXCHANGE_ACTIVE) {
          xt_request_wait(&(exchange->request));
          exchange->state = EXCHANGE_IDLE;
        }

        // ensure that we are in idle state
        YAC_ASSERT_F(
          exchange->state == EXCHANGE_IDLE,
          "ERROR(%s): the \"%s\"-exchange is not in idle state, "
          "are you missing a get or wait?", routine_name, exchange->name);

        // do the actual exchange
        xt_redist_a_exchange(
          exchange->redist, exchange->count,
          (void const **)send_data, exchange->dummy_buffer,
          &(exchange->request));
        exchange->state = EXCHANGE_ACTIVE;
        break;
      }
      // both source and target
      case(3): {

        // if a previous exchange operation has not yet been completed
        if (exchange->state == EXCHANGE_ACTIVE) {
          xt_request_wait(&(exchange->request));
          exchange->state = EXCHANGE_IDLE;
        }

        YAC_ASSERT_F(
          (exchange->state == EXCHANGE_IDLE) ||
          (exchange->state == EXCHANGE_WAIT_GET),
          "ERROR(%s): the \"%s\"-exchange is in an invalid state",
          routine_name, exchange->name);

        switch (exchange->state) {

          default:
          case (EXCHANGE_WAIT_GET): {

            // do the actual exchange
            xt_redist_a_exchange(
              exchange->redist, exchange->count, (void const **)send_data,
              exchange->temp_recv_buffer, &(exchange->request));
            exchange->state = EXCHANGE_ACTIVE;

            break;
          }
          case (EXCHANGE_IDLE): {

            // we have to wait for the receive array before we can start the
            // actual exchange, therefore we just buffer the send array
            for (size_t i = 0; i < exchange->count; ++i)
              exchange->temp_send_buffer[i] =
                (void const *)(send_data[i]);
            exchange->state = EXCHANGE_WAIT_PUT;
          }
        }
        break;
      }
    }
  }
}

static void yac_interpolation_exchange_execute_get_(
  struct yac_interpolation_exchange * exchange, double ** recv_data,
  int is_async, char const * routine_name) {

  // if we have to do an exchange
  if (exchange->redist != NULL) {

    switch ((exchange->is_source << 1) +
            (exchange->is_target)) {

      // neither source nor target
      default:
      case(0): {

        // we should be in idle state
        YAC_ASSERT_F(
          exchange->state == EXCHANGE_IDLE,
          "ERROR(%s): the \"%s\"-exchange is not in idle state, "
          "internal error?", routine_name, exchange->name);
        do_empty_exchange(exchange, 0, routine_name);
        break;
      }
      // only target
      case(1): {

        // if a previous exchange operation has not yet been completed
        if (exchange->state == EXCHANGE_ACTIVE) {
          xt_request_wait(&(exchange->request));
          exchange->state = EXCHANGE_IDLE;
        }

        // we should be in idle state
        YAC_ASSERT_F(
          exchange->state == EXCHANGE_IDLE,
          "ERROR(%s): the \"%s\"-exchange is not in idle state, "
          "internal error?", routine_name, exchange->name);

        // if we have to send data
        void const ** send_data =
          (void const **)(exchange->dummy_buffer);

        // do the actual exchange
        if (is_async) {
          xt_redist_a_exchange(
            exchange->redist, exchange->count, send_data, (void **)recv_data,
            &(exchange->request));
          exchange->state = EXCHANGE_ACTIVE;
        } else {
          xt_redist_s_exchange(
            exchange->redist, exchange->count, send_data, (void **)recv_data);
          exchange->state = EXCHANGE_IDLE;
        }
        break;
      }
      // only source
      case(2): {

        // if this exchange is currently active
        if (exchange->state == EXCHANGE_ACTIVE) {

          if (is_async) {

            yac_interpolation_exchange_get_test(exchange, routine_name);

          } else {

            xt_request_wait(&(exchange->request));
            exchange->state = EXCHANGE_IDLE;
          }
        }
        break;
      }
      // both source and target
      case(3): {

        // if a previous exchange operation has not yet been completed
        if (exchange->state == EXCHANGE_ACTIVE) {
          xt_request_wait(&(exchange->request));
          exchange->state = EXCHANGE_IDLE;
        }

        YAC_ASSERT_F(
          (exchange->state == EXCHANGE_IDLE) ||
          (exchange->state == EXCHANGE_WAIT_PUT),
          "ERROR(%s): the \"%s\"-exchange is in an invalid state",
          routine_name, exchange->name);

        switch (exchange->state) {

          default:
          case (EXCHANGE_WAIT_PUT): {

            void const ** send_data =
                exchange->temp_send_buffer;

            // do the actual exchange
            if (is_async) {
              xt_redist_a_exchange(
                exchange->redist, exchange->count, send_data,
                (void **)recv_data, &(exchange->request));
              exchange->state = EXCHANGE_ACTIVE;
            } else {
              xt_redist_s_exchange(
                exchange->redist, exchange->count, send_data,
                (void **)recv_data);
              exchange->state = EXCHANGE_IDLE;
            }
            break;
          };

          case (EXCHANGE_IDLE): {

            YAC_ASSERT_F(
              is_async,
              "ERROR(%s): the \"%s\"-exchange is in an invalid state, "
              "is there a put missing?", routine_name, exchange->name);

            // we have to wait for the send array before we can start the
            // actual exchange, therefore we just buffer the recv array
            for (size_t i = 0; i < exchange->count; ++i)
              exchange->temp_recv_buffer[i] = (void **)(recv_data[i]);
            exchange->state = EXCHANGE_WAIT_GET;

            break;
          };
        }
        break;
      }
    }
  }
}

void yac_interpolation_exchange_execute_get(
  struct yac_interpolation_exchange * exchange, double ** recv_data,
  char const * routine_name) {
  yac_interpolation_exchange_execute_get_(
    exchange, recv_data, 0, routine_name);
}

void yac_interpolation_exchange_execute_get_async(
  struct yac_interpolation_exchange * exchange, double ** recv_data,
  char const * routine_name) {
  yac_interpolation_exchange_execute_get_(
    exchange, recv_data, 1, routine_name);
}

void yac_interpolation_exchange_execute(
  struct yac_interpolation_exchange * exchange,
  double const ** send_data_, double ** recv_data_,
  char const * routine_name) {

  // if we have to do an exchange
  if (exchange->redist != NULL) {

    // if a previous put has not yet been completed
    if (exchange->state == EXCHANGE_ACTIVE) {
      xt_request_wait(&(exchange->request));
      exchange->state = EXCHANGE_IDLE;
    }

    // ensure that the exchange is in idle state
    YAC_ASSERT_F(
      exchange->state == EXCHANGE_IDLE, "ERROR(%s): "
      "the \"%s\"-exchange is not in idle state, are you missing a get or wait?",
      routine_name, exchange->name);

    void const ** send_data =
      exchange->is_source?
      (void const **)(send_data_):(void const **)(exchange->dummy_buffer);
    void ** recv_data =
      exchange->is_target?(void **)(recv_data_):exchange->dummy_buffer;

    xt_redist_s_exchange(
      exchange->redist, exchange->count, send_data, recv_data);

    exchange->state = EXCHANGE_IDLE;
  } else {

    // ensure that the exchange is in idle state
    YAC_ASSERT_F(
      exchange->state == EXCHANGE_IDLE, "ERROR(%s): "
      "the \"%s\"-exchange is not in idle state",
      routine_name, exchange->name);
  }
}

void yac_interpolation_exchange_delete(
  struct yac_interpolation_exchange * exchange, char const * routine_name) {

  if (exchange->state != EXCHANGE_IDLE)
    yac_interpolation_exchange_wait(exchange, routine_name);

  free(exchange->name);

  if (exchange->redist != NULL) {
    xt_request_wait(&(exchange->request));
    xt_redist_delete(exchange->redist);
  }
  free(exchange->temp_send_buffer);
  free(exchange->temp_recv_buffer);
  free(exchange->dummy_buffer);
  free(exchange);
}
