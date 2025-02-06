#ifndef PPM_EXTENTS_H
#define PPM_EXTENTS_H
/**
 * @file ppm_extents.h
 * @brief declarations for functions on extents
 *
 * @copyright Copyright  (C)  2011-2017  Thomas Jahns <jahns@dkrz.de>
 *
 * @version 1.0
 * @author Thomas Jahns <jahns@dkrz.de>
 */
/*
 * Maintainer: Thomas Jahns <jahns@dkrz.de>
 * URL: https://www.dkrz.de/redmine/projects/scales-ppm
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are  permitted provided that the following conditions are
 * met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the DKRZ GmbH nor the names of its contributors
 * may be used to endorse or promote products derived from this software
 * without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
 * PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
 * OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <inttypes.h>
#include <stdbool.h>
#include <stdlib.h>

#include "ppm_math_extensions.h"

/**
 * integer range description
 */
struct PPM_extent
{
  int32_t first, /**< the range is anchored at \a first */
    size;  /**< the range has size \a size, i.e. the range equals
            * \li [first, first + size - 1] for size > 0
            * \li [first + size + 1, first] for size < 0
            * \li the empty set if size = 0
            */
};

static inline int32_t
PPM_extent_size(struct PPM_extent extent)
{
  return extent.size;
}

static inline int32_t
PPM_extents_size(size_t nextents, const struct PPM_extent extent[])
{
  int32_t size = nextents > 0;
  for (size_t i = 0; i < nextents; ++i)
    size *= PPM_extent_size(extent[i]);
  return size;
}

static inline int64_t
PPM_extents_size64(size_t nextents, const struct PPM_extent extent[])
{
  int64_t size = nextents > 0;
  for (size_t i = 0; i < nextents; ++i)
    size *= PPM_extent_size(extent[i]);
  return size;
}

static inline int32_t
PPM_extent_start(struct PPM_extent extent)
{
  return extent.first;
}

static inline int32_t
PPM_extent_end(struct PPM_extent extent)
{
  return extent.first + extent.size  - ((extent.size >= 0) * 2 - 1);
}

static inline int
PPM_extent_eq(struct PPM_extent a, struct PPM_extent b)
{
  return (a.first == b.first)
    && (a.size == b.size || (abs(a.size) == 1 && a.size == -b.size));
}

static inline int
PPM_extents_eq(size_t nextents, struct PPM_extent a[nextents],
               struct PPM_extent b[nextents])
{
  size_t i;
  for (i = 0; i < nextents && PPM_extent_eq(a[i], b[i]); ++i)
    ;
  return i == nextents;
}

static inline bool
PPM_int32_is_contained_in_extent(int32_t i, struct PPM_extent rng)
{
  return  i >= rng.first && i < rng.first + rng.size;
}

static inline bool
PPM_coord_is_contained_in_extents(size_t dim,
                                  const int32_t coords[],
                                  const struct PPM_extent rng[])
{
  bool is_contained = true;
  for (size_t i = 0; i < dim; ++i)
    is_contained = is_contained && coords[i] >= rng[i].first
      && coords[i] < rng[i].first + rng[i].size;
  return is_contained;
}

/**
 * integer range description
 */
struct PPM_extent64
{
  int64_t first, /**< the range is anchored at \a first */
    size;  /**< the range has size \a size, i.e. the range equals
            * \li [first, first + size - 1] for size > 0
            * \li [first + size + 1, first] for size < 0
            * \li the empty set if size = 0
            */
};

static inline int64_t
PPM_extent64_size(struct PPM_extent64 extent)
{
  return extent.size;
}

static inline int64_t
PPM_extents64_size(size_t nextents, const struct PPM_extent64 extent[])
{
  int64_t size;
  if (nextents)
  {
    size = 1;
    for (size_t i = 0; i < nextents; ++i)
      size *= PPM_extent64_size(extent[i]);
  } else
    size = 0;
  return size;
}

static inline int64_t
PPM_extent64_start(struct PPM_extent64 extent)
{
  return extent.first;
}

static inline int64_t
PPM_extent64_end(struct PPM_extent64 extent)
{
  return extent.first + extent.size  - ((extent.size >= 0) * 2 - 1);
}

static inline int
PPM_extent64_eq(struct PPM_extent64 a, struct PPM_extent64 b)
{
  return (a.first == b.first)
    && (a.size == b.size || (llabs(a.size) == 1 && a.size == -b.size));
}

static inline int
PPM_extents64_eq(size_t nextents, struct PPM_extent64 a[nextents],
               struct PPM_extent64 b[nextents])
{
  size_t i;
  for (i = 0; i < nextents && PPM_extent64_eq(a[i], b[i]); ++i)
    ;
  return i == nextents;
}

static inline bool
PPM_int64_is_contained_in_extent64(int64_t i, struct PPM_extent64 rng)
{
  return  i >= rng.first && i < rng.first + rng.size;
}

static inline bool
PPM_coord_is_contained_in_extents64(size_t dim,
                                    const int64_t coords[],
                                    const struct PPM_extent64 rng[])
{
  bool is_contained = true;
  for (size_t i = 0; i < dim; ++i)
    is_contained = is_contained && coords[i] >= rng[i].first
      && coords[i] < rng[i].first + rng[i].size;
  return is_contained;
}

/**
 * inclusive integer range description
 *
 * describes range [first,last]
 */
struct PPM_iinterval
{
  int32_t first,   /**< the range is anchored at \a first */
    last;      /**< the range contains \a last */
};

static inline int32_t
PPM_iinterval_size(struct PPM_iinterval interval)
{
  int32_t rng_size = interval.last - interval.first;
  return rng_size + ((rng_size >= 0) * 2 - 1);
}

static inline int32_t
PPM_iintervals_size(size_t nintervals, const struct PPM_iinterval intervals[])
{
  int32_t size = nintervals > 0;
  for (size_t i = 0; i < nintervals; ++i)
    size *= PPM_iinterval_size(intervals[i]);
  return size;
}

static inline int32_t
PPM_iinterval_start(struct PPM_iinterval iinterval)
{
  return iinterval.first;
}

static inline int32_t
PPM_iinterval_end(struct PPM_iinterval iinterval)
{
  return iinterval.last;
}

static inline struct PPM_iinterval
PPM_extent2iinterval(const struct PPM_extent src)
{
  struct PPM_iinterval dst =
    { .first = src.first, .last = PPM_extent_end(src) };
  return dst;
}

void
PPM_extents2iintervals(int ndims, struct PPM_iinterval dst[ndims],
                       const struct PPM_extent src[ndims]);

/**
 * inclusive integer range description
 *
 * describes range [first,last]
 */
struct PPM_iinterval64
{
  int64_t first,   /**< the range is anchored at \a first */
    last;      /**< the range contains \a last */
};

static inline int64_t
PPM_iinterval64_size(struct PPM_iinterval64 interval)
{
  int64_t rng_size = interval.last - interval.first;
  return rng_size + ((rng_size >= 0) * 2 - 1);
}

static inline int64_t
PPM_iintervals64_size(size_t nintervals,
                      const struct PPM_iinterval64 intervals[])
{
  int64_t size = nintervals > 0;
  for (size_t i = 0; i < nintervals; ++i)
    size *= PPM_iinterval64_size(intervals[i]);
  return size;
}

static inline int64_t
PPM_iinterval64_start(struct PPM_iinterval64 iinterval)
{
  return iinterval.first;
}

static inline int64_t
PPM_iinterval64_end(struct PPM_iinterval iinterval)
{
  return iinterval.last;
}

static inline struct PPM_iinterval64
PPM_extent2iinterval64(const struct PPM_extent64 src)
{
  struct PPM_iinterval64 dst =
    { .first = src.first, .last = PPM_extent64_end(src) };
  return dst;
}

void
PPM_extents2iintervals64(int ndims, struct PPM_iinterval64 dst[ndims],
                         const struct PPM_extent64 src[ndims]);

/**
 * inclusive float range description
 *
 * describes range [first,last]
 */
struct PPM_iinterval_sp
{
  float first,  /**< the range is anchored at \a first */
    last;       /**< the range contains \a last */
};

/**
 * inclusive double range description
 *
 * describes range [first,last]
 */
struct PPM_iinterval_dp
{
  double first, /**< the range is anchored at \a first */
    last;       /**< the range contains \a last */
};

/**
 * Write string representation of \a extent to \a buf.
 * @param buf must point to array of appropriate size, 26 is
 * guaranteed to work.
 * @param ext range to print
 * @return number of characters written to \a buf
 */
int
PPM_sprint_extent(char buf[], const struct PPM_extent *ext);

/**
 * Write string representation of \a extent to \a buf.
 * @param buf must point to array of appropriate size, 44 is
 * guaranteed to work.
 * @param ext range to print
 * @return number of characters written to \a buf
 */
int
PPM_sprint_extent64(char buf[], const struct PPM_extent64 *ext);

/**
 * Write string representation of \a iinterval to \a buf.
 * @param buf must point to array of appropriate size, 26 is
 * guaranteed to work.
 * @param iinterval range to print
 * @return number of characters written to \a buf
 */
int
PPM_sprint_iinterval(char buf[], const struct PPM_iinterval *iinterval);

/**
 * Write string representation of \a iinterval to \a buf.
 * @param buf must point to array of appropriate size, 44 is
 * guaranteed to work.
 * @param iinterval range to print
 * @return number of characters written to \a buf
 */
int
PPM_sprint_iinterval64(char buf[], const struct PPM_iinterval64 *iinterval);

enum {
/**
 * number of characters needed at most to print a range of type
 *  \a PPM_iinterval_sp
 */
  PPM_IINTERVAL_SP_BUF_MAX = 4 + 2 * PPM_FLT_DECIMAL_WIDTH,
/**
 * number of characters needed at most to print a range of type
 *  \a PPM_iinterval_dp
 */
  PPM_IINTERVAL_DP_BUF_MAX = 4 + 2 * PPM_DBL_DECIMAL_WIDTH,
};


/**
 * Write string representation of \a iinterval to \a buf.
 * @param buf must point to character array large enough,
 * \a PPM_IINTERVAL_SP_BUF_MAX is guaranteed to suffice.
 * @param iinterval range to print
 * @return number of characters written to \a buf
 */
int
PPM_sprint_iinterval_sp(char buf[], const struct PPM_iinterval_sp *iinterval);


/**
 * Write string representation of \a iinterval to \a buf.
 * @param buf must point to a sufficiently large array,
 * \a PPM_IINTERVAL_DP_BUF_MAX is guaranteed to suffice.
 * @param iinterval range to print
 * @return number of characters written to \a buf
 */
int
PPM_sprint_iinterval_dp(char buf[], const struct PPM_iinterval_dp *iinterval);

/*
 * Local Variables:
 * license-project-url: "https://www.dkrz.de/redmine/projects/scales-ppm"
 * license-markup: "doxygen"
 * license-default: "bsd"
 * End:
 */
#endif
