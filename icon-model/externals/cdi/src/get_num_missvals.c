#include "get_num_missvals.h"
#include "cdi_int.h"

size_t
get_num_missvalsSP(size_t size, float *data, float missval)
{
  size_t numMissVals = 0;

  if (DBL_IS_NAN(missval))
    {
      for (size_t i = 0; i < size; i++)
        if (DBL_IS_EQUAL(data[i], missval))
          {
            data[i] = missval;
            numMissVals++;
          }
    }
  else
    {
      for (size_t i = 0; i < size; i++)
        if (IS_EQUAL(data[i], missval))
          {
            data[i] = missval;
            numMissVals++;
          }
    }

  return numMissVals;
}

size_t
get_num_missvalsDP(size_t size, double *data, double missval)
{
  size_t numMissVals = 0;

  if (DBL_IS_NAN(missval))
    {
      for (size_t i = 0; i < size; i++)
        if (DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float) missval))
          {
            data[i] = missval;
            numMissVals++;
          }
    }
  else
    {
      for (size_t i = 0; i < size; i++)
        if (IS_EQUAL(data[i], missval) || IS_EQUAL(data[i], (float) missval))
          {
            data[i] = missval;
            numMissVals++;
          }
    }

  return numMissVals;
}

size_t
get_cplx_num_missvalsSP(size_t size, float *data, float missval)
{
  size_t numMissVals = 0;

  if (DBL_IS_NAN(missval))
    {
      for (size_t i = 0; i < 2 * size; i += 2)
        if (DBL_IS_EQUAL(data[i], missval))
          {
            data[i] = missval;
            numMissVals++;
          }
    }
  else
    {
      for (size_t i = 0; i < 2 * size; i += 2)
        if (IS_EQUAL(data[i], missval))
          {
            data[i] = missval;
            numMissVals++;
          }
    }

  return numMissVals;
}

size_t
get_cplx_num_missvalsDP(size_t size, double *data, double missval)
{
  size_t numMissVals = 0;

  if (DBL_IS_NAN(missval))
    {
      for (size_t i = 0; i < 2 * size; i += 2)
        if (DBL_IS_EQUAL(data[i], missval) || DBL_IS_EQUAL(data[i], (float) missval))
          {
            data[i] = missval;
            numMissVals++;
          }
    }
  else
    {
      for (size_t i = 0; i < 2 * size; i += 2)
        if (IS_EQUAL(data[i], missval) || IS_EQUAL(data[i], (float) missval))
          {
            data[i] = missval;
            numMissVals++;
          }
    }

  return numMissVals;
}
