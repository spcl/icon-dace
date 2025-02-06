#include <assert.h>
#include <stdlib.h>

#include "cdi.h"

/* function called by CDO */
extern int getByteswap(int);

int
main(void)
{
  assert(getByteswap(-1) == -1);
  return EXIT_SUCCESS;
}
