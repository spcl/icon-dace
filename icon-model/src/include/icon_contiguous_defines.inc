! ICON
!
! ---------------------------------------------------------------
! Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
! Contact information: icon-model.org
!
! See AUTHORS.TXT for a list of authors
! See LICENSES/ for license information
! SPDX-License-Identifier: BSD-3-Clause
! ---------------------------------------------------------------

#ifdef HAVE_FC_ATTRIBUTE_CONTIGUOUS
#define CONTIGUOUS_TARGET CONTIGUOUS, TARGET
#define CONTIGUOUS_POINTER CONTIGUOUS, POINTER
#define CONTIGUOUS_ARGUMENT(intent) CONTIGUOUS, INTENT(intent)
#else
#define CONTIGUOUS_TARGET TARGET
#define CONTIGUOUS_POINTER POINTER
#define CONTIGUOUS_ARGUMENT(intent) INTENT(intent)
#endif
