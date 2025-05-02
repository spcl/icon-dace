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

#ifdef _OPENMP
!-----------------------------------------------
! define default omp parameters
#define ICON_OMP_PARALLEL $omp parallel
#define ICON_OMP_END_PARALLEL $omp end parallel
#define ICON_OMP_DO $omp do
#define ICON_OMP_END_DO $omp end do
#define ICON_OMP_END_DO_NOWAIT $omp end do nowait
#define ICON_OMP $omp
#define ICON_OMP_PARALLEL_DO $omp parallel do
#define ICON_OMP_END_PARALLEL_DO $omp end parallel do
#define ICON_OMP_BARRIER $omp barrier
#define ICON_OMP_SINGLE $omp single
#define ICON_OMP_END_SINGLE $omp end single
#define ICON_OMP_MASTER $omp master
#define ICON_OMP_END_MASTER $omp end master
#define ICON_OMP_WORKSHARE $omp workshare
#define ICON_OMP_END_WORKSHARE $omp end workshare
!-----------------------------------------------

#if ( defined(_CRAYFTN) || defined(__INTEL_COMPILER) )
#define ICON_OMP_DEFAULT_SCHEDULE SCHEDULE(guided)
#define ICON_OMP_RUNTIME_SCHEDULE SCHEDULE(runtime)
#define ICON_OMP_GUIDED_SCHEDULE  SCHEDULE(guided)
#else 
! #define ICON_OMP_DEFAULT_SCHEDULE SCHEDULE(runtime)
#define ICON_OMP_DEFAULT_SCHEDULE SCHEDULE(dynamic,1)
#define ICON_OMP_RUNTIME_SCHEDULE SCHEDULE(runtime)
#define ICON_OMP_GUIDED_SCHEDULE  SCHEDULE(guided)
#endif
!-----------------------------------------------

!-----------------------------------------------
! define ICON HAMOCC omp parameters
#define ICON_HAMOCC_OMP_PARALLEL $omp parallel
#define ICON_HAMOCC_OMP_END_PARALLEL $omp end parallel
#define ICON_HAMOCC_OMP_DO $omp do
#define ICON_HAMOCC_OMP_END_DO $omp end do
#define ICON_HAMOCC_OMP_END_DO_NOWAIT $omp end do nowait
#define ICON_HAMOCC_OMP $omp
#define ICON_HAMOCC_OMP_PARALLEL_DO $omp parallel do
#define ICON_HAMOCC_OMP_END_PARALLEL_DO $omp end parallel do
#define ICON_HAMOCC_OMP_BARRIER $omp barrier
#define ICON_HAMOCC_OMP_SINGLE $omp single
#define ICON_HAMOCC_OMP_END_SINGLE $omp end single
#define ICON_HAMOCC_OMP_MASTER $omp master
#define ICON_HAMOCC_OMP_END_MASTER $omp end master
#define ICON_HAMOCC_OMP_WORKSHARE $omp workshare
#define ICON_HAMOCC_OMP_END_WORKSHARE $omp end workshare


!-----------------------------------------------
#endif
