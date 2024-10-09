// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef YAC_LAPACK_INTERFACE_H
#define YAC_LAPACK_INTERFACE_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#ifndef YAC_LAPACK_INTERFACE_ID
#error None of the supported LAPACK interfaces is available
#endif

#if YAC_LAPACK_INTERFACE_ID == 1 // Intel MKL LAPACKE

#include <mkl_lapacke.h>

#elif YAC_LAPACK_INTERFACE_ID == 2 // Netlib LAPACKE

#include <lapacke.h>

#elif YAC_LAPACK_INTERFACE_ID == 3 // ATLAS CLAPACK

#include <clapack.h>

#ifndef ATL_INT
#define ATL_INT int
#endif

#ifndef lapack_int
#define lapack_int ATL_INT
#endif

#define LAPACK_COL_MAJOR CblasColMajor

#if 0
lapack_int LAPACKE_dgels_work( int matrix_layout, char trans, lapack_int m,
                               lapack_int n, lapack_int nrhs, double* a,
                               lapack_int lda, double* b, lapack_int ldb,
                               double* work, lapack_int lwork );
#endif

lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb );

lapack_int LAPACKE_dgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv );

lapack_int LAPACKE_dgetri_work( int matrix_layout, lapack_int n, double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                double* work, lapack_int lwork );

#define YAC_LAPACK_NO_DSYTR
#define YAC_LAPACK_C_INDEXING

#elif YAC_LAPACK_INTERFACE_ID == 4 // Netlib CLAPACK

#include <f2c.h>

#ifndef lapack_int
#define lapack_int integer
#endif

#define LAPACK_COL_MAJOR 102

#if 0
lapack_int LAPACKE_dgels_work( int matrix_layout, char trans, lapack_int m,
                               lapack_int n, lapack_int nrhs, double* a,
                               lapack_int lda, double* b, lapack_int ldb,
                               double* work, lapack_int lwork );
#endif

lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb );

lapack_int LAPACKE_dgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv );

lapack_int LAPACKE_dgetri_work( int matrix_layout, lapack_int n, double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                double* work, lapack_int lwork );

lapack_int LAPACKE_dsytrf_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda, lapack_int* ipiv,
                                double* work, lapack_int lwork );

lapack_int LAPACKE_dsytri_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda,
                                const lapack_int* ipiv, double* work );

#elif YAC_LAPACK_INTERFACE_ID == 5 // Fortran LAPACK

#ifndef lapack_int
#define lapack_int int
#endif

#define LAPACK_COL_MAJOR 102

#if 0
lapack_int LAPACKE_dgels_work( int matrix_layout, char trans, lapack_int m,
                               lapack_int n, lapack_int nrhs, double* a,
                               lapack_int lda, double* b, lapack_int ldb,
                               double* work, lapack_int lwork );
#endif

lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb );

lapack_int LAPACKE_dgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv );

lapack_int LAPACKE_dgetri_work( int matrix_layout, lapack_int n, double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                double* work, lapack_int lwork );

lapack_int LAPACKE_dsytrf_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda, lapack_int* ipiv,
                                double* work, lapack_int lwork );

lapack_int LAPACKE_dsytri_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda,
                                const lapack_int* ipiv, double* work );

#else

#error Unexpected value for YAC_LAPACK_INTERFACE_ID

#endif

#endif // YAC_LAPACK_INTERFACE_H

