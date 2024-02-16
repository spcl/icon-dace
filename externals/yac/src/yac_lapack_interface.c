#include "yac_lapack_interface.h"

#if YAC_LAPACK_INTERFACE_ID == 3 // ATLAS CLAPACK

#include <assert.h>

#if 0
lapack_int LAPACKE_dgels_work( int matrix_layout, char trans, lapack_int m,
                               lapack_int n, lapack_int nrhs, double* a,
                               lapack_int lda, double* b, lapack_int ldb,
                               double* work, lapack_int lwork )
{
  assert(matrix_layout == LAPACK_COL_MAJOR);
  return (lapack_int) clapack_dgels(LAPACK_COL_MAJOR,
                                    trans == 'N' ? CblasNoTrans : CblasTrans,
                                    (ATL_INT) m, (ATL_INT) n, (ATL_INT) nrhs,
                                    a, (ATL_INT) lda, b, (int) ldb);
}
#endif

lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb )
{
  assert(matrix_layout == LAPACK_COL_MAJOR);
  if (sizeof(int) == sizeof(lapack_int))
  {
    return (lapack_int) clapack_dgesv(LAPACK_COL_MAJOR, (int) n, (int) nrhs,
                                      a, (int) lda, (int*) ipiv, b, (int) ldb);
  }
  else
  {
    int i, result, ipiv_size = (int) n;
    int ipiv_[ipiv_size];
    result = clapack_dgesv(LAPACK_COL_MAJOR, (int) n, (int) nrhs,
                           a, (int) lda, ipiv_, b, (int) ldb);
    for (i = 0; i != ipiv_size; ++i) { ipiv[i] = (lapack_int) ipiv_[i]; }
    return (lapack_int) result;
  }
}

lapack_int LAPACKE_dgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv )
{
  assert(matrix_layout == LAPACK_COL_MAJOR);
  if (sizeof(int) == sizeof(lapack_int))
  {
    return (lapack_int) clapack_dgetrf(LAPACK_COL_MAJOR, (int) m, (int) n,
                                       a, (int) lda, (int*) ipiv);
  }
  else
  {
    int i, result, ipiv_size = (int) (n < m ? n : m);
    int ipiv_[ipiv_size];
    result = clapack_dgetrf(LAPACK_COL_MAJOR, (int) m, (int) n,
                            a, (int) lda, ipiv_);
    for (i = 0; i != ipiv_size; ++i) { ipiv[i] = (lapack_int) ipiv_[i]; }
    return (lapack_int) result;
  } 
}

lapack_int LAPACKE_dgetri_work( int matrix_layout, lapack_int n, double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                double* work, lapack_int lwork )
{
  assert(matrix_layout == LAPACK_COL_MAJOR);
  if (sizeof(int) == sizeof(lapack_int))
  {
    return (lapack_int) clapack_dgetri(LAPACK_COL_MAJOR, (int) n, a,
                                     (int) lda, (int*) ipiv);
  }
  else
  {
    int i, ipiv_size = (int) n;
    int ipiv_[n];
    for (i = 0; i != ipiv_size; ++i) { ipiv_[i] = (int) ipiv[i]; }
    return (lapack_int) clapack_dgetri(LAPACK_COL_MAJOR, (int) n, a,
                                       (int) lda, ipiv_);
  }
}

#elif YAC_LAPACK_INTERFACE_ID == 4 // Netlib CLAPACK

#include <assert.h>
#include <clapack.h>

#if 0
lapack_int LAPACKE_dgels_work( int matrix_layout, char trans, lapack_int m,
                               lapack_int n, lapack_int nrhs, double* a,
                               lapack_int lda, double* b, lapack_int ldb,
                               double* work, lapack_int lwork )
{
  integer info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR &&
         sizeof(doublereal) == sizeof(double));
  if (sizeof(integer) == sizeof(lapack_int))
  {
    dgels_(&trans, (integer*) &m, (integer*) &n, (integer*) &nrhs,
           (doublereal*) a, (integer*) &lda,
           (doublereal*) b, (integer*) &ldb,
           (doublereal*) work, (integer*) &lwork,
           &info);
  }
  else
  {
    integer m_ = (integer) m, n_ = (integer) n,
            nrhs_ = (integer) nrhs, lda_ = (integer) lda,
            ldb_ = (integer) ldb, lwork_ = (integer) lwork;
    dgels_(&trans, &m_, &n_, &nrhs_,
           (doublereal*) a, &lda_,
           (doublereal*) b, &ldb_,
           (doublereal*) work, &lwork_,
           &info);
  }
  return (lapack_int) info;
}
#endif

lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb )
{
  integer info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR &&
         sizeof(doublereal) == sizeof(double));
  if (sizeof(integer) == sizeof(lapack_int))
  {
    dgesv_((integer*) &n, (integer*) &nrhs,
           (doublereal*) a, (integer*) &lda , (integer*) ipiv,
           (doublereal*) b, (integer*) &ldb, &info);
  }
  else
  {
    int i, ipiv_size = (int) n;
    integer n_ = (integer) n, nrhs_ = (integer) nrhs,
            lda_ = (integer) lda, ldb_ = (integer) ldb,
            ipiv_[n];
    dgesv_(&n_, &nrhs_,
           (doublereal*) a, &lda_, ipiv_,
           (doublereal*) b, &ldb_, &info);
    for (i = 0; i != ipiv_size; ++i) { ipiv[i] = (lapack_int) ipiv_[i]; }
  }
  return (lapack_int) info;
}

lapack_int LAPACKE_dgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv )
{
  integer info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR &&
         sizeof(doublereal) == sizeof(double));
  if (sizeof(integer) == sizeof(lapack_int))
  {
    dgetrf_((integer*) &m, (integer*) &n,
            (doublereal*) a, (integer*) &lda , (integer*) ipiv,
            &info);
  }
  else
  {
    int i, ipiv_size = (int) (n < m ? n : m);
    integer m_ = (integer) m, n_ = (integer) n,
            lda_ = (integer) lda, ipiv_[ipiv_size];
    dgetrf_(&m_, &n_,
            (doublereal*) a, &lda_, ipiv_,
            &info);
    for (i = 0; i != ipiv_size; ++i) {ipiv[i] = (lapack_int) ipiv_[i]; }
  }
  return (lapack_int) info;
}

lapack_int LAPACKE_dgetri_work( int matrix_layout, lapack_int n, double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                double* work, lapack_int lwork )
{
  integer info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR &&
         sizeof(doublereal) == sizeof(double));
  if (sizeof(integer) == sizeof(lapack_int))
  {
    dgetri_((integer*) &n, (doublereal*) a,
            (integer*) &lda , (integer*) ipiv,
            (doublereal*) work, (integer*) &lwork,
            &info);
  }
  else
  {
    int i, size_ipiv = (int) n;
    integer n_ = (integer) n, lda_ = (integer) lda,
            lwork_ = (integer) lwork, ipiv_[size_ipiv];
    for (i = 0; i != size_ipiv; ++i) { ipiv_[i] = (integer) ipiv[i]; }
    dgetri_(&n_, (doublereal*) a,
            &lda_, ipiv_,
            (doublereal*) work, &lwork_,
            &info);
  }
  return (lapack_int) info; 
}

lapack_int LAPACKE_dsytrf_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda, lapack_int* ipiv,
                                double* work, lapack_int lwork )
{
  integer info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR &&
         sizeof(doublereal) == sizeof(double));
  if (sizeof(integer) == sizeof(lapack_int))
  {
    dsytrf_(&uplo, (integer*) &n,
            (doublereal*) a, (integer*) &lda, (integer*) ipiv,
            (doublereal*) work, (integer*) &lwork, &info);
  }
  else
  {
    int i, size_ipiv = (int) n;
    integer n_ = (integer) n, lda_ = (integer) lda,
            lwork_ = (integer) lwork, ipiv_[size_ipiv];
    dsytrf_(&uplo, &n_,
            (doublereal*) a, &lda_, ipiv_,
            (doublereal*) work, &lwork_, &info);
    for (i = 0; i != size_ipiv; ++i) { ipiv[i] = (lapack_int) ipiv_[i]; }
  }
  return (lapack_int) info;
}

lapack_int LAPACKE_dsytri_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda,
                                const lapack_int* ipiv, double* work )
{
  integer info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR &&
         sizeof(doublereal) == sizeof(double));
  if (sizeof(integer) == sizeof(lapack_int))
  {
    dsytri_(&uplo, (integer*) &n,
            (doublereal*) a, (integer*) &lda, (integer*) ipiv,
            (doublereal*) work, &info);
  }
  else
  {
    int i, size_ipiv = (int) n;
    integer n_ = (integer) n, lda_ = (integer) lda,
            ipiv_[size_ipiv];
    for (i = 0; i != size_ipiv; ++i) { ipiv_[i] = (integer) ipiv[i]; }
    dsytri_(&uplo, &n_,
            (doublereal*) a, &lda_, ipiv_,
            (doublereal*) work, &info);
  }
  return (lapack_int) info; 
}

#elif YAC_LAPACK_INTERFACE_ID == 5 // Fortran LAPACK

#include <assert.h>

#ifdef __cplusplus
extern "C" {
#endif

#define LAPACK_dgels  YAC_FC_GLOBAL(dgels,DGELS)
#define LAPACK_dgesv  YAC_FC_GLOBAL(dgesv,DGESV)
#define LAPACK_dgetrf YAC_FC_GLOBAL(dgetrf,DGETRF)
#define LAPACK_dgetri YAC_FC_GLOBAL(dgetri,DGETRI)
#define LAPACK_dsytrf YAC_FC_GLOBAL(dsytrf,DSYTRF)
#define LAPACK_dsytri YAC_FC_GLOBAL(dsytri,DSYTRI)

#if 0
void LAPACK_dgels( char* trans, lapack_int* m, lapack_int* n, lapack_int* nrhs,
                   double* a, lapack_int* lda, double* b, lapack_int* ldb,
                   double* work, lapack_int* lwork, lapack_int *info );
#endif
void LAPACK_dgesv( lapack_int* n, lapack_int* nrhs, double* a, lapack_int* lda,
                   lapack_int* ipiv, double* b, lapack_int* ldb,
                   lapack_int *info );
void LAPACK_dgetrf( lapack_int* m, lapack_int* n, double* a, lapack_int* lda,
                    lapack_int* ipiv, lapack_int *info );
void LAPACK_dgetri( lapack_int* n, double* a, lapack_int* lda,
                    const lapack_int* ipiv, double* work, lapack_int* lwork,
                    lapack_int *info );
void LAPACK_dsytrf( char* uplo, lapack_int* n, double* a, lapack_int* lda,
                    lapack_int* ipiv, double* work, lapack_int* lwork,
                    lapack_int *info );
void LAPACK_dsytri( char* uplo, lapack_int* n, double* a, lapack_int* lda,
                    const lapack_int* ipiv, double* work, lapack_int *info );

#ifdef __cplusplus
}
#endif


#if 0
lapack_int LAPACKE_dgels_work( int matrix_layout, char trans, lapack_int m,
                               lapack_int n, lapack_int nrhs, double* a,
                               lapack_int lda, double* b, lapack_int ldb,
                               double* work, lapack_int lwork )
{
  lapack_int info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR);
  LAPACK_dgels(&trans, &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
  if( info < 0 ) { info = info - 1; }
  return info;
}
#endif

lapack_int LAPACKE_dgesv( int matrix_layout, lapack_int n, lapack_int nrhs,
                          double* a, lapack_int lda, lapack_int* ipiv,
                          double* b, lapack_int ldb )
{
  lapack_int info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR);
  LAPACK_dgesv( &n, &nrhs, a, &lda, ipiv, b, &ldb, &info );
  if( info < 0 ) { info = info - 1; }
  return info;
}

lapack_int LAPACKE_dgetrf( int matrix_layout, lapack_int m, lapack_int n,
                           double* a, lapack_int lda, lapack_int* ipiv )
{
  lapack_int info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR);
  LAPACK_dgetrf( &m, &n, a, &lda, ipiv, &info );
  if( info < 0 ) { info = info - 1; }
  return info;
}

lapack_int LAPACKE_dgetri_work( int matrix_layout, lapack_int n, double* a,
                                lapack_int lda, const lapack_int* ipiv,
                                double* work, lapack_int lwork )
{
  lapack_int info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR);
  LAPACK_dgetri( &n, a, &lda, ipiv, work, &lwork, &info );
  if( info < 0 ) { info = info - 1; }
  return info;
}

lapack_int LAPACKE_dsytrf_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda, lapack_int* ipiv,
                                double* work, lapack_int lwork )
{
  lapack_int info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR);
  LAPACK_dsytrf( &uplo, &n, a, &lda, ipiv, work, &lwork, &info );
  if( info < 0 ) { info = info - 1; }
  return info;
}

lapack_int LAPACKE_dsytri_work( int matrix_layout, char uplo, lapack_int n,
                                double* a, lapack_int lda,
                                const lapack_int* ipiv, double* work )
{
  lapack_int info = 0;
  assert(matrix_layout == LAPACK_COL_MAJOR);
  LAPACK_dsytri( &uplo, &n, a, &lda, ipiv, work, &info );
  if( info < 0 ) { info = info - 1; }
  return info;
}

#endif
