#ifndef LAPACK_H
#define LAPACK_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef BLAS_ENUM_DEFINED_H
   #define BLAS_ENUM_DEFINED_H
   enum BLAS_TRANSPOSE {BlasNoTrans=111, BlasTrans=112, BlasConjTrans=113};
   enum BLAS_UPLO  {BlasUpper=121, BlasLower=122};
   enum BLAS_DIAG  {BlasNonUnit=131, BlasUnit=132};
   enum BLAS_SIDE  {BlasLeft=141, BlasRight=142};
#endif

#define LAPACK_INDEX int

int lapack_dgeqrf (const int M, const int N, double *A, const int ldA,
                   double *tau, double *work, const int lwork);
int lapack_zgeqrf (const int M, const int N, void *A, const int ldA,
                   void *tau, void *work, const int lwork);

int lapack_dgelqf (const int M, const int N, double *A, const int ldA,
                   double *tau, double *work, const int lwork);
int lapack_zgelqf (const int M, const int N, void *A, const int ldA,
                   void *tau, void *work, const int lwork);

int lapack_dormqr (const enum BLAS_SIDE side, const enum BLAS_TRANSPOSE trans,
                   const int M, const int N, const int K, const double *A,
                   const int ldA, double *tau, double *C, const int ldC,
                   double *Work, const int ldWork);
int lapack_zunmqr (const enum BLAS_SIDE side, const enum BLAS_TRANSPOSE trans,
                   const int M, const int N, const int K, const void *A,
                   const int ldA, void *tau, void *C, const int ldC,
                   void *Work, const int ldWork);

int lapack_dormlq (const enum BLAS_SIDE side, const enum BLAS_TRANSPOSE trans,
                   const int M, const int N, const int K, const double *A,
                   const int ldA, double *tau, double *C, const int ldC,
                   double *Work, const int ldWork);
int lapack_zunmlq (const enum BLAS_SIDE side, const enum BLAS_TRANSPOSE trans,
                   const int M, const int N, const int K, const void *A,
                   const int ldA, void *tau, void *C, const int ldC,
                   void *Work, const int ldWork);

void lapack_dlarfg (const int N, double *alpha, double *X, const int incX, double *tau);
void lapack_zlarfg (const int N, void *alpha, void *X, const int incX, void *tau);


#ifdef __cplusplus
}
#endif

#endif
