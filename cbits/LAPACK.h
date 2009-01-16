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

void lapack_dlarfg (const int N, double *alpha, double *X, const int incX, double *tau);
void lapack_zlarfg (const int N, void *alpha, void *X, const int incX, void *tau);

#ifdef __cplusplus
}
#endif

#endif
