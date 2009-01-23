#include "LAPACK.h"
#include "config.h"

static char *BLAS_TRANS_CODES[] = { "N", "T", "C" };
#define TRANS(x) BLAS_TRANS_CODES[(int) (x) - (int) BlasNoTrans]

static char *BLAS_UPLO_CODES[] = { "U", "L" };
#define UPLO(x) BLAS_UPLO_CODES[(int) (x) - (int) BlasUpper]

static char *BLAS_DIAG_CODES[] = { "N", "U" };
#define DIAG(x) BLAS_DIAG_CODES[(int) (x) - (int) BlasNonUnit]

static char *BLAS_SIDE_CODES[] = { "L", "R" };
#define SIDE(x) BLAS_SIDE_CODES[(int) (x) - (int) BlasLeft]

extern void
F77_FUNC(zlarfg) (const int *N, void *alpha, void *X, const int *incX, void *tau);

void 
lapack_zlarfg (const int N, void *alpha, void *X, const int incX, void *tau)
{
    F77_FUNC(zlarfg) (&N, alpha, X, &incX, tau);
}

extern void
F77_FUNC(zunmqr) (const char *side, const char *trans,
                  const int *M, const int *N, const int *K, const void *A,
                  const int *ldA, void *tau, void *C, const int *ldC,
                  void *Work, const int *ldWork, const int *info);
                   
int 
lapack_zunmqr (const enum BLAS_SIDE side, const enum BLAS_TRANSPOSE trans,
               const int M, const int N, const int K, const void *A,
               const int ldA, void *tau, void *C, const int ldC,
               void *Work, const int ldWork)
{
    int info = 0;
    F77_FUNC(zunmqr) (SIDE(side), TRANS(trans), &M, &N, &K, A, &ldA, tau, C, 
                      &ldC, Work, &ldWork, &info);
    return info;
}

extern void
F77_FUNC(zunmlq) (const char *side, const char *trans,
                  const int *M, const int *N, const int *K, const void *A,
                  const int *ldA, void *tau, void *C, const int *ldC,
                  void *Work, const int *ldWork, const int *info);
                   
int 
lapack_zunmlq (const enum BLAS_SIDE side, const enum BLAS_TRANSPOSE trans,
               const int M, const int N, const int K, const void *A,
               const int ldA, void *tau, void *C, const int ldC,
               void *Work, const int ldWork)
{
    int info = 0;
    F77_FUNC(zunmlq) (SIDE(side), TRANS(trans), &M, &N, &K, A, &ldA, tau, C, 
                      &ldC, Work, &ldWork, &info);
    return info;
}
