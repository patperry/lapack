
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
F77_FUNC(dlarfg) (const int *N, double *alpha, double *X, const int *incX, double *tau);

void 
lapack_dlarfg (const int N, double *alpha, double *X, const int incX, double *tau)
{
    F77_FUNC(dlarfg) (&N, alpha, X, &incX, tau);
}

extern void
F77_FUNC(dormqr) (const char *side, const char *trans,
                  const int *M, const int *N, const int *K, const double *A,
                  const int *ldA, double *tau, double *C, const int *ldC,
                  double *Work, const int *ldWork, const int *info);
                   
int 
lapack_dormqr (const enum BLAS_SIDE side, const enum BLAS_TRANSPOSE trans,
               const int M, const int N, const int K, const double *A,
               const int ldA, double *tau, double *C, const int ldC,
               double *Work, const int ldWork)
{
    int info = 0;
    F77_FUNC(dormqr) (SIDE(side), TRANS(trans), &M, &N, &K, A, &ldA, tau, C,
                      &ldC, Work, &ldWork, &info);
    return info;
}

extern void
F77_FUNC(dormlq) (const char *side, const char *trans,
                  const int *M, const int *N, const int *K, const double *A,
                  const int *ldA, double *tau, double *C, const int *ldC,
                  double *Work, const int *ldWork, const int *info);
                   
int 
lapack_dormlq (const enum BLAS_SIDE side, const enum BLAS_TRANSPOSE trans,
               const int M, const int N, const int K, const double *A,
               const int ldA, double *tau, double *C, const int ldC,
               double *Work, const int ldWork)
{
    int info = 0;
    F77_FUNC(dormlq) (SIDE(side), TRANS(trans), &M, &N, &K, A, &ldA, tau, C,
                      &ldC, Work, &ldWork, &info);
    return info;
}
