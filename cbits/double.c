
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

