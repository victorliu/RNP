#ifndef _RNP_GENERALIZED_EIGENSYSTEMS_H_
#define _RNP_GENERALIZED_EIGENSYSTEMS_H_

namespace RNP{

// GeneralizedEigensystem
// ----------------------
// Computes for a pair of N-by-N complex nonsymmetric matrices (A,B),
// the generalized eigenvalues, and optionally, the left and/or right
// generalized eigenvectors.
//
// A generalized eigenvalue for a pair of matrices (A,B) is a scalar
// lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
// singular. It is usually represented as the pair (alpha,beta), as
// there is a reasonable interpretation for beta=0, and even for both
// being zero.
//
// The right generalized eigenvector v(j) corresponding to the
// generalized eigenvalue lambda(j) of (A,B) satisfies
//              A * v(j) = lambda(j) * B * v(j).
// The left generalized eigenvector u(j) corresponding to the
// generalized eigenvalues lambda(j) of (A,B) satisfies
//              u(j)^H * A = lambda(j) * u(j)^H * B
// where u(j)^H is the conjugate-transpose of u(j).
//
// Arguments
// =========
//
// N       The order of the matrices A, B, VL, and VR.  N >= 0.
// A       (input/output) dimension (LDA, N)
//         On entry, the matrix A in the pair (A,B).
//         On exit, A has been overwritten.
// LDA     The leading dimension of A.  LDA >= max(1,N).
// B       (input/output) dimension (LDB, N)
//         On entry, the matrix B in the pair (A,B).
//         On exit, B has been overwritten.
// LDB     The leading dimension of B.  LDB >= max(1,N).
//
// ALPHA   (output) dimension (N)
// BETA    (output) dimension (N)
//         On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
//         generalized eigenvalues.
//
//         Note: the quotients ALPHA(j)/BETA(j) may easily over- or
//         underflow, and BETA(j) may even be zero.  Thus, the user
//         should avoid naively computing the ratio alpha/beta.
//         However, ALPHA will be always less than and usually
//         comparable with norm(A) in magnitude, and BETA always less
//         than and usually comparable with norm(B).
//
// VL      (output) dimension (LDVL,N)
//         If VL != NULL, the left generalized eigenvectors u(j) are
//         stored one after another in the columns of VL, in the same
//         order as their eigenvalues.
//         Each eigenvector is scaled so the largest component has
//         abs(real part) + abs(imag. part) = 1.
// LDVL    The leading dimension of the matrix VL. LDVL >= 1, and
//         if VL != NULL, LDVL >= N.
//
// VR      (output) dimension (LDVR,N)
//         If VR != NULL, the right generalized eigenvectors v(j) are
//         stored one after another in the columns of VR, in the same
//         order as their eigenvalues.
//         Each eigenvector is scaled so the largest component has
//         abs(real part) + abs(imag. part) = 1.
// LDVR    The leading dimension of the matrix VR. LDVR >= 1, and
//         if VR != NULL, LDVR >= N.
//
// WORK    (workspace) dimension (MAX(1,2*N))
// RWORK   (workspace) dimension (8*N)
//
// return: = 0:  successful exit
//         < 0:  if INFO = -i, the i-th argument had an illegal value.
//         =1,...,N:
//               The QZ iteration failed.  No eigenvectors have been
//               calculated, but ALPHA(j) and BETA(j) should be
//               correct for j=INFO+1,...,N.
//         > N:  =N+1: other then QZ iteration failed in DHGEQZ,
//               =N+2: error return from DTGEVC.
int GeneralizedEigensystem(size_t n, 
	std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb, 
	std::complex<double> *alpha, std::complex<double> *beta,
	std::complex<double> *vl, size_t ldvl, std::complex<double> *vr, size_t ldvr,
	std::complex<double> *work, double *rwork);

// GeneralizedSchurDecomposition
// -----------------------------
// Computes for a pair of N-by-N complex nonsymmetric matrices (A,B),
// the generalized eigenvalues, the generalized complex Schur form
// (S, T), and optionally left and/or right Schur vectors (VSL and
// VSR). This gives the generalized Schur factorization
//         (A,B) = ( (VSL)*S*(VSR)^H, (VSL)*T*(VSR)^H )
// where (VSR)^H is the conjugate-transpose of VSR.
//
// A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
// or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
// usually represented as the pair (alpha,beta), as there is a
// reasonable interpretation for beta=0, and even for both being zero.
//
// A pair of matrices (S,T) is in generalized complex Schur form if S
// and T are upper triangular and, in addition, the diagonal elements
// of T are non-negative real numbers.
//
// Arguments
// =========
// N       The order of the matrices A, B, VSL, and VSR.  N >= 0.
// A       (input/output) dimension (LDA, N)
//         On entry, the first of the pair of matrices.
//         On exit, A has been overwritten by its generalized Schur
//         form S.
// LDA     The leading dimension of A.  LDA >= max(1,N).
// B       On entry, the second of the pair of matrices.
//         On exit, B has been overwritten by its generalized Schur
//         form T.
// LDB     The leading dimension of B.  LDB >= max(1,N).
//
// ALPHA   (output) dimension (N)
// BETA    (output) dimension (N)
//         On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the
//         generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j),
//         j=1,...,N  are the diagonals of the complex Schur form (A,B)
//         of the output. The  BETA(j) will be non-negative real.
//
//         Note: the quotients ALPHA(j)/BETA(j) may easily over- or
//         underflow, and BETA(j) may even be zero.  Thus, the user
//         should avoid naively computing the ratio alpha/beta.
//         However, ALPHA will be always less than and usually
//         comparable with norm(A) in magnitude, and BETA always less
//         than and usually comparable with norm(B).
//
// VSL     (output) dimension (LDVSL,N)
//         If VSL != NULL, VSL will contain the left Schur vectors.
//         Ignored if NULL
// LDVSL   The leading dimension of the matrix VSL. LDVSL >= 1, and
//         if VSL != NULL, LDVSL >= N.
//
// VSR     (output) dimension (LDVSR,N)
//         If VSR != NULL, VSR will contain the right Schur vectors.
// LDVSR   The leading dimension of the matrix VSR. LDVSR >= 1, and
//         if VSR != NULL, LDVSR >= N.
//
// WORK    (workspace) dimension (MAX(1,2*N))
// RWORK   (workspace) dimension (8*N)
//
// return: = 0:  successful exit
//         < 0:  if INFO = -i, the i-th argument had an illegal value.
//         =1,...,N:
//               The QZ iteration failed.  (A,B) are not in Schur
//               form, but ALPHA(j) and BETA(j) should be correct for
//               j=INFO+1,...,N.
//         > N:  =N+1: other than QZ iteration failed in ZHGEQZ
int GeneralizedSchurDecomposition(size_t n,
	std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb,
	std::complex<double> *alpha, std::complex<double> *beta,
	std::complex<double> *vsl, size_t ldvsl, std::complex<double> *vsr, size_t ldvsr,
	std::complex<double> *work, double *rwork);

}; // namespace RNP

#endif // _RNP_GENERALIZED_EIGENSYSTEMS_H_
