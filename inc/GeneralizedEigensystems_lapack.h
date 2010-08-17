#ifndef _RNP_GENERALIZED_EIGENSYSTEMS_LAPACK_H_
#define _RNP_GENERALIZED_EIGENSYSTEMS_LAPACK_H_

#include <cstddef>
#include <complex>

extern "C" void zggev_(const char *jobvl, const char *jobvr, const long int &n,
	std::complex<double> *a, const long int &lda,
	std::complex<double> *b, const long int &ldb,
	std::complex<double> *alpha, std::complex<double> *beta,
	std::complex<double> *vl, const long int &ldvl,
	std::complex<double> *vr, const long int &ldvr,
	std::complex<double> *work, const long int &lwork,
	double *rwork, long int *info);

extern "C" void zgges_(const char *jobvsl, const char *jobvsr, const char *sort, const char *selctg,
	const long int &n,
	std::complex<double> *a, const long int &lda,
	std::complex<double> *b, const long int &ldb,
	long int *sdim,
	std::complex<double> *alpha, std::complex<double> *beta,
	std::complex<double> *vsl, const long int &ldvsl,
	std::complex<double> *vsr, const long int &ldvsr,
	std::complex<double> *work, const long int &lwork,
	double *rwork, bool *bwork, long int *info);

extern "C" void zggevx_(const char *balanc, const char *jobvl, const char *jobvr, const char *sense,
	const long int &n,
	std::complex<double> *a, const long int &lda,
	std::complex<double> *b, const long int &ldb,
	std::complex<double> *alpha, std::complex<double> *beta,
	std::complex<double> *vl, const long int &ldvl,
	std::complex<double> *vr, const long int &ldvr,
	long int *ilo, long int *ihi,
	double *lscale, double *rscale, double *abnrm, double *bbnrm,
	double *rconde, double *rcondv,
	std::complex<double> *work, const long int &lwork,
	double *rwork, long int *iwork, bool *bwork, long int *info);

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
	std::complex<double> *work_, double *rwork_, size_t lwork_ = 0){

	if(n == 0) {
		return 0;
	}
	char jobvl[2] = "N";
	char jobvr[2] = "N";
	if(vl != NULL){ jobvl[0] = 'V'; }
	if(vr != NULL){ jobvr[0] = 'V'; }
	long int info;

	if((size_t)-1 == lwork_){
		zggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work_, lwork_, rwork_, &info);
		return 0;
	}
	
	std::complex<double> *work = work_;
	double *rwork = rwork_;
	if(NULL == rwork_){
		rwork = new double[2*n];
	}
	if(0 == lwork_){ lwork_ = 2*n; }
	long int lwork = lwork_;
	if(NULL == work_ || lwork < (long int)(2*n)){
		lwork = -1;
		std::complex<double> zlen;
		zggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, &zlen, lwork, rwork, &info);
		lwork = (long int)(zlen.real());
		work = new std::complex<double>[lwork];
	}
	zggev_(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, vr, ldvr, work, lwork, rwork, &info);
	
	if(NULL == work_ || lwork < (long int)(2*n)){
		delete [] work;
	}
	if(NULL == rwork_){
		delete [] rwork;
	}
	
	return info;
}

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
	std::complex<double> *work, double *rwork){
}

}; // namespace RNP

#endif // _RNP_GENERALIZED_EIGENSYSTEMS_H_
