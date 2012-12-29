#define RNP_TBLAS_USE_RANDOM
#include "TBLAS.h" // includes cmath with the proper _USE_MATH_DEFINES
#include "TLASupport.h"
#include "IRA.h"
#include <limits>
#include <algorithm>
#include <cmath>
#include "Eigensystems.h"

bool RNP::IRA::LargestMagnitude (const std::complex<double> &a, const std::complex<double> &b, void *data){
	return std::abs(a) < std::abs(b);
}
bool RNP::IRA::SmallestMagnitude(const std::complex<double> &a, const std::complex<double> &b, void *data){
	return std::abs(a) > std::abs(b);
}
bool RNP::IRA::LargestRealPart (const std::complex<double> &a, const std::complex<double> &b, void *data){
	return a.real() < b.real();
}
bool RNP::IRA::SmallestRealPart(const std::complex<double> &a, const std::complex<double> &b, void *data){
	return a.real() > b.real();
}
bool RNP::IRA::LargestImaginaryPart (const std::complex<double> &a, const std::complex<double> &b, void *data){
	return a.imag() < b.imag();
}
bool RNP::IRA::SmallestImaginaryPart(const std::complex<double> &a, const std::complex<double> &b, void *data){
	return a.imag() > b.imag();
}

struct znaupd_data{
	size_t ldh, ldq, nb;
	size_t ih, ritz, bounds, iq, iw, next;
	size_t np, nev0;
	size_t mxiter;
	int mode, iupd;
	int ishift;
};

struct znaup2_data{
	int np0;
	int kplusp;
	size_t nev0;
	size_t iter;
	size_t nconv;
	size_t nevbef;
	double rnorm;
	double eps23;
	bool initv;
	bool update, ushift;
	bool getv0, cnorm;
};

struct znaitr_data{
	size_t j;
	int ipj, irj, ivj;
	int ierr, iter;
	int itry;
	bool orth1, orth2, step3, step4;
	bool rstart;
	double betaj;
	double wnorm;
	double rnorm1;
};

struct zgetv0_data{
	size_t iter;
	bool orth;
	bool first;
	double rnorm0;
};

struct arz_data{
	void *sort_data;
	zgetv0_data s_zgetv0;
	znaitr_data s_znaitr;
	znaup2_data s_znaup2;
	znaupd_data s_znaupd;
};

static void zsortc_(RNP::IRA::SortFunc zsortc_func, size_t n, std::complex<double> *x, std::complex<double> *y, bool reverse, void *data){
	size_t igap = n / 2;

	while(igap != 0){
		for(size_t i = igap; i < n; ++i){
			int j = i - igap;
			while(j >= 0){
				if (j < 0) {
					break;
				}

				if(reverse == zsortc_func(x[j], x[j+igap], data)){
					std::swap(x[j], x[j+igap]);
					if(NULL != y){
						std::swap(y[j], y[j+igap]);
					}
				}else{
					break;
				}
				j -= igap;
			}
		}
		igap /= 2;
	}
}


// \BeginDoc

// \Name: zngets

// \Description:
//  Given the eigenvalues of the upper Hessenberg matrix H,
//  computes the NP shifts AMU that are zeros of the polynomial of
//  degree NP which filters out components of the unwanted eigenvectors
//  corresponding to the AMU's based on some given criteria.

//  NOTE: call this even in the case of user specified shifts in order
//  to sort the eigenvalues, and error bounds of H for later use.

// \Usage:
//  call zngets
//      ( ISHIFT, WHICH, KEV, NP, RITZ, BOUNDS )

// \Arguments
//  ISHIFT  Integer.  (INPUT)
//          Method for selecting the implicit shifts at each iteration.
//          ISHIFT = 0: user specified shifts
//          ISHIFT = 1: exact shift with respect to the matrix H.

//  WHICH   Character*2.  (INPUT)
//          Shift selection criteria.
//          'LM' -> want the KEV eigenvalues of largest magnitude.
//          'SM' -> want the KEV eigenvalues of smallest magnitude.
//          'LR' -> want the KEV eigenvalues of largest REAL part.
//          'SR' -> want the KEV eigenvalues of smallest REAL part.
//          'LI' -> want the KEV eigenvalues of largest imaginary part.
//          'SI' -> want the KEV eigenvalues of smallest imaginary part.

//  KEV     Integer.  (INPUT)
//          The number of desired eigenvalues.

//  NP      Integer.  (INPUT)
//          The number of shifts to compute.

//  RITZ    Complex*16 array of length KEV+NP.  (INPUT/OUTPUT)
//          On INPUT, RITZ contains the the eigenvalues of H.
//          On OUTPUT, RITZ are sorted so that the unwanted
//          eigenvalues are in the first NP locations and the wanted
//          portion is in the last KEV locations.  When exact shifts are
//          selected, the unwanted part corresponds to the shifts to
//          be applied. Also, if ISHIFT .eq. 1, the unwanted eigenvalues
//          are further sorted so that the ones with largest Ritz values
//          are first.

//  BOUNDS  Complex*16 array of length KEV+NP.  (INPUT/OUTPUT)
//          Error bounds corresponding to the ordering in RITZ.



// \EndDoc

// -----------------------------------------------------------------------

// \BeginLib

// \Local variables:
//     xxxxxx  Complex*16

// \Routines called:
//     zsortc  ARPACK sorting routine.

// \Author
//     Danny Sorensen               Phuong Vu
//     Richard Lehoucq              CRPC / Rice University
//     Dept. of Computational &     Houston, Texas
//     Applied Mathematics
//     Rice University
//     Houston, Texas

// \SCCS Information: @(#)
// FILE: ngets.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2

// \Remarks
//     1. This routine does not keep complex conjugate pairs of
//        eigenvalues together.

// \EndLib

// -----------------------------------------------------------------------

void zngets_(int ishift, RNP::IRA::SortFunc sorter, size_t kev, size_t np, std::complex<double> *ritz, std::complex<double> *bounds, void *data){

    zsortc_(sorter, kev + np, ritz, bounds, false, data);

    if (ishift == 1) {
		// Sort the unwanted Ritz values used as shifts so that  | */
		// the ones with largest Ritz estimates are first        | */
		// This will tend to minimize the effects of the         | */
		// forward instability of the iteration when the shifts  | */
		// are applied in subroutine znapps.                     | */
		// Be careful and use 'SM' since we want to sort BOUNDS! | */
		zsortc_(&RNP::IRA::SmallestMagnitude, np, bounds, ritz, false, NULL);
    }
}


// \BeginDoc

// \Name: zneigh

// \Description:
//  Compute the eigenvalues of the current upper Hessenberg matrix
//  and the corresponding Ritz estimates given the current residual norm.

// \Usage:
//  call zneigh
//     ( RNORM, N, H, LDH, RITZ, BOUNDS, Q, LDQ, WORKL, RWORK, IERR )

// \Arguments
//  RNORM   Double precision scalar.  (INPUT)
//          Residual norm corresponding to the current upper Hessenberg
//          matrix H.

//  N       int.  (INPUT)
//          Size of the matrix H.

//  H       Complex*16 N by N array.  (INPUT)
//          H contains the current upper Hessenberg matrix.

//  LDH     int.  (INPUT)
//          Leading dimension of H exactly as declared in the calling
//          program.

//  RITZ    Complex*16 array of length N.  (OUTPUT)
//          On output, RITZ(1:N) contains the eigenvalues of H.

//  BOUNDS  Complex*16 array of length N.  (OUTPUT)
//          On output, BOUNDS contains the Ritz estimates associated with
//          the eigenvalues held in RITZ.  This is equal to RNORM
//          times the last components of the eigenvectors corresponding
//          to the eigenvalues in RITZ.

//  Q       Complex*16 N by N array.  (WORKSPACE)
//          Workspace needed to store the eigenvectors of H.

//  LDQ     int.  (INPUT)
//          Leading dimension of Q exactly as declared in the calling
//          program.

//  WORKL   Complex*16 work array of length N**2 + 3*N.  (WORKSPACE)
//          Private (replicated) array on each PE or array allocated on
//          the front end.  This is needed to keep the full Schur form
//          of H and also in the calculation of the eigenvectors of H.

//  RWORK   Double precision  work array of length N (WORKSPACE)
//          Private (replicated) array on each PE or array allocated on
//          the front end.

//  IERR    int.  (OUTPUT)
//          Error exit flag from zlahqr or ztrevc.

// \EndDoc

// -----------------------------------------------------------------------

// \BeginLib

// \Local variables:
//     xxxxxx  Complex*16

// \Routines called:
//     zlahqr  LAPACK routine to compute the Schur form of an
//             upper Hessenberg matrix.
//     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
//             in upper triangular form


// \Author
//     Danny Sorensen               Phuong Vu
//     Richard Lehoucq              CRPC / Rice University
//     Dept. of Computational &     Houston, Texas
//     Applied Mathematics
//     Rice University
//     Houston, Texas

// \SCCS Information: @(#)
// FILE: neigh.F   SID: 2.2   DATE OF SID: 4/20/96   RELEASE: 2

// \Remarks
//     None

// \EndLib

// -----------------------------------------------------------------------

int zneigh_(const double &rnorm, size_t n, const std::complex<double> *h,
	size_t ldh, std::complex<double> *ritz, std::complex<double> *bounds,
	std::complex<double> *q, size_t ldq, std::complex<double> *workl, double *rwork)
{
	using namespace std;
	int ierr;

	extern int ztrevc_(char howmny, bool *select, size_t n, std::complex<double> *_t, size_t ldt, std::complex<double> *_vl, size_t ldvl, std::complex<double> *_vr, size_t ldvr, size_t mm, size_t *m, std::complex<double> *_work, double *rwork);

	extern size_t zlahqr_(
		bool wantt, bool wantz,
		size_t n, size_t ilo, size_t ihi,
		std::complex<double> *h, int ldh,
		std::complex<double> *w,
		size_t iloz, size_t ihiz,
		std::complex<double> *z, size_t ldz);

	// 1. Compute the eigenvalues, the last components of the
	//    corresponding Schur vectors and the full Schur form T
	//    of the current upper Hessenberg matrix H.
	//    zlahqr returns the full Schur form of H
	//    in WORKL(1:N**2), and the Schur vectors in q.
	RNP::TBLAS::CopyMatrix<'A'>(n, n, h, ldh, workl, n);
	RNP::TBLAS::SetMatrix<'A'>(n, n, 0.0, 1.0, q, ldq);
	ierr = zlahqr_(true, true, n, 1, n, workl, ldh, ritz, 1, n, q, ldq);
	if (ierr != 0) {
		goto L9000;
	}
	RNP::TBLAS::Copy(n, &q[n-2+0*ldq], ldq, bounds, 1);

	// 2. Compute the eigenvectors of the full Schur form T and
	//    apply the Schur vectors to get the corresponding
	//    eigenvectors.
	{size_t dummy_n;
		ierr = ztrevc_('B', NULL, n, workl, n, NULL, n, q, ldq, n, &dummy_n, &workl[n*n], rwork);
	}

	if (ierr != 0) {
		goto L9000;
	}

	// Scale the returning eigenvectors so that their
	// Euclidean norms are all one. LAPACK subroutine
	// ztrevc returns each eigenvector normalized so
	// that the element of largest magnitude has
	// magnitude 1; here the magnitude of a complex
	// number (x,y) is taken to be |x| + |y|.
	for (size_t j = 0; j < n; ++j) {
		double temp = RNP::TBLAS::Norm2(n, &q[0+j*ldq], 1);
		RNP::TBLAS::Scale(n, 1. / temp, &q[0+j*ldq], 1);
	}


	// Compute the Ritz estimates
	RNP::TBLAS::Copy(n, &q[n-1+0*ldq], n, bounds, 1);
	RNP::TBLAS::Scale(n, rnorm, bounds, 1);

L9000:
	return ierr;
}

static inline double _abs1(const std::complex<double> &z){
	return fabs(z.real()) + fabs(z.imag());
}

// \BeginDoc

// \Name: znapps

// \Description:
//  Given the Arnoldi factorization

//     A*V_{k} - V_{k}*H_{k} = r_{k+p}*e_{k+p}^T,

//  apply NP implicit shifts resulting in

//     A*(V_{k}*Q) - (V_{k}*Q)*(Q^T* H_{k}*Q) = r_{k+p}*e_{k+p}^T * Q

//  where Q is an orthogonal matrix which is the product of rotations
//  and reflections resulting from the NP bulge change sweeps.
//  The updated Arnoldi factorization becomes:

//     A*VNEW_{k} - VNEW_{k}*HNEW_{k} = rnew_{k}*e_{k}^T.

// \Usage:
//  call znapps
//     ( N, KEV, NP, SHIFT, V, LDV, H, LDH, RESID, Q, LDQ,
//       WORKL, WORKD )

// \Arguments
//  N       Integer.  (INPUT)
//          Problem size, i.e. size of matrix A.

//  KEV     Integer.  (INPUT/OUTPUT)
//          KEV+NP is the size of the input matrix H.
//          KEV is the size of the updated matrix HNEW.

//  NP      Integer.  (INPUT)
//          Number of implicit shifts to be applied.

//  SHIFT   Complex*16 array of length NP.  (INPUT)
//          The shifts to be applied.

//  V       Complex*16 N by (KEV+NP) array.  (INPUT/OUTPUT)
//          On INPUT, V contains the current KEV+NP Arnoldi vectors.
//          On OUTPUT, V contains the updated KEV Arnoldi vectors
//          in the first KEV columns of V.

//  LDV     Integer.  (INPUT)
//          Leading dimension of V exactly as declared in the calling
//          program.

//  H       Complex*16 (KEV+NP) by (KEV+NP) array.  (INPUT/OUTPUT)
//          On INPUT, H contains the current KEV+NP by KEV+NP upper
//          Hessenberg matrix of the Arnoldi factorization.
//          On OUTPUT, H contains the updated KEV by KEV upper Hessenberg
//          matrix in the KEV leading submatrix.

//  LDH     Integer.  (INPUT)
//          Leading dimension of H exactly as declared in the calling
//          program.

//  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
//          On INPUT, RESID contains the the residual vector r_{k+p}.
//          On OUTPUT, RESID is the update residual vector rnew_{k}
//          in the first KEV locations.

//  Q       Complex*16 KEV+NP by KEV+NP work array.  (WORKSPACE)
//          Work array used to accumulate the rotations and reflections
//          during the bulge chase sweep.

//  LDQ     Integer.  (INPUT)
//          Leading dimension of Q exactly as declared in the calling
//          program.

//  WORKL   Complex*16 work array of length (KEV+NP).  (WORKSPACE)
//          Private (replicated) array on each PE or array allocated on
//          the front end.

//  WORKD   Complex*16 work array of length 2*N.  (WORKSPACE)
//          Distributed array used in the application of the accumulated
//          orthogonal matrix Q.

// \EndDoc

// -----------------------------------------------------------------------

// \BeginLib

// \Local variables:
//     xxxxxx  Complex*16

// \References:
//  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
//     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
//     pp 357-385.

// \Author
//     Danny Sorensen               Phuong Vu
//     Richard Lehoucq              CRPC / Rice University
//     Dept. of Computational &     Houston, Texas
//     Applied Mathematics
//     Rice University
//     Houston, Texas

// \SCCS Information: @(#)
// FILE: napps.F   SID: 2.3   DATE OF SID: 3/28/97   RELEASE: 2

// \Remarks
//  1. In this version, each shift is applied to all the sublocks of
//     the Hessenberg matrix H and not just to the submatrix that it
//     comes from. Deflation as in LAPACK routine zlahqr (QR algorithm
//     for upper Hessenberg matrices ) is used.
//     Upon output, the subdiagonals of H are enforced to be non-negative
//     real numbers.

// \EndLib

// -----------------------------------------------------------------------

void znapps_(size_t n, size_t *kev, size_t np,
	std::complex<double> *shift, std::complex<double> *v, size_t ldv, std::complex<double> *
	h, size_t ldh, std::complex<double> *resid, std::complex<double> *q, size_t
	ldq, std::complex<double> *workl, std::complex<double> *workd)
{
	using namespace std;

	/* Local variables */
	double c;
	std::complex<double> f, g;
	std::complex<double> r, s, t;
	std::complex<double> h11, h21;
	size_t iend;
	std::complex<double> sigma;
	size_t istart, kplusp;

	// Set machine-dependent constants for the
	// stopping criterion. If norm(H) <= sqrt(OVFL),
	// overflow should not occur.
	// REFERENCE: LAPACK subroutine zlahqr
	const double unfl = std::numeric_limits<double>::min();
	//const double ovfl = 1. / unfl;
	const double ulp = 2.0*std::numeric_limits<double>::epsilon();
	const double smlnum = unfl * (n / ulp);

	kplusp = *kev + np;

	// Initialize Q to the identity to accumulate
	// the rotations and reflections
	RNP::TBLAS::SetMatrix<'A'>(kplusp, kplusp, 0.0, 1.0, q, ldq);

	// Quick return if there are no shifts to apply
	if (np == 0) {
		goto L9000;
	}

	// Chase the bulge with the application of each
	// implicit shift. Each shift is applied to the
	// whole matrix including each block.
	for (size_t jj = 0; jj < np; ++jj) {
		sigma = shift[jj];

		istart = 0;
L20:

		for (size_t i = istart; i < kplusp-1; ++i) {
			// Check for splitting and deflation. Use
			// a standard test as in the QR algorithm
			// REFERENCE: LAPACK subroutine zlahqr
			double tst1 = _abs1(h[i+i*ldh]) + _abs1(h[(i+1)+(i+1)*ldh]);
			if (tst1 == 0.) {
				RNP::TLASupport::CheapHessenbergNorm<'1'>(kplusp - jj, h, ldh, &tst1);
			}
			/* Computing MAX */
			if (std::abs(h[(i+1)+i*ldh].real()) <= max(ulp * tst1,smlnum)) {
				iend = i;
				h[(i+1)+i*ldh] = 0;
				goto L40;
			}
		}
		iend = kplusp-1;
L40:
		// No reason to apply a shift to block of order 1
		// or if the current block starts after the point
		// of compression since we'll discard this stuff
		if (istart == iend || istart >= *kev) {
			goto L100;
		}

		h11 = h[istart+istart*ldh];
		h21 = h[(istart+1)+istart*ldh];
		f = h11 - sigma;
		g = h21;

		for (size_t i = istart; i < iend; ++i) {
			// Construct the plane rotation G to zero out the bulge
			RNP::TLASupport::GeneratePlaneRotation(f, g, &c, &s, &r);
			if (i > istart) {
				h[i+(i-1)*ldh] = r;
				h[(i+1)+(i-1)*ldh] = 0;
			}

			// Apply rotation to the left of H;  H <- G'*H
			for (size_t j = i; j < kplusp; ++j) {
				t = c*h[i+j*ldh] + s*h[(i+1)+j*ldh];
				h[(i+1)+j*ldh] = -std::conj(s)*h[i+j*ldh] + c * h[(i+1)+j*ldh];
				h[i+j*ldh] = t;
			}

			// Apply rotation to the right of H;  H <- H*G

			/* Computing MIN */
			{size_t jend = min(i + 2,iend);
			for (size_t j = 0; j <= jend; ++j) {
				t = c*h[j+i*ldh] + std::conj(s)*h[j+(i+1)*ldh];
				h[j+(i+1)*ldh] = -s*h[j+i*ldh] + c*h[j+(i+1)*ldh];
				h[j+i*ldh] = t;
			}}

			// Accumulate the rotation in the matrix Q;  Q <- Q*G'

			/* Computing MIN */
			{size_t jend = min(i + jj+2,kplusp);
			for (size_t j = 0; j < jend; ++j) {
				t = c*q[j+i*ldq] + std::conj(s)*q[j+(i+1)*ldq];
				q[j+(i+1)*ldq] = -s*q[j+i*ldq] + c*q[j+(i+1)*ldq];
				q[j+i*ldq] = t;
			}}

			// Prepare for next rotation
			if (i < iend - 1) {
				f = h[(i+1)+i*ldh];
				g = h[(i+2)+i*ldh];
			}
		}

		// Finished applying the shift.
L100:
		// Apply the same shift to the next block if there is any.
		istart = iend + 1;
		if (iend+1 < kplusp) {
			goto L20;
		}

		// Loop back to the top to get the next shift.
	}

	// Perform a similarity transformation that makes
	// sure that the compressed H will have non-negative
	// real subdiagonal elements.
	{size_t jend = *kev;
	for (size_t j = 0; j < jend; ++j) {
		if (h[(j+1)+j*ldh].real() < 0. || h[(j+1)+j*ldh].imag() != 0.){
			t = h[(j+1)+j*ldh] / std::abs(h[(j+1)+j*ldh]);
			RNP::TBLAS::Scale(kplusp - j, std::conj(t), &h[(j+1)+j*ldh], ldh);
			/* Computing MIN */
			size_t i__2 = min(j+3,kplusp);
			RNP::TBLAS::Scale(i__2, t, &h[0+(j+1)*ldh], 1);
			/* Computing MIN */
			i__2 = min(j+1 + np + 1,kplusp);
			RNP::TBLAS::Scale(i__2, t, &q[0+(j+1)*ldq], 1);
			h[(j+1)+j*ldh] = h[(j+1)+j*ldh].real();
		}
	}}

	{size_t ilst = *kev;
	for (size_t i = 0; i < ilst; ++i) {
		// Final check for splitting and deflation.
		// Use a standard test as in the QR algorithm
		// REFERENCE: LAPACK subroutine zlahqr.
		// Note: Since the subdiagonals of the
		// compressed H are nonnegative real numbers,
		// we take advantage of this.
		double tst1 = _abs1(h[i+i*ldh]) + _abs1(h[(i+1)+(i+1)*ldh]);
		if (tst1 == 0.) {
			RNP::TLASupport::CheapHessenbergNorm<'1'>(*kev, h, ldh, &tst1);
		}
		/* Computing MAX */
		if (h[(i+1)+i*ldh].real() <= max(ulp * tst1,smlnum)) {
			h[(i+1)+i*ldh] = 0;
		}
	}}

	// Compute the (kev+1)-st column of (V*Q) and
	// temporarily store the result in WORKD(N+1:2*N).
	// This is needed in the residual update since we
	// cannot GUARANTEE that the corresponding entry
	// of H would be zero as in exact arithmetic.
	if (h[*kev+(*kev-1)*ldh].real() > 0.) {
		RNP::TBLAS::MultMV<'N'>(n, kplusp, 1.0, v, ldv, &q[0+(*kev)*ldq], 1, 0.0, &workd[n], 1);
	}

	// Compute column 1 to kev of (V*Q) in backward order
	// taking advantage of the upper Hessenberg structure of Q.
	{size_t ilst = *kev;
	for (size_t i = 0; i < ilst; ++i) {
		RNP::TBLAS::MultMV<'N'>(n, kplusp - i, 1.0, v, ldv, &q[0+(*kev - i-1)*ldq], 1, 0.0, workd, 1);
		RNP::TBLAS::Copy(n, workd, 1, &v[0+(kplusp - i-1)*ldv], 1);
	}}

	//  Move v(:,kplusp-kev+1:kplusp) into v(:,1:kev).
	RNP::TBLAS::CopyMatrix<'A'>(n, *kev, &v[0+(kplusp - *kev)*ldv], ldv, v, ldv);

	// Copy the (kev+1)-st column of (V*Q) in the appropriate place
	if (h[*kev+(*kev-1)*ldh].real() > 0.) {
		RNP::TBLAS::Copy(n, &workd[n], 1, &v[0+(*kev)*ldv], 1);
	}

	// Update the residual vector:
	//    r <- sigmak*r + betak*v(:,kev+1)
	// where
	//    sigmak = (e_{kev+p}'*Q)*e_{kev}
	//    betak = e_{kev+1}'*H*e_{kev}
	RNP::TBLAS::Scale(n, q[(kplusp-1)+(*kev-1)*ldq], resid, 1);
	if (h[*kev+(*kev-1)*ldh].real() > 0.) {
		RNP::TBLAS::Axpy(n, h[*kev+(*kev-1)*ldh], &v[0+(*kev)*ldv], 1, resid, 1);
	}

L9000:
	return;
}


// \Name: zgetv0

// \Description:
//  Generate a random initial residual vector for the Arnoldi process.
//  Force the residual vector to be in the range of the operator OP.

// \Usage:
//  call zgetv0
//     ( IDO, BMAT, ITRY, INITV, N, J, V, LDV, RESID, RNORM,
//       IPNTR, WORKD, IERR )

// \Arguments
//  IDO     int.  (INPUT/OUTPUT)
//          Reverse communication flag.  IDO must be zero on the first
//          call to zgetv0.
//          -------------------------------------------------------------
//          IDO =  0: first call to the reverse communication interface
//          IDO = -1: compute  Y = OP * X  where
//                    IPNTR(1) is the pointer into WORKD for X,
//                    IPNTR(2) is the pointer into WORKD for Y.
//                    This is for the initialization phase to force the
//                    starting vector into the range of OP.
//          IDO =  2: compute  Y = B * X  where
//                    IPNTR(1) is the pointer into WORKD for X,
//                    IPNTR(2) is the pointer into WORKD for Y.
//          IDO = 99: done
//          -------------------------------------------------------------

//  BMAT    Character*1.  (INPUT)
//          BMAT specifies the type of the matrix B in the (generalized)
//          eigenvalue problem A*x = lambda*B*x.
//          B = 'I' -> standard eigenvalue problem A*x = lambda*x
//          B = 'G' -> generalized eigenvalue problem A*x = lambda*B*x

//  ITRY    int.  (INPUT)
//          ITRY counts the number of times that zgetv0 is called.
//          It should be set to 1 on the initial call to zgetv0.

//  INITV   bool variable.  (INPUT)
//          .TRUE.  => the initial residual vector is given in RESID.
//          .FALSE. => generate a random initial residual vector.

//  N       int.  (INPUT)
//          Dimension of the problem.

//  J       int.  (INPUT)
//          Index of the residual vector to be generated, with respect to
//          the Arnoldi process.  J > 1 in case of a "restart".

//  V       Complex*16 N by J array.  (INPUT)
//          The first J-1 columns of V contain the current Arnoldi basis
//          if this is a "restart".

//  LDV     int.  (INPUT)
//          Leading dimension of V exactly as declared in the calling
//          program.

//  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
//          Initial residual vector to be generated.  If RESID is
//          provided, force RESID into the range of the operator OP.

//  RNORM   Double precision scalar.  (OUTPUT)
//          B-norm of the generated residual.

//  IPNTR   int array of length 3.  (OUTPUT)

//  WORKD   Complex*16 work array of length 2*N.  (REVERSE COMMUNICATION).
//          On exit, WORK(1:N) = B*RESID to be used in SSAITR.

//  IERR    int.  (OUTPUT)
//          =  0: Normal exit.
//          = -1: Cannot generate a nontrivial restarted residual vector
//                in the range of the operator OP.


// -----------------------------------------------------------------------

// \References:
//  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
//     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
//     pp 357-385.


// \Author
//     Danny Sorensen               Phuong Vu
//     Richard Lehoucq              CRPC / Rice University
//     Dept. of Computational &     Houston, Texas
//     Applied Mathematics
//     Rice University
//     Houston, Texas

// -----------------------------------------------------------------------

int zgetv0_(arz_data *data, int *ido, char bmat, size_t itry, bool
	initv, size_t n, size_t j, const std::complex<double> *v, size_t ldv,
	std::complex<double> *resid, double *rnorm, int *ipntr,
	std::complex<double> *workd)
{
	if(*ido == 0){
		data->s_zgetv0.iter = 0;
		data->s_zgetv0.first = false;
		data->s_zgetv0.orth = false;

		// Possibly generate a random starting vector in RESID
		// Use a LAPACK random number generator used by the
		// matrix generation routines.
		//    idist = 1: uniform (0,1)  distribution;
		//    idist = 2: uniform (-1,1) distribution;
		//    idist = 3: normal  (0,1)  distribution;
		if(!initv){
			RNP::TBLAS::RandomVector<2>(n, resid);
		}

		// Force the starting vector into the range of OP to handle
		// the generalized problem when B is possibly (singular).
		if(bmat == 'G'){
			ipntr[0] = 0;
			ipntr[1] = n + ipntr[0];
			RNP::TBLAS::Copy(n, resid, 1, workd, 1);
			*ido = -1;
			return 0;
		}
    }

	// Back from computing B*(initial-vector)
	if(data->s_zgetv0.first){
		goto L20;
	}

	// Back from computing B*(orthogonalized-vector)
	if(data->s_zgetv0.orth){
		goto L40;
	}


	// Starting vector is now in the range of OP; r = OP*r;
	// Compute B-norm of starting vector.
	data->s_zgetv0.first = true;
	if(bmat == 'G'){
		RNP::TBLAS::Copy(n, &workd[n], 1, resid, 1);
		ipntr[1] = 0;
		ipntr[0] = n + ipntr[1];
		*ido = 2;
		return 0;
	}else if(bmat == 'I'){
		RNP::TBLAS::Copy(n, resid, 1, workd, 1);
	}

L20:
	data->s_zgetv0.first = false;
	if(bmat == 'G'){
		std::complex<double> cnorm = RNP::TBLAS::ConjugateDot(n, resid, 1, workd, 1);
		data->s_zgetv0.rnorm0 = sqrt(std::abs(cnorm));
	}else if(bmat == 'I'){
		data->s_zgetv0.rnorm0 = RNP::TBLAS::Norm2(n, resid, 1);
	}
	*rnorm = data->s_zgetv0.rnorm0;

	// Exit if this is the very first Arnoldi step
	if(j == 0){
		goto L50;
	}

	// Otherwise need to B-orthogonalize the starting vector against
	// the current Arnoldi basis using Gram-Schmidt with iter. ref.
	// This is the case where an invariant subspace is encountered
	// in the middle of the Arnoldi factorization.
	//
	//       s = V^{T}*B*r;   r = r - V*s;
	//
	// Stopping criteria used for iter. ref. is discussed in
	// Parlett's book, page 107 and in Gragg & Reichel TOMS paper.
	data->s_zgetv0.orth = true;
L30:
	RNP::TBLAS::MultMV<'C'>(n, j, 1.0, v, ldv, workd, 1, 0.0, &workd[n], 1);
	RNP::TBLAS::MultMV<'N'>(n, j, -1.0, v, ldv, &workd[n], 1, 1.0, resid, 1);

	// Compute the B-norm of the orthogonalized starting vector
	if(bmat == 'G'){
		RNP::TBLAS::Copy(n, resid, 1, &workd[n], 1);
		ipntr[1] = 0;
		ipntr[0] = n + ipntr[1];
		*ido = 2;
		return 0;
	}else if(bmat == 'I'){
		RNP::TBLAS::Copy(n, resid, 1, workd, 1);
	}

L40:

	if(bmat == 'G'){
		std::complex<double> cnorm = RNP::TBLAS::ConjugateDot(n, resid, 1, workd, 1);
		*rnorm = sqrt(std::abs(cnorm));
	}else if(bmat == 'I'){
		*rnorm = RNP::TBLAS::Norm2(n, resid, 1);
	}

	// Check for further orthogonalization.
	if(*rnorm > data->s_zgetv0.rnorm0 * .717f){
		goto L50;
	}

	++data->s_zgetv0.iter;
	if(data->s_zgetv0.iter <= 1){
		// Perform iterative refinement step
		data->s_zgetv0.rnorm0 = *rnorm;
		goto L30;
	}else{
		// Iterative refinement step "failed"
		for(size_t jj = 0; jj < n; ++jj){
			resid[jj] = 0;
		}
		*rnorm = 0.;
		*ido = 99;
		return -1;
	}

L50:
	*ido = 99;

	return 0;
}

// \BeginDoc

// \Name: znaitr

// \Description:
//  Reverse communication interface for applying NP additional steps to
//  a K step nonsymmetric Arnoldi factorization.

//  Input:  OP*V_{k}  -  V_{k}*H = r_{k}*e_{k}^T

//          with (V_{k}^T)*B*V_{k} = I, (V_{k}^T)*B*r_{k} = 0.

//  Output: OP*V_{k+p}  -  V_{k+p}*H = r_{k+p}*e_{k+p}^T

//          with (V_{k+p}^T)*B*V_{k+p} = I, (V_{k+p}^T)*B*r_{k+p} = 0.

//  where OP and B are as in znaupd.  The B-norm of r_{k+p} is also
//  computed and returned.

// \Usage:
//  call znaitr
//     ( IDO, BMAT, N, K, NP, NB, RESID, RNORM, V, LDV, H, LDH,
//       IPNTR, WORKD, INFO )

// \Arguments
//  IDO     int.  (INPUT/OUTPUT)
//          Reverse communication flag.
//          -------------------------------------------------------------
//          IDO =  0: first call to the reverse communication interface
//          IDO = -1: compute  Y = OP * X  where
//                    IPNTR(1) is the pointer into WORK for X,
//                    IPNTR(2) is the pointer into WORK for Y.
//                    This is for the restart phase to force the new
//                    starting vector into the range of OP.
//          IDO =  1: compute  Y = OP * X  where
//                    IPNTR(1) is the pointer into WORK for X,
//                    IPNTR(2) is the pointer into WORK for Y,
//                    IPNTR(3) is the pointer into WORK for B * X.
//          IDO =  2: compute  Y = B * X  where
//                    IPNTR(1) is the pointer into WORK for X,
//                    IPNTR(2) is the pointer into WORK for Y.
//          IDO = 99: done
//          -------------------------------------------------------------
//          When the routine is used in the "shift-and-invert" mode, the
//          vector B * Q is already available and do not need to be
//          recomputed in forming OP * Q.

//  BMAT    Character*1.  (INPUT)
//          BMAT specifies the type of the matrix B that defines the
//          semi-inner product for the operator OP.  See znaupd.
//          B = 'I' -> standard eigenvalue problem A*x = lambda*x
//          B = 'G' -> generalized eigenvalue problem A*x = lambda*M**x

//  N       int.  (INPUT)
//          Dimension of the eigenproblem.

//  K       int.  (INPUT)
//          Current size of V and H.

//  NP      int.  (INPUT)
//          Number of additional Arnoldi steps to take.

//  NB      int.  (INPUT)
//          Blocksize to be used in the recurrence.
//          Only work for NB = 1 right now.  The goal is to have a
//          program that implement both the block and non-block method.

//  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
//          On INPUT:  RESID contains the residual vector r_{k}.
//          On OUTPUT: RESID contains the residual vector r_{k+p}.

//  RNORM   Double precision scalar.  (INPUT/OUTPUT)
//          B-norm of the starting residual on input.
//          B-norm of the updated residual r_{k+p} on output.

//  V       Complex*16 N by K+NP array.  (INPUT/OUTPUT)
//          On INPUT:  V contains the Arnoldi vectors in the first K
//          columns.
//          On OUTPUT: V contains the new NP Arnoldi vectors in the next
//          NP columns.  The first K columns are unchanged.

//  LDV     int.  (INPUT)
//          Leading dimension of V exactly as declared in the calling
//          program.

//  H       Complex*16 (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
//          H is used to store the generated upper Hessenberg matrix.

//  LDH     int.  (INPUT)
//          Leading dimension of H exactly as declared in the calling
//          program.

//  IPNTR   int array of length 3.  (OUTPUT)
//          Pointer to mark the starting locations in the WORK for
//          vectors used by the Arnoldi iteration.
//          -------------------------------------------------------------
//          IPNTR(1): pointer to the current operand vector X.
//          IPNTR(2): pointer to the current result vector Y.
//          IPNTR(3): pointer to the vector B * X when used in the
//                    shift-and-invert mode.  X is the current operand.
//          -------------------------------------------------------------

//  WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
//          Distributed array to be used in the basic Arnoldi iteration
//          for reverse communication.  The calling program should not
//          use WORKD as temporary workspace during the iteration !!!!!!
//          On input, WORKD(1:N) = B*RESID and is used to save some
//          computation at the first step.

//  INFO    int.  (OUTPUT)
//          = 0: Normal exit.
//          > 0: Size of the spanning invariant subspace of OP found.

// \EndDoc

// -----------------------------------------------------------------------

// \BeginLib

// \Local variables:
//     xxxxxx  Complex*16

// \References:
//  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
//     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
//     pp 357-385.
//  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
//     Restarted Arnoldi Iteration", Rice University Technical Report
//     TR95-13, Department of Computational and Applied Mathematics.

// \Routines called:
//     zgetv0  ARPACK routine to generate the initial vector.

// \Author
//     Danny Sorensen               Phuong Vu
//     Richard Lehoucq              CRPC / Rice University
//     Dept. of Computational &     Houston, Texas
//     Applied Mathematics
//     Rice University
//     Houston, Texas

// \SCCS Information: @(#)
// FILE: naitr.F   SID: 2.3   DATE OF SID: 8/27/96   RELEASE: 2

// \Remarks
//  The algorithm implemented is:

//  restart = .false.
//  Given V_{k} = [v_{1}, ..., v_{k}], r_{k};
//  r_{k} contains the initial residual vector even for k = 0;
//  Also assume that rnorm = || B*r_{k} || and B*r_{k} are already
//  computed by the calling program.

//  betaj = rnorm ; p_{k+1} = B*r_{k} ;
//  For  j = k+1, ..., k+np  Do
//     1) if ( betaj < tol ) stop or restart depending on j.
//        ( At present tol is zero )
//        if ( restart ) generate a new starting vector.
//     2) v_{j} = r(j-1)/betaj;  V_{j} = [V_{j-1}, v_{j}];
//        p_{j} = p_{j}/betaj
//     3) r_{j} = OP*v_{j} where OP is defined as in znaupd
//        For shift-invert mode p_{j} = B*v_{j} is already available.
//        wnorm = || OP*v_{j} ||
//     4) Compute the j-th step residual vector.
//        w_{j} =  V_{j}^T * B * OP * v_{j}
//        r_{j} =  OP*v_{j} - V_{j} * w_{j}
//        H(:,j) = w_{j};
//        H(j,j-1) = rnorm
//        rnorm = || r_(j) ||
//        If (rnorm > 0.717*wnorm) accept step and go back to 1)
//     5) Re-orthogonalization step:
//        s = V_{j}'*B*r_{j}
//        r_{j} = r_{j} - V_{j}*s;  rnorm1 = || r_{j} ||
//        alphaj = alphaj + s_{j};
//     6) Iterative refinement step:
//        If (rnorm1 > 0.717*rnorm) then
//           rnorm = rnorm1
//           accept step and go back to 1)
//        Else
//           rnorm = rnorm1
//           If this is the first time in step 6), go to 5)
//           Else r_{j} lies in the span of V_{j} numerically.
//              Set r_{j} = 0 and rnorm = 0; go to 1)
//        EndIf
//  End Do

// \EndLib

// -----------------------------------------------------------------------

int znaitr_(arz_data *data, int *ido, char bmat, size_t n, size_t k,
	 size_t np, size_t nb, std::complex<double> *resid, double *rnorm,
	std::complex<double> *v, size_t ldv, std::complex<double> *h, size_t ldh,
	int *ipntr, std::complex<double> *workd)
{
	using namespace std;

	int info = 0;

	extern int zgetv0_(arz_data *data, int *ido, char bmat, size_t itry, bool
		initv, size_t n, size_t j, const std::complex<double> *v, size_t ldv,
		std::complex<double> *resid, double *rnorm, int *ipntr,
		std::complex<double> *workd);


	// Set machine-dependent constants for the
	// the splitting and deflation criterion.
	// If norm(H) <= sqrt(OVFL),
	// overflow should not occur.
	// REFERENCE: LAPACK subroutine zlahqr
	static const double unfl = std::numeric_limits<double>::min();
	//const double ovfl = 1. / unfl;
	static const double ulp = 2.*std::numeric_limits<double>::epsilon();
	const double smlnum = unfl * (n / ulp);

	if (*ido == 0) {
		// Initial call to this routine

		info = 0;
		data->s_znaitr.step3 = false;
		data->s_znaitr.step4 = false;
		data->s_znaitr.rstart = false;
		data->s_znaitr.orth1 = false;
		data->s_znaitr.orth2 = false;
		data->s_znaitr.j = k;
		data->s_znaitr.ipj = 0;
		data->s_znaitr.irj = data->s_znaitr.ipj + n;
		data->s_znaitr.ivj = data->s_znaitr.irj + n;
	}
	
	std::complex<double> *work_pj = &workd[data->s_znaitr.ipj];
	std::complex<double> *work_rj = &workd[data->s_znaitr.irj];
	std::complex<double> *work_vj = &workd[data->s_znaitr.ivj];

	// When in reverse communication mode one of:
	// STEP3, STEP4, ORTH1, ORTH2, RSTART
	// will be .true. when ....
	// STEP3: return from computing OP*v_{j}.
	// STEP4: return from computing B-norm of OP*v_{j}
	// ORTH1: return from computing B-norm of r_{j+1}
	// ORTH2: return from computing B-norm of
	//        correction to the residual vector.
	// RSTART: return from OP computations needed by
	//         zgetv0.

	if (data->s_znaitr.step3) {
		goto L50;
	}
	if (data->s_znaitr.step4) {
		goto L60;
	}
	if (data->s_znaitr.orth1) {
		goto L70;
	}
	if (data->s_znaitr.orth2) {
		goto L90;
	}
	if (data->s_znaitr.rstart) {
		goto L30;
	}

	// Else this is the first step

	//        A R N O L D I     I T E R A T I O N     L O O P
	//
	// Note:  B*r_{j-1} is already in WORKD(1:N)=WORKD(IPJ:IPJ+N-1)
L1000:

	// STEP 1: Check if the B norm of j-th residual
	// vector is zero. Equivalent to determine whether
	// an exact j-step Arnoldi factorization is present.
	data->s_znaitr.betaj = *rnorm;
	if (*rnorm > 0.) {
		goto L40;
	}

	// Invariant subspace found, generate a new starting
	// vector which is orthogonal to the current Arnoldi
	// basis and continue the iteration.

	// ITRY is the loop variable that controls the
	// maximum amount of times that a restart is
	// attempted. NRSTRT is used by stat.h

	data->s_znaitr.betaj = 0.;
	data->s_znaitr.itry = 1;
L20:
	data->s_znaitr.rstart = true;
	*ido = 0;
L30:
	// If in reverse communication mode and
	// RSTART = .true. flow returns here.
	data->s_znaitr.ierr = zgetv0_(data, ido, bmat, data->s_znaitr.itry, false, n, data->s_znaitr.j, v, ldv, resid, rnorm, ipntr, workd);
	if (*ido != 99) {
		goto L9000;
	}
	if (data->s_znaitr.ierr < 0) {
		++data->s_znaitr.itry;
		if (data->s_znaitr.itry <= 3) {
			goto L20;
		}

		// Give up after several restart attempts.
		// Set INFO to the size of the invariant subspace
		// which spans OP and exit.
		info = (int)data->s_znaitr.j;
		*ido = 99;
		goto L9000;
	}

L40:

	// STEP 2:  v_{j} = r_{j-1}/rnorm and p_{j} = p_{j}/rnorm
	// Note that p_{j} = B*r_{j-1}. In order to avoid overflow
	// when reciprocating a small RNORM, test against lower
	// machine bound.

	RNP::TBLAS::Copy(n, resid, 1, &v[0+data->s_znaitr.j*ldv], 1);
	if (*rnorm >= unfl) {
		double temp1 = 1. / *rnorm;
		RNP::TBLAS::Scale(n, temp1, &v[0+data->s_znaitr.j*ldv], 1);
		RNP::TBLAS::Scale(n, temp1, work_pj, 1);
	} else {
		// To scale both v_{j} and p_{j} carefully
		// use LAPACK routine zlascl
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, *rnorm, 1.0, n, 1, &v[0+data->s_znaitr.j*ldv], n);
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, *rnorm, 1.0, n, 1, work_pj, n);
	}

	// STEP 3:  r_{j} = OP*v_{j}; Note that p_{j} = B*v_{j}
	// Note that this is not quite yet r_{j}. See STEP 4

	data->s_znaitr.step3 = true;
	RNP::TBLAS::Copy(n, &v[0+data->s_znaitr.j*ldv], 1, work_vj, 1);
	ipntr[0] = data->s_znaitr.ivj;
	ipntr[1] = data->s_znaitr.irj;
	ipntr[2] = data->s_znaitr.ipj;
	*ido = 1;

	// Exit in order to compute OP*v_{j}
	goto L9000;
L50:

	// Back from reverse communication;
	// WORKD(IRJ:IRJ+N-1) := OP*v_{j}
	// if step3 = .true.
	data->s_znaitr.step3 = false;

	// Put another copy of OP*v_{j} into RESID.
	RNP::TBLAS::Copy(n, work_rj, 1, resid, 1);

	// STEP 4:  Finish extending the Arnoldi
	//          factorization to length j.
	if (bmat == 'G') {
		data->s_znaitr.step4 = true;
		ipntr[0] = data->s_znaitr.irj;
		ipntr[1] = data->s_znaitr.ipj;
		*ido = 2;

		// Exit in order to compute B*OP*v_{j}
		goto L9000;
	} else if (bmat == 'I') {
		RNP::TBLAS::Copy(n, resid, 1, work_pj, 1);
	}
L60:

	// Back from reverse communication;
	// WORKD(IPJ:IPJ+N-1) := B*OP*v_{j}
	// if step4 = .true.
	data->s_znaitr.step4 = false;

	// The following is needed for STEP 5.
	// Compute the B-norm of OP*v_{j}.

	if (bmat == 'G') {
		std::complex<double> cnorm = RNP::TBLAS::ConjugateDot(n, resid, 1, work_pj, 1);
		data->s_znaitr.wnorm = sqrt(std::abs(cnorm));
	} else if (bmat == 'I') {
		data->s_znaitr.wnorm = RNP::TBLAS::Norm2(n, resid, 1);
	}

	// Compute the j-th residual corresponding
	// to the j step factorization.
	// Use Classical Gram Schmidt and compute:
	// w_{j} <-  V_{j}^T * B * OP * v_{j}
	// r_{j} <-  OP*v_{j} - V_{j} * w_{j}

	// Compute the j Fourier coefficients w_{j}
	// WORKD(IPJ:IPJ+N-1) contains B*OP*v_{j}.
	RNP::TBLAS::MultMV<'C'>(n, data->s_znaitr.j+1, 1.0, v, ldv, work_pj, 1, 0.0, &h[0+data->s_znaitr.j*ldh], 1);

	// Orthogonalize r_{j} against V_{j}.
	// RESID contains OP*v_{j}. See STEP 3.
	RNP::TBLAS::MultMV<'N'>(n, data->s_znaitr.j+1, -1.0, v, ldv, &h[0+data->s_znaitr.j*ldh], 1, 1.0, resid, 1);

	if (data->s_znaitr.j > 0) {
		h[data->s_znaitr.j+(data->s_znaitr.j-1)*ldh] = data->s_znaitr.betaj;
	}

	data->s_znaitr.orth1 = true;

	if (bmat == 'G') {
		RNP::TBLAS::Copy(n, resid, 1, work_rj, 1);
		ipntr[0] = data->s_znaitr.irj;
		ipntr[1] = data->s_znaitr.ipj;
		*ido = 2;

		// Exit in order to compute B*r_{j}
		goto L9000;
	} else if (bmat == 'I') {
		RNP::TBLAS::Copy(n, resid, 1, work_pj, 1);
	}
L70:
	// Back from reverse communication if ORTH1 = .true.
	// WORKD(IPJ:IPJ+N-1) := B*r_{j}.
	data->s_znaitr.orth1 = false;

	// Compute the B-norm of r_{j}.
	if (bmat == 'G') {
		std::complex<double> cnorm = RNP::TBLAS::ConjugateDot(n, resid, 1, work_pj, 1);
		*rnorm = sqrt(std::abs(cnorm));
	} else if (bmat == 'I') {
		*rnorm = RNP::TBLAS::Norm2(n, resid, 1);
	}

	// STEP 5: Re-orthogonalization / Iterative refinement phase
	// Maximum NITER_ITREF tries.
	//
	//          s      = V_{j}^T * B * r_{j}
	//          r_{j}  = r_{j} - V_{j}*s
	//          alphaj = alphaj + s_{j}
	//
	// The stopping criteria used for iterative refinement is
	// discussed in Parlett's book SEP, page 107 and in Gragg &
	// Reichel ACM TOMS paper; Algorithm 686, Dec. 1990.
	// Determine if we need to correct the residual. The goal is
	// to enforce ||v(:,1:j)^T * r_{j}|| .le. eps * || r_{j} ||
	// The following test determines whether the sine of the
	// angle between  OP*x and the computed residual is less
	// than or equal to 0.717.
	if (*rnorm > data->s_znaitr.wnorm * .717f) {
		goto L100;
	}

	data->s_znaitr.iter = 0;

	// Enter the Iterative refinement phase. If further
	// refinement is necessary, loop back here. The loop
	// variable is ITER. Perform a step of Classical
	// Gram-Schmidt using all the Arnoldi vectors V_{j}
L80:
	// Compute V_{j}^T * B * r_{j}.
	// WORKD(IRJ:IRJ+J-1) = v(:,1:J)'*WORKD(IPJ:IPJ+N-1).

	RNP::TBLAS::MultMV<'C'>(n, data->s_znaitr.j+1, 1.0, v, ldv, work_pj, 1, 0.0, work_rj, 1);

	// Compute the correction to the residual:
	// r_{j} = r_{j} - V_{j} * WORKD(IRJ:IRJ+J-1).
	// The correction to H is v(:,1:J)*H(1:J,1:J)
	// + v(:,1:J)*WORKD(IRJ:IRJ+J-1)*e'_j.

	RNP::TBLAS::MultMV<'N'>(n, data->s_znaitr.j+1, -1.0, v, ldv, work_rj, 1, 1.0, resid, 1);
	RNP::TBLAS::Axpy(data->s_znaitr.j+1, 1.0, work_rj, 1, &h[0+data->s_znaitr.j*ldh], 1);

	data->s_znaitr.orth2 = true;
	if (bmat == 'G') {
		RNP::TBLAS::Copy(n, resid, 1, work_rj, 1);
		ipntr[0] = data->s_znaitr.irj;
		ipntr[1] = data->s_znaitr.ipj;
		*ido = 2;

		// Exit in order to compute B*r_{j}.
		// r_{j} is the corrected residual.
		goto L9000;
	} else if (bmat == 'I') {
		RNP::TBLAS::Copy(n, resid, 1, work_pj, 1);
	}
L90:
	// Back from reverse communication if ORTH2 = .true.

	// Compute the B-norm of the corrected residual r_{j}.
	if (bmat == 'G') {
		std::complex<double> cnorm = RNP::TBLAS::ConjugateDot(n, resid, 1, work_pj, 1);
		data->s_znaitr.rnorm1 = sqrt(std::abs(cnorm));
	} else if (bmat == 'I') {
		data->s_znaitr.rnorm1 = RNP::TBLAS::Norm2(n, resid, 1);
	}

	// Determine if we need to perform another
	// step of re-orthogonalization.
	if (data->s_znaitr.rnorm1 > *rnorm * .717f) {
		// No need for further refinement.
		// The cosine of the angle between the
		// corrected residual vector and the old
		// residual vector is greater than 0.717
		// In other words the corrected residual
		// and the old residual vector share an
		// angle of less than arcCOS(0.717)
		*rnorm = data->s_znaitr.rnorm1;
	} else {
		// Another step of iterative refinement step
		// is required. NITREF is used by stat.h
		*rnorm = data->s_znaitr.rnorm1;
		++data->s_znaitr.iter;
		if (data->s_znaitr.iter <= 1) {
			goto L80;
		}

		// Otherwise RESID is numerically in the span of V
		for(size_t jj = 0; jj < n; ++jj){
			resid[jj] = 0;
		}
		*rnorm = 0.;
	}

	// Branch here directly if iterative refinement
	// wasn't necessary or after at most NITER_REF
	// steps of iterative refinement.
L100:
	data->s_znaitr.rstart = false;
	data->s_znaitr.orth2 = false;

	// STEP 6: Update  j = j+1;  Continue
	++data->s_znaitr.j;
	if (data->s_znaitr.j >= k + np) {
		*ido = 99;
		for (size_t i = max((size_t)1,k); i < k + np; ++i) {
			// Check for splitting and deflation.
			// Use a standard test as in the QR algorithm
			// REFERENCE: LAPACK subroutine zlahqr
			double tst1 = std::abs(h[(i-1)+(i-1)*ldh]) + std::abs(h[i+i*ldh]);
			if (tst1 == 0.) {
				RNP::TLASupport::CheapHessenbergNorm<'1'>(k + np, h, ldh, &tst1);
			}
			/* Computing MAX */
			if (std::abs(h[i+i*ldh]) <= max(ulp * tst1,smlnum)) {
				h[i+(i-1)*ldh] = 0;
			}
		}

		goto L9000;
	}

	// Loop back to extend the factorization by another step.
	goto L1000;

	//  E N D     O F     M A I N     I T E R A T I O N     L O O P

L9000:
	return info;
}

// \BeginDoc

// \Name: znaup2

// \Description:
//  Intermediate level interface called by znaupd .

// \Usage:
//  call znaup2
//     ( IDO, BMAT, N, WHICH, NEV, NP, TOL, RESID, MODE, IUPD,
//       ISHIFT, MXITER, V, LDV, H, LDH, RITZ, BOUNDS,
//       Q, LDQ, WORKL, IPNTR, WORKD, RWORK, INFO )

// \Arguments

//  IDO, BMAT, N, WHICH, NEV, TOL, RESID: same as defined in znaupd .
//  MODE, ISHIFT, MXITER: see the definition of IPARAM in znaupd .

//  NP      int.  (INPUT/OUTPUT)
//          Contains the number of implicit shifts to apply during
//          each Arnoldi iteration.
//          If ISHIFT=1, NP is adjusted dynamically at each iteration
//          to accelerate convergence and prevent stagnation.
//          This is also roughly equal to the number of matrix-vector
//          products (involving the operator OP) per Arnoldi iteration.
//          The logic for adjusting is contained within the current
//          subroutine.
//          If ISHIFT=0, NP is the number of shifts the user needs
//          to provide via reverse comunication. 0 < NP < NCV-NEV.
//          NP may be less than NCV-NEV since a leading block of the current
//          upper Hessenberg matrix has split off and contains "unwanted"
//          Ritz values.
//          Upon termination of the IRA iteration, NP contains the number
//          of "converged" wanted Ritz values.

//  IUPD    int.  (INPUT)
//          IUPD .EQ. 0: use explicit restart instead implicit update.
//          IUPD .NE. 0: use implicit update.

//  V       Complex*16  N by (NEV+NP) array.  (INPUT/OUTPUT)
//          The Arnoldi basis vectors are returned in the first NEV
//          columns of V.

//  LDV     int.  (INPUT)
//          Leading dimension of V exactly as declared in the calling
//          program.

//  H       Complex*16  (NEV+NP) by (NEV+NP) array.  (OUTPUT)
//          H is used to store the generated upper Hessenberg matrix

//  LDH     int.  (INPUT)
//          Leading dimension of H exactly as declared in the calling
//          program.

//  RITZ    Complex*16  array of length NEV+NP.  (OUTPUT)
//          RITZ(1:NEV)  contains the computed Ritz values of OP.

//  BOUNDS  Complex*16  array of length NEV+NP.  (OUTPUT)
//          BOUNDS(1:NEV) contain the error bounds corresponding to
//          the computed Ritz values.

//  Q       Complex*16  (NEV+NP) by (NEV+NP) array.  (WORKSPACE)
//          Private (replicated) work array used to accumulate the
//          rotation in the shift application step.

//  LDQ     int.  (INPUT)
//          Leading dimension of Q exactly as declared in the calling
//          program.

//  WORKL   Complex*16  work array of length at least
//          (NEV+NP)**2 + 3*(NEV+NP).  (WORKSPACE)
//          Private (replicated) array on each PE or array allocated on
//          the front end.  It is used in shifts calculation, shifts
//          application and convergence checking.


//  IPNTR   int array of length 3.  (OUTPUT)
//          Pointer to mark the starting locations in the WORKD for
//          vectors used by the Arnoldi iteration.
//          -------------------------------------------------------------
//          IPNTR(1): pointer to the current operand vector X.
//          IPNTR(2): pointer to the current result vector Y.
//          IPNTR(3): pointer to the vector B * X when used in the
//                    shift-and-invert mode.  X is the current operand.
//          -------------------------------------------------------------

//  WORKD   Complex*16  work array of length 3*N.  (WORKSPACE)
//          Distributed array to be used in the basic Arnoldi iteration
//          for reverse communication.  The user should not use WORKD
//          as temporary workspace during the iteration !!!!!!!!!!
//          See Data Distribution Note in ZNAUPD .

//  RWORK   Double precision    work array of length  NEV+NP ( WORKSPACE)
//          Private (replicated) array on each PE or array allocated on
//          the front end.

//  INFO    int.  (INPUT/OUTPUT)
//          If INFO .EQ. 0, a randomly initial residual vector is used.
//          If INFO .NE. 0, RESID contains the initial residual vector,
//                          possibly from a previous run.
//          Error flag on output.
//          =     0: Normal return.
//          =     1: Maximum number of iterations taken.
//                   All possible eigenvalues of OP has been found.
//                   NP returns the number of converged Ritz values.
//          =     2: No shifts could be applied.
//          =    -8: Error return from LAPACK eigenvalue calculation;
//                   This should never happen.
//          =    -9: Starting vector is zero.
//          = -9999: Could not build an Arnoldi factorization.
//                   Size that was built in returned in NP.

// \EndDoc

// -----------------------------------------------------------------------

// \BeginLib

// \Local variables:
//     xxxxxx  Complex*16

// \References:
//  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
//     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
//     pp 357-385.
//  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
//     Restarted Arnoldi Iteration", Rice University Technical Report
//     TR95-13, Department of Computational and Applied Mathematics.

// \Routines called:
//     zgetv0   ARPACK initial vector generation routine.
//     znaitr   ARPACK Arnoldi factorization routine.
//     znapps   ARPACK application of implicit shifts routine.
//     zneigh   ARPACK compute Ritz values and error bounds routine.
//     zngets   ARPACK reorder Ritz values and error bounds routine.
//     zsortc   ARPACK sorting routine.

// \Author
//     Danny Sorensen               Phuong Vu
//     Richard Lehoucq              CRPC / Rice Universitya
//     Chao Yang                    Houston, Texas
//     Dept. of Computational &
//     Applied Mathematics
//     Rice University
//     Houston, Texas

// \SCCS Information: @(#)
// FILE: naup2.F   SID: 2.6   DATE OF SID: 06/01/00   RELEASE: 2

// \Remarks
//     1. None

// \EndLib

// -----------------------------------------------------------------------

void znaup2_(arz_data *data, int *ido, char bmat, size_t n, RNP::IRA::SortFunc sorter,
	size_t *nev, size_t *np, double tol, std::complex<double> *
	resid, int mode, int ishift, size_t *mxiter,
	 std::complex<double> *v, size_t ldv, std::complex<double> *h, size_t ldh,
	std::complex<double> *ritz, std::complex<double> *bounds, std::complex<double> *q, size_t
	ldq, std::complex<double> *workl, int *ipntr, std::complex<double> *workd,
	double *rwork, int *info)
{
	using namespace std;

	/* Function Body */
	if (*ido == 0) {
		data->s_znaup2.nev0 = *nev;
		data->s_znaup2.np0 = *np;

		// kplusp is the bound on the largest
		//        Lanczos factorization built.
		// nconv is the current number of
		//        "converged" eigenvalues.
		// iter is the counter on the current
		//      iteration step.
		data->s_znaup2.kplusp = *nev + *np;
		data->s_znaup2.nconv = 0;
		data->s_znaup2.iter = 0;

		// Get machine dependent constant.
		data->s_znaup2.eps23 = pow(std::numeric_limits<double>::epsilon(), (double)2/(double)3);

		// Set flags for computing the first NEV
		// steps of the Arnoldi factorization.
		data->s_znaup2.getv0 = true;
		data->s_znaup2.update = false;
		data->s_znaup2.ushift = false;
		data->s_znaup2.cnorm = false;

		if (*info != 0) {
			// User provides the initial residual vector.
			data->s_znaup2.initv = true;
			*info = 0;
		} else {
			data->s_znaup2.initv = false;
		}
	}

	// Get a possibly random starting vector and
	// force it into the range of the operator OP.
	if (data->s_znaup2.getv0) {
		*info = zgetv0_(data, ido, bmat, 1, data->s_znaup2.initv, n, 0, v, ldv, resid, &data->s_znaup2.rnorm, ipntr, workd);
		if (*ido != 99) {
			goto L9000;
		}

		if (data->s_znaup2.rnorm == 0.) {
			// The initial vector is zero. Error exit.
			*info = -9;
			goto L1100;
		}
		data->s_znaup2.getv0 = false;
		*ido = 0;
	}

	// Back from reverse communication :
	// continue with update step
	if (data->s_znaup2.update) {
		goto L20;
	}

	// Back from computing user specified shifts
	if (data->s_znaup2.ushift) {
		goto L50;
	}

	// Back from computing residual norm
	// at the end of the current iteration
	if (data->s_znaup2.cnorm) {
		goto L100;
	}

	// Compute the first NEV steps of the Arnoldi factorization
	*info = znaitr_(data, ido, bmat, n, 0, *nev, mode, resid, &data->s_znaup2.rnorm, v, ldv, h, ldh, ipntr, workd);

	if (*ido != 99) {
		goto L9000;
	}

	if (*info > 0) {
		*np = *info;
		*mxiter = data->s_znaup2.iter;
		*info = -9999;
		goto L1200;
	}

	//           M A I N  ARNOLDI  I T E R A T I O N  L O O P
	//           Each iteration implicitly restarts the Arnoldi
	//           factorization in place.
L1000:
	++data->s_znaup2.iter;

	// Compute NP additional steps of the Arnoldi factorization.
	// Adjust NP since NEV might have been updated by last call
	// to the shift application routine znapps .
	*np = data->s_znaup2.kplusp - *nev;

	// Compute NP additional steps of the Arnoldi factorization.
	*ido = 0;
L20:
	data->s_znaup2.update = true;

	*info = znaitr_(data, ido, bmat, n, *nev, *np, mode, resid, &data->s_znaup2.rnorm, v, ldv, h, ldh, ipntr, workd);

	if (*ido != 99) {
		goto L9000;
	}

	if (*info > 0) {
		*np = *info;
		*mxiter = data->s_znaup2.iter;
		*info = -9999;
		goto L1200;
	}
	data->s_znaup2.update = false;

	// Compute the eigenvalues and corresponding error bounds
	// of the current upper Hessenberg matrix.
	{
		int ierr = zneigh_(data->s_znaup2.rnorm, data->s_znaup2.kplusp, h, ldh, ritz, bounds, q, ldq, workl, rwork);
		if (ierr != 0) {
			*info = -8;
			goto L1200;
		}
	}
	// Select the wanted Ritz values and their bounds
	// to be used in the convergence test.
	// The wanted part of the spectrum and corresponding
	// error bounds are in the last NEV loc. of RITZ,
	// and BOUNDS respectively.
	*nev = data->s_znaup2.nev0;
	*np = data->s_znaup2.np0;

	// Make a copy of Ritz values and the corresponding
	// Ritz estimates obtained from zneigh .
	RNP::TBLAS::Copy(data->s_znaup2.kplusp, ritz, 1, &workl[data->s_znaup2.kplusp * data->s_znaup2.kplusp], 1);
	RNP::TBLAS::Copy(data->s_znaup2.kplusp, bounds, 1, &workl[data->s_znaup2.kplusp * data->s_znaup2.kplusp + data->s_znaup2.kplusp], 1);

	// Select the wanted Ritz values and their bounds
	// to be used in the convergence test.
	// The wanted part of the spectrum and corresponding
	// bounds are in the last NEV loc. of RITZ
	// BOUNDS respectively.
	zngets_(ishift, sorter, *nev, *np, ritz, bounds, data->sort_data);

	// Convergence test: currently we use the following criteria.
	// The relative accuracy of a Ritz value is considered
	// acceptable if:
	//
	// error_bounds(i) .le. tol*max(eps23, magnitude_of_ritz(i)).
	//
	data->s_znaup2.nconv = 0;

	{size_t ilst = *nev;
	for (size_t i = 0; i < ilst; ++i) {
		/* Computing MAX */
		double rtemp = max(data->s_znaup2.eps23,std::abs(ritz[*np + i]));
		if (std::abs(bounds[*np + i]) <= tol * rtemp) {
			++data->s_znaup2.nconv;
		}
	}}

	// Count the number of unwanted Ritz values that have zero
	// Ritz estimates. If any Ritz estimates are equal to zero
	// then a leading block of H of order equal to at least
	// the number of Ritz values with zero Ritz estimates has
	// split off. None of these Ritz values may be removed by
	// shifting. Decrease NP the number of shifts to apply. If
	// no shifts may be applied, then prepare to exit
	{size_t nptemp = *np;
	for (size_t j = 0; j < nptemp; ++j) {
		if (bounds[j] == 0.) {
			--(*np);
			++(*nev);
		}
	}}

	if (data->s_znaup2.nconv >= data->s_znaup2.nev0 || data->s_znaup2.iter > *mxiter || *np == 0) {
		// Prepare to exit. Put the converged Ritz values
		// and corresponding bounds in RITZ(1:NCONV) and
		// BOUNDS(1:NCONV) respectively. Then sort. Be
		// careful when NCONV > NP

		//  Use h( 3,1 ) as storage to communicate
		//  rnorm to zneupd  if needed
		h[ldh + 3] = data->s_znaup2.rnorm;

		// Sort Ritz values so that converged Ritz
		// values appear within the first NEV locations
		// of ritz and bounds, and the most desired one
		// appears at the front.
		zsortc_(sorter, data->s_znaup2.kplusp, ritz, bounds, true, data->sort_data);

		// Scale the Ritz estimate of each Ritz value
		// by 1 / max(eps23, magnitude of the Ritz value).
		for (size_t j = 0; j < data->s_znaup2.nev0; ++j) {
			/* Computing MAX */
			double rtemp = max(data->s_znaup2.eps23,std::abs(ritz[j]));
			bounds[j] /= rtemp;
		}

		// Sort the Ritz values according to the scaled Ritz
		// estimates.  This will push all the converged ones
		// towards the front of ritz, bounds (in the case
		// when NCONV < NEV.)
		zsortc_(&RNP::IRA::LargestMagnitude, data->s_znaup2.nev0, bounds, ritz, false, NULL);

		// Scale the Ritz estimate back to its original
		// value.
		for (size_t j = 0; j < data->s_znaup2.nev0; ++j) {
			double rtemp = max(data->s_znaup2.eps23,std::abs(ritz[j]));
			bounds[j] *= rtemp;
		}

		// Sort the converged Ritz values again so that
		// the "threshold" value appears at the front of
		// ritz and bound.
		zsortc_(sorter, data->s_znaup2.nconv, ritz, bounds, false, data->sort_data);

		// Max iterations have been exceeded.
		if (data->s_znaup2.iter > *mxiter && data->s_znaup2.nconv < data->s_znaup2.nev0) {
			*info = 1;
		}

		// No shifts to apply.
		if (*np == 0 && data->s_znaup2.nconv < data->s_znaup2.nev0) {
			*info = 2;
		}

		*np = data->s_znaup2.nconv;
		goto L1100;
	}else if (data->s_znaup2.nconv < data->s_znaup2.nev0 && ishift == 1){
		// Do not have all the requested eigenvalues yet.
		// To prevent possible stagnation, adjust the size
		// of NEV.
		data->s_znaup2.nevbef = *nev;
		/* Computing MIN */
		*nev += min(data->s_znaup2.nconv, *np / 2);
		if (*nev == 1 && data->s_znaup2.kplusp >= 6) {
			*nev = data->s_znaup2.kplusp / 2;
		} else if (*nev == 1 && data->s_znaup2.kplusp > 3) {
			*nev = 2;
		}
		*np = data->s_znaup2.kplusp - *nev;

		// If the size of NEV was just increased
		// resort the eigenvalues.
		if (data->s_znaup2.nevbef < *nev) {
			zngets_(ishift, sorter, *nev, *np, ritz, bounds, data->sort_data);
		}
	}

	if (ishift == 0) {
		// User specified shifts: pop back out to get the shifts
		// and return them in the first 2*NP locations of WORKL.

		data->s_znaup2.ushift = true;
		*ido = 3;
		goto L9000;
	}
L50:
	data->s_znaup2.ushift = false;

	if (ishift != 1) {
		// Move the NP shifts from WORKL to
		// RITZ, to free up WORKL
		// for non-exact shift case.
		RNP::TBLAS::Copy(*np, workl, 1, ritz, 1);
	}

	// Apply the NP implicit shifts by QR bulge chasing.
	// Each shift is applied to the whole upper Hessenberg
	// matrix H.
	// The first 2*N locations of WORKD are used as workspace.
	znapps_(n, nev, *np, ritz, v, ldv, h, ldh, resid, q, ldq, workl, workd);

	// Compute the B-norm of the updated residual.
	// Keep B*RESID in WORKD(1:N) to be used in
	// the first step of the next call to znaitr .
	data->s_znaup2.cnorm = true;
	if (bmat == 'G') {
		RNP::TBLAS::Copy(n, resid, 1, &workd[n], 1);
		ipntr[1] = 0;
		ipntr[0] = n + ipntr[1];
		*ido = 2;

		// Exit in order to compute B*RESID
		goto L9000;
	} else if (bmat == 'I') {
		RNP::TBLAS::Copy(n, resid, 1, workd, 1);
	}

L100:
	// Back from reverse communication;
	// WORKD(1:N) := B*RESID
	if (bmat == 'G') {
		std::complex<double> cmpnorm = RNP::TBLAS::ConjugateDot(n, resid, 1, workd, 1);
		data->s_znaup2.rnorm = sqrt(std::abs(cmpnorm));
	} else if (bmat == 'I') {
		data->s_znaup2.rnorm = RNP::TBLAS::Norm2(n, resid, 1);
	}
	data->s_znaup2.cnorm = false;

	goto L1000;

	//  E N D     O F     M A I N     I T E R A T I O N     L O O P

L1100:
	*mxiter = data->s_znaup2.iter;
	*nev = data->s_znaup2.nconv;

L1200:
	*ido = 99;

	// Error Exit
L9000:
	return;
}

// \BeginDoc

// \Name: znaupd

// \Description:
//  Reverse communication interface for the Implicitly Restarted Arnoldi
//  iteration. This is intended to be used to find a few eigenpairs of a
//  complex linear operator OP with respect to a semi-inner product defined
//  by a hermitian positive semi-definite real matrix B. B may be the identity
//  matrix.  NOTE: if both OP and B are real, then dsaupd or dnaupd should
//  be used.


//  The computed approximate eigenvalues are called Ritz values and
//  the corresponding approximate eigenvectors are called Ritz vectors.

//  znaupd is usually called iteratively to solve one of the
//  following problems:

//  Mode 1:  A*x = lambda*x.
//           ===> OP = A  and  B = I.

//  Mode 2:  A*x = lambda*M*x, M hermitian positive definite
//           ===> OP = inv[M]*A  and  B = M.
//           ===> (If M can be factored see remark 3 below)

//  Mode 3:  A*x = lambda*M*x, M hermitian semi-definite
//           ===> OP =  inv[A - sigma*M]*M   and  B = M.
//           ===> shift-and-invert mode
//           If OP*x = amu*x, then lambda = sigma + 1/amu.


//  NOTE: The action of w <- inv[A - sigma*M]*v or w <- inv[M]*v
//        should be accomplished either by a direct method
//        using a sparse matrix factorization and solving

//           [A - sigma*M]*w = v  or M*w = v,

//        or through an iterative method for solving these
//        systems.  If an iterative method is used, the
//        convergence test must be more stringent than
//        the accuracy requirements for the eigenvalue
//        approximations.

// \Usage:
//  call znaupd
//     ( IDO, BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM,
//       IPNTR, WORKD, WORKL, LWORKL, RWORK, INFO )

// \Arguments
//  IDO     int.  (INPUT/OUTPUT)
//          Reverse communication flag.  IDO must be zero on the first
//          call to znaupd.  IDO will be set internally to
//          indicate the type of operation to be performed.  Control is
//          then given back to the calling routine which has the
//          responsibility to carry out the requested operation and call
//          znaupd with the result.  The operand is given in
//          WORKD(IPNTR(1)), the result must be put in WORKD(IPNTR(2)).
//          -------------------------------------------------------------
//          IDO =  0: first call to the reverse communication interface
//          IDO = -1: compute  Y = OP * X  where
//                    IPNTR(1) is the pointer into WORKD for X,
//                    IPNTR(2) is the pointer into WORKD for Y.
//                    This is for the initialization phase to force the
//                    starting vector into the range of OP.
//          IDO =  1: compute  Y = OP * X  where
//                    IPNTR(1) is the pointer into WORKD for X,
//                    IPNTR(2) is the pointer into WORKD for Y.
//                    In mode 3, the vector B * X is already
//                    available in WORKD(ipntr(3)).  It does not
//                    need to be recomputed in forming OP * X.
//          IDO =  2: compute  Y = M * X  where
//                    IPNTR(1) is the pointer into WORKD for X,
//                    IPNTR(2) is the pointer into WORKD for Y.
//          IDO =  3: compute and return the shifts in the first
//                    NP locations of WORKL.
//          IDO = 99: done
//          -------------------------------------------------------------
//          After the initialization phase, when the routine is used in
//          the "shift-and-invert" mode, the vector M * X is already
//          available and does not need to be recomputed in forming OP*X.

//  BMAT    Character*1.  (INPUT)
//          BMAT specifies the type of the matrix B that defines the
//          semi-inner product for the operator OP.
//          BMAT = 'I' -> standard eigenvalue problem A*x = lambda*x
//          BMAT = 'G' -> generalized eigenvalue problem A*x = lambda*M*x

//  N       int.  (INPUT)
//          Dimension of the eigenproblem.

//  WHICH   Character*2.  (INPUT)
//          'LM' -> want the NEV eigenvalues of largest magnitude.
//          'SM' -> want the NEV eigenvalues of smallest magnitude.
//          'LR' -> want the NEV eigenvalues of largest real part.
//          'SR' -> want the NEV eigenvalues of smallest real part.
//          'LI' -> want the NEV eigenvalues of largest imaginary part.
//          'SI' -> want the NEV eigenvalues of smallest imaginary part.

//  NEV     int.  (INPUT)
//          Number of eigenvalues of OP to be computed. 0 < NEV < N-1.

//  TOL     Double precision  scalar.  (INPUT)
//          Stopping criteria: the relative accuracy of the Ritz value
//          is considered acceptable if BOUNDS(I) .LE. TOL*ABS(RITZ(I))
//          where ABS(RITZ(I)) is the magnitude when RITZ(I) is complex.
//          DEFAULT = dlamch('EPS')  (machine precision as computed
//                    by the LAPACK auxiliary subroutine dlamch).

//  RESID   Complex*16 array of length N.  (INPUT/OUTPUT)
//          On INPUT:
//          If INFO .EQ. 0, a random initial residual vector is used.
//          If INFO .NE. 0, RESID contains the initial residual vector,
//                          possibly from a previous run.
//          On OUTPUT:
//          RESID contains the final residual vector.

//  NCV     int.  (INPUT)
//          Number of columns of the matrix V. NCV must satisfy the two
//          inequalities 1 <= NCV-NEV and NCV <= N.
//          This will indicate how many Arnoldi vectors are generated
//          at each iteration.  After the startup phase in which NEV
//          Arnoldi vectors are generated, the algorithm generates
//          approximately NCV-NEV Arnoldi vectors at each subsequent update
//          iteration. Most of the cost in generating each Arnoldi vector is
//          in the matrix-vector operation OP*x. (See remark 4 below.)

//  V       Complex*16 array N by NCV.  (OUTPUT)
//          Contains the final set of Arnoldi basis vectors.

//  LDV     int.  (INPUT)
//          Leading dimension of V exactly as declared in the calling program.

//  IPARAM  int array of length 11.  (INPUT/OUTPUT)
//          IPARAM(1) = ISHIFT: method for selecting the implicit shifts.
//          The shifts selected at each iteration are used to filter out
//          the components of the unwanted eigenvector.
//          -------------------------------------------------------------
//          ISHIFT = 0: the shifts are to be provided by the user via
//                      reverse communication.  The NCV eigenvalues of
//                      the Hessenberg matrix H are returned in the part
//                      of WORKL array corresponding to RITZ.
//          ISHIFT = 1: exact shifts with respect to the current
//                      Hessenberg matrix H.  This is equivalent to
//                      restarting the iteration from the beginning
//                      after updating the starting vector with a linear
//                      combination of Ritz vectors associated with the
//                      "wanted" eigenvalues.
//          ISHIFT = 2: other choice of internal shift to be defined.
//          -------------------------------------------------------------

//          IPARAM(2) = No longer referenced

//          IPARAM(3) = MXITER
//          On INPUT:  maximum number of Arnoldi update iterations allowed.
//          On OUTPUT: actual number of Arnoldi update iterations taken.

//          IPARAM(4) = NB: blocksize to be used in the recurrence.
//          The code currently works only for NB = 1.

//          IPARAM(5) = NCONV: number of "converged" Ritz values.
//          This represents the number of Ritz values that satisfy
//          the convergence criterion.

//          IPARAM(6) = IUPD
//          No longer referenced. Implicit restarting is ALWAYS used.

//          IPARAM(7) = MODE
//          On INPUT determines what type of eigenproblem is being solved.
//          Must be 1,2,3; See under \Description of znaupd for the
//          four modes available.

//          IPARAM(8) = NP
//          When ido = 3 and the user provides shifts through reverse
//          communication (IPARAM(1)=0), _naupd returns NP, the number
//          of shifts the user is to provide. 0 < NP < NCV-NEV.

//          IPARAM(9) = NUMOP, IPARAM(10) = NUMOPB, IPARAM(11) = NUMREO,
//          OUTPUT: NUMOP  = total number of OP*x operations,
//                  NUMOPB = total number of B*x operations if BMAT='G',
//                  NUMREO = total number of steps of re-orthogonalization.

//  IPNTR   int array of length 14.  (OUTPUT)
//          Pointer to mark the starting locations in the WORKD and WORKL
//          arrays for matrices/vectors used by the Arnoldi iteration.
//          -------------------------------------------------------------
//          IPNTR(1): pointer to the current operand vector X in WORKD.
//          IPNTR(2): pointer to the current result vector Y in WORKD.
//          IPNTR(3): pointer to the vector B * X in WORKD when used in
//                    the shift-and-invert mode.
//          IPNTR(4): pointer to the next available location in WORKL
//                    that is untouched by the program.
//          IPNTR(5): pointer to the NCV by NCV upper Hessenberg
//                    matrix H in WORKL.
//          IPNTR(6): pointer to the  ritz value array  RITZ
//          IPNTR(7): pointer to the (projected) ritz vector array Q
//          IPNTR(8): pointer to the error BOUNDS array in WORKL.
//          IPNTR(14): pointer to the NP shifts in WORKL. See Remark 5 below.

//          Note: IPNTR(9:13) is only referenced by zneupd. See Remark 2 below.

//          IPNTR(9): pointer to the NCV RITZ values of the
//                    original system.
//          IPNTR(10): Not Used
//          IPNTR(11): pointer to the NCV corresponding error bounds.
//          IPNTR(12): pointer to the NCV by NCV upper triangular
//                     Schur matrix for H.
//          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
//                     of the upper Hessenberg matrix H. Only referenced by
//                     zneupd if RVEC = .TRUE. See Remark 2 below.

//          -------------------------------------------------------------

//  WORKD   Complex*16 work array of length 3*N.  (REVERSE COMMUNICATION)
//          Distributed array to be used in the basic Arnoldi iteration
//          for reverse communication.  The user should not use WORKD
//          as temporary workspace during the iteration !!!!!!!!!!
//          See Data Distribution Note below.

//  WORKL   Complex*16 work array of length LWORKL.  (OUTPUT/WORKSPACE)
//          Private (replicated) array on each PE or array allocated on
//          the front end.  See Data Distribution Note below.

//  LWORKL  int.  (INPUT)
//          LWORKL must be at least 3*NCV**2 + 5*NCV.

//  RWORK   Double precision  work array of length NCV (WORKSPACE)
//          Private (replicated) array on each PE or array allocated on
//          the front end.


//  INFO    int.  (INPUT/OUTPUT)
//          If INFO .EQ. 0, a randomly initial residual vector is used.
//          If INFO .NE. 0, RESID contains the initial residual vector,
//                          possibly from a previous run.
//          Error flag on output.
//          =  0: Normal exit.
//          =  1: Maximum number of iterations taken.
//                All possible eigenvalues of OP has been found. IPARAM(5)
//                returns the number of wanted converged Ritz values.
//          =  2: No longer an informational error. Deprecated starting
//                with release 2 of ARPACK.
//          =  3: No shifts could be applied during a cycle of the
//                Implicitly restarted Arnoldi iteration. One possibility
//                is to increase the size of NCV relative to NEV.
//                See remark 4 below.
//          = -1: N must be positive.
//          = -2: NEV must be positive.
//          = -3: NCV-NEV >= 1 and less than or equal to N.
//          = -4: The maximum number of Arnoldi update iteration
//                must be greater than zero.
//          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
//          = -6: BMAT must be one of 'I' or 'G'.
//          = -7: Length of private work array is not sufficient.
//          = -8: Error return from LAPACK eigenvalue calculation;
//          = -9: Starting vector is zero.
//          = -10: IPARAM(7) must be 1,2,3.
//          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
//          = -12: IPARAM(1) must be equal to 0 or 1.
//          = -9999: Could not build an Arnoldi factorization.
//                   User input error highly likely.  Please
//                   check actual array dimensions and layout.
//                   IPARAM(5) returns the size of the current Arnoldi
//                   factorization.

// \Remarks
//  1. The computed Ritz values are approximate eigenvalues of OP. The
//     selection of WHICH should be made with this in mind when using
//     Mode = 3.  When operating in Mode = 3 setting WHICH = 'LM' will
//     compute the NEV eigenvalues of the original problem that are
//     closest to the shift SIGMA . After convergence, approximate eigenvalues
//     of the original problem may be obtained with the ARPACK subroutine zneupd.

//  2. If a basis for the invariant subspace corresponding to the converged Ritz
//     values is needed, the user must call zneupd immediately following
//     completion of znaupd. This is new starting with release 2 of ARPACK.

//  3. If M can be factored into a Cholesky factorization M = LL`
//     then Mode = 2 should not be selected.  Instead one should use
//     Mode = 1 with  OP = inv(L)*A*inv(L`).  Appropriate triangular
//     linear systems should be solved with L and L` rather
//     than computing inverses.  After convergence, an approximate
//     eigenvector z of the original problem is recovered by solving
//     L`z = x  where x is a Ritz vector of OP.

//  4. At present there is no a-priori analysis to guide the selection
//     of NCV relative to NEV.  The only formal requirement is that NCV > NEV + 1.
//     However, it is recommended that NCV .ge. 2*NEV.  If many problems of
//     the same type are to be solved, one should experiment with increasing
//     NCV while keeping NEV fixed for a given test problem.  This will
//     usually decrease the required number of OP*x operations but it
//     also increases the work and storage required to maintain the orthogonal
//     basis vectors.  The optimal "cross-over" with respect to CPU time
//     is problem dependent and must be determined empirically.
//     See Chapter 8 of Reference 2 for further information.

//  5. When IPARAM(1) = 0, and IDO = 3, the user needs to provide the
//     NP = IPARAM(8) complex shifts in locations
//     WORKL(IPNTR(14)), WORKL(IPNTR(14)+1), ... , WORKL(IPNTR(14)+NP).
//     Eigenvalues of the current upper Hessenberg matrix are located in
//     WORKL(IPNTR(6)) through WORKL(IPNTR(6)+NCV-1). They are ordered
//     according to the order defined by WHICH.  The associated Ritz estimates
//     are located in WORKL(IPNTR(8)), WORKL(IPNTR(8)+1), ... ,
//     WORKL(IPNTR(8)+NCV-1).

// -----------------------------------------------------------------------

// \Data Distribution Note:

//  Fortran-D syntax:
//  ================
//  Complex*16 resid(n), v(ldv,ncv), workd(3*n), workl(lworkl)
//  decompose  d1(n), d2(n,ncv)
//  align      resid(i) with d1(i)
//  align      v(i,j)   with d2(i,j)
//  align      workd(i) with d1(i)     range (1:n)
//  align      workd(i) with d1(i-n)   range (n+1:2*n)
//  align      workd(i) with d1(i-2*n) range (2*n+1:3*n)
//  distribute d1(block), d2(block,:)
//  replicated workl(lworkl)

//  Cray MPP syntax:
//  ===============
//  Complex*16 resid(n), v(ldv,ncv), workd(n,3), workl(lworkl)
//  shared     resid(block), v(block,:), workd(block,:)
//  replicated workl(lworkl)

//  CM2/CM5 syntax:
//  ==============

// -----------------------------------------------------------------------

//     include   'ex-nonsym.doc'

// -----------------------------------------------------------------------

// \BeginLib

// \Local variables:
//     xxxxxx  Complex*16

// \References:
//  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
//     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
//     pp 357-385.
//  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
//     Restarted Arnoldi Iteration", Rice University Technical Report
//     TR95-13, Department of Computational and Applied Mathematics.
//  3. B.N. Parlett & Y. Saad, "_Complex_ Shift and Invert Strategies for
//     Double precision Matrices", Linear Algebra and its Applications, vol 88/89,
//     pp 575-595, (1987).

// \Routines called:
//     znaup2  ARPACK routine that implements the Implicitly Restarted
//             Arnoldi Iteration.

// \Author
//     Danny Sorensen               Phuong Vu
//     Richard Lehoucq              CRPC / Rice University
//     Dept. of Computational &     Houston, Texas
//     Applied Mathematics
//     Rice University
//     Houston, Texas

// \SCCS Information: @(#)
// FILE: naupd.F   SID: 2.9   DATE OF SID: 07/21/02   RELEASE: 2

// \Remarks

// \EndLib

// -----------------------------------------------------------------------

void znaupd_(arz_data *data, int *ido, char bmat, size_t n, RNP::IRA::SortFunc sorter,
	size_t nev, double tol, std::complex<double> *resid, size_t
	ncv, std::complex<double> *v, size_t ldv, int *iparam, int *ipntr,
	std::complex<double> *workd, std::complex<double> *workl, size_t lworkl,
	double *rwork, int *info)
{
	if (*ido == 0) {
		// Error checking
		int ierr = 0;
		data->s_znaupd.ishift = iparam[0];
		// levec  = iparam(1)
		data->s_znaupd.mxiter = iparam[2];
		// nb     = iparam(3)
		data->s_znaupd.nb = 1;

		// Revision 2 performs only implicit restart.

		data->s_znaupd.iupd = 1;
		data->s_znaupd.mode = iparam[6];

		if (n <= 0) {
			ierr = -1;
		} else if (nev <= 0) {
			ierr = -2;
		} else if (ncv <= nev || ncv > n) {
			ierr = -3;
		} else if (data->s_znaupd.mxiter <= 0) {
			ierr = -4;
		} else if (bmat != 'I' && bmat != 'G') {
			ierr = -6;
		} else if (lworkl < ncv * ncv * 3 + ncv * 5) {
			ierr = -7;
		} else if (data->s_znaupd.mode < 1 || data->s_znaupd.mode > 3) {
			ierr = -10;
		} else if (data->s_znaupd.mode == 1 && bmat == 'G') {
			ierr = -11;
		}

		// Error Exit
		if (ierr != 0) {
			*info = ierr;
			*ido = 99;
			return;
		}

		// Set default parameters
		if (data->s_znaupd.nb <= 0) {
			data->s_znaupd.nb = 1;
		}
		if (tol <= 0.) {
			tol = std::numeric_limits<double>::epsilon();
		}
		if (data->s_znaupd.ishift != 0 && data->s_znaupd.ishift != 1 && data->s_znaupd.ishift != 2) {
			data->s_znaupd.ishift = 1;
		}

		// NP is the number of additional steps to
		// extend the length NEV Lanczos factorization.
		// NEV0 is the local variable designating the
		// size of the invariant subspace desired.
		data->s_znaupd.np = ncv - nev;
		data->s_znaupd.nev0 = nev;

		// Zero out internal workspace
		{
			size_t wksize = ncv * ncv * 3 + ncv * 5;
			for (size_t j = 0; j < wksize; ++j) {
				workl[j] = 0;
			}
		}

		// Pointer into WORKL for address of H, RITZ, BOUNDS, Q
		// etc... and the remaining workspace.
		// Also update pointer to be used on output.
		// Memory is laid out as follows:
		// workl(1:ncv*ncv) := generated Hessenberg matrix
		// workl(ncv*ncv+1:ncv*ncv+ncv) := the ritz values
		// workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv)   := error bounds
		// workl(ncv*ncv+2*ncv+1:2*ncv*ncv+2*ncv) := rotation matrix Q
		// workl(2*ncv*ncv+2*ncv+1:3*ncv*ncv+5*ncv) := workspace
		// The final workspace is needed by subroutine zneigh called
		// by znaup2. Subroutine zneigh calls LAPACK routines for
		// calculating eigenvalues and the last row of the eigenvector
		// matrix.

		data->s_znaupd.ldh = ncv;
		data->s_znaupd.ldq = ncv;
		data->s_znaupd.ih = 0;
		data->s_znaupd.ritz = data->s_znaupd.ih + data->s_znaupd.ldh * ncv;
		data->s_znaupd.bounds = data->s_znaupd.ritz + ncv;
		data->s_znaupd.iq = data->s_znaupd.bounds + ncv;
		data->s_znaupd.iw = data->s_znaupd.iq + data->s_znaupd.ldq * ncv;
		data->s_znaupd.next = data->s_znaupd.iw + ncv*ncv + ncv * 3;

		ipntr[3] = data->s_znaupd.next;
		ipntr[4] = data->s_znaupd.ih;
		ipntr[5] = data->s_znaupd.ritz;
		ipntr[6] = data->s_znaupd.iq;
		ipntr[7] = data->s_znaupd.bounds;
		ipntr[13] = data->s_znaupd.iw;
	}

	// Carry out the Implicitly restarted Arnoldi Iteration.
	znaup2_(data, ido, bmat, n, sorter, &data->s_znaupd.nev0, &data->s_znaupd.np, tol, resid, data->s_znaupd.mode, data->s_znaupd.ishift, &data->s_znaupd.mxiter, v, ldv, &workl[data->s_znaupd.ih], data->s_znaupd.ldh, &workl[data->s_znaupd.ritz], &workl[data->s_znaupd.bounds], &workl[data->s_znaupd.iq], data->s_znaupd.ldq, &workl[data->s_znaupd.iw], ipntr, workd, rwork, info);

	// ido .ne. 99 implies use of reverse communication
	// to compute operations involving OP.

	if (*ido == 3) {
		iparam[7] = data->s_znaupd.np;
	}
	if (*ido != 99) {
		return;
	}

	iparam[2] = data->s_znaupd.mxiter;
	iparam[4] = data->s_znaupd.np;
	//iparam[8] = timing_1.nopx;
	//iparam[9] = timing_1.nbx;
	//iparam[10] = timing_1.nrorth;

	// Exit if there was an informational
	// error within znaup2.
	if (*info < 0) {
		return;
	}
	if (*info == 2) {
		*info = 3;
	}
}

// \BeginDoc

// \Name: zneupd

// \Description:
//  This subroutine returns the converged approximations to eigenvalues
//  of A*z = lambda*B*z and (optionally):

//      (1) The corresponding approximate eigenvectors;

//      (2) An orthonormal basis for the associated approximate
//          invariant subspace;

//      (3) Both.

//  There is negligible additional cost to obtain eigenvectors.  An orthonormal
//  basis is always computed.  There is an additional storage cost of n*nev
//  if both are requested (in this case a separate array Z must be supplied).

//  The approximate eigenvalues and eigenvectors of  A*z = lambda*B*z
//  are derived from approximate eigenvalues and eigenvectors of
//  of the linear operator OP prescribed by the MODE selection in the
//  call to ZNAUPD.  ZNAUPD must be called before this routine is called.
//  These approximate eigenvalues and vectors are commonly called Ritz
//  values and Ritz vectors respectively.  They are referred to as such
//  in the comments that follow.   The computed orthonormal basis for the
//  invariant subspace corresponding to these Ritz values is referred to as a
//  Schur basis.

//  The definition of OP as well as other terms and the relation of computed
//  Ritz values and vectors of OP with respect to the given problem
//  A*z = lambda*B*z may be found in the header of ZNAUPD.  For a brief
//  description, see definitions of IPARAM(7), MODE and WHICH in the
//  documentation of ZNAUPD.

// \Usage:
//  call zneupd
//     ( RVEC, HOWMNY, SELECT, D, Z, LDZ, SIGMA, WORKEV, BMAT,
//       N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR, WORKD,
//       WORKL, LWORKL, RWORK, INFO )

// \Arguments:
//  RVEC    LOGICAL  (INPUT)
//          Specifies whether a basis for the invariant subspace corresponding
//          to the converged Ritz value approximations for the eigenproblem
//          A*z = lambda*B*z is computed.

//             RVEC = .FALSE.     Compute Ritz values only.

//             RVEC = .TRUE.      Compute Ritz vectors or Schur vectors.
//                                See Remarks below.

//  HOWMNY  Character*1  (INPUT)
//          Specifies the form of the basis for the invariant subspace
//          corresponding to the converged Ritz values that is to be computed.

//          = 'A': Compute NEV Ritz vectors;
//          = 'P': Compute NEV Schur vectors;
//          = 'S': compute some of the Ritz vectors, specified
//                 by the logical array SELECT.

//  SELECT  Logical array of dimension NCV.  (INPUT)
//          If HOWMNY = 'S', SELECT specifies the Ritz vectors to be
//          computed. To select the  Ritz vector corresponding to a
//          Ritz value D(j), SELECT(j) must be set to .TRUE..
//          If HOWMNY = 'A' or 'P', SELECT need not be initialized
//          but it is used as internal workspace.

//  D       Complex*16 array of dimension NEV+1.  (OUTPUT)
//          On exit, D contains the  Ritz  approximations
//          to the eigenvalues lambda for A*z = lambda*B*z.

//  Z       Complex*16 N by NEV array    (OUTPUT)
//          On exit, if RVEC = .TRUE. and HOWMNY = 'A', then the columns of
//          Z represents approximate eigenvectors (Ritz vectors) corresponding
//          to the NCONV=IPARAM(5) Ritz values for eigensystem
//          A*z = lambda*B*z.

//          If RVEC = .FALSE. or HOWMNY = 'P', then Z is NOT REFERENCED.

//          NOTE: If if RVEC = .TRUE. and a Schur basis is not required,
//          the array Z may be set equal to first NEV+1 columns of the Arnoldi
//          basis array V computed by ZNAUPD.  In this case the Arnoldi basis
//          will be destroyed and overwritten with the eigenvector basis.

//  LDZ     int.  (INPUT)
//          The leading dimension of the array Z.  If Ritz vectors are
//          desired, then  LDZ .ge.  max( 1, N ) is required.
//          In any case,  LDZ .ge. 1 is required.

//  SIGMA   Complex*16  (INPUT)
//          If IPARAM(7) = 3 then SIGMA represents the shift.
//          Not referenced if IPARAM(7) = 1 or 2.

//  WORKEV  Complex*16 work array of dimension 2*NCV.  (WORKSPACE)

//  **** The remaining arguments MUST be the same as for the   ****
//  **** call to ZNAUPD that was just completed.               ****

//  NOTE: The remaining arguments

//           BMAT, N, WHICH, NEV, TOL, RESID, NCV, V, LDV, IPARAM, IPNTR,
//           WORKD, WORKL, LWORKL, RWORK, INFO

//         must be passed directly to ZNEUPD following the last call
//         to ZNAUPD.  These arguments MUST NOT BE MODIFIED between
//         the the last call to ZNAUPD and the call to ZNEUPD.

//  Three of these parameters (V, WORKL and INFO) are also output parameters:

//  V       Complex*16 N by NCV array.  (INPUT/OUTPUT)

//          Upon INPUT: the NCV columns of V contain the Arnoldi basis
//                      vectors for OP as constructed by ZNAUPD .

//          Upon OUTPUT: If RVEC = .TRUE. the first NCONV=IPARAM(5) columns
//                       contain approximate Schur vectors that span the
//                       desired invariant subspace.

//          NOTE: If the array Z has been set equal to first NEV+1 columns
//          of the array V and RVEC=.TRUE. and HOWMNY= 'A', then the
//          Arnoldi basis held by V has been overwritten by the desired
//          Ritz vectors.  If a separate array Z has been passed then
//          the first NCONV=IPARAM(5) columns of V will contain approximate
//          Schur vectors that span the desired invariant subspace.

//  WORKL   Double precision work array of length LWORKL.  (OUTPUT/WORKSPACE)
//          WORKL(1:ncv*ncv+2*ncv) contains information obtained in
//          znaupd.  They are not changed by zneupd.
//          WORKL(ncv*ncv+2*ncv+1:3*ncv*ncv+4*ncv) holds the
//          untransformed Ritz values, the untransformed error estimates of
//          the Ritz values, the upper triangular matrix for H, and the
//          associated matrix representation of the invariant subspace for H.

//          Note: IPNTR(9:13) contains the pointer into WORKL for addresses
//          of the above information computed by zneupd.
//          -------------------------------------------------------------
//          IPNTR(9):  pointer to the NCV RITZ values of the
//                     original system.
//          IPNTR(10): Not used
//          IPNTR(11): pointer to the NCV corresponding error estimates.
//          IPNTR(12): pointer to the NCV by NCV upper triangular
//                     Schur matrix for H.
//          IPNTR(13): pointer to the NCV by NCV matrix of eigenvectors
//                     of the upper Hessenberg matrix H. Only referenced by
//                     zneupd if RVEC = .TRUE. See Remark 2 below.
//          -------------------------------------------------------------

//  INFO    int.  (OUTPUT)
//          Error flag on output.
//          =  0: Normal exit.

//          =  1: The Schur form computed by LAPACK routine csheqr
//                could not be reordered by LAPACK routine ztrsen.
//                Re-enter subroutine zneupd with IPARAM(5)=NCV and
//                increase the size of the array D to have
//                dimension at least dimension NCV and allocate at least NCV
//                columns for Z. NOTE: Not necessary if Z and V share
//                the same space. Please notify the authors if this error
//                occurs.

//          = -1: N must be positive.
//          = -2: NEV must be positive.
//          = -3: NCV-NEV >= 1 and less than or equal to N.
//          = -5: WHICH must be one of 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
//          = -6: BMAT must be one of 'I' or 'G'.
//          = -7: Length of private work WORKL array is not sufficient.
//          = -8: Error return from LAPACK eigenvalue calculation.
//                This should never happened.
//          = -9: Error return from calculation of eigenvectors.
//                Informational error from LAPACK routine ztrevc.
//          = -10: IPARAM(7) must be 1,2,3
//          = -11: IPARAM(7) = 1 and BMAT = 'G' are incompatible.
//          = -12: HOWMNY = 'S' not yet implemented
//          = -13: HOWMNY must be one of 'A' or 'P' if RVEC = .true.
//          = -14: ZNAUPD did not find any eigenvalues to sufficient
//                 accuracy.
//          = -15: ZNEUPD got a different count of the number of converged
//                 Ritz values than ZNAUPD got.  This indicates the user
//                 probably made an error in passing data from ZNAUPD to
//                 ZNEUPD or that the data was modified before entering
//                 ZNEUPD

// \BeginLib

// \References:
//  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
//     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
//     pp 357-385.
//  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly
//     Restarted Arnoldi Iteration", Rice University Technical Report
//     TR95-13, Department of Computational and Applied Mathematics.
//  3. B. Nour-Omid, B. N. Parlett, T. Ericsson and P. S. Jensen,
//     "How to Implement the Spectral Transformation", Math Comp.,
//     Vol. 48, No. 178, April, 1987 pp. 664-673.

// \Routines called:
//     ztrevc  LAPACK routine to compute the eigenvectors of a matrix
//             in upper triangular form.
//     ztrsen  LAPACK routine that re-orders the Schur form.

// \Remarks

//  1. Currently only HOWMNY = 'A' and 'P' are implemented.

//  2. Schur vectors are an orthogonal representation for the basis of
//     Ritz vectors. Thus, their numerical properties are often superior.
//     If RVEC = .true. then the relationship
//             A * V(:,1:IPARAM(5)) = V(:,1:IPARAM(5)) * T, and
//       transpose( V(:,1:IPARAM(5)) ) * V(:,1:IPARAM(5)) = I
//     are approximately satisfied.
//     Here T is the leading submatrix of order IPARAM(5) of the
//     upper triangular matrix stored workl(ipntr(12)).

// \Authors
//     Danny Sorensen               Phuong Vu
//     Richard Lehoucq              CRPC / Rice University
//     Chao Yang                    Houston, Texas
//     Dept. of Computational &
//     Applied Mathematics
//     Rice University
//     Houston, Texas

// \SCCS Information: @(#)
// FILE: neupd.F   SID: 2.8   DATE OF SID: 07/21/02   RELEASE: 2

// \EndLib

// -----------------------------------------------------------------------
void zneupd_(bool rvec, char howmny, bool *select,
	std::complex<double> *d, std::complex<double> *z, size_t ldz, const std::complex<double> *
	sigma, std::complex<double> *workev, char bmat, size_t n, RNP::IRA::SortFunc sorter, void *sort_data,
	size_t nev, double tol, std::complex<double> *resid, size_t ncv,
	std::complex<double> *v, size_t ldv, int *iparam, int *ipntr,
	std::complex<double> *workd, std::complex<double> *workl, size_t lworkl,
	double *rwork, int *info)
{
	using namespace std;

	/* Local variables */
	int mode;
	double eps23;
	int ierr;
	std::complex<double> temp;
	enum shift_type{
		REGULAR,
		SHIFTI
	} type__;
	bool reord;
	size_t nconv;
	double rtemp;
	std::complex<double> rnorm;
	size_t numcnv;
	extern int ztrevc_(char howmny, bool *select, size_t n, std::complex<double> *_t, size_t ldt, std::complex<double> *_vl, size_t ldvl, std::complex<double> *_vr, size_t ldvr, size_t mm, size_t *m, std::complex<double> *_work, double *rwork);
	extern void ztrexc_(int n, std::complex<double> *t, int ldt, std::complex<double> *q, int ldq, int ifst, int ilst);
	extern size_t zlahqr_(
		bool wantt, bool wantz,
		size_t n, size_t ilo, size_t ihi,
		std::complex<double> *h, int ldh,
		std::complex<double> *w,
		size_t iloz, size_t ihiz,
		std::complex<double> *z, size_t ldz);

	/* Function Body */
	mode = iparam[6];
	nconv = iparam[4];
	*info = 0;


	// Get machine dependent constant.
	eps23 = pow(std::numeric_limits<double>::epsilon(), (double)2/(double)3);

	// Quick return
	// Check for incompatible input

	ierr = 0;

	if (nconv <= 0) {
		ierr = -14;
	} else if (n <= 0) {
		ierr = -1;
	} else if (nev <= 0) {
		ierr = -2;
	} else if (ncv <= nev || ncv > n) {
		ierr = -3;
	} else if (bmat != 'I' && bmat != 'G'){
		ierr = -6;
	} else if (lworkl < ncv * ncv * 3 + 4*ncv) {
		ierr = -7;
	} else if (howmny != 'A' && howmny != 'P' && howmny != 'S' && rvec) {
		ierr = -13;
	} else if (howmny == 'S') {
		ierr = -12;
	}

	if (mode == 1 || mode == 2) {
		type__ = REGULAR;
	} else if (mode == 3) {
		type__ = SHIFTI;
	} else {
		ierr = -10;
	}
	if (mode == 1 && bmat == 'G') {
		ierr = -11;
	}

	// Error Exit

	if (ierr != 0) {
		*info = ierr;
		return;
	}

	// Pointer into WORKL for address of H, RITZ, WORKEV, Q
	// etc... and the remaining workspace.
	// Also update pointer to be used on output.
	// Memory is laid out as follows:
	// workl(1:ncv*ncv) := generated Hessenberg matrix
	// workl(ncv*ncv+1:ncv*ncv+ncv) := ritz values
	// workl(ncv*ncv+ncv+1:ncv*ncv+2*ncv) := error bounds

	// The following is used and set by ZNEUPD.
	// workl(ncv*ncv+2*ncv+1:ncv*ncv+3*ncv) := The untransformed
	//                                      Ritz values.
	// workl(ncv*ncv+3*ncv+1:ncv*ncv+4*ncv) := The untransformed
	//                                      error bounds of
	//                                      the Ritz values
	// workl(ncv*ncv+4*ncv+1:2*ncv*ncv+4*ncv) := Holds the upper
	//                                      triangular matrix
	//                                      for H.
	// workl(2*ncv*ncv+4*ncv+1: 3*ncv*ncv+4*ncv) := Holds the
	//                                      associated matrix
	//                                      representation of
	//                                      the invariant
	//                                      subspace for H.
	// GRAND total of NCV * ( 3 * NCV + 4 ) locations.

	const size_t ih = ipntr[4];
	const size_t ritz = ipntr[5];
	//const size_t iq = ipntr[6];
	const size_t bounds = ipntr[7];
	const size_t ldh = ncv;
	const size_t ldq = ncv;
	const size_t iheig = bounds + ldh;
	const size_t ihbds = iheig + ldh;
	const size_t iuptri = ihbds + ldh;
	const size_t invsub = iuptri + ldh * ncv;
	ipntr[8] = iheig;
	ipntr[10] = ihbds;
	ipntr[11] = iuptri;
	ipntr[12] = invsub;
	//const size_t wr = 0;
	//const size_t iwev = wr + ncv;

	// irz points to the Ritz values computed
	//     by _neigh before exiting _naup2.
	// ibd points to the Ritz estimates
	//     computed by _neigh before exiting
	//     _naup2.

	const size_t irz = ipntr[13] + ncv * ncv;
	const size_t ibd = irz + ncv;
	
	std::complex<double> *work_h = &workl[ih];
	std::complex<double> *work_ritz = &workl[ritz];
	std::complex<double> *work_bounds = &workl[bounds];
	std::complex<double> *work_heig = &workl[iheig];
	std::complex<double> *work_hbds = &workl[ihbds];
	std::complex<double> *work_uptri = &workl[iuptri];
	std::complex<double> *work_invsub = &workl[invsub];
	std::complex<double> *work_rz = &workl[irz];
	std::complex<double> *work_bd = &workl[ibd];
	
	// RNORM is B-norm of the RESID(1:N).
	rnorm = work_h[2];
	work_h[2] = 0;

	if (rvec) {
		reord = false;

		// Use the temporary bounds array to store indices
		// These will be used to mark the select array later
		for (size_t j = 0; j < ncv; ++j) {
			work_bounds[j] = (double)j;
			select[j] = false;
		}

		// Select the wanted Ritz values.
		// Sort the Ritz values so that the
		// wanted ones appear at the tailing
		// NEV positions of workl(irr) and
		// workl(iri).  Move the corresponding
		// error estimates in workl(ibd)
		// accordingly.
		zngets_(0, sorter, nev, ncv - nev, work_rz, work_bounds, sort_data);

		// Record indices of the converged wanted Ritz values
		// Mark the select array for possible reordering
		numcnv = 0;
		for (size_t j = 0; j < ncv; ++j) {
			/* Computing MAX */
			rtemp = std::max(eps23,std::abs(work_rz[ncv - j-1]));
			size_t jj = (int)work_bounds[ncv - j-1].real();
			if (numcnv < nconv && std::abs(work_bd[jj]) <= tol * rtemp) {
				select[jj] = true;
				++numcnv;
				if (jj >= nev) {
					reord = true;
				}
			}
		}

		// Check the count (numcnv) of converged Ritz values with
		// the number (nconv) reported by dnaupd.  If these two
		// are different then there has probably been an error
		// caused by incorrect passing of the dnaupd data.
		if (numcnv != nconv) {
			*info = -15;
			goto L9000;
		}

		// Call LAPACK routine zlahqr to compute the Schur form
		// of the upper Hessenberg matrix returned by ZNAUPD.
		// Make a copy of the upper Hessenberg matrix.
		// Initialize the Schur vector matrix Q to the identity.
		RNP::TBLAS::Copy(ldh * ncv, work_h, 1, work_uptri, 1);
		RNP::TBLAS::SetMatrix<'A'>(ncv, ncv, 0.0, 1.0, work_invsub, ldq);
		ierr = zlahqr_(true, true, ncv, 1, ncv, work_uptri, ldh, work_heig, 1, ncv, work_invsub, ldq);
		RNP::TBLAS::Copy(ncv, &work_invsub[ncv - 1], ldq, work_hbds, 1);

		if (ierr != 0) {
			*info = -8;
			goto L9000;
		}

		if (reord) {
			// Reorder the computed upper triangular matrix.
			//ztrsen_(&select[1], ncv, &workl[iuptri], ldh, &workl[invsub], ldq, &workl[iheig], &nconv, NULL, NULL);
			{
				nconv = 0;
				for (size_t k = 0; k < ncv; ++k) {
					if (select[k]) {
						++nconv;
					}
				}

				if(!(nconv == ncv || nconv == 0)){
					// Collect the selected eigenvalues at the top left corner of T.
					size_t ks = 0;
					for (size_t k = 0; k < ncv; ++k) {
						if (select[k]) {
							// Swap the K-th eigenvalue to position KS.
							if (k != ks) {
								ztrexc_(ncv, work_uptri, ldh, work_invsub, ldq, k, ks);
							}
							++ks;
						}
					}
				}
				// Copy reordered eigenvalues to W.
				for (size_t k = 0; k < ncv; ++k) {
					work_heig[k] = work_uptri[k+k*ldh];
				}
			}
			if (ierr == 1) {
				*info = 1;
				goto L9000;
			}
		}

		// Copy the last row of the Schur basis matrix
		// to workl(ihbds).  This vector will be used
		// to compute the Ritz estimates of converged
		// Ritz values.
		RNP::TBLAS::Copy(ncv, &work_invsub[ncv - 1], ldq, work_hbds, 1);

		// Place the computed eigenvalues of H into D
		// if a spectral transformation was not used.
		if (type__ == REGULAR) {
			RNP::TBLAS::Copy(nconv, work_heig, 1, d, 1);
		}

		// Compute the QR factorization of the matrix representing
		// the wanted invariant subspace located in the first NCONV
		// columns of workl(invsub,ldq).
		RNP::TLASupport::QRFactorization(ncv, nconv, work_invsub, ldq, workev, &workev[ncv]);

		// * Postmultiply V by Q using zunm2r.
		// * Copy the first NCONV columns of VQ into Z.
		// * Postmultiply Z by R.
		// The N by NCONV matrix Z is now a matrix representation
		// of the approximate invariant subspace associated with
		// the Ritz values in workl(iheig). The first NCONV
		// columns of V are now approximate Schur vectors
		// associated with the upper triangular matrix of order
		// NCONV in workl(iuptri).
		RNP::TLASupport::ApplyOrthognalMatrixFromElementaryReflectors<'R','N'>(n, ncv, nconv, work_invsub, ldq, workev, v, ldv, &workd[n]);
		RNP::TBLAS::CopyMatrix<'A'>(n, nconv, v, ldv, z, ldz);

		for (size_t j = 0; j < nconv; ++j) {
			// Perform both a column and row scaling if the
			// diagonal element of workl(invsub,ldq) is negative
			// I'm lazy and don't take advantage of the upper
			// triangular form of workl(iuptri,ldq).
			// Note that since Q is orthogonal, R is a diagonal
			// matrix consisting of plus or minus ones.
			if (work_invsub[j+j*ldq].real() < 0.) {
				RNP::TBLAS::Scale(nconv, -1.0, &work_uptri[j], ldq);
				RNP::TBLAS::Scale(nconv, -1.0, &work_uptri[j*ldq], 1);
			}
		}

		if (howmny == 'A') {
			// Compute the NCONV wanted eigenvectors of T
			// located in workl(iuptri,ldq).
			for (size_t j = 0; j < ncv; ++j) {
				if (j < nconv) {
					select[j] = true;
				} else {
					select[j] = false;
				}
			}

			{
				size_t dummy_ncv;
				ierr = ztrevc_('S', select, ncv, work_uptri, ldq, NULL, 1, work_invsub, ldq, ncv, &dummy_ncv, workev, rwork);
			}
			if (ierr != 0) {
				*info = -9;
				goto L9000;
			}

			// Scale the returning eigenvectors so that their
			// Euclidean norms are all one. LAPACK subroutine
			// ztrevc returns each eigenvector normalized so
			// that the element of largest magnitude has
			// magnitude 1.

			for (size_t j = 0; j < nconv; ++j) {
				rtemp = RNP::TBLAS::Norm2(ncv, &work_invsub[j*ldq], 1);
				RNP::TBLAS::Scale(ncv, 1. / rtemp, &work_invsub[j*ldq], 1);

				// Ritz estimates can be obtained by taking
				// the inner product of the last row of the
				// Schur basis of H with eigenvectors of T.
				// Note that the eigenvector matrix of T is
				// upper triangular, thus the length of the
				// inner product can be set to j.
				workev[j] = RNP::TBLAS::ConjugateDot(j+1, work_hbds, 1, &work_invsub[j*ldq], 1);
			}


			// Copy Ritz estimates into workl(ihbds)
			RNP::TBLAS::Copy(nconv, workev, 1, work_hbds, 1);

			// The eigenvector matrix Q of T is triangular.
			// Form Z*Q.
			RNP::TBLAS::MultTrM<'R','U','N','N'>(n, nconv, 1.0, work_invsub, ldq, z, ldz);
		}
	}else{
		// An approximate invariant subspace is not needed.
		// Place the Ritz values computed ZNAUPD into D.

		RNP::TBLAS::Copy(nconv, work_ritz, 1, d, 1);
		RNP::TBLAS::Copy(nconv, work_ritz, 1, work_heig, 1);
		RNP::TBLAS::Copy(nconv, work_bounds, 1, work_hbds, 1);
	}

	// Transform the Ritz values and possibly vectors
	// and corresponding error bounds of OP to those
	// of A*x = lambda*B*x.
	if (type__ == REGULAR) {
		if (rvec) {
			RNP::TBLAS::Scale(ncv, rnorm, work_hbds, 1);
		}
	}else{
		//   A spectral transformation was used.
		// * Determine the Ritz estimates of the
		//   Ritz values in the original system.
		if (rvec) {
			RNP::TBLAS::Scale(ncv, rnorm, work_hbds, 1);
		}
		for (size_t k = 0; k < ncv; ++k) {
			std::complex<double> temp = work_heig[k];
			work_hbds[k] /= temp;
			work_hbds[k] /= temp;
		}
	}

	// *  Transform the Ritz values back to the original system.
	//    For TYPE = 'SHIFTI' the transformation is
	//             lambda = 1/theta + sigma
	// NOTES:
	// *The Ritz vectors are not affected by the transformation.
	if (type__ == SHIFTI) {
		for (size_t k = 0; k < nconv; ++k) {
			d[k] = 1.0/work_heig[k] + *sigma;
		}
	}

	// Eigenvector Purification step. Formally perform
	// one of inverse subspace iteration. Only used
	// for MODE = 3. See reference 3.
	if (rvec && howmny == 'A' && type__ == SHIFTI) {
		// Purify the computed Ritz vectors by adding a
		// little bit of the residual vector:
		//                      T
		//          resid(:)*( e    s ) / theta
		//                      NCV
		// where H s = s theta.
		for (size_t j = 0; j < nconv; ++j) {
			if (work_heig[j] != 0.) {
				workev[j] = work_invsub[j*ldq + ncv - 1] / work_heig[j];
			}
		}
		// Perform a rank one update to Z and
		// purify all the Ritz vectors together.
		RNP::TBLAS::Rank1Update(n, nconv, 1.0, resid, 1, workev, 1, z, ldz);
	}

L9000:
	return;
}

bool RNP::IRA::AllocateWorkspace(size_t n, size_t n_wanted, size_t n_arnoldi, Workspace *workspace){
	if(NULL != workspace){
		
		workspace->resid = new std::complex<double>[
			n                           // resid
			+ 3*n                       // workd
			+ (3*n_arnoldi+5)*n_arnoldi // workl
			+ 2*n_arnoldi               // workev
		];
		workspace->workd  = workspace->resid + n;
		workspace->workl  = workspace->workd + 3*n;
		workspace->workev = workspace->workl + (3*n_arnoldi+5)*n_arnoldi;
		/*
		workspace->resid  = new std::complex<double>[n];
		workspace->workd  = new std::complex<double>[3*n];
		workspace->workl  = new std::complex<double>[(3*n_arnoldi+5)*n_arnoldi];
		workspace->workev = new std::complex<double>[2*n_arnoldi];
		*/
		workspace->rwork  = new double[n_arnoldi];
		workspace->bwork  = new bool[n_arnoldi];
		return true;
	}
	return false;
}
void RNP::IRA::FreeWorkspace(size_t n, size_t n_wanted, size_t n_arnoldi, Workspace *workspace){
	if(NULL != workspace){
		delete [] workspace->resid;
		/*
		delete [] workspace->workd;
		delete [] workspace->workl;
		delete [] workspace->workev;
		*/
		delete [] workspace->rwork;
		delete [] workspace->bwork;
	}
}

int RNP::IRA::Regular(
	size_t n, RNP::IRA::Areg Aop, RNP::IRA::Bfunc Bop,
	size_t n_wanted, size_t n_arnoldi,   // n_wanted+1 <= n_arnoldi <= n
	RNP::IRA::SortFunc sorter,
	std::complex<double> *w,             // eigenvalues (length n_wanted)
	std::complex<double> *v, size_t ldv, // eigenvectors (size ldv-by-n_wanted), ldv >= n
	RNP::IRA::Params *params,
	RNP::IRA::Workspace *workspace,
	void *Adata, void *Bdata, void *SortData
){
	RNP::IRA::Params *u_params = params;
	if(NULL == params){
		u_params = new Params();
	}
	RNP::IRA::Workspace *u_workspace = workspace;
	if(NULL == workspace){
		u_workspace = new RNP::IRA::Workspace();
		RNP::IRA::AllocateWorkspace(n, n_wanted, n_arnoldi, u_workspace);
	}
	
	if(u_params->tol <= 0){
		u_params->tol = std::numeric_limits<double>::epsilon();
	}
	
	int iparam[11], ipntr[14];
	int ido = 0;
	int info = 0;
	iparam[0] = 1; // ishift
	iparam[2] = u_params->max_iterations;
	iparam[6] = (NULL == Bop) ? 1 : 2;
	char bmat = (NULL == Bop) ? 'I' : 'G';
	arz_data static_data;
	static_data.sort_data = SortData;
	size_t lworkl = (3*n_arnoldi+5)*n_arnoldi;
	
	do{
		znaupd_(
			&static_data, &ido, bmat,
			n, sorter, n_wanted, u_params->tol,
			u_workspace->resid,
			n_arnoldi, v, ldv,
			iparam, ipntr,
			u_workspace->workd, u_workspace->workl, lworkl, u_workspace->rwork, &info);

		if(ido == -1 || ido == 1) {
			Aop(n, &u_workspace->workd[ipntr[0]], &u_workspace->workd[ipntr[1]], Adata);
			continue;
		}else if(ido == 2){
			Bop(n, &u_workspace->workd[ipntr[0]], &u_workspace->workd[ipntr[1]], Bdata);
			continue;
		}
		break;
	}while(1);
	if(info >= 0){
		int ierr;
		zneupd_(
			!(u_params->no_eigenvectors), 'A', u_workspace->bwork,
			w, v, ldv,
			NULL, u_workspace->workev,
			bmat, n, sorter, SortData, n_wanted,
			u_params->tol, u_workspace->resid, n_arnoldi, v, ldv,
			iparam, ipntr,
			u_workspace->workd, u_workspace->workl, lworkl,
			u_workspace->rwork, &ierr);
		if(ierr == 0){
			info = iparam[4]; // nconv
		}else if(ierr > 0){
			ierr = -100*ierr;
		}
	}
	
	if(NULL == params){
		delete u_params;
	}
	if(NULL == workspace){
		RNP::IRA::FreeWorkspace(n, n_wanted, n_arnoldi, u_workspace);
		delete u_workspace;
	}
	return info;
}



int RNP::IRA::ShiftInvert(
	size_t n, const std::complex<double> &shift, RNP::IRA::Ashifted Aop, RNP::IRA::Bfunc Bop,
	size_t n_wanted, size_t n_arnoldi,   // n_wanted+1 <= n_arnoldi <= n
	RNP::IRA::SortFunc sorter,
	std::complex<double> *w,             // eigenvalues (length n_wanted)
	std::complex<double> *v, size_t ldv, // eigenvectors (size ldv-by-n_wanted), ldv >= n
	RNP::IRA::Params *params,
	RNP::IRA::Workspace *workspace,
	void *Adata, void *Bdata, void *SortData
){
	RNP::IRA::Params *u_params = params;
	if(NULL == params){
		u_params = new Params();
	}
	RNP::IRA::Workspace *u_workspace = workspace;
	if(NULL == workspace){
		u_workspace = new RNP::IRA::Workspace();
		RNP::IRA::AllocateWorkspace(n, n_wanted, n_arnoldi, u_workspace);
	}
	
	if(u_params->tol <= 0){
		u_params->tol = std::numeric_limits<double>::epsilon();
	}
	
	int iparam[11], ipntr[14];
	int ido = 0;
	int info = 0;
	iparam[0] = 1; // ishift
	iparam[2] = u_params->max_iterations;
	iparam[6] = 3;
	char bmat = (NULL == Bop) ? 'I' : 'G';
	arz_data static_data;
	static_data.sort_data = SortData;
	size_t lworkl = (3*n_arnoldi+5)*n_arnoldi;
	
	do{
		znaupd_(
			&static_data, &ido, bmat,
			n, sorter, n_wanted, u_params->tol,
			u_workspace->resid,
			n_arnoldi, v, ldv,
			iparam, ipntr,
			u_workspace->workd, u_workspace->workl, lworkl, u_workspace->rwork, &info);

		if(ido == -1){
			Bop(n, &u_workspace->workd[ipntr[0]], &u_workspace->workd[ipntr[1]], Bdata);
			RNP::TBLAS::Copy(n, &u_workspace->workd[ipntr[1]], 1, &u_workspace->workd[ipntr[0]], 1);
			Aop(n, shift, &u_workspace->workd[ipntr[0]], &u_workspace->workd[ipntr[1]], Adata);
		}else if(ido == 1) {
			Aop(n, shift, &u_workspace->workd[ipntr[2]], &u_workspace->workd[ipntr[1]], Adata);
		}else if(ido == 2){
			Bop(n, &u_workspace->workd[ipntr[0]], &u_workspace->workd[ipntr[1]], Bdata);
		}else{
			break;
		}
	}while(1);
	if(info >= 0){
		int ierr;
		zneupd_(
			!(u_params->no_eigenvectors), 'A', u_workspace->bwork,
			w, v, ldv,
			&shift, u_workspace->workev,
			bmat, n, sorter, SortData, n_wanted,
			u_params->tol, u_workspace->resid, n_arnoldi, v, ldv,
			iparam, ipntr,
			u_workspace->workd, u_workspace->workl, lworkl,
			u_workspace->rwork, &ierr);
		if(ierr == 0){
			info = iparam[4]; // nconv
		}else if(ierr > 0){
			info = -100*ierr;
		}
	}
	
	if(NULL == params){
		delete u_params;
	}
	if(NULL == workspace){
		RNP::IRA::FreeWorkspace(n, n_wanted, n_arnoldi, u_workspace);
		delete u_workspace;
	}
	return info;
}

