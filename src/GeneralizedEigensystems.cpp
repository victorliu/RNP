#ifdef DEBUG_ZGG
# include <iostream>
# include "IO.h"
#endif

#include <limits>
#include "TBLAS.h"
#include "TLASupport.h"
#include "GeneralizedEigensystems.h"
#include <cmath>
#include <complex>
#include <algorithm>

static inline double _abs1(const std::complex<double> &z){
	return fabs(z.real()) + fabs(z.imag());
}

/*
Copyright (c) 1992-2010 The University of Tennessee.  All rights reserved.

Additional copyrights may follow

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

- Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer. 
  
- Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer listed
  in this license in the documentation and/or other materials
  provided with the distribution.
  
- Neither the name of the copyright holders nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.
  
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT  
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT 
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT  
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

*/

static void zggbal_(const char *job, size_t n, std::complex<double> *a, size_t 
		lda, std::complex<double> *b, size_t ldb, int *ilo, int *ihi, 
		double *lscale, double *rscale, double *work)
{
	using namespace std;

	int i, j, k, l;
	int ip1, jp1;

/*  ZGGBAL balances a pair of general complex matrices (A,B).  This */
/*  involves, first, permuting A and B by similarity transformations to */
/*  isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N */
/*  elements on the diagonal; and second, applying a diagonal similarity */
/*  transformation to rows and columns ILO to IHI to make the rows */
/*  and columns as close in norm as possible. Both steps are optional. */

/*  Balancing may reduce the 1-norm of the matrices, and improve the */
/*  accuracy of the computed eigenvalues and/or eigenvectors in the */
/*  generalized eigenvalue problem A*x = lambda*B*x. */

/*  Arguments */
/*  ========= */

/*  JOB     (input) CHARACTER*1 */
/*          Specifies the operations to be performed on A and B: */
/*          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0 */
/*                  and RSCALE(I) = 1.0 for i=1,...,N; */
/*          = 'P':  permute only; */
/*          = 'S':  scale only; */
/*          = 'B':  both permute and scale. */

/*  N       (input) INTEGER */
/*          The order of the matrices A and B.  N >= 0. */

/*  A       (input/output) COMPLEX*16 array, dimension (LDA,N) */
/*          On entry, the input matrix A. */
/*          On exit, A is overwritten by the balanced matrix. */
/*          If JOB = 'N', A is not referenced. */

/*  LDA     (input) INTEGER */
/*          The leading dimension of the array A. LDA >= max(1,N). */

/*  B       (input/output) COMPLEX*16 array, dimension (LDB,N) */
/*          On entry, the input matrix B. */
/*          On exit, B is overwritten by the balanced matrix. */
/*          If JOB = 'N', B is not referenced. */

/*  LDB     (input) INTEGER */
/*          The leading dimension of the array B. LDB >= max(1,N). */

/*  ILO     (output) INTEGER */
/*  IHI     (output) INTEGER */
/*          ILO and IHI are set to integers such that on exit */
/*          A(i,j) = 0 and B(i,j) = 0 if i > j and */
/*          j = 1,...,ILO-1 or i = IHI+1,...,N. */
/*          If JOB = 'N' or 'S', ILO = 1 and IHI = N. */

/*  LSCALE  (output) DOUBLE PRECISION array, dimension (N) */
/*          Details of the permutations and scaling factors applied */
/*          to the left side of A and B.  If P(j) is the index of the */
/*          row interchanged with row j, and D(j) is the scaling factor */
/*          applied to row j, then */
/*            LSCALE(j) = P(j)    for J = 1,...,ILO-1 */
/*                      = D(j)    for J = ILO,...,IHI */
/*                      = P(j)    for J = IHI+1,...,N. */
/*          The order in which the interchanges are made is N to IHI+1, */
/*          then 1 to ILO-1. */

/*  RSCALE  (output) DOUBLE PRECISION array, dimension (N) */
/*          Details of the permutations and scaling factors applied */
/*          to the right side of A and B.  If P(j) is the index of the */
/*          column interchanged with column j, and D(j) is the scaling */
/*          factor applied to column j, then */
/*            RSCALE(j) = P(j)    for J = 1,...,ILO-1 */
/*                      = D(j)    for J = ILO,...,IHI */
/*                      = P(j)    for J = IHI+1,...,N. */
/*          The order in which the interchanges are made is N to IHI+1, */
/*          then 1 to ILO-1. */

/*  WORK    (workspace) REAL array, dimension (lwork) */
/*          lwork must be at least max(1,6n) when JOB = 'S' or 'B', and */
/*          at least 1 when JOB = 'N' or 'P'. */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

/*  Further Details */
/*  =============== */

/*  See R.C. WARD, Balancing the generalized eigenvalue problem, */
/*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152. */

/*  ===================================================================== */

/*     Test the input parameters */

	/* Parameter adjustments */
	size_t a_offset = 1 + lda;
	a -= a_offset;
	size_t b_offset = 1 + ldb;
	b -= b_offset;
	--lscale;
	--rscale;
	--work;

/*     Quick return if possible */

	if (n == 0) {
		*ilo = 1;
		*ihi = n;
		return;
	}

	if (n == 1) {
		*ilo = 1;
		*ihi = n;
		lscale[1] = 1.;
		rscale[1] = 1.;
		return;
	}

	if (job[0] == 'N'){
		*ilo = 1;
		*ihi = n;
		for (i = 1; i <= (int)n; ++i) {
			lscale[i] = 1.;
			rscale[i] = 1.;
		}
		return;
	}

	l = n;
	if (job[0] != 'S') {
		// Permute the matrices A and B to isolate the eigenvalues.

		// Find row with one nonzero in columns 1 through L
		bool found_row_perm;
		do{
			found_row_perm = false;
			for (i = l; i >= 1; --i) {
				bool found = false;
				for (j = 1; j <= l-1; ++j) {
					jp1 = j + 1;
					if ((0. != a[i+j*lda]) || (0. != b[i+j*ldb])) {
						// We found a nonzero in column j, to the left of (i,l)
						found = true;
						break;
					}
				}
				if(!found){
					j = l;
				}else{
					bool found2 = false;
					for (j = jp1; j <= l; ++j) {
						if ((0. != a[i+j*lda]) || (0. != b[i+j*ldb])) {
							// We found another nonzero in column j, between the first nonzero and (i,l+1)
							found2 = true;
							break;
						}
					}
					if(found2){
						// We found more than 1 nonzero, so let's try the next row up
						continue;
					}else{
						j = jp1 - 1;
					}
				}
				int m = l;
				// Permute rows M and I
				lscale[m] = (double) i;
				if (i != m) {
					RNP::TBLAS::Swap(n, &a[i + 1 * lda], lda, &a[m + 1 * lda], lda);
					RNP::TBLAS::Swap(n, &b[i + 1 * ldb], ldb, &b[m + 1 * ldb], ldb);
				}
				// Permute columns M and J
				rscale[m] = (double) j;
				if (j != m) {
					RNP::TBLAS::Swap(l, &a[j * lda + 1], 1, &a[m * lda + 1], 1);
					RNP::TBLAS::Swap(l, &b[j * ldb + 1], 1, &b[m * ldb + 1], 1);
				}

				--l;
				if (l == 1) {
					// We have completely deflated the matrix from bottom up
					rscale[1] = 1.;
					lscale[1] = 1.;
				}else{
					found_row_perm = true;
				}
				break;
			}
		}while(found_row_perm);

		// Find column with one nonzero in rows K through N */
		k = 1;
		bool found_col_perm;

		do{
			found_col_perm = false;
			for (j = k; j <= l; ++j) {
				bool found = false;
				for (i = k; i <= l-1; ++i) {
					ip1 = i + 1;
					if ((0. != a[i+j*lda]) || (0. != b[i+j*ldb])) {
						found = true;
						break;
					}
				}
				if(!found){
					i = l;
				}else{
					bool found2 = false;
					for (i = ip1; i <= l; ++i) {
						if ((0. != a[i+j*lda]) || (0. != b[i+j*ldb])) {
							found2 = true;
							break;
						}
					}
					if(found2){
						continue;
					}else{
						i = ip1 - 1;
					}
				}
				int m = k;
				// Permute rows M and I
				lscale[m] = (double) i;
				if (i != m) {
					RNP::TBLAS::Swap(n - k + 1, &a[i + k * lda], lda, &a[m + k * lda], lda);
					RNP::TBLAS::Swap(n - k + 1, &b[i + k * ldb], ldb, &b[m + k * ldb], ldb);
				}
				// Permute columns M and J
				rscale[m] = (double) j;
				if (j != m) {
					RNP::TBLAS::Swap(l, &a[j * lda + 1], 1, &a[m * lda + 1], 1);
					RNP::TBLAS::Swap(l, &b[j * ldb + 1], 1, &b[m * ldb + 1], 1);
				}

				++k;
				found_col_perm = true;
				break;
			}
		}while(found_col_perm);
	}

	// End of permutations
	*ilo = k;
	*ihi = l;

	if (job[0] == 'P') { // if permutations was all that was requested, end here
		for (i = *ilo; i <= *ihi; ++i) {
			lscale[i] = 1.;
			rscale[i] = 1.;
		}
		return;
	}

	if (*ilo == *ihi) {
		return;
	}

	// Balance the submatrix in rows ILO to IHI.

	const size_t nr = *ihi - *ilo + 1;
	for (i = *ilo; i <= *ihi; ++i) {
		rscale[i] = 0.;
		lscale[i] = 0.;

		work[i] = 0.;
		work[i + n] = 0.;
		work[i + 2*n] = 0.;
		work[i + n * 3] = 0.;
		work[i + 2*n] = 0.;
		work[i + n * 5] = 0.;
	}

	// Compute right side vector in resulting linear equations
	for (i = *ilo; i <= *ihi; ++i) {
		for (j = *ilo; j <= *ihi; ++j) {
			double ta;
			if (0. == a[i+j*lda]) {
				ta = 0.;
			}else{
				ta = log10(_abs1(a[i+j*lda]));
			}
			double tb;
			if (0. == b[i+j*ldb]) {
				tb = 0.;
			}else{
				tb = log10(_abs1(b[i+j*ldb]));
			}
			work[i + 2*n] = work[i + 2*n] - ta - tb;
			work[j + n * 5] = work[j + n * 5] - ta - tb;
		}
	}

	const double coef = 1. / (double) (2*nr);
	const double coef2 = coef * coef;
	const double coef5 = coef2 * .5;
	double beta = 0.;
	
	const size_t nrp2 = nr + 2; // iteration limit
	size_t it = 1;
	double gamma_prev = 1;
	do{ // Start generalized conjugate gradient iteration
		double gamma = RNP::TBLAS::Dot(nr, &work[*ilo + 2*n], 1, &work[*ilo + 2*n], 1)
			+ RNP::TBLAS::Dot(nr, &work[*ilo + n * 5], 1, &work[*ilo + n * 5], 1);

		double ew = 0.;
		double ewc = 0.;
		for (i = *ilo; i <= *ihi; ++i) {
			ew += work[i + 2*n];
			ewc += work[i + n * 5];
		}

		double diff = ew - ewc;
		gamma = coef * gamma - coef2 * (ew*ew + ewc*ewc) - coef5 * (diff*diff);
		if (gamma == 0.) {
			break;
		}
		if (it > 1) {
			beta = gamma / gamma_prev;
		}
		double t = coef5 * (ewc - ew * 3.);
		double tc = coef5 * (ew - ewc * 3.);

		RNP::TBLAS::Scale(nr, beta, &work[*ilo], 1);
		RNP::TBLAS::Scale(nr, beta, &work[*ilo + n], 1);

		RNP::TBLAS::Axpy(nr, coef, &work[*ilo + 2*n], 1, &work[*ilo + n], 1);
		RNP::TBLAS::Axpy(nr, coef, &work[*ilo + n * 5], 1, &work[*ilo], 1);

		for (i = *ilo; i <= *ihi; ++i) {
			work[i] += tc;
			work[i + n] += t;
		}

		// Apply matrix to vector
		for (i = *ilo; i <= *ihi; ++i) {
			int kount = 0;
			double sum = 0.;
			for (j = *ilo; j <= *ihi; ++j) {
				if (0. != a[i+j*lda]) {
					++kount;
					sum += work[j];
				}
				if (0. != b[i+j*ldb]) {
					++kount;
					sum += work[j];
				}
			}
			work[i + 2*n] = (double) kount * work[i + n] + sum;
		}

		for (j = *ilo; j <= *ihi; ++j) {
			int kount = 0;
			double sum = 0.;
			for (i = *ilo; i <= *ihi; ++i) {
				if (0. != a[i+j*lda]) {
					++kount;
					sum += work[i + n];
				}
				if (0. != b[i+j*ldb]) {
					++kount;
					sum += work[i + n];
				}
			}
			work[j + 3*n] = (double) kount * work[j] + sum;
		}

		double sum = RNP::TBLAS::Dot(nr, &work[*ilo + n], 1, &work[*ilo + 2*n], 1) 
				+ RNP::TBLAS::Dot(nr, &work[*ilo], 1, &work[*ilo + 3*n], 1);
		double alpha = gamma / sum;

		// Determine correction to current iteration
		double cmax = 0.;
		for (i = *ilo; i <= *ihi; ++i) {
			double cor = alpha * work[i + n];
			if (abs(cor) > cmax) {
				cmax = abs(cor);
			}
			lscale[i] += cor;
			cor = alpha * work[i];
			if (abs(cor) > cmax) {
				cmax = abs(cor);
			}
			rscale[i] += cor;
		}
		if (cmax < .5) {
			break;
		}

		RNP::TBLAS::Axpy(nr, -alpha, &work[*ilo + 2*n], 1, &work[*ilo + 2*n], 1);
		RNP::TBLAS::Axpy(nr, -alpha, &work[*ilo + 3*n], 1, &work[*ilo + 5*n], 1);

		gamma_prev = gamma;
		++it;
	}while(it <= nrp2); // End generalized conjugate gradient iteration

	const double sfmin = std::numeric_limits<double>::min();
	const double sfmax = 1. / sfmin;
	const int lsfmin = (int) (log10(sfmin) + 1.);
	const int lsfmax = (int) (log10(sfmax));
	for (i = *ilo; i <= *ihi; ++i) {
		double da;
		
		size_t irab = 1+RNP::TBLAS::MaximumIndex(n - *ilo + 1, (std::complex<double>*)&a[i + *ilo * lda], lda);
		double rab = abs(a[i + (irab + *ilo - 1) * lda]);
		irab = 1+RNP::TBLAS::MaximumIndex(n - *ilo + 1, (std::complex<double>*)&b[i + *ilo * ldb], ldb);
		da = abs(b[i + (irab + *ilo - 1) * ldb]);
		if(da > rab){ rab = da; }
		int lrab = (int)(log10(rab + sfmin) + 1.);
		int ir;
		if(lscale[i] >= 0){
			ir = (int) (lscale[i] + 0.5);	
		}else{
			ir = (int) (lscale[i] - 0.5);
		}
		if(lsfmin > ir){ ir = lsfmin; }
		if(lsfmax < ir){ ir = lsfmax; }
		if(lsfmax-lrab < ir){ ir = lsfmax-lrab; }
		lscale[i] = pow(10., ir);
		
		size_t icab = 1+RNP::TBLAS::MaximumIndex(*ihi, (std::complex<double>*)&a[i * lda + 1], 1);
		double cab = abs(a[icab + i * lda]);
		icab = 1+RNP::TBLAS::MaximumIndex(*ihi, (std::complex<double>*)&b[i * ldb + 1], 1);
		da = abs(b[icab + i * ldb]);
		if(da > cab){ cab = da; }
		int lcab = (int) (log10(cab + sfmin) + 1.);
		int jc;
		if(rscale[i] >= 0){
			jc = (int) (rscale[i] + 0.5);
		}else{
			jc = (int) (rscale[i] - 0.5);
		}
		if(lsfmin > jc){ jc = lsfmin; }
		if(lsfmax < jc){ jc = lsfmax; }
		if(lsfmax-lcab < jc){ jc = lsfmax-lcab; }
		rscale[i] = pow(10., jc);
	}

	// Row scaling of matrices A and B
	for (i = *ilo; i <= *ihi; ++i) {
		RNP::TBLAS::Scale(n - *ilo + 1, lscale[i], (std::complex<double>*)&a[i + *ilo * lda], lda);
		RNP::TBLAS::Scale(n - *ilo + 1, lscale[i], (std::complex<double>*)&b[i + *ilo * ldb], ldb);
	}

	// Column scaling of matrices A and B
	for (j = *ilo; j <= *ihi; ++j) {
		RNP::TBLAS::Scale(*ihi, rscale[j], (std::complex<double>*)&a[j * lda + 1], 1);
		RNP::TBLAS::Scale(*ihi, rscale[j], (std::complex<double>*)&b[j * ldb + 1], 1);
	}
}


void zggbak_(char *job, char *side, int n, int ilo, 
		int ihi, const double *lscale, const double *rscale, int m, 
		std::complex<double> *v, int ldv)
{
	using namespace std;
	

// Purpose
// =======

// ZGGBAK forms the right or left eigenvectors of a complex generalized
// eigenvalue problem A*x = lambda*B*x, by backward transformation on
// the computed eigenvectors of the balanced pair of matrices output by
// ZGGBAL.

// Arguments
// =========

// JOB     (input) CHARACTER*1
//          Specifies the type of backward transformation required:
//          = 'N':  do nothing, return immediately;
//          = 'P':  do backward transformation for permutation only;
//          = 'S':  do backward transformation for scaling only;
//          = 'B':  do backward transformations for both permutation and scaling.
//          JOB must be the same as the argument JOB supplied to ZGGBAL.

// SIDE    (input) CHARACTER*1
//          = 'R':  V contains right eigenvectors;
//          = 'L':  V contains left eigenvectors.

// N       (input) int
//          The number of rows of the matrix V.  N >= 0.

// ILO     (input) int
// IHI     (input) int
//          The ints ILO and IHI determined by ZGGBAL.
//          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.

// LSCALE  (input) DOUBLE PRECISION array, dimension (N)
//          Details of the permutations and/or scaling factors applied
//          to the left side of A and B, as returned by ZGGBAL.

// RSCALE  (input) DOUBLE PRECISION array, dimension (N)
//          Details of the permutations and/or scaling factors applied
//          to the right side of A and B, as returned by ZGGBAL.

// M       (input) int
//          The number of columns of the matrix V.  M >= 0.

// V       (input/output) COMPLEX*16 array, dimension (LDV,M)
//          On entry, the matrix of right or left eigenvectors to be
//          transformed, as returned by ZTGEVC.
//          On exit, V is overwritten by the transformed eigenvectors.

// LDV     (input) int
//          The leading dimension of the matrix V. LDV >= max(1,N).

// INFO    (output) int
//          = 0:  successful exit.
//          < 0:  if INFO = -i, the i-th argument had an illegal value.

// Further Details
// ===============

// See R.C. Ward, Balancing the generalized eigenvalue problem,
//                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.

	int v_offset = 1 + ldv;
	v -= v_offset;

	const bool rightv = (side[0] == 'R');
	const bool leftv = (side[0] == 'L');

	if (n == 0 || m == 0 || job[0] == 'N') {
		return;
	}

	if (ilo != ihi) {
		// Backward balance
		if ((job[0] == 'S') || (job[0] == 'B')) {
			// Backward transformation on right eigenvectors
			if (rightv) {
				for (int i = ilo; i <= ihi; ++i) {
					RNP::TBLAS::Scale(m, rscale[i-1], &v[i + ldv], ldv);
				}
			}

			// Backward transformation on left eigenvectors
			if (leftv) {
				for (int i = ilo; i <= ihi; ++i) {
					RNP::TBLAS::Scale(m, lscale[i-1], &v[i + ldv], ldv);
				}
			}
		}
	}

	// Backward permutation
	if ((job[0] == 'P') || (job[0] == 'B')) {
		// Backward permutation on right eigenvectors

		if (rightv) {
			if (ilo != 1) {
				for (int i = ilo - 1; i >= 1; --i) {
					int k = (int) rscale[i-1];
					if (k != i) {
						RNP::TBLAS::Swap(m, &v[i + ldv], ldv, &v[k + ldv], ldv);
					}
				}
			}
			if (ihi != n) {
				for (int i = ihi + 1; i <= n; ++i) {
					int k = (int) rscale[i-1];
					if (k != i) {
						RNP::TBLAS::Swap(m, &v[i + ldv], ldv, &v[k + ldv], ldv);
					}
				}
			}
		}
		// Backward permutation on left eigenvectors
		if (leftv) {
			if (ilo != 1) {
				for (int i = ilo - 1; i >= 1; --i) {
					int k = (int) lscale[i-1];
					if (k != i) {
						RNP::TBLAS::Swap(m, &v[i + ldv], ldv, &v[k + ldv], ldv);
					}
				}
			}
			if (ihi != n) {
				for (int i = ihi + 1; i <= n; ++i) {
					int k = (int) lscale[i-1];
					if (k != i) {
						RNP::TBLAS::Swap(m, &v[i + ldv], ldv, &v[k + ldv], ldv);
					}
				}
			}
		}
	}
}

static int zhgeqz_(char *job, char *compq, char *compz, size_t n, 
	size_t ilo, size_t ihi, std::complex<double> *h, size_t ldh, 
	std::complex<double> *t, size_t ldt, std::complex<double> *alpha, std::complex<double> *
	beta, std::complex<double> *q, size_t ldq, std::complex<double> *z, size_t 
	ldz, std::complex<double> *work, size_t lwork, double *rwork)
{
	using namespace std;
	// System generated locals
	size_t h_offset, q_offset, t_offset, z_offset;

	// Local variables
	double c;
	int j;
	std::complex<double> s;
	int jc;
	int jr;
	int jch;
	bool ilq, ilz;
	
	std::complex<double> ctemp;
	bool ilschr;
	int icompq, ilastm;
	int ischur;
	int icompz, ifirst;
	int ifrstm;

	// ZHGEQZ computes the eigenvalues of a complex matrix pair (H,T),
	// where H is an upper Hessenberg matrix and T is upper triangular,
	// using the single-shift QZ method.
	// Matrix pairs of this type are produced by the reduction to
	// generalized upper Hessenberg form of a complex matrix pair (A,B):

	//    A = Q1*H*Z1**H,  B = Q1*T*Z1**H,

	// as computed by ZGGHRD.

	// If JOB='S', then the Hessenberg-triangular pair (H,T) is
	// also reduced to generalized Schur form,

	//    H = Q*S*Z**H,  T = Q*P*Z**H,

	// where Q and Z are unitary matrices and S and P are upper triangular.

	// Optionally, the unitary matrix Q from the generalized Schur
	// factorization may be postmultiplied into an input matrix Q1, and the
	// unitary matrix Z may be postmultiplied into an input matrix Z1.
	// If Q1 and Z1 are the unitary matrices from ZGGHRD that reduced
	// the matrix pair (A,B) to generalized Hessenberg form, then the output
	// matrices Q1*Q and Z1*Z are the unitary factors from the generalized
	// Schur factorization of (A,B):

	//    A = (Q1*Q)*S*(Z1*Z)**H,  B = (Q1*Q)*P*(Z1*Z)**H.

	// To avoid overflow, eigenvalues of the matrix pair (H,T)
	// (equivalently, of (A,B)) are computed as a pair of complex values
	// (alpha,beta).  If beta is nonzero, lambda = alpha / beta is an
	// eigenvalue of the generalized nonsymmetric eigenvalue problem (GNEP)
	//    A*x = lambda*B*x
	// and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
	// alternate form of the GNEP
	//    mu*A*y = B*y.
	// The values of alpha and beta for the i-th eigenvalue can be read
	// directly from the generalized Schur form:  alpha = S(i,i),
	// beta = P(i,i).

	// Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
	//      Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
	//      pp. 241--256.

	// Arguments
	// =========

	// JOB     (input) CHARACTER*1
	//         = 'E': Compute eigenvalues only;
	//         = 'S': Computer eigenvalues and the Schur form.

	// COMPQ   (input) CHARACTER*1
	//         = 'N': Left Schur vectors (Q) are not computed;
	//         = 'I': Q is initialized to the unit matrix and the matrix Q
	//                of left Schur vectors of (H,T) is returned;
	//         = 'V': Q must contain a unitary matrix Q1 on entry and
	//                the product Q1*Q is returned.

	// COMPZ   (input) CHARACTER*1
	//         = 'N': Right Schur vectors (Z) are not computed;
	//         = 'I': Q is initialized to the unit matrix and the matrix Z
	//                of right Schur vectors of (H,T) is returned;
	//         = 'V': Z must contain a unitary matrix Z1 on entry and
	//                the product Z1*Z is returned.

	// N       (input) INTEGER
	//         The order of the matrices H, T, Q, and Z.  N >= 0.

	// ILO     (input) INTEGER
	// IHI     (input) INTEGER
	//         ILO and IHI mark the rows and columns of H which are in
	//         Hessenberg form.  It is assumed that A is already upper
	//         triangular in rows and columns 1:ILO-1 and IHI+1:N.
	//         If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.

	// H       (input/output) COMPLEX*16 array, dimension (LDH, N)
	//         On entry, the N-by-N upper Hessenberg matrix H.
	//         On exit, if JOB = 'S', H contains the upper triangular
	//         matrix S from the generalized Schur factorization.
	//         If JOB = 'E', the diagonal of H matches that of S, but
	//         the rest of H is unspecified.

	// LDH     (input) INTEGER
	//         The leading dimension of the array H.  LDH >= max( 1, N ).

	// T       (input/output) COMPLEX*16 array, dimension (LDT, N)
	//         On entry, the N-by-N upper triangular matrix T.
	//         On exit, if JOB = 'S', T contains the upper triangular
	//         matrix P from the generalized Schur factorization.
	//         If JOB = 'E', the diagonal of T matches that of P, but
	//         the rest of T is unspecified.

	// LDT     (input) INTEGER
	//         The leading dimension of the array T.  LDT >= max( 1, N ).

	// ALPHA   (output) COMPLEX*16 array, dimension (N)
	//         The complex scalars alpha that define the eigenvalues of
	//         GNEP.  ALPHA(i) = S(i,i) in the generalized Schur
	//         factorization.

	// BETA    (output) COMPLEX*16 array, dimension (N)
	//         The real non-negative scalars beta that define the
	//         eigenvalues of GNEP.  BETA(i) = P(i,i) in the generalized
	//         Schur factorization.

	//         Together, the quantities alpha = ALPHA(j) and beta = BETA(j)
	//         represent the j-th eigenvalue of the matrix pair (A,B), in
	//         one of the forms lambda = alpha/beta or mu = beta/alpha.
	//         Since either lambda or mu may overflow, they should not,
	//         in general, be computed.

	// Q       (input/output) COMPLEX*16 array, dimension (LDQ, N)
	//         On entry, if COMPZ = 'V', the unitary matrix Q1 used in the
	//         reduction of (A,B) to generalized Hessenberg form.
	//         On exit, if COMPZ = 'I', the unitary matrix of left Schur
	//         vectors of (H,T), and if COMPZ = 'V', the unitary matrix of
	//         left Schur vectors of (A,B).
	//         Not referenced if COMPZ = 'N'.

	// LDQ     (input) INTEGER
	//         The leading dimension of the array Q.  LDQ >= 1.
	//         If COMPQ='V' or 'I', then LDQ >= N.

	// Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
	//         On entry, if COMPZ = 'V', the unitary matrix Z1 used in the
	//         reduction of (A,B) to generalized Hessenberg form.
	//         On exit, if COMPZ = 'I', the unitary matrix of right Schur
	//         vectors of (H,T), and if COMPZ = 'V', the unitary matrix of
	//         right Schur vectors of (A,B).
	//         Not referenced if COMPZ = 'N'.

	// LDZ     (input) INTEGER
	//         The leading dimension of the array Z.  LDZ >= 1.
	//         If COMPZ='V' or 'I', then LDZ >= N.

	// WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
	//         On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.

	// LWORK   (input) INTEGER
	//         The dimension of the array WORK.  LWORK >= max(1,N).

	//         If LWORK = -1, then a workspace query is assumed; the routine
	//         only calculates the optimal size of the WORK array, returns
	//         this value as the first entry of the WORK array, and no error
	//         message related to LWORK is issued by XERBLA.

	// RWORK   (workspace) DOUBLE PRECISION array, dimension (N)

	// INFO    (output) INTEGER
	//         = 0: successful exit
	//         < 0: if INFO = -i, the i-th argument had an illegal value
	//         = 1,...,N: the QZ iteration did not converge.  (H,T) is not
	//                    in Schur form, but ALPHA(i) and BETA(i),
	//                    i=INFO+1,...,N should be correct.
	//         = N+1,...,2n: the shift calculation failed.  (H,T) is not
	//                    in Schur form, but ALPHA(i) and BETA(i),
	//                    i=INFO-N+1,...,N should be correct.

	// Further Details
	// ===============

	// We assume that complex ABS works as long as its value is less than
	// overflow.

	h_offset = 1 + ldh;
	h -= h_offset;
	t_offset = 1 + ldt;
	t -= t_offset;
	--alpha;
	--beta;
	q_offset = 1 + ldq;
	q -= q_offset;
	z_offset = 1 + ldz;
	z -= z_offset;
	--work;
	--rwork;


	if(job[0] == 'E'){
		ilschr = false;
		ischur = 1;
	}else if(job[0] == 'S'){
		ilschr = true;
		ischur = 2;
	}else{
		ischur = 0;
	}

	if(compq[0] == 'N'){
		ilq = false;
		icompq = 1;
	}else if(compq[0] == 'V'){
		ilq = true;
		icompq = 2;
	}else if(compq[0] == 'I'){
		ilq = true;
		icompq = 3;
	}else{
		icompq = 0;
	}

	if(compz[0] == 'N'){
		ilz = false;
		icompz = 1;
	}else if(compz[0] == 'V'){
		ilz = true;
		icompz = 2;
	}else if(compz[0] == 'I'){
		ilz = true;
		icompz = 3;
	}else{
		icompz = 0;
	}

	// Check Argument Values
	int info = 0;
	if(ischur == 0){
		info = -1;
	}else if(icompq == 0){
		info = -2;
	}else if(icompz == 0){
		info = -3;
	}else if(n < 0){
		info = -4;
	}else if(ilo < 1){
		info = -5;
	}else if(ihi > n || ihi < ilo - 1){
		info = -6;
	}else if(ldh < n){
		info = -8;
	}else if(ldt < n){
		info = -10;
	}else if(ldq < 1 || ilq && ldq < n){
		info = -14;
	}else if(ldz < 1 || ilz && ldz < n){
		info = -16;
	}
	if(info != 0){
		return info;
	}
	if(n <= 0){
		return info;
	}

	// Initialize Q and Z
	if(icompq == 3){
		RNP::TBLAS::SetMatrix<'F'>(n, n, std::complex<double>(0), std::complex<double>(1), &q[q_offset], ldq);
	}
	if(icompz == 3){
		RNP::TBLAS::SetMatrix<'F'>(n, n, std::complex<double>(0), std::complex<double>(1), &z[z_offset], ldz);
	}

#ifdef DEBUG_ZGG
	std::cout << "q upon entering QZ iteration:" << std::endl;
	RNP::IO::PrintMatrix(n,n,&q[q_offset],ldq) << std::endl;
#endif
	// Machine Constants
	const int in = ihi+1 - ilo;
	const double safmin = std::numeric_limits<double>::min();
	const double ulp = std::numeric_limits<double>::epsilon() * 2;
	double anorm, bnorm;
	RNP::TLASupport::CheapHessenbergNorm<'F'>(in, &h[ilo+ilo* ldh], ldh, &anorm);
	RNP::TLASupport::CheapHessenbergNorm<'F'>(in, &t[ilo+ilo* ldt], ldt, &bnorm);

	// Computing MAX
	const double atol = max(safmin,ulp*anorm);
	// Computing MAX
	const double btol = max(safmin,ulp*bnorm);
	const double ascale = 1. / max(safmin,anorm);
	const double bscale = 1. / max(safmin,bnorm);


	// Set Eigenvalues IHI+1:N
	for(j = ihi+1; j <= (int)n; ++j){
		double absb = abs(t[j+j*ldt]);
		if(absb > safmin){
			std::complex<double> signbc = std::conj(t[j+j*ldt]/absb);
			t[j+j*ldt] = absb;
			if(ilschr){
				RNP::TBLAS::Scale(j-1, signbc, &t[j*ldt+1], 1);
				RNP::TBLAS::Scale(j, signbc, &h[j*ldh+1], 1);
			}else{
				h[j+j*ldh] *= signbc;
			}
			if(ilz){
				RNP::TBLAS::Scale(n, signbc, &z[j*ldz+1], 1);
			}
		}else{
			t[j+j*ldt] = 0;
		}
		alpha[j] = h[j+j*ldh];
		beta[j] = t[j+j*ldt];
	}

	if(ilo <= ihi){ // If IHI < ILO, skip QZ steps
		// MAIN QZ ITERATION LOOP

		// Initialize dynamic indices

		// Eigenvalues ILAST+1:N have been found.
		// Column operations modify rows IFRSTM:whatever
		// Row operations modify columns whatever:ILASTM

		// If only eigenvalues are being computed, then
		// IFRSTM is the row of the last splitting row above row ILAST;
		// this is always at least ILO.
		// IITER counts iterations since the last eigenvalue was found,
		// to tell when to use an extraordinary shift.
		// MAXIT is the maximum number of QZ sweeps allowed.

		int ilast = ihi;
		if(ilschr){
			ifrstm = 1;
			ilastm = n;
		}else{
			ifrstm = ilo;
			ilastm = ihi;
		}
		int iiter = 0;
		std::complex<double> eshift = 0;
		const int maxit = (ihi - ilo+1) * 30;

		bool converged = false;
		for(int jiter = 1; jiter <= maxit; ++jiter){
			// Split the matrix if possible.
			// Two tests:
			// 1: H(j,j-1)=0  or  j=ILO
			// 2: T(j,j)=0
			
			bool ilazro;
			bool ilazr2;
			// Special case: j=ILAST
			if(ilast == (int)ilo){
				goto L60_zhgeqz;
			}else if(_abs1(h[ilast+(ilast-1)*ldh]) <= atol){
				h[ilast+(ilast-1)*ldh] = 0;
				goto L60_zhgeqz;
			}

			if(abs(t[ilast+ilast*ldt]) <= btol){
				t[ilast+ilast*ldt] = 0;
				goto L50_zhgeqz;
			}

			// General case: j<ILAST
			for(j = ilast - 1; j >= (int)ilo; --j){
				// Test 1: for H(j,j-1)=0 or j=ILO
				if(j == (int)ilo){
					ilazro = true;
				}else{
					if(_abs1(h[j+(j-1)*ldh]) <= atol){
						h[j+(j-1)*ldh] = 0;
						ilazro = true;
					}else{
						ilazro = false;
					}
				}

				// Test 2: for T(j,j)=0
				if(abs(t[j+j*ldt]) < btol){
					t[j+j*ldt] = 0;

					// Test 1a: Check for 2 consecutive small subdiagonals in A
					ilazr2 = false;
					if(! ilazro){
						if( _abs1(h[j+(j-1)*ldh])*( ascale*_abs1(h[j+1+j*ldh])) <= _abs1(h[j+j*ldh])*( ascale*atol ) ){
							ilazr2 = true;
						}
					}

					// If both tests pass (1 & 2), i.e., the leading diagonal
					// element of B in the block is zero, split a 1x1 block off
					// at the top. (I.e., at the J-th row/column) The leading
					// diagonal element of the remainder can also be zero, so
					// this may have to be done repeatedly.
					if(ilazro || ilazr2){
						for(jch = j; jch <= ilast-1; ++jch){
							ctemp = h[jch+jch*ldh];
							RNP::TLASupport::GeneratePlaneRotation(ctemp, h[jch+1+jch*ldh], &c, &s, &h[jch + jch*ldh]);
							h[jch+1+jch*ldh] = 0;
							RNP::TLASupport::ApplyPlaneRotation(ilastm-jch, &h[jch+(jch+1)*ldh], ldh, &h[jch+1+(jch+1)*ldh], ldh, c, s);
							RNP::TLASupport::ApplyPlaneRotation(ilastm-jch, &t[jch+(jch+1)*ldt], ldt, &t[jch+1+(jch+1)*ldt], ldt, c, s);
							if(ilq){
								RNP::TLASupport::ApplyPlaneRotation(n, &q[jch*ldq+1], 1, &q[(jch+1)*ldq+1], 1, c, std::conj(s));
							}
							if(ilazr2){
								h[jch+(jch-1)*ldh] *= c;
							}
							ilazr2 = false;
							if(_abs1(t[jch+1+(jch+1)*ldt]) >= btol){
								if(jch+1 >= ilast){
									goto L60_zhgeqz;
								}else{
									ifirst = jch+1;
									goto L70_zhgeqz;
								}
							}
							t[jch+1+(jch+1)*ldt] = 0;
						}
						goto L50_zhgeqz;
					}else{
						// Only test 2 passed -- chase the zero to T(ILAST,ILAST)
						// Then process as in the case T(ILAST,ILAST)=0
						for(jch = j; jch <= ilast-1; ++jch){
							ctemp = t[jch+(jch+1)*ldt];
							RNP::TLASupport::GeneratePlaneRotation(ctemp, t[jch+1+(jch+1)*ldt], &c, &s, &t[jch+(jch+1)*ldt]);
							t[jch+1+(jch+1)*ldt] = 0;
							if(jch < ilastm - 1){
								RNP::TLASupport::ApplyPlaneRotation(ilastm-jch-1, &t[jch+(jch+2)*ldt], ldt, &t[jch+1+(jch+2)*ldt], ldt, c, s);
							}
							RNP::TLASupport::ApplyPlaneRotation(ilastm-jch+2, &h[jch+(jch-1)*ldh], ldh, &h[jch+1+(jch-1)*ldh], ldh, c, s);
							if(ilq){
								RNP::TLASupport::ApplyPlaneRotation(n, &q[jch*ldq+1], 1, &q[(jch+1)*ldq+1], 1, c, std::conj(s));
							}
							ctemp = h[jch+1+jch*ldh];
							RNP::TLASupport::GeneratePlaneRotation(ctemp, h[jch+1+(jch-1)*ldh], &c, &s, &h[jch+1+jch*ldh]);
							h[jch+1+(jch-1)*ldh] = 0;
							RNP::TLASupport::ApplyPlaneRotation(jch+1-ifrstm, &h[ifrstm + jch*ldh], 1, &h[ifrstm + (jch-1)*ldh], 1, c, s);
							RNP::TLASupport::ApplyPlaneRotation(jch-ifrstm, &t[ifrstm + jch*ldt], 1, &t[ifrstm + (jch-1)*ldt], 1, c, s);
							if(ilz){
								RNP::TLASupport::ApplyPlaneRotation(n, &z[jch*ldz+1], 1, &z[(jch-1)*ldz+1], 1, c, s);
							}
						}
						goto L50_zhgeqz;
					}
				}else if(ilazro){
					// Only test 1 passed -- work on J:ILAST
					ifirst = j;
					goto L70_zhgeqz;
				}
				// Neither test passed -- try next J
			}

			// (Drop-through is "impossible")
			return 2*n+1;

			// T(ILAST,ILAST)=0 -- clear H(ILAST,ILAST-1) to split off a 1x1 block.
L50_zhgeqz:
			ctemp = h[ilast+ilast*ldh];
			RNP::TLASupport::GeneratePlaneRotation(ctemp, h[ilast + (ilast-1)*ldh], &c, &s, &h[ilast + ilast*ldh]);
			h[ilast+(ilast-1)*ldh] = 0;
			RNP::TLASupport::ApplyPlaneRotation(ilast-ifrstm, &h[ifrstm+ilast*ldh], 1, &h[ifrstm + (ilast-1)*ldh], 1, c, s);
			RNP::TLASupport::ApplyPlaneRotation(ilast-ifrstm, &t[ifrstm+ilast*ldt], 1, &t[ifrstm + (ilast-1)*ldt], 1, c, s);
			if(ilz){
				RNP::TLASupport::ApplyPlaneRotation(n, &z[ilast*ldz+1], 1, &z[(ilast-1)*ldz+1], 1, c, s);
			}

			// H(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA

L60_zhgeqz:
			double absb = abs(t[ilast + ilast*ldt]);
			if(absb > safmin){
				std::complex<double> signbc = std::conj(t[ilast+ilast*ldt]/absb);
				t[ilast+ilast*ldt] = absb;
				if(ilschr){
					RNP::TBLAS::Scale(ilast-ifrstm, signbc, &t[ifrstm + ilast*ldt], 1);
					RNP::TBLAS::Scale(ilast+1-ifrstm, signbc, &h[ifrstm + ilast*ldh], 1);
				}else{
					h[ilast+ilast*ldh] *= signbc;
				}
				if(ilz){
					RNP::TBLAS::Scale(n, signbc, &z[ilast*ldz+1], 1);
				}
			}else{
				t[ilast+ilast*ldt] = 0;
			}
			alpha[ilast] = h[ilast+ilast*ldh];
			beta[ilast] = t[ilast+ilast*ldt];

			// Go to next block -- exit if finished.
			--ilast;
			if(ilast < (int)ilo){
				converged = true;
				break;
			}

			// Reset counters
			iiter = 0;
			eshift = 0;
			if(! ilschr){
				ilastm = ilast;
				if(ifrstm > ilast){
					ifrstm = ilo;
				}
			}
			continue;

			// QZ step

			// This iteration only involves rows/columns IFIRST:ILAST.  We
			// assume IFIRST < ILAST, and that the diagonal of B is non-zero.

L70_zhgeqz:
			++iiter;
			if(! ilschr){
				ifrstm = ifirst;
			}

			// Compute the Shift.

			// At this point, IFIRST < ILAST, and the diagonal elements of
			// T(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
			// magnitude)
			std::complex<double> shift;
			if(iiter / 10 * 10 != iiter){
				// The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
				// the bottom-right 2x2 block of A inv(B) which is nearest to
				// the bottom-right element.

				// We factor B as U*D, where U has unit diagonals, and
				// compute (A*inv(D))*inv(U).

				std::complex<double> u12 = (bscale*t[ilast - 1 + ilast*ldt]) / (bscale * t[ilast + ilast*ldt]);
				std::complex<double> ad11 = (ascale*h[ilast - 1+(ilast-1)*ldh]) / (bscale*t[ilast - 1+(ilast-1)*ldt]);
				std::complex<double> ad21 = (ascale*h[ilast + (ilast-1)*ldh]) / (bscale*t[ilast - 1+(ilast-1)*ldt]);
				std::complex<double> ad12 = (ascale*h[ilast - 1 + ilast*ldh]) / (bscale*t[ilast + ilast*ldt]);
				std::complex<double> ad22 = (ascale*h[ilast + ilast*ldh]) / (bscale*t[ilast + ilast*ldt]);
				std::complex<double> abi22 = ad22-u12*ad21;
				
				std::complex<double> t1 = (ad11+abi22) / double(2);
				std::complex<double> rtdisc = sqrt(t1*t1 + ad12*ad21 - ad11*ad22);
				if(((t1-abi22)*std::conj(rtdisc)).real() <= 0.){
					shift = t1 + rtdisc;
				}else{
					shift = t1 - rtdisc;
				}
			}else{
				// Exceptional shift.  Chosen for no particularly good reason.
				eshift += std::conj((ascale*h[ilast - 1 + ilast*ldh]) / (bscale*t[ilast - 1+(ilast-1)*ldt]));
				shift = eshift;
			}

			// Now check for two consecutive small subdiagonals.
			bool found = false;
			int istart;
			for(j = ilast - 1; j >= ifirst+1; --j){
				istart = j;
				ctemp = ascale*h[j+j*ldh] - shift*(bscale*t[j+j*ldt]);
				double abs1ctemp = _abs1(ctemp);
				double abs1next = ascale * _abs1(h[j+1+j*ldh]);
				double maxabs1 = max(abs1ctemp,abs1next);
				if(maxabs1 < 1. && maxabs1 != 0.){
					abs1ctemp /= maxabs1;
					abs1next /= maxabs1;
				}
				if(_abs1(h[j+(j-1)*ldh]) * abs1next <= abs1ctemp * atol){
					found = true;
					break;
				}
			}
			if(!found){
				istart = ifirst;
				ctemp = ascale*h[ifirst+ifirst*ldh] - shift*(bscale*t[ifirst+ifirst*ldt]);
			}

			// Do an implicit-shift QZ sweep.

			// Initial Q
			{
				std::complex<double> ctemp3;
				RNP::TLASupport::GeneratePlaneRotation(ctemp, ascale * h[istart+1+istart*ldh], &c, &s, &ctemp3);
			}

			// Sweep
			for(j = istart; j <= ilast-1; ++j){
				if(j > istart){
					ctemp = h[j+(j-1)*ldh];
					RNP::TLASupport::GeneratePlaneRotation(ctemp, h[j+1+(j-1)*ldh], &c, &s, &h[j+(j-1)*ldh]);
					h[j+1+(j-1)*ldh] = 0;
				}

				for(jc = j; jc <= ilastm; ++jc){
					ctemp = c*h[j+jc*ldh] + s*h[j+1+jc*ldh];
					h[j+1+jc*ldh] = -std::conj(s)*h[j+jc*ldh] + c*h[j+1+jc*ldh];
					h[j+jc*ldh] = ctemp;
					ctemp = c*t[j+jc*ldt] + s*t[j+1+jc*ldt];
					t[j+1+jc*ldt] = -std::conj(s)*t[j+jc*ldt] + c*t[j+1+jc*ldt];
					t[j+jc*ldt] = ctemp;
				}
				if(ilq){
					for(jr = 1; jr <= (int)n; ++jr){
						ctemp = c*q[jr+j*ldq] + std::conj(s)*q[jr+(j+1)*ldq];
						q[jr+(j+1)*ldq] = -s*q[jr + j*ldq] + c*q[jr+(j+1)*ldq];
						q[jr+j*ldq] = ctemp;
					}
				}

				ctemp = t[j+1+(j+1)*ldt];
				RNP::TLASupport::GeneratePlaneRotation(ctemp, t[j+1+j*ldt], &c, &s, &t[j+1+(j+1)*ldt]);
				t[j+1+j*ldt] = 0;

				// Computing MIN
				for(jr = ifrstm; jr <= min(j+2,ilast); ++jr){
					ctemp = c*h[jr+(j+1)*ldh] + s*h[jr+j*ldh];
					h[jr+j*ldh] = -std::conj(s)*h[jr+(j+1)*ldh] + c*h[jr+j*ldh];
					h[jr+(j+1)*ldh] = ctemp;
				}
				for(jr = ifrstm; jr <= j; ++jr){
					ctemp = c*t[jr+(j+1)*ldt] + s*t[jr+j*ldt];
					t[jr+j*ldt] = -std::conj(s)*t[jr+(j+1)*ldt] + c*t[jr+j*ldt];
					t[jr+(j+1)*ldt] = ctemp;
				}
				if(ilz){
					for(jr = 1; jr <= (int)n; ++jr){
						ctemp = c*z[jr+(j+1)*ldz] + s*z[jr+j*ldz];
						z[jr+j*ldz] = -std::conj(s)*z[jr+(j+1)*ldz] + c*z[jr+j*ldz];
						z[jr+(j+1)*ldz] = ctemp;
					}
				}
			}
		}
		if(!converged){
			return ilast;
		}
		// Successful completion of all QZ steps
	}

	// Set Eigenvalues 1:ILO-1
	for(j = 1; j <= (int)ilo-1; ++j){
		double absb = abs(t[j+j*ldt]);
		if(absb > safmin){
			std::complex<double> signbc = std::conj(t[j+j*ldt]/absb);
			t[j+j*ldt] = absb;
			if(ilschr){
				RNP::TBLAS::Scale(j-1, signbc, &t[j*ldt+1], 1);
				RNP::TBLAS::Scale(j, signbc, &h[j*ldh+1], 1);
			}else{
				h[j+j*ldh] *= signbc;
			}
			if(ilz){
				RNP::TBLAS::Scale(n, signbc, &z[j*ldz+1], 1);
			}
		}else{
			t[j+j*ldt] = 0;
		}
		alpha[j] = h[j+j*ldh];
		beta[j] = t[j+j*ldt];
	}
#ifdef DEBUG_ZGG
	std::cout << "q at end of QZ iteration:" << std::endl;
	RNP::IO::PrintMatrix(n,n,&q[q_offset],ldq) << std::endl;
#endif
	return 0;
}

static int ztgevc_(const char *howmny, bool *select, 
		size_t n, std::complex<double> *s, size_t lds, std::complex<double> *p, size_t 
		ldp, std::complex<double> *vl, size_t ldvl, std::complex<double> *vr, size_t 
		ldvr, size_t mm, size_t *m, std::complex<double> *work, double *rwork)
{
	using namespace std;
	
	/* System generated locals */
	size_t p_offset, s_offset, vl_offset, 
			vr_offset;
	double d__1, d__2, d__3, d__4, d__5, d__6;

	/* Local variables */
	std::complex<double> d;
	int i, j;
	int je, im, jr;
	double big;
	bool lsa, lsb;
	double ulp;
	int ibeg, ieig, iend;
	double dmin__;
	int isrc;
	double temp;
	double xmax, scale;
	bool ilall;
	int iside;
	double sbeta;

	double small;
	bool bcoml;
	double anorm, bnorm;
	bool compr;
	bool ilbbad;
	double acoefa, bcoefa, acoeff;
	std::complex<double> bcoeff;
	bool ilback;
	double ascale, bscale;
	std::complex<double> salpha;
	double safmin;
	double bignum;
	bool ilcomp;
	int ihwmny;

/*  Purpose */
/*  ======= */

/*  ZTGEVC computes some or all of the right and/or left eigenvectors of */
/*  a pair of complex matrices (S,P), where S and P are upper triangular. */
/*  Matrix pairs of this type are produced by the generalized Schur */
/*  factorization of a complex matrix pair (A,B): */

/*     A = Q*S*Z**H,  B = Q*P*Z**H */

/*  as computed by ZGGHRD + ZHGEQZ. */

/*  The right eigenvector x and the left eigenvector y of (S,P) */
/*  corresponding to an eigenvalue w are defined by: */

/*     S*x = w*P*x,  (y**H)*S = w*(y**H)*P, */

/*  where y**H denotes the conjugate tranpose of y. */
/*  The eigenvalues are not input to this routine, but are computed */
/*  directly from the diagonal elements of S and P. */

/*  This routine returns the matrices X and/or Y of right and left */
/*  eigenvectors of (S,P), or the products Z*X and/or Q*Y, */
/*  where Z and Q are input matrices. */
/*  If Q and Z are the unitary factors from the generalized Schur */
/*  factorization of a matrix pair (A,B), then Z*X and Q*Y */
/*  are the matrices of right and left eigenvectors of (A,B). */

/*  Arguments */
/*  ========= */

/*  SIDE    (input) CHARACTER*1 */
/*          = 'R': compute right eigenvectors only; */
/*          = 'L': compute left eigenvectors only; */
/*          = 'B': compute both right and left eigenvectors. */

/*  HOWMNY  (input) CHARACTER*1 */
/*          = 'A': compute all right and/or left eigenvectors; */
/*          = 'B': compute all right and/or left eigenvectors, */
/*                 backtransformed by the matrices in VR and/or VL; */
/*          = 'S': compute selected right and/or left eigenvectors, */
/*                 specified by the bool array SELECT. */

/*  SELECT  (input) bool array, dimension (N) */
/*          If HOWMNY='S', SELECT specifies the eigenvectors to be */
/*          computed.  The eigenvector corresponding to the j-th */
/*          eigenvalue is computed if SELECT(j) = .TRUE.. */
/*          Not referenced if HOWMNY = 'A' or 'B'. */

/*  N       (input) INTEGER */
/*          The order of the matrices S and P.  N >= 0. */

/*  S       (input) COMPLEX*16 array, dimension (LDS,N) */
/*          The upper triangular matrix S from a generalized Schur */
/*          factorization, as computed by ZHGEQZ. */

/*  LDS     (input) INTEGER */
/*          The leading dimension of array S.  LDS >= max(1,N). */

/*  P       (input) COMPLEX*16 array, dimension (LDP,N) */
/*          The upper triangular matrix P from a generalized Schur */
/*          factorization, as computed by ZHGEQZ.  P must have real */
/*          diagonal elements. */

/*  LDP     (input) INTEGER */
/*          The leading dimension of array P.  LDP >= max(1,N). */

/*  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM) */
/*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must */
/*          contain an N-by-N matrix Q (usually the unitary matrix Q */
/*          of left Schur vectors returned by ZHGEQZ). */
/*          On exit, if SIDE = 'L' or 'B', VL contains: */
/*          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P); */
/*          if HOWMNY = 'B', the matrix Q*Y; */
/*          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by */
/*                      SELECT, stored consecutively in the columns of */
/*                      VL, in the same order as their eigenvalues. */
/*          Not referenced if SIDE = 'R'. */

/*  LDVL    (input) INTEGER */
/*          The leading dimension of array VL.  LDVL >= 1, and if */
/*          SIDE = 'L' or 'l' or 'B' or 'b', LDVL >= N. */

/*  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM) */
/*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must */
/*          contain an N-by-N matrix Q (usually the unitary matrix Z */
/*          of right Schur vectors returned by ZHGEQZ). */
/*          On exit, if SIDE = 'R' or 'B', VR contains: */
/*          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P); */
/*          if HOWMNY = 'B', the matrix Z*X; */
/*          if HOWMNY = 'S', the right eigenvectors of (S,P) specified by */
/*                      SELECT, stored consecutively in the columns of */
/*                      VR, in the same order as their eigenvalues. */
/*          Not referenced if SIDE = 'L'. */

/*  LDVR    (input) INTEGER */
/*          The leading dimension of the array VR.  LDVR >= 1, and if */
/*          SIDE = 'R' or 'B', LDVR >= N. */

/*  MM      (input) INTEGER */
/*          The number of columns in the arrays VL and/or VR. MM >= M. */

/*  M       (output) INTEGER */
/*          The number of columns in the arrays VL and/or VR actually */
/*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M */
/*          is set to N.  Each selected eigenvector occupies one column. */

/*  WORK    (workspace) COMPLEX*16 array, dimension (2*N) */

/*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N) */

/*  INFO    (output) INTEGER */
/*          = 0:  successful exit. */
/*          < 0:  if INFO = -i, the i-th argument had an illegal value. */

/*  ===================================================================== */


	if(NULL != vl){
		if(NULL != vr){
			iside = 3;
			bcoml = true;
			compr = true;
		}else{
			iside = 2;
			bcoml = true;
			compr = false;
		}
	}else if(NULL != vr){
		iside = 1;
		bcoml = false;
		compr = true;
	}else{
		iside = -1;
	}
	/* Parameter adjustments */
	--select;
	s_offset = 1 + lds;
	s -= s_offset;
	p_offset = 1 + ldp;
	p -= p_offset;
	vl_offset = 1 + ldvl;
	vl -= vl_offset;
	vr_offset = 1 + ldvr;
	vr -= vr_offset;
	--work;
	--rwork;

	/* Function Body */
	if (howmny[0] == 'A') {
		ihwmny = 1;
		ilall = true;
		ilback = false;
	} else if (howmny[0] == 'S') {
		ihwmny = 2;
		ilall = false;
		ilback = false;
	} else if (howmny[0] == 'B') {
		ihwmny = 3;
		ilall = true;
		ilback = true;
	} else {
		ihwmny = -1;
	}

	int info = 0;
	if (iside < 0) {
		info = -1;
	} else if (ihwmny < 0) {
		info = -2;
	} else if (n < 0) {
		info = -4;
	} else if (lds < max((size_t)1,n)) {
		info = -6;
	} else if (ldp < max((size_t)1,n)) {
		info = -8;
	}
	if (info != 0) {
		return info;
	}

/*     Count the number of eigenvectors */

	if (! ilall) {
		im = 0;
		for (j = 1; j <= (int)n; ++j) {
			if (select[j]) {
				++im;
			}
		}
	} else {
		im = n;
	}

/*     Check diagonal of B */

	ilbbad = false;
	for (j = 1; j <= (int)n; ++j) {
		if (p[j+j*ldp].imag() != 0.) {
			ilbbad = true;
		}
	}

	if (ilbbad) {
		info = -7;
	} else if (bcoml && ldvl < n || ldvl < 1) {
		info = -10;
	} else if (compr && ldvr < n || ldvr < 1) {
		info = -12;
	} else if ((int)mm < im) {
		info = -13;
	}
	if (info != 0) {
		return info;
	}

/*     Quick return if possible */

	*m = im;
	if (n == 0) {
		return 0;
	}

/*     Machine Constants */

	safmin = std::numeric_limits<double>::min();
	ulp = std::numeric_limits<double>::epsilon() * 2;
	small = safmin * n / ulp;
	big = 1. / small;
	bignum = 1. / (safmin * n);

/*     Compute the 1-norm of each column of the strictly upper triangular */
/*     part of A and B to check for possible overflow in the triangular */
/*     solver. */

	anorm = _abs1(s[1+1*lds]);
	bnorm = _abs1(p[1+1*ldp]);
	rwork[1] = 0.;
	rwork[n + 1] = 0.;
	for (j = 2; j <= (int)n; ++j) {
		rwork[j] = 0.;
		rwork[n + j] = 0.;
		for (i = 1; i <= j-1; ++i) {
			rwork[j] += _abs1(s[i+j*lds]);
			rwork[n + j] += _abs1(p[i+j*ldp]);
		}
/* Computing MAX */
		anorm = max(anorm,rwork[j] + _abs1(s[j+j*lds]));
/* Computing MAX */
		bnorm = max(bnorm,rwork[n+j] + _abs1(p[j+j*ldp]));
	}

	ascale = 1. / max(anorm,safmin);
	bscale = 1. / max(bnorm,safmin);

/*     Left eigenvectors */

	if (bcoml) {
		ieig = 0;

/*        Main loop over eigenvalues */

		for (je = 1; je <= (int)n; ++je) {
			if (ilall) {
				ilcomp = true;
			} else {
				ilcomp = select[je];
			}
			if (ilcomp) {
				++ieig;

				if (_abs1(s[je+je*lds]) <= safmin && abs(p[je+je*ldp].real()) <= safmin) {
/*                 Singular matrix pencil -- return unit eigenvector */
					for (jr = 1; jr <= (int)n; ++jr) {
						vl[jr+ieig*ldvl] = 0;
					}
					vl[ieig+ieig*ldvl] = 1;
					continue;
				}

/*              Non-singular eigenvalue: */
/*              Compute coefficients  a  and  b  in */
/*                   H */
/*                 y  ( a A - b B ) = 0 */

/* Computing MAX */
				d__4 = _abs1(s[je+je*lds]) * ascale, d__5 = abs(p[je+je*ldp].real()) * bscale, d__4 = max(d__4,d__5);
				temp = 1. / max(d__4,safmin);
				salpha = (temp*s[je+je*lds]) * ascale;
				sbeta = temp * p[je+je*ldp].real() * bscale;
				acoeff = sbeta * ascale;
				bcoeff = bscale * salpha;

/*              Scale to avoid underflow */

				lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
				lsb = _abs1(salpha) >= safmin && _abs1(bcoeff) < small;

				scale = 1.;
				if (lsa) {
					scale = small / abs(sbeta) * min(anorm,big);
				}
				if (lsb) {
/* Computing MAX */
					d__4 = small / _abs1(salpha) * min(bnorm,big);
					scale = max(scale,d__4);
				}
				if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
					d__5 = 1., d__6 = abs(acoeff), d__5 = max(d__5,d__6), 
							d__6 = _abs1(bcoeff);
					d__3 = scale, d__4 = 1. / (safmin * max(d__5,d__6));
					scale = min(d__3,d__4);
					if (lsa) {
						acoeff = ascale * (scale * sbeta);
					} else {
						acoeff = scale * acoeff;
					}
					if (lsb) {
						bcoeff = bscale * (scale*salpha);
					} else {
						bcoeff *= scale;
					}
				}

				acoefa = abs(acoeff);
				bcoefa = _abs1(bcoeff);
				xmax = 1.;
				for (jr = 1; jr <= (int)n; ++jr) {
					work[jr] = 0;
				}
				work[je] = 1;
/* Computing MAX */
				d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, 
						d__1 = max(d__1,d__2);
				dmin__ = max(d__1,safmin);

/*                                              H */
/*              Triangular solve of  (a A - b B)  y = 0 */

/*                                      H */
/*              (rowwise in  (a A - b B) , or columnwise in a A - b B) */

				for (j = je + 1; j <= (int)n; ++j) {

/*                 Compute */
/*                       j-1 */
/*                 SUM = sum  conjg( a*S(k,j) - b*P(k,j) )*x(k) */
/*                       k=je */
/*                 (Scale if necessary) */

					temp = 1. / xmax;
					if (acoefa * rwork[j] + bcoefa * rwork[n + j] > bignum * 
							temp) {
						for (jr = je; jr <= j-1; ++jr) {
							work[jr] *= temp;
						}
						xmax = 1.;
					}
					std::complex<double> suma(0), sumb(0);

					for (jr = je; jr <= j-1; ++jr) {
						suma += std::conj(((std::complex<double>*)(s))[jr+j*lds]) * work[jr];
						sumb += std::conj(((std::complex<double>*)(p))[jr+j*ldp]) * work[jr];
					}
					std::complex<double> sum = acoeff*suma - std::conj(bcoeff)*sumb;


/*                 Form x(j) = - SUM / conjg( a*S(j,j) - b*P(j,j) ) */

/*                 with scaling and perturbation of the denominator */
					d = std::conj(acoeff*s[j+j*lds] - bcoeff*p[j+j*ldp]);
					if (_abs1(d) <= dmin__) {
						d = dmin__;
					}

					double abs1d = _abs1(d);
					if (abs1d < 1.) {
						double abs1sum = RNP::TBLAS::_RealOrComplexChooser<std::complex<double> >::_abs1(sum);
						if (abs1sum >= bignum * abs1d) {
							temp = 1. / abs1sum;
							for (jr = je; jr <= j-1; ++jr) {
								work[jr] *= temp;
							}
							xmax *= temp;
							sum *= temp;
						}
					}
					work[j] = -sum/d;
/* Computing MAX */
					xmax = max(xmax,RNP::TBLAS::_RealOrComplexChooser<std::complex<double> >::_abs1(work[j]));
				}

/*              Back transform eigenvector if HOWMNY='B'. */

				if (ilback) {
					RNP::TBLAS::MultMV<'N'>(n, n + 1 - je, double(1), &vl[1+je*ldvl], ldvl, &work[je], 1, double(0), &work[n+1], 1);
					isrc = 2;
					ibeg = 1;
				} else {
					isrc = 1;
					ibeg = je;
				}

/*              Copy and scale eigenvector into column of VL */

				xmax = 0.;
				for (jr = ibeg; jr <= (int)n; ++jr) {
					xmax = max(xmax,RNP::TBLAS::_RealOrComplexChooser<std::complex<double> >::_abs1(work[(isrc-1)*n+jr]));
				}

				if (xmax > safmin) {
					temp = 1. / xmax;
					for (jr = ibeg; jr <= (int)n; ++jr) {
						vl[jr+ieig*ldvl] = temp * work[(isrc-1)*n+jr];
					}
				} else {
					ibeg = n + 1;
				}

				for (jr = 1; jr <= ibeg - 1; ++jr) {
					vl[jr+ieig*ldvl] = 0;
				}

			}
		}
	}

/*     Right eigenvectors */

	if (compr) {
		ieig = im + 1;

/*        Main loop over eigenvalues */

		for (je = n; je >= 1; --je) {
			if (ilall) {
				ilcomp = true;
			} else {
				ilcomp = select[je];
			}
			if (ilcomp) {
				--ieig;

				if (_abs1(s[je+je*lds]) <= safmin && abs(p[je+je*ldp].real()) <= safmin) {

/*                 Singular matrix pencil -- return unit eigenvector */

					for (jr = 1; jr <= (int)n; ++jr) {
						vr[jr+ieig*ldvr] = 0;
					}
					vr[ieig+ieig*ldvr] = 1;
					continue;
				}

/*              Non-singular eigenvalue: */
/*              Compute coefficients  a  and  b  in */

/*              ( a A - b B ) x  = 0 */

/* Computing MAX */
				d__4 = _abs1(s[je+je*lds]) * ascale, d__5 = _abs1(
						p[je+je*ldp].real()) * bscale, d__4 = max(d__4,d__5);
				temp = 1. / max(d__4,safmin);
				salpha = (temp*s[je+je*lds]) * ascale;
				sbeta = (temp * p[je+je*ldp].real()) * bscale;
				acoeff = sbeta * ascale;
				bcoeff = bscale * salpha;

/*              Scale to avoid underflow */

				lsa = abs(sbeta) >= safmin && abs(acoeff) < small;
				lsb = _abs1(salpha) >= safmin && _abs1(bcoeff) < small;

				scale = 1.;
				if (lsa) {
					scale = small / abs(sbeta) * min(anorm,big);
				}
				if (lsb) {
/* Computing MAX */
					d__3 = scale, d__4 = small / _abs1(salpha) * min(bnorm,big);
					scale = max(d__3,d__4);
				}
				if (lsa || lsb) {
/* Computing MIN */
/* Computing MAX */
					d__5 = 1., d__6 = abs(acoeff), d__5 = max(d__5,d__6), 
							d__6 = _abs1(bcoeff);
					d__3 = scale, d__4 = 1. / (safmin * max(d__5,d__6));
					scale = min(d__3,d__4);
					if (lsa) {
						acoeff = ascale * (scale * sbeta);
					} else {
						acoeff = scale * acoeff;
					}
					if (lsb) {
						bcoeff = bscale * (scale*salpha);
					} else {
						bcoeff *= scale;
					}
				}

				acoefa = abs(acoeff);
				bcoefa = _abs1(bcoeff);
				xmax = 1.;
				for (jr = 1; jr <= (int)n; ++jr) {
					work[jr] = 0;
				}
				work[je] = 1;
/* Computing MAX */
				d__1 = ulp * acoefa * anorm, d__2 = ulp * bcoefa * bnorm, 
						d__1 = max(d__1,d__2);
				dmin__ = max(d__1,safmin);

/*              Triangular solve of  (a A - b B) x = 0  (columnwise) */

/*              WORK(1:j-1) contains sums w, */
/*              WORK(j+1:JE) contains x */

				for (jr = 1; jr <= je - 1; ++jr) {
					work[jr] = acoeff*((std::complex<double>*)s)[jr+je*lds] - bcoeff*((std::complex<double>*)p)[jr+je*ldp];
				}
				work[je] = 1;

				for (j = je - 1; j >= 1; --j) {

/*                 Form x(j) := - w(j) / d */
/*                 with scaling and perturbation of the denominator */
					d = acoeff*s[j+j*lds] - bcoeff*p[j+j*ldp];
					if (_abs1(d) <= dmin__) {
						d = dmin__;
					}

					double _abs1d = _abs1(d);
					if (_abs1d < 1.) {
						double abs1workj = RNP::TBLAS::_RealOrComplexChooser<std::complex<double> >::_abs1(work[j]);
						if (abs1workj >= bignum * _abs1d) {
							temp = 1. / abs1workj;
							for (jr = 1; jr <= je; ++jr) {
								work[jr] *= temp;
							}
						}
					}
					work[j] = -work[j]/d;

					if (j > 1) {

/*                    w = w + x(j)*(a S(*,j) - b P(*,j) ) with scaling */

						double abs1workj = RNP::TBLAS::_RealOrComplexChooser<std::complex<double> >::_abs1(work[j]);
						if (abs1workj > 1.) {
							temp = 1. / abs1workj;
							if (acoefa * rwork[j] + bcoefa * rwork[n + j] >= 
									bignum * temp) {
								for (jr = 1; jr <= je; ++jr) {
									work[jr] *= temp;
								}
							}
						}

						std::complex<double> ca = acoeff*work[j];
						std::complex<double> cb = bcoeff*work[j];
						for (jr = 1; jr <= j - 1; ++jr) {
							work[jr] += ca*((std::complex<double>*)s)[jr+j*lds] - cb*((std::complex<double>*)p)[jr+j*ldp];
						}
					}
				}

/*              Back transform eigenvector if HOWMNY='B'. */

				if (ilback) {
					RNP::TBLAS::MultMV<'N'>(n, je, double(1), &vr[vr_offset], ldvr, &work[1], 1, double(0), &work[n+1], 1);
					isrc = 2;
					iend = n;
				} else {
					isrc = 1;
					iend = je;
				}

/*              Copy and scale eigenvector into column of VR */
				xmax = 0.;
				for (jr = 1; jr <= iend; ++jr) {
					xmax = max(xmax,RNP::TBLAS::_RealOrComplexChooser<std::complex<double> >::_abs1(work[(isrc-1)*n+jr]));
				}

				if (xmax > safmin) {
					temp = 1. / xmax;
					for (jr = 1; jr <= iend; ++jr) {
						vr[jr+ieig*ldvr] = temp * work[(isrc-1)*n+jr];
					}
				} else {
					iend = 0;
				}

				for (jr = iend + 1; jr <= (int)n; ++jr) {
					vr[jr+ieig*ldvr] = 0;
				}

			}
		}
	}
	return 0;
}


// Computes for a pair of N-by-N complex nonsymmetric matrices (A,B),
// the generalized eigenvalues, and optionally, the left and/or right
// generalized eigenvectors.

// A generalized eigenvalue for a pair of matrices (A,B) is a scalar
// lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
// singular. It is usually represented as the pair (alpha,beta), as
// there is a reasonable interpretation for beta=0, and even for both
// being zero.

// The right generalized eigenvector v(j) corresponding to the
// generalized eigenvalue lambda(j) of (A,B) satisfies
//              A * v(j) = lambda(j) * B * v(j).
// The left generalized eigenvector u(j) corresponding to the
// generalized eigenvalues lambda(j) of (A,B) satisfies
//              u(j)^H * A = lambda(j) * u(j)^H * B
// where u(j)^H is the conjugate-transpose of u(j).

// Arguments
// =========

// N       The order of the matrices A, B, VL, and VR.  N >= 0.
// A       (input/output) COMPLEX*16 array, dimension (LDA, N)
//         On entry, the matrix A in the pair (A,B).
//         On exit, A has been overwritten.
// LDA     The leading dimension of A.  LDA >= max(1,N).
// B       (input/output) COMPLEX*16 array, dimension (LDB, N)
//         On entry, the matrix B in the pair (A,B).
//         On exit, B has been overwritten.
// LDB     The leading dimension of B.  LDB >= max(1,N).

// ALPHA   (output) COMPLEX*16 array, dimension (N)
// BETA    (output) COMPLEX*16 array, dimension (N)
//         On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
//         generalized eigenvalues.

//         Note: the quotients ALPHA(j)/BETA(j) may easily over- or
//         underflow, and BETA(j) may even be zero.  Thus, the user
//         should avoid naively computing the ratio alpha/beta.
//         However, ALPHA will be always less than and usually
//         comparable with norm(A) in magnitude, and BETA always less
//         than and usually comparable with norm(B).

// VL      (output) COMPLEX*16 array, dimension (LDVL,N)
//         If VL != NULL, the left generalized eigenvectors u(j) are
//         stored one after another in the columns of VL, in the same
//         order as their eigenvalues.
//         Each eigenvector is scaled so the largest component has
//         abs(real part) + abs(imag. part) = 1.
// LDVL    The leading dimension of the matrix VL. LDVL >= 1, and
//         if VL != NULL, LDVL >= N.

// VR      (output) COMPLEX*16 array, dimension (LDVR,N)
//         If VR != NULL, the right generalized eigenvectors v(j) are
//         stored one after another in the columns of VR, in the same
//         order as their eigenvalues.
//         Each eigenvector is scaled so the largest component has
//         abs(real part) + abs(imag. part) = 1.
// LDVR    The leading dimension of the matrix VR. LDVR >= 1, and
//         if VR != NULL, LDVR >= N.

// WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,2*N))
// RWORK   (workspace/output) DOUBLE PRECISION array, dimension (8*N)

// return: = 0:  successful exit
//         < 0:  if INFO = -i, the i-th argument had an illegal value.
//         =1,...,N:
//               The QZ iteration failed.  No eigenvectors have been
//               calculated, but ALPHA(j) and BETA(j) should be
//               correct for j=INFO+1,...,N.
//         > N:  =N+1: other then QZ iteration failed in DHGEQZ,
//               =N+2: error return from DTGEVC.
int RNP::GeneralizedEigensystem(size_t n, 
		std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb, 
		std::complex<double> *alpha, std::complex<double> *beta,
		std::complex<double> *vl, size_t ldvl, std::complex<double> *vr, size_t ldvr,
		std::complex<double> *work, double *rwork)
{
	/* System generated locals */
	size_t a_offset, b_offset, vl_offset,
			vr_offset;

	using namespace std;

	/* Local variables */
	int ihi, ilo;
	double anrm, bnrm;
	int itau;
	double temp;
	int iwrk;
	int ileft, icols, irwrk, irows;
	bool ilascl, ilbscl;
	bool ldumma[1];
	double anrmto;
	double bnrmto;

	const bool ilvl = (NULL != vl);
	const bool ilvr = (NULL != vr);
	const bool ilv = (ilvl || ilvr);

	/* Parameter adjustments */
	a_offset = 1 + lda;
	a -= a_offset;
	b_offset = 1 + ldb;
	b -= b_offset;
	--alpha;
	--beta;
	vl_offset = 1 + ldvl;
	vl -= vl_offset;
	vr_offset = 1 + ldvr;
	vr -= vr_offset;
	--work;
	--rwork;


	// Test the input arguments

	if (n < 0) {
		return -3;
	} else if (lda < max((size_t)1,n)) {
		return -5;
	} else if (ldb < max((size_t)1,n)) {
		return -7;
	} else if (ldvl < 1 || ilvl && ldvl < n) {
		return -11;
	} else if (ldvr < 1 || ilvr && ldvr < n) {
		return -13;
	}

	if (n == 0) {
		return 0;
	}

	const double eps = std::numeric_limits<double>::epsilon() * 2;
	const double smlnum = sqrt(std::numeric_limits<double>::min()) / eps;
	const double bignum = 1. / smlnum;

	// Scale A if max element outside range [SMLNUM,BIGNUM]
	RNP::TLASupport::CheapMatrixNorm<'M'>(n, n, &a[a_offset], lda, &anrm);
	ilascl = false;
	if (anrm > 0. && anrm < smlnum) {
		anrmto = smlnum;
		ilascl = true;
	} else if (anrm > bignum) {
		anrmto = bignum;
		ilascl = true;
	}
	if (ilascl) {
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, anrm, anrmto, n, n, &a[a_offset], lda);
	}

	// Scale B if max element outside range [SMLNUM,BIGNUM]
	RNP::TLASupport::CheapMatrixNorm<'M'>(n, n, &b[b_offset], ldb, &bnrm);
	ilbscl = false;
	if (bnrm > 0. && bnrm < smlnum) {
		bnrmto = smlnum;
		ilbscl = true;
	} else if (bnrm > bignum) {
		bnrmto = bignum;
		ilbscl = true;
	}
	if (ilbscl) {
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, bnrm, bnrmto, n, n, &b[b_offset], ldb);
	}

	// Permute the matrices A, B to isolate eigenvalues if possible
	// (Real Workspace: need 6*N)
	ileft = 1;
	int iright = n + 1;
	irwrk = iright + n;
	zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[ileft], &rwork[iright], &rwork[irwrk]);

	// Reduce B to triangular form (QR decomposition of B)
	// (Complex Workspace: need N)

	irows = ihi + 1 - ilo;
	if (ilv) {
		icols = n + 1 - ilo;
	} else {
		icols = irows;
	}
	itau = 1;
	iwrk = itau + irows;
	RNP::TLASupport::QRFactorization(irows, icols, &b[ilo + ilo * ldb], ldb, &work[itau], &work[iwrk]);

	// Apply the orthogonal transformation to matrix A
	// (Complex Workspace: need N)
	RNP::TLASupport::ApplyOrthognalMatrixFromElementaryReflectors<'L','C'>(
		irows, icols, irows, &b[ilo + ilo * ldb], ldb, &
		work[itau], &a[ilo + ilo * lda], lda, &work[iwrk]);

	// Initialize VL
	// (Complex Workspace: need N)
	if (ilvl) {
		RNP::TBLAS::SetMatrix<'F'>(n, n, std::complex<double>(0), std::complex<double>(1), &vl[vl_offset], ldvl);
		if (irows > 1) {
			RNP::TBLAS::CopyMatrix<'L'>(irows-1, irows-1, &b[ilo + 1 + ilo * ldb], ldb, &vl[ilo + 1 + ilo * ldvl], ldvl);
		}
		RNP::TLASupport::GenerateOrthognalMatrixFromElementaryReflectors(irows, irows, irows, &vl[ilo + ilo * ldvl], ldvl, &work[
				itau], &work[iwrk]);
	}

	// Initialize VR
	if (ilvr) {
		RNP::TBLAS::SetMatrix<'F'>(n, n, std::complex<double>(0), std::complex<double>(1), &vr[vr_offset], ldvr);
	}

	// Reduce to generalized Hessenberg form
	if(ilvl){
		if(ilvr){
			RNP::TLASupport::GeneralizedHessenbergReduction<'V','V'>(n, ilo-1, ihi-1, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr);
		}else{
			RNP::TLASupport::GeneralizedHessenbergReduction<'V','N'>(n, ilo-1, ihi-1, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr);
		}
	}else{
		if(ilvr){
			RNP::TLASupport::GeneralizedHessenbergReduction<'N','V'>(n, ilo-1, ihi-1, &a[a_offset], lda, &b[b_offset], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr);
		}else{
			RNP::TLASupport::GeneralizedHessenbergReduction<'N','N'>(irows, 0, irows-1, &a[ilo + ilo * lda], lda, &b[ilo + ilo * ldb], ldb, &vl[vl_offset], ldvl, &vr[vr_offset], ldvr);
		}
	}

	// Perform QZ algorithm (Compute eigenvalues, and optionally, the
	// Schur form and Schur vectors)
	// (Complex Workspace: need N)
	// (Real Workspace: need N)

	iwrk = itau;
	char chtemp[1];
	if (ilv) {
		chtemp[0] = 'S';
	} else {
		chtemp[0] = 'E';
	}
	char jobvl[1], jobvr[1];
	if(ilvl){
		jobvl[0] = 'V';
	}else{
		jobvl[0] = 'N';
	}
	if(ilvr){
		jobvr[0] = 'V';
	}else{
		jobvr[0] = 'N';
	}
	int ierr = zhgeqz_(chtemp, jobvl, jobvr, n, ilo, ihi, &a[a_offset], lda, &b[
			b_offset], ldb, &alpha[1], &beta[1], &vl[vl_offset], ldvl, &vr[
			vr_offset], ldvr, &work[iwrk], 2*n + 1 - iwrk, &rwork[irwrk]);
	if (ierr != 0) {
		if (ierr > 0 && ierr <= (int)n) {
			return ierr;
		} else if (ierr > (int)n && ierr <= (int)n << 1) {
			return ierr - n;
		} else {
			return n + 1;
		}
		return ierr;
	}

	// Compute Eigenvectors
	// (Real Workspace: need 2*N)
	// (Complex Workspace: need 2*N)

	if (ilv) {
		size_t in;
		ierr = ztgevc_("B", ldumma, n, &a[a_offset], lda, &b[b_offset], ldb, 
				&vl[vl_offset], ldvl, &vr[vr_offset], ldvr, n, &in, &work[
				iwrk], &rwork[irwrk]);
		if (ierr != 0) {
			return n+2;
		}

		// Undo balancing on VL and VR and normalization
		// (Workspace: none needed)

		if (ilvl) {
			zggbak_("P", "L", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vl[vl_offset], ldvl);
			for (size_t jc = 1; jc <= n; ++jc) {
				temp = 0.;
				for (size_t jr = 1; jr <= n; ++jr) {
					double abs1vl = RNP::TBLAS::_RealOrComplexChooser<std::complex<double> >::_abs1(vl[jr+jc*ldvl]);
					if(abs1vl > temp){ temp = abs1vl; }
				}
				if (temp >= smlnum) {
					temp = 1. / temp;
					for (size_t jr = 1; jr <= n; ++jr) {
						vl[jr+jc*ldvl] *= temp;
					}
				}
			}
		}
		if (ilvr) {
			zggbak_("P", "R", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vr[vr_offset], ldvr);
			for (size_t jc = 1; jc <= n; ++jc) {
				temp = 0.;
				for (size_t jr = 1; jr <= n; ++jr) {
					double abs1vr = RNP::TBLAS::_RealOrComplexChooser<std::complex<double> >::_abs1(vr[jr+jc*ldvr]);
					if(abs1vr > temp){ temp = abs1vr; }
				}
				if (temp >= smlnum) {
					temp = 1. / temp;
					for (size_t jr = 1; jr <= n; ++jr) {
						vr[jr+jc*ldvr] *= temp;
					}
				}
			}
		}
	}

	// Undo scaling if necessary
	if (ilascl) {
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, anrmto, anrm, n, 1, &alpha[1], n);
	}
	if (ilbscl) {
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, bnrmto, bnrm, n, 1, &beta[1], n);
	}

	return 0;
}







// Computes for a pair of N-by-N complex nonsymmetric matrices (A,B),
// the generalized eigenvalues, the generalized complex Schur form
// (S, T), and optionally left and/or right Schur vectors (VSL and
// VSR). This gives the generalized Schur factorization
//         (A,B) = ( (VSL)*S*(VSR)^H, (VSL)*T*(VSR)^H )
// where (VSR)^H is the conjugate-transpose of VSR.

// A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
// or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
// usually represented as the pair (alpha,beta), as there is a
// reasonable interpretation for beta=0, and even for both being zero.

// A pair of matrices (S,T) is in generalized complex Schur form if S
// and T are upper triangular and, in addition, the diagonal elements
// of T are non-negative real numbers.

// Arguments
// =========
// N       The order of the matrices A, B, VSL, and VSR.  N >= 0.
// A       (input/output) COMPLEX*16 array, dimension (LDA, N)
//         On entry, the first of the pair of matrices.
//         On exit, A has been overwritten by its generalized Schur
//         form S.
// LDA     The leading dimension of A.  LDA >= max(1,N).
// B       On entry, the second of the pair of matrices.
//         On exit, B has been overwritten by its generalized Schur
//         form T.
// LDB     The leading dimension of B.  LDB >= max(1,N).

// ALPHA   (output) COMPLEX*16 array, dimension (N)
// BETA    (output) COMPLEX*16 array, dimension (N)
//         On exit,  ALPHA(j)/BETA(j), j=1,...,N, will be the
//         generalized eigenvalues.  ALPHA(j), j=1,...,N  and  BETA(j),
//         j=1,...,N  are the diagonals of the complex Schur form (A,B)
//         of the output. The  BETA(j) will be non-negative real.

//         Note: the quotients ALPHA(j)/BETA(j) may easily over- or
//         underflow, and BETA(j) may even be zero.  Thus, the user
//         should avoid naively computing the ratio alpha/beta.
//         However, ALPHA will be always less than and usually
//         comparable with norm(A) in magnitude, and BETA always less
//         than and usually comparable with norm(B).

// VSL     (output) COMPLEX*16 array, dimension (LDVSL,N)
//         If VSL != NULL, VSL will contain the left Schur vectors.
//         Ignored if NULL
// LDVSL   The leading dimension of the matrix VSL. LDVSL >= 1, and
//         if VSL != NULL, LDVSL >= N.

// VSR     (output) COMPLEX*16 array, dimension (LDVSR,N)
//         If VSR != NULL, VSR will contain the right Schur vectors.
// LDVSR   The leading dimension of the matrix VSR. LDVSR >= 1, and
//         if VSR != NULL, LDVSR >= N.

// WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,2*N))
// RWORK   (workspace) DOUBLE PRECISION array, dimension (8*N)

// return: = 0:  successful exit
//         < 0:  if INFO = -i, the i-th argument had an illegal value.
//         =1,...,N:
//               The QZ iteration failed.  (A,B) are not in Schur
//               form, but ALPHA(j) and BETA(j) should be correct for
//               j=INFO+1,...,N.
//         > N:  =N+1: other than QZ iteration failed in ZHGEQZ
int RNP::GeneralizedSchurDecomposition(size_t n,
	std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb,
	std::complex<double> *alpha, std::complex<double> *beta,
	std::complex<double> *vsl, size_t ldvsl, std::complex<double> *vsr, size_t ldvsr,
	std::complex<double> *work, double *rwork)
{
	using namespace std;
    /* System generated locals */
    size_t a_offset, b_offset, vsl_offset, vsr_offset;

    /* Local variables */
    int ihi, ilo;
    double anrm, bnrm;
    int itau, iwrk;
    int ileft, icols;
    int irwrk, irows;
    
    bool ilascl, ilbscl;
    int iright;
   
    double anrmto;
    double bnrmto;

	const bool ilvsl = (NULL != vsl);
	const bool ilvsr = (NULL != vsr);

    /* Parameter adjustments */
    a_offset = 1 + lda;
    a -= a_offset;
    b_offset = 1 + ldb;
    b -= b_offset;
    --alpha;
    --beta;
    vsl_offset = 1 + ldvsl;
    vsl -= vsl_offset;
    vsr_offset = 1 + ldvsr;
    vsr -= vsr_offset;
    --work;
    --rwork;

    if (n < 0) {
		return -1;
    } else if (lda < max((size_t)1,n)) {
		return -3;
    } else if (ldb < max((size_t)1,n)) {
		return -5;
    } else if (ldvsl < 1 || ilvsl && ldvsl < n) {
		return -9;
    } else if (ldvsr < 1 || ilvsr && ldvsr < n) {
		return -11;
    }

    if (n == 0) {
		return 0;
    }

    const double eps = 2*std::numeric_limits<double>::epsilon();
    const double smlnum = sqrt(std::numeric_limits<double>::min()) / eps;
    const double bignum = 1. / smlnum;

	// Scale A if max element outside range [SMLNUM,BIGNUM]
	RNP::TLASupport::CheapMatrixNorm<'M'>(n, n, &a[a_offset], lda, &anrm);
    ilascl = false;
    if (anrm > 0. && anrm < smlnum) {
		anrmto = smlnum;
		ilascl = true;
    } else if (anrm > bignum) {
		anrmto = bignum;
		ilascl = true;
    }

    if (ilascl) {
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, anrm, anrmto, n, n, &a[a_offset], lda);
    }

	// Scale B if max element outside range [SMLNUM,BIGNUM]
	RNP::TLASupport::CheapMatrixNorm<'M'>(n, n, &b[b_offset], ldb, &bnrm);
    ilbscl = false;
    if (bnrm > 0. && bnrm < smlnum) {
		bnrmto = smlnum;
		ilbscl = true;
    } else if (bnrm > bignum) {
		bnrmto = bignum;
		ilbscl = true;
    }

    if (ilbscl) {
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, bnrm, bnrmto, n, n, &b[b_offset], ldb);
    }

	// Permute the matrix to make it more nearly triangular
	// (Real Workspace: need 6*N)
    ileft = 1;
    iright = n + 1;
    irwrk = iright + n;
    zggbal_("P", n, &a[a_offset], lda, &b[b_offset], ldb, &ilo, &ihi, &rwork[ileft], &rwork[iright], &rwork[irwrk]);

	// Reduce B to triangular form (QR decomposition of B)
	// (Complex Workspace: need N, prefer N*NB)

    irows = ihi + 1 - ilo;
    icols = n + 1 - ilo;
    itau = 1;
    iwrk = itau + irows;
	RNP::TLASupport::QRFactorization(irows, icols, &b[ilo + ilo * ldb], ldb, &work[itau], &work[iwrk]);

	// Apply the orthogonal transformation to matrix A
	// (Complex Workspace: need N, prefer N*NB)

	RNP::TLASupport::ApplyOrthognalMatrixFromElementaryReflectors<'L','C'>(
		irows, icols, irows, &b[ilo + ilo * ldb], ldb, &
		work[itau], &a[ilo + ilo * lda], lda, &work[iwrk]);

	// Initialize VSL
	// (Complex Workspace: need N, prefer N*NB)
    if (ilvsl) {
		RNP::TBLAS::SetMatrix<'F'>(n, n, std::complex<double>(0), std::complex<double>(1), &vsl[vsl_offset], ldvsl);
		if (irows > 1) {
			RNP::TBLAS::CopyMatrix<'L'>(irows-1, irows-1, &b[ilo + 1 + ilo * ldb], ldb, &vsl[ilo + 1 + ilo * ldvsl], ldvsl);
		}
		RNP::TLASupport::GenerateOrthognalMatrixFromElementaryReflectors(irows, irows, irows, &vsl[ilo + ilo * ldvsl], ldvsl, &work[itau], &work[iwrk]);
    }

    if (ilvsr) { // Initialize VSR
		RNP::TBLAS::SetMatrix<'F'>(n, n, std::complex<double>(0), std::complex<double>(1), &vsr[vsr_offset], ldvsr);
    }

	// Reduce to generalized Hessenberg form
	// (Workspace: none needed)
	if(ilvsl){
		if(ilvsr){
			RNP::TLASupport::GeneralizedHessenbergReduction<'V','V'>(n, ilo-1, ihi-1, &a[a_offset], lda, &b[b_offset], ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr);
		}else{
			RNP::TLASupport::GeneralizedHessenbergReduction<'V','N'>(n, ilo-1, ihi-1, &a[a_offset], lda, &b[b_offset], ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr);
		}
	}else{
		if(ilvsr){
			RNP::TLASupport::GeneralizedHessenbergReduction<'N','V'>(n, ilo-1, ihi-1, &a[a_offset], lda, &b[b_offset], ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr);
		}else{
			RNP::TLASupport::GeneralizedHessenbergReduction<'N','N'>(n, ilo-1, ihi-1, &a[a_offset], lda, &b[b_offset], ldb, &vsl[vsl_offset], ldvsl, &vsr[vsr_offset], ldvsr);
		}
	}

	// Perform QZ algorithm, computing Schur vectors if desired
	// (Complex Workspace: need N)
	// (Real Workspace: need N)
    iwrk = itau;
	char jobvl[1], jobvr[1];
	if(ilvsl){
		jobvl[0] = 'V';
	}else{
		jobvl[0] = 'N';
	}
	if(ilvsr){
		jobvr[0] = 'V';
	}else{
		jobvr[0] = 'N';
	}
	int ierr = zhgeqz_("S", jobvl, jobvr, n, ilo, ihi, &a[a_offset], lda, &b[
			b_offset], ldb, &alpha[1], &beta[1], &vsl[vsl_offset], ldvsl, &vsr[
			vsr_offset], ldvsr, &work[iwrk], 2*n + 1 - iwrk, &rwork[irwrk]);

    if (ierr != 0) {
		if (ierr > 0 && ierr <= (int)n) {
			return ierr;
		} else if (ierr > (int)n && ierr <= 2*(int)n) {
			return ierr - n;
		} else {
			return n + 1;
		}
		return ierr;
    }
    
	// Apply back-permutation to VSL and VSR
	// (Workspace: none needed)
    if (ilvsl) {
		zggbak_("P", "L", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vsl[vsl_offset], ldvsl);
    }
    if (ilvsr) {
		zggbak_("P", "R", n, ilo, ihi, &rwork[ileft], &rwork[iright], n, &vsr[vsr_offset], ldvsr);
    }

	// Undo scaling
    if (ilascl) {
		RNP::TLASupport::RescaleMatrix<'U'>(0, 0, anrmto, anrm, n, n, &a[a_offset], lda);
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, anrmto, anrm, n, 1, &alpha[1], n);
    }

    if (ilbscl) {
		RNP::TLASupport::RescaleMatrix<'U'>(0, 0, bnrmto, bnrm, n, n, &b[b_offset], ldb);
		RNP::TLASupport::RescaleMatrix<'G'>(0, 0, bnrmto, bnrm, n, 1, &beta[1], n);
    }

    return 0;

}
