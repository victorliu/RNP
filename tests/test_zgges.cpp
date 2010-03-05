#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <ctime>
#include <GeneralizedEigensystems.h>
#include <IO.h>
#include <TBLAS.h>
#include <TLASupport.h>

int GeneralizedSchurDecomposition(size_t n, 
		std::complex<double> *a, size_t lda, std::complex<double> *b, size_t ldb, 
		std::complex<double> *alpha, std::complex<double> *beta, std::complex<double> *vl, size_t 
		ldvl, std::complex<double> *vr, size_t ldvr, std::complex<double> *work, double *rwork);

double frand(){
	return (double)rand()/(double)RAND_MAX;
}
inline std::complex<double> crand(){
	return std::complex<double>(frand(),frand());
}

void print(size_t rows, size_t lda, std::complex<double> *A){
	std::cout << "{";
	for(size_t i = 0; i < rows; ++i){
		std::cout << "{";
		for(size_t j = 0; j < rows; ++j){
			std::cout << A[i+j*lda];
			if(j != rows-1){
				std::cout << ", ";
			}
		}
		std::cout << "}";
		if(i != rows-1){
			std::cout << ",";
		}
		std::cout << std::endl;
	}
	std::cout << "}";
}
void print(size_t n, std::complex<double> *V, size_t incv = 1){
	std::cout << "{";
	for(size_t i = 0; i < n; ++i){
		std::cout << V[i*incv];
		if(i != n-1){
			std::cout << ", ";
		}
	}
	std::cout << "}";
}

void makeid(size_t rows, size_t cols, std::complex<double> *A){
	for(size_t i = 0; i < rows; ++i){
		for(size_t j = 0; j < cols; ++j){
			if(i == j){
				A[i+j*rows] = 1;
			}else{
				A[i+j*rows] = 0;
			}
		}
	}
}

double CheckGeneralizedSchurDecomposition(size_t n,
	const std::complex<double> *a, size_t lda,
	const std::complex<double> *b, size_t ldb,
	const std::complex<double> *s, size_t lds,
	const std::complex<double> *t, size_t ldt,
	const std::complex<double> *u, size_t ldu,
	const std::complex<double> *v, size_t ldv)
{
	std::complex<double> *p = new std::complex<double>[n*n];
	std::complex<double> *q = new std::complex<double>[n*n];
	double maxval1, maxval2;
	
	RNP::TBLAS::MultMM<'N','N'>(n,n,n,1.,u,n,s,lds,0.,p,n);
	RNP::TBLAS::MultMM<'N','C'>(n,n,n,1.,p,n,v,ldv,0.,q,n);
	RNP::TBLAS::Axpy(n*n,-1.,a,1,q,1);
	RNP::TLASupport::CheapMatrixNorm<'M'>(n,n,q,n,&maxval1);
	
	RNP::TBLAS::MultMM<'N','N'>(n,n,n,1.,u,n,t,ldt,0.,p,n);
	RNP::TBLAS::MultMM<'N','C'>(n,n,n,1.,p,n,v,ldv,0.,q,n);
	RNP::TBLAS::Axpy(n*n,-1.,b,1,q,1);
	RNP::TLASupport::CheapMatrixNorm<'M'>(n,n,q,n,&maxval2);

	delete [] q;
	delete [] p;

	return (maxval1 > maxval2 ? maxval1 : maxval2);
}

int main(){
	time_t seed = time(0);
	//std::cout << "seed = " << seed << std::endl;
	srand((unsigned int)seed);
	size_t n = 4;
	std::complex<double> *A = new std::complex<double>[n*n];
	std::complex<double> *B = new std::complex<double>[n*n];
	std::complex<double> *Acopy = new std::complex<double>[n*n];
	std::complex<double> *Bcopy = new std::complex<double>[n*n];
	std::complex<double> *alpha = new std::complex<double>[n];
	std::complex<double> *beta  = new std::complex<double>[n];
	std::complex<double> *U = new std::complex<double>[n*n];
	std::complex<double> *V = new std::complex<double>[n*n];
	
	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			A[i+j*n] = crand();
		}
	}
	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			B[i+j*n] = crand();
		}
	}
	
	//makeid(n,n,A);
	//makeid(n,n,B);
	/*
	A[1+0*n] = frand();
	A[2+1*n] = frand();
	A[2+0*n] = frand();
	makesym(n,n,A);
	*/
	
	//print(n,n,A); std::cout << std::endl;
	
	
	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			Acopy[i+j*n] = A[i+j*n];
			Bcopy[i+j*n] = B[i+j*n];
		}
	}
	
	size_t lwork = 2*n;
	std::complex<double> *work = new std::complex<double>[lwork];
	size_t lrwork = 8*n;
	double *rwork = new double[lrwork];
	int info;
	info = RNP::GeneralizedSchurDecomposition(n,A,n,B,n,alpha,beta,U,n,V,n, work,rwork);
	delete [] work;
	delete [] rwork;
	
	std::cout << "A:" << std::endl;
	RNP::IO::PrintMatrix(n,n,A,n) << std::endl;
	std::cout << "B:" << std::endl;
	RNP::IO::PrintMatrix(n,n,B,n) << std::endl;
	std::cout << "U:" << std::endl;
	RNP::IO::PrintMatrix(n,n,U,n) << std::endl;
	std::cout << "V:" << std::endl;
	RNP::IO::PrintMatrix(n,n,V,n) << std::endl;
	std::cout << "alpha:" << std::endl;
	RNP::IO::PrintVector(n,alpha) << std::endl;
	std::cout << "beta:" << std::endl;
	RNP::IO::PrintVector(n,beta) << std::endl;
	
	std::cout << "Info = " << info << std::endl;

	std::cout << "max err = " << CheckGeneralizedSchurDecomposition(n, Acopy, n, Bcopy, n, A, n, B, n, U, n, V, n) << std::endl;
/*
	QZSort(n, A,n, B,n, U,n, V,n, &QZSorter_MaxReal, NULL);
	std::cout << "max err = " << CheckGeneralizedSchurDecomposition(n, Acopy, n, Bcopy, n, A, n, B, n, U, n, V, n) << std::endl;
	*/
	//print(n,evals); std::cout << std::endl;
	//print(n,n,evecs); std::cout << std::endl;
	
	delete [] A;
	delete [] B;
	delete [] Acopy;
	delete [] Bcopy;
	delete [] alpha;
	delete [] beta;
	delete [] U;
	delete [] V;
	return 0;
}
