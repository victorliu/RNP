#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <ctime>
#include "GeneralizedEigensystems.h"

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
void print(size_t n, std::complex<double> *V){
	std::cout << "{";
	for(size_t i = 0; i < n; ++i){
		std::cout << V[i];
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
int main(){
	srand(time(0));
	size_t n = 400;
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
	info = RNP::GeneralizedEigensystem(n,A,n,B,n,alpha,beta,U,n,V,n, work,rwork);
	delete [] work;
	delete [] rwork;
	
	//std::cout << "alpha = "; print(n,alpha); std::cout << std::endl;
	//std::cout << "beta  = "; print(n,beta ); std::cout << std::endl;
	
	std::cout << "Info = " << info << std::endl;
	
	// Verify
	std::complex<double> *temp = new std::complex<double>[n*n];
	double maxerr = 0;
	
	for(size_t i = 0; i < n; ++i){
		// Compute beta*A*V - alpha*B*V
		for(size_t j = 0; j < n; ++j){
			std::complex<double> asum(0), bsum(0);
			for(size_t k = 0; k < n; ++k){
				asum += Acopy[i+k*n]*V[k+j*n];
				bsum += Bcopy[i+k*n]*V[k+j*n];
			}
			temp[i+j*n] = beta[j]*asum-alpha[j]*bsum;
			double curerr = std::abs(temp[i+j*n]);
			if(curerr > maxerr){ maxerr = curerr; }
		}
	}
	//print(n,n,temp); std::cout << std::endl;
	std::cout << "Max err = " << maxerr << std::endl;
	
	maxerr = 0;
	for(size_t i = 0; i < n; ++i){
		// Compute beta*U^H*A - alpha*U^H*B
		for(size_t j = 0; j < n; ++j){
			std::complex<double> asum(0), bsum(0);
			for(size_t k = 0; k < n; ++k){
				asum += Acopy[k+j*n]*std::conj(U[k+i*n]);
				bsum += Bcopy[k+j*n]*std::conj(U[k+i*n]);
			}
			temp[i+j*n] = beta[i]*asum-alpha[i]*bsum;
			double curerr = std::abs(temp[i+j*n]);
			if(curerr > maxerr){ maxerr = curerr; }
		}
	}
	//print(n,n,temp); std::cout << std::endl;
	std::cout << "Max err = " << maxerr << std::endl;
	
	delete [] temp;
	
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
