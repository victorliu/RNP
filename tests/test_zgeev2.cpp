#include <cstdlib>
#include <cmath>
#include <complex>
#include <limits>
#include <iostream>
#include <ctime>
#include "TBLAS.h"
#include "TLASupport.h"
#include "Eigensystems.h"

double frand(){
	return (double)rand()/(double)RAND_MAX;
}
inline std::complex<double> crand(){
	double r = frand();
	double i = frand();
	return std::complex<double>(r,i);
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
	size_t n = 100;
	std::complex<double> *A = new std::complex<double>[n*n];
	std::complex<double> *Acopy = new std::complex<double>[n*n];
	std::complex<double> *w = new std::complex<double>[n];
	std::complex<double> *U = new std::complex<double>[n*n];
	std::complex<double> *V = new std::complex<double>[n*n];
	
	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			A[i+j*n] = crand();
		}
	}
	
	//makeid(n,n,A);
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
		}
	}
	
	size_t lwork = 2*n;
	std::complex<double> *work = new std::complex<double>[lwork];
	size_t lrwork = 8*n;
	double *rwork = new double[lrwork];
	int info;
	info = RNP::Eigensystem(n,A,n,w,U,n,V,n, work,rwork);
	delete [] work;
	delete [] rwork;
	
	//std::cout << "alpha = "; print(n,alpha); std::cout << std::endl;
	//std::cout << "beta  = "; print(n,beta ); std::cout << std::endl;
	
	std::cout << "Info = " << info << std::endl;
	
	// Verify
	std::complex<double> *temp = new std::complex<double>[n*n];
	double maxerr;
	
	
	// Verify A^H*U = conj(Lambda)*U
	
	for(size_t i = 0; i < n; ++i){
		RNP::TBLAS::MultMV<'C'>(n,n, 1.,Acopy,n, &U[0+i*n],1, 0.,&temp[0+i*n], 1);
		RNP::TBLAS::Axpy(n, -std::conj(w[i]), &U[0+i*n],1, &temp[0+i*n], 1);
	}
	
	RNP::TLASupport::CheapMatrixNorm<'M'>(n,n,temp,n, &maxerr);
	std::cout << "Max err = " << maxerr << std::endl;
	
	// Verify A*V = V*Lambda
	RNP::TBLAS::MultMM<'N','N'>(n,n,n, 1.,Acopy,n, V,n, 0.,temp,n);
	for(size_t i = 0; i < n; ++i){
		RNP::TBLAS::Axpy(n, -w[i], &V[0+i*n], 1, &temp[0+i*n], 1);
	}
	RNP::TLASupport::CheapMatrixNorm<'M'>(n,n,temp,n, &maxerr);
	std::cout << "Max err = " << maxerr << std::endl;
	
	delete [] temp;
	
	/*
	for(size_t i = 0; i < n; ++i){
		std::cout << "w[" << i << "] = " << w[i] << std::endl;
	}
	*/
	
	//print(n,evals); std::cout << std::endl;
	//print(n,n,evecs); std::cout << std::endl;
	
	delete [] A;
	delete [] Acopy;
	delete [] w;
	delete [] U;
	delete [] V;
	return 0;
}
