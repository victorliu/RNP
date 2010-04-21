#include "Optimization.h"
#include "Random.h"
#define RNP_TBLAS_USE_RANDOM
#include "TBLAS.h"
#include "IO.h"
#include <complex>
#include <cstdlib>

typedef double real_t;
typedef std::complex<real_t> complex_t;

void generate_random_LP(
	size_t m, size_t n,
	real_t *A, real_t *b, real_t *c,
	real_t *bl, real_t *bu,
	real_t &d1, real_t &d2
){
	RNP::TBLAS::RandomVector<2>(m*n, A);
	RNP::TBLAS::RandomVector<1>(n, c); // use as dummy first
	RNP::TBLAS::MultMV<'N'>(m, n, 10., A, m, c, 1, 0., b, 1); // generate b to make problem feasible
	RNP::TBLAS::RandomVector<2>(n, c);
	RNP::TBLAS::Fill(n, 0., bl, 1);
	RNP::TBLAS::Fill(n, 10., bu, 1);
	d1 = 1e-3;
	d2 = 1e-3;
}

struct data_LP{
	real_t *A, *c;
};

void phi_LP(size_t n, const double *x, double *phi, double *dphi, double *ddphi, void *data_){
	data_LP *data = (data_LP*)data_;
	RNP::TBLAS::Fill(n, 0., dphi, 1);
	RNP::TBLAS::Fill(n, 0., ddphi, 1);
	*phi = RNP::TBLAS::Dot(n, data->c, 1, x, 1);
	std::cout << "phi = " << *phi << std::endl;
}
void A_LP(bool trans, size_t m, size_t n, double *x, double *y, void *data_){
	data_LP *data = (data_LP*)data_;
	if(trans){
		RNP::TBLAS::MultMV<'N'>(n, m, 1., data->A, m, y, 1, 0., x, 1);
	}else{
		RNP::TBLAS::MultMV<'N'>(m, n, 1., data->A, m, x, 1, 0., y, 1);
	}
}

int main(){
	const size_t n = 10;
	const size_t m = 20;
	const size_t lda = m;
	
	data_LP data;
	data.A = new real_t[n*m];
	real_t *b = new real_t[m];
	data.c = new real_t[n];
	real_t *bl = new real_t[n];
	real_t *bu = new real_t[n];
	real_t d1, d2;
	real_t *x = new real_t[n];
	real_t *y = new real_t[m];
	real_t *z = new real_t[n];
	
	generate_random_LP(m, n, data.A, b, data.c, bl, bu, d1, d2);
	
	RNP::Optimization::PDCO::Params params;
	RNP::Optimization::PDCO::Minimize(
		m, n,
		phi_LP, A_LP,
		b, bl, bu, &d1, &d2,
		x, y, z,
		&params,
		NULL,
		(void*)&data);
	
	delete [] data.A;
	delete [] b;
	delete [] data.c;
	delete [] bl;
	delete [] bu;
	delete [] x;
	delete [] y;
	delete [] z;
	return EXIT_SUCCESS;
}


