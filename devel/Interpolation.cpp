#include "Interpolation.h"
#include <cmath>
#include <iostream>

struct RNP::Interpolation::CubicSpline::Spline{
	double *data; // length 3*n; 3 fields are: x, y, M
	size_t n;
	EndpointType endtype;
};

void RNP::Interpolation::CubicSpline::Generate(size_t npts, const Point2 *pt, Spline *spline, EndpointType type){
	if(NULL == spline){ return; }
	spline->n = npts;
	spline->endtype = type;
	spline->data = new double[3*npts];
	
	// The spline is a_i(x-x_i)^3 + b_i(x-x_i)^2 + c_i(x-x_i) + d_i
	// for x in [x_i, x_{i+1}]
	// Working through the continuity relations, we obtain
	// M_{i+2}(x_{i+2}-x_{i+1}) + 2M_{i+1}(x_{i+2}-x_i) + M_i(x_{i+1}-x_i)
	// = 6 [ (y_{i+2}-y_{i+1})/(x_{i+2}-x_{i+1}) - (y_{i+1}-y_i)/(x_{i+1}-x_i) ]
	// where M_i = 2b_i, and we can recover the other coefficients:
	// a_i = (M_{i+1}-M_i)/(6(x_{i+1}-x_i))
	// d_i = y_i
	// c_i = (y_{i+1}-y_i)/(x_{i+1}-x_i) - (M_{i+1}-M_i) (x_{i+1}-x_i) / 6
	
	// we need to solve a symmetric positive definite tridiagonal system
	// for the M variables.
	double *Adiag = spline->data; // length n-2
	double *Adiag1 = &spline->data[npts]; // length n-3
	double *rhs = &spline->data[2*npts+1]; // length n-2
	
	// fill in RHS
	for(size_t i = 0; i < npts-2; ++i){
		rhs[i] = 6*( (pt[i+2].y-pt[i+1].y)/(pt[i+2].x-pt[i+1].x) - (pt[i+1].y-pt[i+0].y)/(pt[i+1].x-pt[i+0].x) );
		std::cout << "rhs[" << i << "] = " << rhs[i] << std::endl;
	}
	// fill in matrix assuming linear endpoints
	for(size_t i = 0; i < npts-2; ++i){
		Adiag[i] = 2*(pt[i+2].x - pt[i+0].x);
		Adiag1[i]= pt[i+2].x - pt[i+1].x;
		std::cout << "Adiag[" << i << "] = " << Adiag[i] << std::endl;
		std::cout << "Adiag1[" << i << "] = " << Adiag1[i] << std::endl;
	}
	switch(type){
	case LINEAR:
		break;
	case QUADRATIC:
		Adiag[0] = 5;
		Adiag[npts-3] = 5;
		break;
	default:
		break;
	}
	// Cholesky factor the A matrix (generate the lower triangular)
	for(size_t i = 0; i < npts-2; ++i){ // dpbtf2
		double ajj = Adiag[i];
		if(ajj <= 0){
			// not positive definite
			break;
		}
		ajj = sqrt(ajj);
		Adiag[i] = ajj;
		if(i < npts-3){
			Adiag1[i] /= ajj;
			Adiag[i+1] -= Adiag1[i]*Adiag1[i];
		}
	}
	// Solve for M
	for(size_t i = 0; i < npts-2; ++i){ // dtbsv
		if(0 != rhs[i]){
			rhs[i] /= Adiag[i];
			if(i < npts-3){
				rhs[i+1] -= rhs[i]*Adiag1[i];
			}
		}
	}
	{size_t i = npts-2; // dtbsv
	while(i --> 0){
		rhs[i] /= Adiag[i];
		if(i > 0){
			rhs[i-1] -= Adiag1[i-1]*rhs[i];
		}
	}}
	for(size_t i = 0; i < npts-2; ++i){
		std::cout << "M[" << i << "] = " << rhs[i] << std::endl;
	}
	switch(type){
	case LINEAR:
		rhs[-1] = 0;
		rhs[npts-2] = 0;
		break;
	case QUADRATIC:
		rhs[-1] = rhs[0];
		rhs[npts-2] = rhs[npts-3];
		break;
	default:
		break;
	}
	for(size_t i = 0; i < npts; ++i){
		spline->data[0+3*i] = pt[i].x;
		spline->data[1+3*i] = pt[i].y;
		spline->data[2+3*i] = rhs[i-1];
		std::cout << "x[" << i << "] = " << spline->data[0+3*i];
		std::cout << "y[" << i << "] = " << spline->data[1+3*i];
		std::cout << "M[" << i << "] = " << spline->data[2+3*i] << std::endl;
	}
}

void RNP::Interpolation::CubicSpline::Destroy(Spline *spline){
	if(NULL == spline){ return; }
	if(NULL != spline->data){
		delete [] spline->data;
	}
}

void RNP::Interpolation::CubicSpline::Evaluate(const Spline *spline, double x, double *y, double *dy, double *ddy){
	if(NULL != y){ *y = 0; }
	if(NULL != dy){ *dy = 0; }
	if(NULL != ddy){ *ddy = 0; }
	if(NULL == spline){ return; }
	size_t i = 0;
	double dx = 0;
	double ai = 0, bi = 0, ci = 0, di = 0;
	if(x < spline->data[0+3*0]){
		i = 0;
		dx = x - spline->data[0+3*0];
		switch(spline->endtype){
		case LINEAR:
			{
				ci = (spline->data[1+3*(i+1)] - spline->data[1+3*(i+0)]) / (spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]) - (spline->data[2+3*(i+1)] +2*spline->data[2+3*(i+0)])/6.*(spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]);
				di = spline->data[1+3*(i+0)];
			}
			break;
		case QUADRATIC:
			{
				bi = spline->data[2+3*(i+0)] / 2;
				ci = (spline->data[1+3*(i+1)] - spline->data[1+3*(i+0)]) / (spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]) - (spline->data[2+3*(i+1)] +2*spline->data[2+3*(i+0)])/6.*(spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]);
				di = spline->data[1+3*(i+0)];
			}
			break;
		default:
			break;
		}
	}else if(x > spline->data[0+3*(spline->n-1)]){
		i = spline->n-2;
		dx = x - spline->data[0+3*(i+1)];
		switch(spline->endtype){
		case LINEAR:
			{
				ci = (spline->data[1+3*(i+1)] - spline->data[1+3*(i+0)]) / (spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]) + (2*spline->data[2+3*(i+1)] + spline->data[2+3*(i+0)])/6.*(spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]);
				di = spline->data[1+3*(i+1)];
			}
			break;
		case QUADRATIC:
			{
				dx = x - spline->data[0+3*(i+0)];
				bi = spline->data[2+3*(i+0)] / 2;
				ci = (spline->data[1+3*(i+1)] - spline->data[1+3*(i+0)]) / (spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]) - (spline->data[2+3*(i+1)] +2*spline->data[2+3*(i+0)])/6.*(spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]);
				di = spline->data[1+3*(i+0)];
			}
			break;
		default:
			break;
		}
	}else{
		// binary search to find the proper interval
		// for now, linear search
		for(i = 1; i < spline->n; ++i){
			if(x < spline->data[0+3*i]){
				--i;
				dx = x - spline->data[0+3*i];
				break;
			}
		}
		ai = (spline->data[2+3*(i+1)] - spline->data[2+3*(i+0)]) / (6*(spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]));
		bi = spline->data[2+3*(i+0)] / 2;
		ci = (spline->data[1+3*(i+1)] - spline->data[1+3*(i+0)]) / (spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]) - (spline->data[2+3*(i+1)] +2*spline->data[2+3*(i+0)])/6.*(spline->data[0+3*(i+1)] - spline->data[0+3*(i+0)]);
		di = spline->data[1+3*(i+0)];
	}
	if(NULL != y){
		*y = ai * dx*dx*dx + bi*dx*dx + ci*dx + di;
	}
	if(NULL != dy){
		*dy = 3*ai * dx*dx + 2*bi*dx + ci;
	}
	if(NULL != ddy){
		*ddy = 6*ai * dx + 2*bi;
	}
}
