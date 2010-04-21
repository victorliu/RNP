#include "RootFinding.h"
#include <limits>
#include <cstdlib>
#include <iostream>

// from http://www.google.com/codesearch/p?hl=en#4FSOSMZ6Pxc/distfiles/camfr-1.2.tgz|f5UKxYPUedM/camfr_1.2/camfr/math/calculus/croot/mueller.cpp&q=mueller%20method%20complex
// Deflation consists of canceling out a previously converged zero.
// This deflation strategy assumes simple zeros. For higher order zeros,
// this will not perform a deflation.
// In case of duplicate zeros in the output, one could post process by deflating within clusters.
void RNP::RootFinding::Mueller::FindRoots(RNP::RootFinding::Mueller::ComplexFunction f, size_t npts, RNP::RootFinding::Mueller::ComplexRoot *pt, RNP::RootFinding::Mueller::Params *params, void *fdata){
	typedef RNP::RootFinding::Mueller::complex_t complex_t;
	if(npts < 1){ return; }
	
	const double eps = (0 == params->tolerance) ? std::numeric_limits<double>::epsilon() : params->tolerance;
	const double root_eps = sqrt(eps);
	const size_t max_eval = (0 == params->max_evaluations) ? 1000 : params->max_evaluations;
	
	// find extents of initial points
	double min_re = std::numeric_limits<double>::max();
	double max_re = -std::numeric_limits<double>::max();
	double min_im = std::numeric_limits<double>::max();
	double max_im = -std::numeric_limits<double>::max();
	for(size_t i = 0; i < npts; ++i){
		if(pt[i].z.real() > max_re){ max_re = pt[i].z.real(); }
		if(pt[i].z.real() < min_re){ min_re = pt[i].z.real(); }
		if(pt[i].z.imag() > max_im){ max_im = pt[i].z.imag(); }
		if(pt[i].z.imag() < min_im){ min_im = pt[i].z.imag(); }
		pt[i].converged = false;
	}
	double dim_re = max_re - min_re;
	double dim_im = max_re - min_im;
	double dim_max = dim_re; if(dim_im > dim_max){ dim_max = dim_im; }
	double characteristic_length = 0.1*(1./(double)npts)*dim_max;
	if(npts < 2){
		characteristic_length = eps;
	}
	
	// initialize starting directions
	for(size_t i = 0; i < npts; ++i){
		if(0 != std::abs(pt[i].dz)){ continue; }
		pt[i].dz = std::complex<double>(
			(double)rand()/(double)RAND_MAX - 0.5,
			(double)rand()/(double)RAND_MAX - 0.5);
		pt[i].dz /= std::abs(pt[i].dz);
		pt[i].dz *= characteristic_length;
	}
	
	size_t max_evals_used = 0;
	
	// Perform Mueller's method on each starting point
	for(size_t i = 0; i < npts; ++i){
		complex_t z3 = pt[i].z;
		complex_t z2 = z3 + pt[i].dz;
		complex_t z1 = z2 + pt[i].dz;
		complex_t f3 = f(z3, fdata);
		complex_t f2 = f(z2, fdata);
		for(size_t j = 0; j < i; ++j){ // deflation
			if(!pt[j].converged){ continue; }
			f3 /= (z3 - pt[j].z);
			f2 /= (z2 - pt[j].z);
		}
		complex_t dz23 = (f2-f3) / (z2-z3);
		
		for(size_t iter = 2; iter < max_eval; ++iter){
			complex_t f1 = f(z1, fdata);
			for(size_t j = 0; j < i; ++j){ // deflation
				if(!pt[j].converged){ continue; }
				f1 /= (z1 - pt[j].z);
			}
			if(std::abs(f1) > 1e4*(std::abs(f2) + std::abs(f3))){
				break;
			}
			if(std::abs(z1-z2) < eps || std::abs(z1-z3) < eps){
				pt[i].converged = true;
				pt[i].z = z1;
				if(iter > max_evals_used){ max_evals_used = iter; }
				break;
			}
			// generate new estimate
			complex_t dz12 = (f1-f2) / (z1-z2);
			complex_t dz123 = (dz12 - dz23) / (z1 - z3);
			complex_t W = (z1 - z2)*dz123 + dz12;
			complex_t w = sqrt(W*W - 4.*f1*dz123);
			complex_t corr = (std::abs(W+w) > std::abs(W-w)) ? 2.*f1/(W+w) : 2.*f1/(W-w);
			if(params->limit_step_size){
				double az12 = std::abs(z1-z2);
				double az23 = std::abs(z2-z3);
				double acorr = std::abs(corr);
				double amax = az12;
				if(az23 > amax){ amax = az23; }
				amax *= 2.;
				if(acorr > amax){
					acorr = amax;
				}
				corr /= std::abs(corr);
				corr *= acorr;
			}
			z3 = z2;
			z2 = z1;
			z1 -= corr;
			f3 = f2;
			f2 = f1;
			dz23 = dz12;
			// convergence check
			double az1 = std::abs(z1);
			if(root_eps*std::abs(corr) > std::abs(z2) || !(az1 == az1)){
				break;
			}
			if(std::abs(corr) < eps*az1){
				pt[i].converged = true;
				pt[i].z = z1;
				if(iter > max_evals_used){ max_evals_used = iter; }
				break;
			}
		}
	}
	params->max_evaluations = max_evals_used;
}

