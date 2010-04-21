#include "Interpolation.h"
#include <iostream>

int main(){
	
	static const size_t n = 20;
	RNP::Interpolation::Point2 pt[n];
	pt[0].x = 0;
	pt[0].y = (double)rand()/(double)RAND_MAX;
	for(size_t i = 1; i < n; ++i){
		pt[i].x = pt[i-1].x + (double)rand()/(double)RAND_MAX + 0.1;
		pt[i].y = pt[i-1].y + (double)rand()/(double)RAND_MAX - 0.5;
	}
	
	/*
	static const size_t n = 4;
	RNP::Interpolation::Point2 pt[n];
	pt[0].x = 0;
	pt[0].y = 0;
	pt[1].x = 1;
	pt[1].y = 1;
	pt[2].x = 2;
	pt[2].y = 1;
	pt[3].x = 3;
	pt[3].y = 0;
	*/
	
	
	
	for(size_t i = 0; i < n; ++i){
		std::cerr << pt[i].x << '\t' << pt[i].y << std::endl;
	}
	RNP::Interpolation::CubicSpline::Spline s;
	RNP::Interpolation::CubicSpline::Generate(n, pt, &s, RNP::Interpolation::CubicSpline::QUADRATIC);
	for(double x = -1.0; x < pt[n-1].x+1.0; x += 0.05){
		double y, dy, ddy;
		RNP::Interpolation::CubicSpline::Evaluate(&s, x, &y, &dy, &ddy);
		std::cout << x << '\t' << y << '\t' << dy << '\t' << ddy << std::endl;
	}
	RNP::Interpolation::CubicSpline::Destroy(&s);
	return 0;
}

