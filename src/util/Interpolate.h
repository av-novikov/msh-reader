#ifndef INTERPOLATE_H_
#define INTERPOLATE_H_

#include <math.h>
#include "adolc/adouble.h"

class Interpolate
{
	public:
	Interpolate();
	Interpolate(const double *ptx, const double *pty, const size_t N);
	Interpolate(const double *ptx, const double *pty, const double *dpty, const size_t N);
	Interpolate(const double *ptx, const double *pty, const double *dpty,const double *d2pty, const size_t N);
	~Interpolate();
	
	double Solve(const double arg);
	adouble Solve(adouble arg);
	double DSolve(const double arg);
	double D2Solve(const double arg);
	
	double *x;
	double *y;
	double *dy;
	double *d2y;
	int *Flag;
	int N;
	int Nf;
	double xmin, xmax, temp;
	bool IsSetDiff;
	bool IsSetDiff2;
};

#endif /* INTERPOLATE_H_ */