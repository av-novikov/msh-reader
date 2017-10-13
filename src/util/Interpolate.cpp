#include "src/util/Interpolate.h"
#include <iostream>

Interpolate::Interpolate()
{
}
Interpolate::~Interpolate()
{
}
Interpolate::Interpolate(const double *ptx, const double *pty, const size_t _N)
{
	Nf = 1000000;
	Flag = new int[Nf+1];
	N  = _N;
	x = new double[N];
	y = new double[N];

	for (size_t i = 0; i < N; i++)
	{
		x[i] = ptx[i];
		y[i] = pty[i];
	}
	
	xmin = x[0];
	xmax = x[N-1];
	int j;
	temp = (xmax - xmin)/Nf;
	
	for (int i = 0; i < N-1; i++)
		for (j = floor((x[i] - xmin)/temp); j <= floor((x[i+1] - xmin)/temp); j++)
			Flag[j] = i;
	
	IsSetDiff = false;

}
Interpolate::Interpolate(const double *ptx, const double *pty, const double *dpty, const size_t N)
{
	Nf = 10000;
	Flag = new int[Nf+1];
	x = new double[N];
	y = new double[N];
	dy = new double[N];
	for (int i= 0; i < N; i++){
		x[i] = ptx[i];
		y[i] = pty[i];
		dy[i] = dpty[i];
	}
	xmin = x[0];
	xmax = x[N-1];
	temp = (xmax - xmin)/Nf;
	for (int i = 0; i < N-1; i++){
		for (int j = floor((x[i] - xmin)/temp); j <= floor((x[i+1] - xmin)/temp); j++)
		Flag[j] = i;
	}
	IsSetDiff = true;

}
Interpolate::Interpolate(const double *ptx, const double *pty, const double *dpty, const double *d2pty, const size_t N)
{
	Nf = 10000;
	Flag = new int[Nf+1];
	x = new double[N];
	y = new double[N];
	dy = new double[N];
	d2y = new double[N];

	for (int i= 0; i < N; i++){
		x[i] = ptx[i];
		y[i] = pty[i];
		dy[i] = dpty[i];
		d2y[i]  = d2pty[i];
	}
	xmin = x[0];
	xmax = x[N-1];
	temp = (xmax - xmin)/Nf;
	for (int i = 0; i < N-1; i++){
		for (int j = floor((x[i] - xmin)/temp); j <= floor((x[i+1] - xmin)/temp); j++)
		Flag[j] = i;
	}
	IsSetDiff = true;
	IsSetDiff2 = true;
}
double Interpolate::Solve(double arg)
{
	if (arg < xmin ) 
		arg = xmin;
	else if(arg > xmax) 
		arg = xmax;

	int index = floor( (arg - xmin) / temp );
	double dx = x[Flag[index]+1] - x[Flag[index]];
	double dy = y[Flag[index]+1] - y[Flag[index]];
	
	return (arg - x[Flag[index]])/dx*dy + y[Flag[index]];
}
adouble Interpolate::Solve(adouble arg)
{
	double tmp = arg.value();
	if (tmp < xmin)
		tmp = xmin;
	else if (tmp > xmax)
		tmp = xmax;
		
	int i = Flag[int((tmp - xmin) / temp)];
	return (adouble)(DSolve(tmp)) * (arg - (adouble)(x[i])) + (adouble)y[i];
}
double Interpolate::DSolve(const double arg)
{
	double result;
	double dx;
	double dcy;
	
	if(arg < xmin)
		return 0;
	else if(arg > xmax)
		return 0;

	int i = floor((arg - xmin)/temp);
	int index = Flag[i];
	
	if (IsSetDiff)
	{
		 dx = x[index+1] - x[index];
		 dcy = dy[index+1] - dy[index];
		 result = (arg - x[index])/dx*dcy + dy[index];
	} else {
		dx = (x[index] - x[index + 1]);
		dcy = (y[index] - y[index + 1]);
		result = dcy/dx;
	}

	return result;
}
double Interpolate::D2Solve(const double arg)
{
	double result;
	double dx;
	double dcy;
	
	int i = floor((arg - xmin)/temp);
	if (i > Nf)
		i= Nf;
	if (i < 0)
		i = 0;
    int index = Flag[i];
	
	if (IsSetDiff2)
	{
		 dx = x[index+1] - x[index];
		 dcy = d2y[index+1] - d2y[index];
		 result = (arg - x[index])/dx*dcy + d2y[index];
	} else 
		result = 0;
	
	return result;
}