#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

#include <valarray>
#include <array>
#include "adolc/adouble.h"

namespace var
{
	namespace containers
	{
		struct Var1phase
		{
			static const int size = 1;
			double& p;

			Var1phase(double* data) : p(data[0]) {};
			Var1phase(const double* data) : p(const_cast<double&>(data[0])) {};
		};
		struct TapeVar1Phase
		{
			static const int size = 1;
			adouble p;
		};
		struct AcidVar
		{
			static const int size = 5;
			double& m;
			double& p;
			double& s;
			double& xa;
			double& xw;

			AcidVar(double* data) : m(data[0]), p(data[1]), s(data[2]), xa(data[3]), xw(data[4]) {};
			AcidVar(const double* data) : m(const_cast<double&>(data[0])), p(const_cast<double&>(data[1])), 
											s(const_cast<double&>(data[2])), xa(const_cast<double&>(data[3])), 
											xw(const_cast<double&>(data[4])) {};
		};
		struct TapeAcidVar
		{
			static const int size = 5;
			adouble m;
			adouble p;
			adouble s;
			adouble xa;
			adouble xw;
		};
	};

	template <typename TVariable>
	struct BaseVarWrapper
	{
		TVariable u_prev, u_iter, u_next;
	};
	template <typename TVariable>
	struct BasicVariables
	{
		static const int size = TVariable::size;
		typedef BaseVarWrapper<TVariable> Wrap;
		std::valarray<double> u_prev, u_iter, u_next;

		Wrap operator[](const size_t idx)
		{
			return{ TVariable(&u_prev[idx * size]), TVariable(&u_iter[idx * size]), TVariable(&u_next[idx * size]) };
		};
		Wrap operator[](const size_t idx) const
		{
			return{ TVariable(&u_prev[idx * size]), TVariable(&u_iter[idx * size]), TVariable(&u_next[idx * size]) };
		};
	};
}

#endif /* VARIABLES_HPP_ */