#ifndef OIL_PROPERTIES_HPP_
#define OIL_PROPERTIES_HPP_

#include <vector>
#include <utility>

#include "adolc/adouble.h"
#include "adolc/taping.h"

namespace oil
{
	struct Skeleton_Props
	{
		double m;
		double beta;
		double dens_stc;

		double kx, ky, kz;

		double height;
		double p_init;
		double p_out;

		inline adouble getPoro(adouble p) const
		{
			return (adouble)(m)* ((adouble)(1.0) + (adouble)(beta)* (p - p_init));
		};
		inline adouble getDensity(adouble p) const
		{
			return (adouble)(dens_stc);
		};
	};
	struct Oil_Props
	{
		double visc;
		double dens_stc;
		double beta;

		double p_ref;
		inline adouble getB(adouble p) const
		{
			return exp(-(adouble)beta * (p - p_ref));
		};
		inline adouble getDensity(adouble p) const
		{
			return dens_stc / getB(p);
		};
		inline adouble getViscosity(const adouble p) const
		{
			return (adouble)(visc);
		};
	};
	struct Properties
	{
		std::vector<double> timePeriods;
		std::vector<double> rates;
		std::vector<double> pwf;

		bool leftBoundIsRate;
		bool rightBoundIsPres;

		std::vector<std::pair<int, int> > perfIntervals;
		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;

		double ht;
		double ht_min;
		double ht_max;

		double alpha;

		double r_w;
		double r_e;

		double R_dim;
	};
};

#endif /* OIL_PROPERTIES_HPP_ */
