#ifndef OILSOLVER_HPP_
#define OILSOLVER_HPP_

#include "src/models/AbstractSolver.hpp"
#include "src/models/oil/Oil.hpp"
#include "src/util/ParalutionInterface.h"
#include <fstream>
#include <iomanip>

namespace oil
{
	class OilSolver : public AbstractSolver<Oil>
	{
	protected:
		void control();
		void solveStep();
		void writeData();

		std::ofstream plot_P, plot_Q;
		ParSolver solver;

		void computeJac();
		void fill();
		void copySolution(const paralution::LocalVector<double>& sol);
	public:
		OilSolver(Model* _model);
		~OilSolver();

		void start();
	};
};

#endif /* OILSOLVER_HPP_ */
