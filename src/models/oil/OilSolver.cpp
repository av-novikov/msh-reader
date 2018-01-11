#include "src/models/oil/OilSolver.hpp"
#include <iostream>

#include "adolc/sparse/sparsedrivers.h"
#include "adolc/drivers/drivers.h"

using namespace oil;
using std::vector;
using std::ofstream;
using std::map;
using std::endl;
using std::setprecision;

OilSolver::OilSolver(Model* _model) : AbstractSolver<Model>(_model)
{
	const int strNum = var_size * model->cellsNum;

	y = new double[strNum];
	ind_i = new int[grid::stencil * var_size * strNum];
	ind_j = new int[grid::stencil * var_size * strNum];
	//cols = new int[strNum];
	a = new double[grid::stencil * var_size * strNum];
	ind_rhs = new int[strNum];
	rhs = new double[strNum];

	//options[0] = 0;          /* sparsity pattern by index domains (default) */
	//options[1] = 0;          /*                         safe mode (default) */
	//options[2] = 0;          /*              not required if options[0] = 0 */
	//options[3] = 0;          /*                column compression (default) */

	plot_P.open("snaps/P.dat", ofstream::out);
	plot_Q.open("snaps/Q.dat", ofstream::out);
};
OilSolver::~OilSolver()
{
	delete[] y;
	//for (int i = 0; i < Model::var_size * size; i++)
	//	delete[] jac[i];
	//delete[] jac;

	delete[] ind_i, ind_j, ind_rhs;
	//delete[] cols;
	delete[] a, rhs;

	plot_P.close();
	plot_Q.close();
};
void OilSolver::writeData()
{
	double p = 0.0, q = 0;

	plot_Q << cur_t * t_dim / 3600.0;

	map<size_t, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->u_next[it->first * Model::var_size] * model->P_dim;
		if (model->leftBoundIsRate)
			plot_Q << "\t" << it->second * model->Q_dim * 86400.0;
		else
		{
			plot_Q << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
			q += model->getRate(it->first);
		}
	}
	plot_P << cur_t * t_dim / 3600.0 << "\t" << p / (double)(model->Qcell.size()) << endl;

	if (model->leftBoundIsRate)
		plot_Q << "\t" << model->Q_sum * model->Q_dim * 86400.0 << endl;
	else
		plot_Q << "\t" << q * model->Q_dim * 86400.0 << endl;
}
void OilSolver::control()
{
	writeData();

	if (cur_t >= model->period[curTimePeriod])
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
	}

	if (model->ht <= model->ht_max && iterations < 6)
		model->ht = model->ht * 1.5;
	else if (iterations > 6 && model->ht > model->ht_min)
		model->ht = model->ht / 1.5;

	if (cur_t + model->ht > model->period[curTimePeriod])
		model->ht = model->period[curTimePeriod] - cur_t;

	cur_t += model->ht;
}
void OilSolver::start()
{
	int counter = 0;
	iterations = 8;

	fillIndices();
	solver.Init(Model::var_size * model->cellsNum, 1.e-15, 1.e-15);

	model->setPeriod(curTimePeriod);
	while (cur_t < Tt)
	{
		control();
		model->snapshot_all(counter++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << endl;
		cout << setprecision(6);
		cout << "time = " << cur_t << endl;
	}
	model->snapshot_all(counter++);
	writeData();
}
void OilSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPrev = averValue(0), aver, dAver = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 /*&& (dAverSat > 1.e-9 || dAverPres > 1.e-7)*/ && iterations < 20)
	{
		copyIterLayer();

		computeJac();
		fill();
		solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
 		solver.Solve(PRECOND::ILU_SERIOUS);
		copySolution(solver.getSolution());

		//if (repeat == 0)
		//	repeat = 1;

		err_newton = convergance(cellIdx, varIdx);
		aver = averValue(0);		dAver = fabs(aver - averPrev);		averPrev = aver;
		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}
void OilSolver::copySolution(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < size; i++)
	{
		auto& var = (*model)[i].u_next;
		var.p += sol[Model::var_size * i];
	}
}

void OilSolver::computeJac()
{
	trace_on(0);

	for (size_t i = 0; i < size; i++)
		model->x[i].p <<= model->u_next[i * var_size];

	for (int i = 0; i < mesh->inner_size; i++)
	{
		const auto& cell = mesh->cells[i];
		/*adouble isFracCell = (cell.type == elem::FRAC_HEX) ? true : false;
		adouble isBorderCell = (cell.type == elem::BORDER_HEX) ? true : false;
		adouble isInnerCell = (cell.type == elem::HEX || cell.type == elem::PRISM) ? true : false;
		condassign(model->h[i], isFracCell, model->solveFrac(cell));
		condassign(model->h[i], isBorderCell, model->solveBorder(cell));
		condassign(model->h[i], isInnerCell, model->solveInner(cell));*/
		if (cell.type == elem::FRAC_HEX)
			model->h[i] = cell.V * model->solveFrac(cell);
		else if (cell.type == elem::BORDER_HEX)
			model->h[i] = model->solveBorder(cell);
		else if (cell.type == elem::HEX || cell.type == elem::PRISM)
			model->h[i] = cell.V * model->solveInner(cell);
	}

	for (int i = 0; i < Model::var_size * size; i++)
		model->h[i] >>= y[i];

	trace_off();
}
void OilSolver::fill()
{
	sparse_jac(0, Model::var_size * model->cellsNum, Model::var_size * model->cellsNum, repeat,
		&model->u_next[0], &elemNum, (unsigned int**)(&ind_i), (unsigned int**)(&ind_j), &a, options);

	int counter = 0;
	for (int j = 0; j < mesh->inner_size; j++)
	{
		const auto& cell = mesh->cells[j];
		//getMatrixStencil(cell);
		for (int i = 0; i < Model::var_size; i++)
		{
			const int str_idx = Model::var_size * cell.id + i;
			/*for (const int idx : stencil_idx)
			{
				for (int j = 0; j < Model::var_size; j++)
					a[counter++] = jac[str_idx][Model::var_size * idx + j];
			}*/

			rhs[str_idx] = -y[str_idx];
		}
		//stencil_idx.clear();
	}
}