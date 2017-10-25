#include "src/models/AbstractSolver.hpp"
#include "src/util/utils.h"
#include <iomanip>
#include <iterator>

#include "src/models/oil/Oil.hpp"

using namespace std;

template <class modelType>
AbstractSolver<modelType>::AbstractSolver(modelType* _model) : model(_model), mesh(_model->getMesh()), size(_model->cellsNum), Tt(model->period[model->period.size() - 1])
{
	NEWTON_STEP = 1.0;
	cur_t = cur_t_log = 0.0;
	curTimePeriod = 0;

	idx1 = int((model->perfIntervals[0].first + model->perfIntervals[0].second) / 2);
	//idx2 = idx1 + model->cellsNum_z + 1;

	t_dim = model->t_dim;
	repeat = 0;
}
template <class modelType>
AbstractSolver<modelType>::~AbstractSolver()
{
}
template <class modelType>
void AbstractSolver<modelType>::start()
{
	int counter = 0;
	iterations = 8;

	model->setPeriod(curTimePeriod);
	while(cur_t < Tt)
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
template <class modelType>
void AbstractSolver<modelType>::doNextStep()
{
	solveStep();
}
template <class modelType>
void AbstractSolver<modelType>::fill()
{
}
template <class modelType>
void AbstractSolver<modelType>::copyIterLayer()
{
	model->u_iter = model->u_next;
}
template <class modelType>
void AbstractSolver<modelType>::revertIterLayer()
{
	model->u_next = model->u_iter;
}
template <class modelType>
void AbstractSolver<modelType>::copyTimeLayer()
{
	model->u_prev = model->u_iter = model->u_next;
}

template <class modelType>
double AbstractSolver<modelType>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	varInd = 0;
	auto diff = std::abs((model->u_next - model->u_iter) / model->u_next);
	auto max_iter = std::max_element(std::begin(diff), std::end(diff));
	ind = std::distance(std::begin(diff), max_iter);
			
	return *max_iter;
}
template <class modelType>
double AbstractSolver<modelType>::averValue(const int varInd)
{
	double tmp = 0.0;
	for (size_t i = 0; i < model->cellsNum; i++)
	{
		auto cell = (*model)[i];
		tmp += cell.u_next.p * mesh->cells[i].V;
	}
	return tmp / mesh->Volume;
}
template <class modelType>
void AbstractSolver<modelType>::averValue(std::array<double, var_size>& aver)
{
	std::fill(aver.begin(), aver.end(), 0.0);

	for (int i = 0; i < var_size; i++)
	{
		const auto var = static_cast<std::valarray<double>>(model->u_next[std::slice(i, model->cellsNum, var_size)]);
		int cell_idx = 0;
		for (const auto& cell : mesh->cells)
			aver[i] += var[cell_idx++] * cell.V;
	}
	for(auto& val : aver)
		val /= mesh->Volume;
}

template <class modelType>
void AbstractSolver<modelType>::checkStability()
{
}

template class AbstractSolver<oil::Oil>;