#include "src/util/ParalutionInterface.h"

#include <fstream>
#include <iostream>

using namespace paralution;
using std::ifstream;
using std::cout;
using std::endl;

ParSolver::ParSolver() : resHistoryFile("snaps/resHistory.dat")
{
	isAssembled = false;
	isPrecondBuilt = false;
	gmres.Init(1.E-17, 1.E-12, 1E+12, 500);
	bicgstab.Init(1.E-17, 1.E-12, 1E+12, 500);
}
ParSolver::~ParSolver()
{
}
void ParSolver::Init(const int vecSize, const double relTol, const double dropTol)
{
	matSize = vecSize;
	x.Allocate("x", vecSize);
}
void ParSolver::Assemble(const int* ind_i, const int* ind_j, const double* a, const int counter, const int* ind_rhs, const double* rhs)
{
	Mat.Zeros();
	Rhs.Zeros();
	x.Zeros();

	/*if (isAssembled)
	{
		Mat.AssembleUpdate(a);
		Rhs.Assemble(ind_rhs, rhs, matSize, "rhs");
	} else {*/
		Mat.Assemble(ind_i, ind_j, a, counter, "A", matSize, matSize);
		Rhs.Assemble(ind_rhs, rhs, matSize, "rhs");

		Mat.MoveToAccelerator();
		Rhs.MoveToAccelerator();
		x.MoveToAccelerator();
	//}	
}
void ParSolver::Solve()
{
	//SolveGMRES();
	SolveBiCGStab();
		
	x.MoveToHost();
}
void ParSolver::Solve(const PRECOND key)
{
	if (key == PRECOND::ILU_SERIOUS)
		SolveBiCGStab();
	else if (key == PRECOND::ILU_SIMPLE)
		SolveBiCGStab_Simple();
	else if (key == PRECOND::ILU_GMRES)
		SolveGMRES();
	else if (key == PRECOND::ILUT)
		SolveBiCGStab_ILUT();

	x.MoveToHost();
}
void ParSolver::SolveBiCGStab_ILUT()
{
	bicgstab.SetOperator(Mat);
	p_ilut.Set(1.E-20, 100);
	bicgstab.SetPreconditioner(p_ilut);
	bicgstab.Build();
	isAssembled = true;

	bicgstab.Init(1.E-20, 1.E-15, 1E+12, 1000);
	Mat.info();

	//bicgstab.RecordResidualHistory();
	bicgstab.Solve(Rhs, &x);
	status = static_cast<RETURN_TYPE>(bicgstab.GetSolverStatus());
	//if(status == RETURN_TYPE::DIV_CRITERIA || status == RETURN_TYPE::MAX_ITER)
	//bicgstab.RecordHistory(resHistoryFile);
	writeSystem();

	bicgstab.Clear();
}
void ParSolver::SolveBiCGStab()
{
	bicgstab.SetOperator(Mat);
	//p.Set(1.E-15, 100);
	p.Set(0);
	bicgstab.SetPreconditioner(p);
	bicgstab.Build();
	isAssembled = true;

	bicgstab.Init(1.E-25, 1.E-9, 1E+12, 1000);
	Mat.info();

	//bicgstab.RecordResidualHistory();
	bicgstab.Solve(Rhs, &x);
	status = static_cast<RETURN_TYPE>(bicgstab.GetSolverStatus());
	//if(status == RETURN_TYPE::DIV_CRITERIA || status == RETURN_TYPE::MAX_ITER)
	//bicgstab.RecordHistory(resHistoryFile);
	writeSystem();

	//getResiduals();
	//cout << "Initial residual: " << initRes << endl;
	//cout << "Final residual: " << finalRes << endl;
	//cout << "Number of iterations: " << iterNum << endl << endl;

	bicgstab.Clear();
}
void ParSolver::SolveBiCGStab_Simple()
{
	bicgstab.SetOperator(Mat);
	//p.Set(1.E-15, 100);
	p.Set(0);
	bicgstab.SetPreconditioner(p);
	bicgstab.Build();
	isAssembled = true;

	bicgstab.Init(1.E-30, 1.E-12, 1E+12, 1000);
	Mat.info();

	//bicgstab.RecordResidualHistory();
	bicgstab.Solve(Rhs, &x);
	status = static_cast<RETURN_TYPE>(bicgstab.GetSolverStatus());
	//if(status == RETURN_TYPE::DIV_CRITERIA || status == RETURN_TYPE::MAX_ITER)
	//bicgstab.RecordHistory(resHistoryFile);
	writeSystem();

	//getResiduals();
	//cout << "Initial residual: " << initRes << endl;
	//cout << "Final residual: " << finalRes << endl;
	//cout << "Number of iterations: " << iterNum << endl << endl;

	bicgstab.Clear();
}
void ParSolver::SolveGMRES()
{
	gmres.SetOperator(Mat);
	p_ilut.Set(1.E-20, 100);
	//p.Set(3);
	gmres.SetPreconditioner(p_ilut);
	gmres.Build();
	isAssembled = true;

	gmres.Init(1.E-17, 1.E-12, 1E+12, 500);
	Mat.info();

	//gmres.RecordResidualHistory();
	gmres.Solve(Rhs, &x);
	status = static_cast<RETURN_TYPE>(bicgstab.GetSolverStatus());
	//gmres.RecordHistory(resHistoryFile);
	//writeSystem();


	//getResiduals();
	//cout << "Initial residual: " << initRes << endl;
	//cout << "Final residual: " << finalRes << endl;
	//cout << "Number of iterations: " << iterNum << endl << endl;

	gmres.Clear();
}
void ParSolver::getResiduals()
{
	double tmp;
	int i = 0;

	ifstream file;
	file.open(resHistoryFile, ifstream::in);

	file >> initRes;
	while ( !file.eof() )
	{
		file >> tmp;
		i++;
	}
	finalRes = tmp;
	iterNum = i;

	file.close();
}