#ifndef PARALUTIONINTERFACE_H_
#define PARALUTIONINTERFACE_H_

#include <string>

#include "paralution.hpp"

enum class PRECOND {ILU_SIMPLE, ILU_SERIOUS, ILUT, ILU_GMRES};

class ParSolver
{
	enum class RETURN_TYPE { NO_CRITERIA, ABS_CRITERION, REL_CRITERION, DIV_CRITERIA, MAX_ITER };
public:
	typedef paralution::LocalMatrix<double> Matrix;
	typedef paralution::LocalVector<double> Vector;
protected:
	Vector x, Rhs;
	Matrix Mat;
	paralution::BiCGStab<Matrix,Vector,double> bicgstab;
	void SolveBiCGStab();
	void SolveBiCGStab_ILUT();
	void SolveBiCGStab_Simple();
	paralution::GMRES<Matrix,Vector,double> gmres;
	void SolveGMRES();
	paralution::ILU<Matrix,Vector,double> p;
	paralution::ILUT<Matrix, Vector, double> p_ilut;

	bool isAssembled;
	bool isPrecondBuilt;
	int matSize;
	RETURN_TYPE status;

	inline void writeSystem()
	{
		Mat.WriteFileMTX("snaps/mat.mtx");
		Rhs.WriteFileASCII("snaps/rhs.dat");
		x.WriteFileASCII("snaps/x.dat");
	};

	double initRes, finalRes;
	int iterNum;
	const std::string resHistoryFile;
	void getResiduals();
public:
	void Init(const int vecSize, const double relTol, const double dropTol);
	void Assemble(const int* ind_i, const int* ind_j, const double* a, const int counter, const int* ind_rhs, const double* rhs);
	void Solve();
	void Solve(const PRECOND key);

	const Vector& getSolution() { return x; };

	ParSolver();
	~ParSolver();
};

#endif /* PARALUTIONINTERFACE_H_ */