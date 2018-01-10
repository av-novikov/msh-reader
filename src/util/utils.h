#ifndef UTILS_H_
#define UTILS_H_

#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>

#include "src/util/Interpolate.h"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

#define BAR_TO_PA 1.E5
#define P_ATM 1.0

using std::vector;
using std::pair;
using std::string;
using std::sort;
using std::make_pair;
using std::ifstream;
using std::function;

inline constexpr double delta(const int i, const int j)
{
	return (i == j);
};
inline void setDataFromFile(vector< pair<double,double> >& vec, string fileName)
{
	ifstream file;
	file.open(fileName.c_str(), ifstream::in);
	
	double temp1, temp2;
	while( !file.eof() )
	{
		file >> temp1;
		if( file.eof() )
			break;
		file >> temp2;
		vec.push_back(make_pair(temp1, temp2));
	}

	file.close();
};
inline bool IsNan(double a)
{
	if (a!=a)  return true;
	return false;
};
template <typename TData>
double sign(TData a)
{
	if (a > (TData)0) return 1.0;
	else if (a < (TData)0) return -1.0;
	else return 0.0;
};

inline double MilliDarcyToM2(double perm)
{
	return perm * 0.986923 * 1.E-15;
};
inline double M2toMilliDarcy(double perm)
{
	return perm * 1.E15 / 0.986923;
};
inline double cPToPaSec(double visc)
{
	return visc / 1000.0;
};
inline double gramToKg(double grams)
{
	return grams / 1000;
};
inline double PaSec2cP(double visc)
{
	return visc * 1000.0;
};

template <typename T>
inline std::string to_string(T value)
{
	std::ostringstream os ;
	os << value;
	return os.str();
};

struct sort_pair_first {
    bool operator() (const std::pair<double,double> &left, const std::pair<double,double> &right) 
	{
        return left.first < right.first;
    }
};
struct sort_pair_second {
    bool operator() (const std::pair<double,double> &left, const std::pair<double,double> &right) 
	{
        return left.second < right.second;
    }
};

inline Interpolate* setDataset(vector< pair<double,double> >& vec, const double xDim, const double yDim)
{
	sort(vec.begin(), vec.end(), sort_pair_first());

	const size_t N = vec.size();
	double* x = new double [N];
	double* y  = new double [N];

	for (size_t i = 0; i < N; i++)
	{
		x[i] = vec[i].first / xDim;
		y[i] = vec[i].second / yDim;
	}

	return new Interpolate(x, y, N);
};
inline Interpolate* setInvDataset(vector< pair<double,double> >& vec, const double xDim, const double yDim)
{
	sort(vec.begin(), vec.end(), sort_pair_second());

	const size_t N = vec.size();
	double* x = new double [N];
	double* y  = new double [N];

	for (size_t i = 0; i < N; i++)
	{
		x[i] = vec[i].second / xDim;
		y[i] = vec[i].first / yDim;
	}

	return new Interpolate(x, y, N);
};


inline adouble** multiply(const adouble** a, const adouble** b, const int N)
{
	adouble** c = new adouble*[N];
	for (int i = 0; i < N; i++)
		c[i] = new adouble[N];
	adouble sum;
	for (int i = 0; i < N; i++)
		for (int j = 0; j < N; j++)
		{
			sum = 0.0;
			for (int k = 0; k < N; k++)
				sum += a[i][k] * b[k][j];
			c[i][j] = sum;
		}
	return c;
};

template <typename ValType>
class QR_Decomp
{
protected:
	void QR_Decomposition()
	{
		int i = 0, j = 0, k = 0;
		// Main loop.
		for (k = 0; k < n; k++) {
			// Compute 2-norm of k-th column without under/overflow.
			adouble nrm = 0;
			for (i = k; i < n; i++)
				nrm = sqrt(nrm * nrm + mat[i][k] * mat[i][k]);

			if (nrm.value() != 0.0) {
				// Form k-th Householder vector.
				adouble tmp = nrm;
				adouble cond = mat[k][k].value() < 0 ? true : false;
				condassign(tmp, cond, -nrm);
				nrm = tmp;
				/*if (mat[k][k].value() < 0) {
				nrm = -nrm;
				}*/
				for (i = k; i < n; i++) {
					mat[i][k] /= nrm;
				}
				mat[k][k] += 1.0;

				// Apply transformation to remaining columns.
				for (j = k + 1; j < n; j++) {
					adouble s = 0.0;
					for (i = k; i < n; i++) {
						s += mat[i][k] * mat[i][j];
					}
					s = -s / mat[k][k];
					for (i = k; i < n; i++) {
						mat[i][j] += s*mat[i][k];
					}
				}
			}

			R_Diag_Vector[k] = -nrm;
		}

	}
	ValType** mat_alloc()
	{
		ValType** b = new ValType*[n];
		for (int i = 0; i < n; i++)
			b[i] = new ValType[n];

		return b;
	};
	ValType* vec_alloc()
	{
		ValType* b = new ValType[n];
		return b;
	};

	ValType* R_Diag_Vector;
	ValType** mat;
	int n;
public:
	QR_Decomp() {};
	~QR_Decomp()
	{
		mat_free(mat);
		vec_free(R_Diag_Vector);
	};

	void mat_free(ValType** b)
	{
		for (int i = 0; i < n; i++)
			delete[] b[i];
		delete[] b;
	}
	void vec_free(ValType* b)
	{
		delete[] b;
	}

	ValType** Invert(ValType** ImpactMatrix, const int _n)
	{
		n = _n;
		mat = mat_alloc();
		ValType** X = mat_alloc();
		R_Diag_Vector = vec_alloc();
		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < n; j++)
			{
				mat[i][j] = ImpactMatrix[i][j];
				X[i][j] = 0.0;
			}
			X[i][i] = 1.0;
		}

		int i = 0, j = 0, k = 0;

		QR_Decomposition();
		// Compute Y = transpose(Q)*B
		for (k = 0; k < n; k++) {
			for (j = 0; j < n; j++) {
				adouble s = 0.0;
				for (i = k; i < n; i++) {
					s += mat[i][k] * X[i][j];
				}
				s = -s / mat[k][k];
				for (i = k; i < n; i++) {
					X[i][j] += s*mat[i][k];
				}
			}
		}

		// Solve R*X = Y;
		for (k = n - 1; k >= 0; k--) {
			for (j = 0; j < n; j++) {
				X[k][j] /= R_Diag_Vector[k];
			}
			for (i = 0; i < k; i++) {
				for (j = 0; j < n; j++) {
					X[i][j] -= X[k][j] * mat[i][k];
				}
			}
		}

		return X;
	}
};


#endif /* UTILS_H_ */