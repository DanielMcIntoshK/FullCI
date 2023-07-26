#ifndef __HARTREEFOCKSOLVER__H_
#define __HARTREEFOCKSOLVER__H_

#include <vector>
#include <array>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <libint2.hpp>

#include "ModelParams.h"

using namespace libint2;


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

typedef std::vector<std::vector<std::vector<std::vector<double>>>> twobodylist;

typedef std::vector<std::vector<std::vector<std::vector<bool>>>> checklist;

class HartreeFockSolver{
public:
	struct HFParams{

	};

	struct HFResults{
		double eelec;
		double enuc;
		Matrix C;
		Matrix D;
	};
public:
	HartreeFockSolver();

	HartreeFockSolver::HFResults RestrictedHF(ModelParams & param, BasisSet & bs);
private:
	Matrix compute1eints(ModelParams & param, BasisSet & bs, libint2::Operator obtype);
	void compute2eints(BasisSet &bs);
	void compute2eints_crappy(BasisSet& bs);
	Matrix computeGMatrix(Matrix & D);
	
	Matrix initialDGuess(ModelParams & param,Matrix & H, Matrix &S); 

	std::array<int,4> get2bodyintcord(int a, int b, int c, int d);

	twobodylist twobodyints;

	int nelectron;
};

#endif

