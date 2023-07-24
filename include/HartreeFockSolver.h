#ifndef __HARTREEFOCKSOLVER__H_
#define __HARTREEFOCKSOLVER__H_

#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <libint2.hpp>

using namespace libint2;


class HartreeFockSolver{
public:
	struct HFParams{

	};

	struct HFResults{
		char result[100];
	};
public:
	HartreeFockSolver();

	HartreeFockSolver::HFResults RestrictedHF(std::vector<Atom> ats, BasisSet bs);
private:
};

#endif

