#ifndef __HARTREEFOCKSOLVER__H_
#define __HARTREEFOCKSOLVER__H_

#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <libint2.hpp>

class HartreeFockSolver{
public:
	HartreeFockSolver();

	HFResults RestrictedHF(std::vector<Atom> ats, BasisSet bs);
public:
	struct HFResults{
		char result[100];
	};
private:
};

#endif

