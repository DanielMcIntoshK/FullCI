#ifndef MBPT__H_
#define MBPT__H_
#include <vector>
#include <list>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Integrals.h"
#include "HartreeFockSolver.h"

class MBPT_Solver{
	MBPT_Solver(){}

	std::vector<Matrix> SolveWavefunction(Matrix & fullH,HartreeFockSolver::HFResults & hfr);
};

#endif

