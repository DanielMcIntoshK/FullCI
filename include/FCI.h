#ifndef FCI__H_
#define FCI__H_

#include <vector>
#include <array>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "ModelParams.h"
#include "Integrals.h"
#include "HartreeFockSolver.h"

class FullCISolver{
	public:
		struct FCIResults{
			Matrix eigenvectors;
		};
	public:
		FullCISolver(){}

		FullCISolver::FCIResults fci(ModelParams & mp, IntegralChugger & ic, HartreeFockSolver::HFResults &hf);
	
		void computeHamiltonian();
		double matrixEl(int x, int y);

		double matel1e(std::vector<int> &diff);
		double matel2e(std::vector<int> &diff,bool verbose=false);

		void cleanup();
	private:
		IntegralChugger * ints;
		
		int strcnt, cisize;

		Matrix CIMat;
};

#endif

