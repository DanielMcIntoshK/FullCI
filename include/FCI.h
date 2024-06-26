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
#include "SlaterDet.h"

class FullCISolver{
	public:
		struct FCIResults{
			Matrix eigenvalues;
			Matrix eigenvectors;
			Matrix * fullH;
		};
	public:
		FullCISolver(){}

		FullCISolver::FCIResults fci(ModelParams & mp, IntegralChugger & ic, HartreeFockSolver::HFResults &hf, double lambda=1.0);
	
		void computeHamiltonian();
		double matrixEl(int x, int y);

		double secondQuantMatel1e(SlaterDet & s1, SlaterDet & s2,SlaterCompare & sc);
		double secondQuantMatel2e(SlaterDet & s1, SlaterDet & s2,SlaterCompare & sc);

		void cleanup();
	private:
		IntegralChugger * ints;
		
		int strcnt, cisize, bssize;

		Matrix CIMat;

		bool lambdaDeriv;
		double lambda;
		double fockTotal;
};

#endif

