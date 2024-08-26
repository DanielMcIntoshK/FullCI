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

struct PHOp{
	PHOp():i{0},j{0}{}
	PHOp(int a, int b):i{a},j{b}{}
	int i,j;

	PHOp adjoint(){return PHOp(j,i);}
};

class FullCISolver{
	public:
		struct FCIResults{
			Matrix eigenvalues;
			Matrix eigenvectors;
			Matrix * fullH;
		};
		struct MBPTResults{
			std::vector<Matrix> wavefunctions;
			std::vector<double> energies;
		};
	public:
		FullCISolver(){}


		FullCISolver::FCIResults fci(IntegralChugger & ic, HartreeFockSolver::HFResults &hf, double lambda=1.0);
		FullCISolver::MBPTResults mbpt(HartreeFockSolver::HFResults &hf,int order);
		std::vector<Matrix> recursivegreen(int order, double E, HartreeFockSolver::HFResults & hf, FullCISolver::MBPTResults & mbptr); 

		void computeHamiltonian();
		double matrixEl(int x, int y);

		double secondQuantMatel1e(SlaterDet & s1, SlaterDet & s2,SlaterCompare & sc);
		double secondQuantMatel2e(SlaterDet & s1, SlaterDet & s2,SlaterCompare & sc);

		void cleanup();
	private:
		int opOnSlater(PHOp op, int det,bool verbose=false);

		IntegralChugger * ints;
		
		int strcnt, cisize, bssize;

		Matrix CIMat;
		Matrix H0;

		bool lambdaDeriv;
		double lambda;
		double fockTotal;

		bool fciSuccess=false;
};

#endif

