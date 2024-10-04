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


		void buildOperators(std::vector<double> avocc, HartreeFockSolver::HFResults & hf);
	
		void cleanup();
		
		std::vector<PHOp> operators;
	private:
		void computeHamiltonian();
		double matrixEl(int x, int y);

		double secondQuantMatel1e(SlaterDet & s1, SlaterDet & s2,SlaterCompare & sc);
		double secondQuantMatel2e(SlaterDet & s1, SlaterDet & s2,SlaterCompare & sc);

	private:
		double  opOnSlater(PHOp op, unsigned char * alpha, unsigned char * beta);

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

