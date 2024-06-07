#ifndef GREENS__H_
#define GREENS__H_

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "FCI.h"
#include "Integrals.h"
#include "HartreeFockSolver.h"
#include "SlaterDet.h"

struct PHOpp{
	PHOpp():i{0},j{0}{}
	PHOpp(int a, int b):i{a},j{b}{}
	int i, j;

	PHOpp adjoint(){return PHOpp(j,i);}
};

struct chslater{
	int index;
	double sign, c;

	bool operator<(const chslater &chs)const{index<chs.index;}
};

class GreensCalculator{
public:
	GreensCalculator();

	void computeExcitationSpectra(ModelParams & mp, IntegralChugger & ic, HartreeFockSolver::HFResults &hf, FullCISolver::FCIResults &fcir);

private:
	std::vector<double> computeAvOc_Pure(Matrix state);
	Matrix computeLMatrix(std::vector<double> & ocavs);
	Matrix computeAMatrix(Matrix state,Matrix ev,Matrix*fullH);
       	Matrix computeBMatrix(Matrix state,Matrix ev,Matrix*fullH);

	double qqH(PHOpp q1, PHOpp q2,Matrix state, double e0);
	double qHq(PHOpp q1, PHOpp q2,Matrix state, Matrix * fullH);
private:
	IntegralChugger * ints;
	std::vector<PHOpp> manifold;
};

double applyPHOpp(PHOpp opp, unsigned char * sd);

#endif 

