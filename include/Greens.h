#ifndef GREENS__H_
#define GREENS__H_

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "FCI.h"
#include "Integrals.h"
#include "HartreeFockSolver.h"
#include "SlaterDet.h"


struct chslater{
	int index;
	double sign, c;

	bool operator<(const chslater &chs)const{return index<chs.index;}
	bool operator==(const chslater &chs)const{return index==chs.index;}
};

class GreensCalculator{
public:
	GreensCalculator();

	void buildManifold(HartreeFockSolver::HFResults & hf);

	Matrix ComputeGreens(double E, IntegralChugger & ic, HartreeFockSolver::HFResults &hf, FullCISolver::FCIResults & fcir);
	Matrix ComputeSelfEnergy(double E, IntegralChugger & ic, HartreeFockSolver::HFResults &hf, FullCISolver::FCIResults &fcir);

public:
	Matrix buildProp(double E, FullCISolver::FCIResults & fcir);
	double calcPropEl(double E,PHOp q1, PHOp q2,Matrix &evals, Matrix &evecs);
	Matrix buildProp_Smart(double E, FullCISolver::FCIResults & fcir);
	double calcPropEl_Smart(int q1, int q2,Matrix & OM, Matrix & MO, std::vector<double> & Ep, std::vector<double> & Em);

	std::vector<double> computeAvOc_Pure(Matrix state);
	Matrix computeOverlap(Matrix state);
	Matrix computeLMatrix(std::vector<double> & ocavs);
	Matrix computeAMatrix(Matrix state,Matrix ev,Matrix*fullH);
       	Matrix computeBMatrix(Matrix state,Matrix ev,Matrix*fullH);

	double qqH(PHOp q1, PHOp q2,Matrix state, double e0);
	double qHq(PHOp q1, PHOp q2,Matrix state, Matrix * fullH);
	double qq(PHOp q1, PHOp q2, Matrix state,Matrix * fullH=nullptr);
	double NqM(PHOp q, Matrix N, Matrix M);
public:
	IntegralChugger * ints;
	std::vector<PHOp> manifold;
};

double applyPHOp(PHOp op, unsigned char * alpha, unsigned char * beta);

#endif 

