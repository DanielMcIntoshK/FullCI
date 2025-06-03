#ifndef POLESEARCH__H_
#define POLESEARCH__H_
#include "FCI.h"
#include "Greens.h"
#include "HartreeFockSolver.h"
#include "Integrals.h"
#include <string>

class PoleSearch{
public:
	PoleSearch(FullCISolver * fci_, GreensCalculator * gr_, IntegralChugger * ic_,HartreeFockSolver::HFResults hfr_,FullCISolver::FCIResults fcir_, FullCISolver::MBPTResults mbptr_);

	std::vector<double> refinePoints(std::vector<double> EE, double threshold);
	std::vector<double> refinePointsOrder(std::vector<double> EE, double threshold, int order);

	std::vector<double> scan(double E_s, double dE, int steps, std::string filename="");
	void mapSelfEnergy(double E_s, double dE, int steps, std::string filename);

	double refinePoint(double point,double threshold);
private:
	double detval(double E);

	FullCISolver * fci;
	GreensCalculator * gr;
	IntegralChugger * ic;

	FullCISolver::FCIResults fcir;
	FullCISolver::MBPTResults mbptr;
	HartreeFockSolver::HFResults hfr;

	Matrix orbEdiff;
};

#endif

