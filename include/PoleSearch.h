#ifndef POLESEARCH__H_
#define POLESEARCH__H_
#include "FCI.h"
#include "Greens.h"
#include "HartreeFockSolver.h"
#include "Integrals.h"
#include <string>

class PoleSearch{
public:
	PoleSearch(FullCISolver * fci_, GreensCalculator * gr_, IntegralChugger * ic_);

	std::vector<double> scan(double E_s, double dE, int steps, FullCISolver::FCIResults fcir, 
			FullCISolver::MBPTResults mbptr,HartreeFockSolver::HFResults hfr, std::string filename="");
private:
	FullCISolver * fci;
	GreensCalculator * gr;
	IntegralChugger * ic;
};

#endif

