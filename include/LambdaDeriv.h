#ifndef LAMBDADERIV__H_
#define LAMBDADERIV__H_
#include "Greens.h"
#include "HartreeFockSolver.h"

class LambdaDeriv{
public:
	LambdaDeriv(){}

	Matrix ComputeSelfEnergyOrder(int M, double x0, std::vector<double> grid, HartreeFockSolver::HFResults hfr,IntegralChugger & ic);	

	std::vector<Matrix> ComputeGreensNumerical(int M, double x0, std::vector<double> grid, HartreeFockSolver::HFResults hfr, IntegralChugger & ic,double E);
private:
	std::vector<Matrix> generateDeltas(int M, int N, double x0, std::vector<double> &grid);
};

#endif

