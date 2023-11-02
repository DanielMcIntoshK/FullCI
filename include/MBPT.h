#ifndef MBPT__H_
#define MBPT__H_
#include <vector>
#include <list>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include "Integrals.h"
#include "HartreeFockSolver.h"

struct diagram{
	int size;
	int nelec;
	int norb;

	std::list<int> structure;
	int nhole, npart;
	std::vector<int> hpcount;
	int chole, cpart;

	diagram(int s, int e);
	bool permuteDiagram();	

	bool incrementhp();
};

class ArbitraryMBPT{
public:
	struct MBPTResults{
		int degree;
		std::vector<double> vd; 
	};
public:
	ArbitraryMBPT();

	ArbitraryMBPT::MBPTResults computeMBPT(int degree, int nelec, hfresults &hfr, IntegralChugger & ic);

	double calcDiagram(diagram & d);
private:

};

#endif 

