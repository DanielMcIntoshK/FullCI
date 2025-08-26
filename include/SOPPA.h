#ifndef SOPPA__H_
#include "SlaterDet.h"
#include "HartreeFockSolver.h"
#include "Integrals.h"
#include <vector>


class SOPPASolver{
public:
	enum ABlabel{
		ABLAB_A0,
		ABLAB_A1,
		ABLAB_A2,
		ABLAB_Astar,
		ABLAB_B1,
		ABLAB_B2
	};
	enum Clabel{
		CLAB_C1
	};
	enum Dlabel{
		DLAB_D0,
		DLAB_D1
	};
	SOPPASolver(hfresults results, IntegralChugger * ints);

	Matrix FOPPA(double w);
	Matrix SOPPA(double w);

	Matrix ssmat(int i);
	Matrix ddmat(int i);
	Matrix dsmat(int i);

	Matrix C(int i);


	double matelA0(PHOp r, PHOp c);
	double matelAstar(PHOp r, PHOp c);
	double matelD0(PHOp r1, PHOp r2, PHOp c1, PHOp c2);

	double matelC1(PHOp r1, PHOp r2, PHOp c);
	double matelD1(PHOp r1, PHOp r2, PHOp c1, PHOp c2);
	double matelA1(PHOp r, PHOp c);
	double matelB1(PHOp r, PHOp c);
	
	double matelA2(PHOp r, PHOp c);
	double matelB2(PHOp r, PHOp c);

	double f(int m, int p,PHOp r1, PHOp r2, PHOp c1, PHOp c2);

	double E;
	hfresults hfr;
	IntegralChugger * ic;
};


#endif
