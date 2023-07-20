#include <iostream>

#include "HartreeFockSolver.h"


int maint(int argc, char** argv){
	HartreeFockSolver hfs;
	HFResults hfr = hfs.RestrictedHF(std::vector<Atom>(), BasisSet());
	std::cout << hfr <<std::endl;

	return 0;
}
