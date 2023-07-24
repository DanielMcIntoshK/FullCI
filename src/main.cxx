#include <iostream>
#include <fstream>

#include "HartreeFockSolver.h"
#include "ModelParams.h"

int main(int argc, char** argv){
	libint2::initialize();
	
	std::ifstream inFile("../inputs/BH_sto3g");
	if(!inFile) {std::cout << "FILE DOES NOT EXIST\n"; return -1;}
	
	ModelParams mp;
	mp.ReadInputFile(inFile);

	inFile.close();

	HartreeFockSolver hfs;
	//HartreeFockSolver::HFResults hfr = hfs.RestrictedHF(std::vector<Atom>(), BasisSet());
	//std::cout << hfr <<std::endl;

	std::cout << "HELLO WORLD!\n";
	libint2::finalize();

	return 0;
}
