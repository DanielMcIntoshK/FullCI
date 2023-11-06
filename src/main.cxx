#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <cmath>

/*
#include "HartreeFockSolver.h"
#include "ModelParams.h"
#include "BasisBuilder.h"
#include "Integrals.h"
*/
#include "SlaterDet.h"

int main(int argc, char** argv){
	SlaterDet::buildStrings(10,5,252);

	std::cout << SlaterDet::codes.strs.size() << std::endl;
	SlaterDet::printStrings();

	SlaterDet s1(0,0),s2(0,251);
	auto diff(s1|s2);
	std::cout << diff.size() << std::endl;
	for(auto a: diff){std::cout << a << " ";}
	std::cout << std::endl;
	/*libint2::initialize();
	
	if(argc < 2) {std::cout << "NEEDS AN INPUT FILE\n";return -1;}

	std::string inputfilename=argv[1];

	std::ifstream inFile(inputfilename);
	//std::ifstream inFile("../inputs/Benzene_sto3g");
	if(!inFile) {std::cout << "FILE DOES NOT EXIST\n"; return -1;}
	
	std::cout << "LOADING PARAMETERS FROM INPUT FILE\n";
	ModelParams params;
	params.ReadInputFile(inFile);
	inFile.close();

	std::cout << "MOLECULAR GEOMETRY\n";
	for(auto at: params.atoms) {
		std::cout << at.atomic_number << ": " << at.x << " " << at.y << " " << at.z << std::endl;
	}

	std::cout << "LOADING BASIS SET: " << params.cparams["BASISSET"] << std::endl;
	BasisBuilder bb;
	std::string basissetdir = "../basisset/";
	bb.Build(basissetdir+params.cparams["BASISSET"], params);
	//std::cout << "Orbital Basis Set Rank = " << bb.bs.nbf() << std::endl;
	inFile.close();

	/BasisSet bstest("sto-3g", params.atoms);
	std::cout << bstest.nbf() << std::endl;
	for(int i = 0; i < bstest.size();i++){
		libint2::Shell sh = bstest[i];
		libint2::Shell sh2= bb.bs[i];
		std::cout << "SHELL" << i << std::endl;
		for(int j=0; j < sh.alpha.size(); j++){
			std::cout << sh.alpha[j] << " " << sh.contr[0].coeff[j]<< " " << sh2.alpha[j] << " " << sh2.contr[0].coeff[j] <<std::endl;
		}
		std::cout << std::endl;
	}/

	std::cout << "PRELOADING INTEGRALS" << std::endl;
	IntegralChugger ic(bb.bs,params);
	ic.compute(IntegralChugger::ALL);

	std::cout << "RUNNING RESTRICTED HF" << std::endl;
	HartreeFockSolver hfs;
	HartreeFockSolver::HFResults hfr = hfs.RestrictedHF(params, bb.bs,ic);
	//HartreeFockSolver::HFResults hfr = hfs.RestrictedHF(params,bstest);

	std::cout << "HARTREE FOCK ENERGY: " << hfr.eelec << std::endl;
	std::cout << "NUCELAR REPULSE:     " << hfr.enuc << std::endl;
	std::cout << "TOTAL ENERGY:        " << hfr.eelec+hfr.enuc << std::endl;


	std::list<int> testvec{1,2,3,4,5};
	do{
		for(auto a = testvec.begin(); a != testvec.end(); a++){
			std::cout << *a << " ";
		}	
		std::cout << std::endl;
	}while(std::next_permutation(++testvec.begin(), testvec.end()));

	libint2::finalize();

	*/
	return 0;
}
