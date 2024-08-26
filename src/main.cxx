#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <cmath>
#include <bitset>

#include "HartreeFockSolver.h"
#include "ModelParams.h"
#include "BasisBuilder.h"
#include "Integrals.h"
#include "FCI.h" 
#include "SlaterDet.h"
#include "Greens.h"
#include "LambdaDeriv.h"

int main(int argc, char** argv){
	libint2::initialize();
	
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
	inFile.close();

	std::cout << "PRELOADING INTEGRALS" << std::endl;
	IntegralChugger ic(bb.bs,params);
	ic.compute(IntegralChugger::ALL);

	std::cout << "RUNNING RESTRICTED HF" << std::endl;
	HartreeFockSolver hfs;
	HartreeFockSolver::HFResults hfr = hfs.RestrictedHF(params, bb.bs,ic);

	std::cout << "HARTREE FOCK ENERGY: " << hfr.eelec << std::endl;
	std::cout << "NUCELAR REPULSE:     " << hfr.enuc << std::endl;
	std::cout << "TOTAL ENERGY:        " << hfr.eelec+hfr.enuc << std::endl;

	ic.TransformInts(hfr.C);
	ic.TransformFock(hfr.C,hfr.G);

	int N = hfr.C.rows();
	SlaterDet::buildStrings(N,params.nelec/2);
	std::cout << "DECODE TEST:\n";
	for(int i =0; i < SlaterDet::codes.strs.size();i++){
		std::cout << i << ": " << SlaterDet::decode(SlaterDet::codes.strs[i]) << std::endl;
	}

	FullCISolver fcis;
	FullCISolver::FCIResults fcir =fcis.fci(ic,hfr,1.0);

	std::cout << "FCI E: " << fcir.eigenvalues(0,0) +hfr.enuc << std::endl;

	int order=14;
	FullCISolver::MBPTResults mbptr=fcis.mbpt(hfr,order);

	double EFC=0.0;
	std::cout << mbptr.wavefunctions.size() << " " << mbptr.energies.size() <<std::endl;
	for(int i = 0; i <= order;i++){
		//std::cout << i << "th Order Wavefunction\n";
		//std::cout << mbptr.wavefunctions[i] << std::endl << std::endl;
	}
	for(int i = 0; i <= order;i++){
		std::cout << i << "th Order Energy: ";
		std::cout << mbptr.energies[i] << std::endl;
		EFC+=mbptr.energies[i];
	}
	std::cout << EFC+hfr.enuc << std::endl;

	double E=-0.02;
	
	std::vector<Matrix> greens=fcis.recursivegreen(order,E,hfr,mbptr);
	Matrix rg=greens[0];
	for(int i = 1; i <=order; i++){
		rg+=greens[i];
		std::cout << i << "th order:\n";
		std::cout << greens[i] << std::endl<<std::endl;
	}
	
	//std::cout << std::setprecision(3) << greens[0] <<std::endl;
	//std::cout << std::setprecision(2)<< greens[1] << std::endl << std::endl;

	GreensCalculator gr;
	gr.buildManifold(hfr);
	Matrix greensfull=gr.ComputeGreens(E,ic,hfr,fcir);

	std::cout << "RECURSIVEGREEN\n";
	std::cout << rg << std::endl;
	
	//for(int i = 0; i < rg.rows();i++){
	//	std::cout << rg(i,i) << std::endl;
	//}

	std::cout << std::endl <<std::endl << std::endl << "FULLGREEN\n";
	for(int i = 0; i < greensfull.rows();i++){
	for(int j = 0; j < greensfull.cols();j++){
		if(greensfull(i,j)<0.000001)greensfull(i,j)=0.0;
	}
	}
	std::cout << greensfull << std::endl;

	/*

	LambdaDeriv ld;
	std::vector<double> grid;
	double gridspace=0.01;
	grid.push_back(0.0);
	for(int i = 1; i <=3; i++){
		grid.push_back(gridspace*i);
		grid.push_back(-gridspace*i);
	}
	Matrix se=ld.ComputeSelfEnergyOrder(0, 0, grid, hfr,ic);
	for(int i = 0; i < se.rows();i++){
		for(int j = 0; j<se.cols();j++){
			if(se(i,j)<0.0000001) se(i,j)=0.0;
		}
	}

	std::cout << se << std::endl;
	*/

	fcis.cleanup();
	SlaterDet::cleanStrings();

	libint2::finalize();


	return 0;
}
