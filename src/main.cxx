#include <iostream>
#include <fstream>
#include <algorithm>
#include <list>
#include <stack>
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
#include "PoleSearch.h"

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
	TDHFSolver tdhfs;
	TDHFSolver::TDHFResults tdhfr=tdhfs.TDHFCalc(hfr,ic);

	std::cout << "TDHF EXCITATION ENERGIES:\n";
	for(int i = 0; i < tdhfr.EE.rows();i++){
		std::cout << tdhfr.EE(i,0)<<std::endl;
	}

	return 0;

	FullCISolver fcis;
	FullCISolver::FCIResults fcir =fcis.fci(ic,hfr,1.0);

	std::cout << "FCI E: " << fcir.eigenvalues(0,0) +hfr.enuc << std::endl;
	
	double E=.2;

	std::cout << "ATEMPTING TO BUILD FULL GREEN'S FUNCTION\n";
	GreensCalculator gr;
	std::vector<double> avocc=gr.computeAvOc_Pure(fcir.eigenvectors.col(0));
	std::vector<double> avocc1=gr.computeAvOc_Pure(fcir.eigenvectors.col(1));

	Matrix greensfull=gr.ComputeGreens(E,ic,hfr,fcir);

	std::cout << "AVERAGE OCCUPATION\n";
	for(int i = 0; i < avocc.size(); i++) std::cout << avocc[i] << " ";
	std::cout << std::endl;
	for(int i = 0; i < avocc.size(); i++) std::cout << avocc1[i] << " ";
	std::cout << std::endl;
	
	int order=10;
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

	
	
	//CALCULATE GREENS FUNCTION RECURSIVELY
	std::vector<Matrix> greens=fcis.recursivegreen(order,E,hfr,mbptr);
	Matrix rg=greens[0];
	std::cout << "0th order: \n";
	std::cout << greens[0] << std::endl<< std::endl << " " << rg(1,2) << std::endl;
	for(int i = 1; i <=order; i++){
		rg+=greens[i];
		//std::cout << i << "th order:\n";
		//std::cout << greens[i] << std::endl<<std::endl;
		double max=0.0;
		int ri=0,ci=0;
		Matrix sdiff=rg-greensfull;
		for(int r = 0; r < sdiff.rows(); r++){
			for(int c = 0; c < sdiff.cols();c++){
				if(std::abs(sdiff(r,c))>max){
					max=std::abs(sdiff(r,c));
					ri=r;
					ci=c;
				}
			}
		}
		std::cout << i << " " << max << " " << ri << " " << ci << " " << rg(ri,ci)<<" " << greensfull(ri,ci)<<std::endl;
	}

	
	std::cout << "RECURSIVE GREEN'S COMPLETE\n";
	
	std::cout << "RECURSIVEGREEN\n";
	std::cout << rg << std::endl<< std::endl;
	

	std::cout << std::endl <<std::endl << std::endl << "FULLGREEN\n";
	std::cout << greensfull << std::endl;

	std::cout << "G0, G0i test\n";
	std::cout << "RECURSEIVE G0\n" << greens[0]<<std::endl<<std::endl<<
		"GENERATED G0\n" << hfr.getG0(E) << std::endl << std::endl<<
		"G0i TEST\n" << greens[0]*hfr.getG0i(E)<<std::endl<<std::endl;

	int p=2, q=1;
	
	std::cout << "EXACT: " <<greensfull(p,q) << std::endl;
	double val = 0;
	for(int i = 0; i <= order; i++){
		val+= greens[i](p,q);
		std::cout << i << ": " << val <<  " " << greens[i](p,q) << " " << greens[i](q,p) <<  std::endl;
	}

	//CALCULATE GREENS FUNCTION NUMERICALLY
	LambdaDeriv ld;
	std::vector<double> grid;
	double gridspace=0.1;
	grid.push_back(0.0);
	int lambdaorder = 11;
	for(int i = 1; i <=lambdaorder/2; i++){
		grid.push_back(gridspace*i);
		grid.push_back(-gridspace*i);
	}
	std::vector<Matrix> greensnum=ld.ComputeGreensNumerical(lambdaorder,0.0,grid,hfr,ic,avocc,E);
	
	std::cout << "GREENS NUMERICAL\n";
	Matrix greensNum=greensnum[0];
	std::cout << "0th:\n" << greensnum[0]<<std::endl<<std::endl;
	for(int i = 1; i <= lambdaorder; i++){
		std::cout << i<<"th:\n" << greensnum[i]<<std::endl<<std::endl;
		greensNum+=greensnum[i];
	}
	std::cout << "SUM\n" << greensNum << std::endl << std::endl;

	std::cout << "COMPARE\n";
	for(int i = 0; i < std::min(greensnum.size(),greens.size());i++){
		Matrix mdiff=greensnum[i]-greens[i];
		double max=0.0;
		int ri=0,ci=0;
		for(int r = 0; r < mdiff.rows();r++){
			for(int c = 0; c < mdiff.cols();c++){
				if(std::abs(mdiff(r,c))>max){
					max=std::abs(mdiff(r,c));
				}
			}
		}
		std::cout << i << " " << max << std::endl;
	}
	for(int i = 0; i < std::min(greensnum.size(),greens.size());i++){
		Matrix mdiff=greensnum[i]-greens[i];
		double max=0.0;
		int ri=0,ci=0;
		for(int r = 0; r < mdiff.rows();r++){
			for(int c = 0; c < mdiff.cols();c++){
				if(std::abs(mdiff(r,c))>max){
					max=std::abs(mdiff(r,c));
				}
			}
		}
		std::cout << i << " " << max << std::endl;
	}

	//SELF ENERGY CALCULATIONS
	Matrix Mfull=gr.ComputeSelfEnergy(E,hfr,greensfull);

	
	//std::vector<Matrix> senum=gr.ComputeSelfEnergies(E,lambdaorder,hfr,greensnum);
	std::vector<Matrix> senum=ld.ComputeSelfEnergyNumerical(lambdaorder,0.0,grid,hfr,ic,avocc,E);
	std::cout << "SELF ENERGY NUMERICAL\n";
	Matrix seSum=Matrix::Zero(senum[1].rows(),senum[1].cols());;
	for(int i = 1; i < senum.size();i++){
		//std::cout << i<<"th:\n"<<senum[i]<<std::endl<<std::endl;
		seSum+=senum[i];
		double max = 0.0;
		int ri=0,ci=0;
		//Matrix sumdiff=serecursive[i]-senum[i];
		Matrix sumdiff=seSum-Mfull;
		for(int r = 0; r < sumdiff.rows();r++){
			for(int c = 0; c < sumdiff.cols();c++){
				if(std::abs(sumdiff(r,c))>max){
					max=std::abs(sumdiff(r,c));
					ri=r;
					ci=c;
				}
			}
		}
		std::cout << i << " " << max << " "<<ri<<" " << ci<<" "<<Mfull(ri,ci) << " " << seSum(ri,ci)<< std::endl;
	}
	std::cout << "SUM\n" << seSum<<std::endl<<std::endl;

	std::cout << greens[1].rows() << " " << greens[1].cols() << " " << hfr.operators.size()<<std::endl;;
	std::cout << "SELF ENERGY RECURSIVE\n"<<hfr.operators.size()<<std::endl;;
	std::vector<Matrix> serecursive=gr.ComputeSelfEnergies(E,order,hfr,greens);
	seSum=Matrix::Zero(seSum.rows(),seSum.cols());
	for(int i = 1; i < std::min(serecursive.size(),senum.size()-1);i++){
		//std::cout << i<<"th:\n"<<serecursive[i]<<std::endl<<std::endl;
		seSum+=serecursive[i];
		double max = 0.0;
		int ri=0,ci=0;
		//Matrix sumdiff=serecursive[i]-senum[i];
		Matrix sumdiff=serecursive[i]-senum[i];
		for(int r = 0; r < sumdiff.rows();r++){
			for(int c = 0; c < sumdiff.cols();c++){
				if(std::abs(sumdiff(r,c))>max){
					max=std::abs(sumdiff(r,c));
					ri=r;
					ci=c;
				}
			}
		}
		std::cout << i << " " << max << " "<<ri<<" " << ci<<" "<<serecursive[i](ri,ci) << " " << senum[i](ri,ci)<< std::endl;
	}
	std::cout << "SUM\n"<<seSum<<std::endl<<std::endl;

	std::cout << "Full M\n" << Mfull <<std::endl<<std::endl;

	std::cout << "HERMITIAN TEST\n";
	for(int n = 0; n < serecursive.size();n++){
		bool hermitian=true;
		for(int i = 0; i < serecursive[n].rows(); i++){
			for(int j = 0; j < serecursive[n].cols();j++){
				if(std::abs(serecursive[n](i,j)-serecursive[n](j,i))>0.000001){
					hermitian=false;
				}
			}
		}
		std::cout << n << " ";
		if(hermitian) std::cout << "HERMITIAN\n";
		else std::cout << "NONHERMITIAN\n";
	}

	return 0;
	//std::vector<double> excitations;
	std::vector<std::vector<double>> excitations;
	excitations.push_back(std::vector<double>());

	excitations[0].resize(tdhfr.EE.rows());
	std::cout << "TDHF VALS:\n"; 
	for(int i = 0; i < excitations[0].size(); i++){
		excitations[0][i]=tdhfr.EE(i,0);
		std::cout << i<< " " << excitations[0][i] << std::endl;
	}
	PoleSearch ps(&fcis,&gr,&ic,hfr,fcir,mbptr);
	for(int i = 1; i <= order; i++){
		
		std::vector<double> rexcitations=ps.refinePointsOrder(excitations[i-1],0.000001,i);
	
		excitations.push_back(std::vector<double>());

		std::cout << "REFINED "<<i<<"st ORDER\n";
		for(int j = 0; j < rexcitations.size();j++){
			if(j==0 || std::abs(rexcitations[j]-rexcitations[j-1])>0.0001){
				excitations[i].push_back(rexcitations[j]);
				std::cout << j << " " <<  rexcitations[j]<<std::endl;
			}
		}
	}
	excitations.push_back(ps.refinePoints(excitations[excitations.size()-1],0.000001));
	std::cout << "REFINED FULL\n";


	for(int i = 0; i < excitations[excitations.size()-1].size();i++){
		std::cout << i << " " <<  excitations[excitations.size()-1][i]<<std::endl;
	}

	std::ofstream eefile("Output/eevals");
	std::stack<double> line;
	for(int i = 0; i < excitations[excitations.size()-1].size();i++){
		line.push(excitations[excitations.size()-1][i]);
		for(int j =excitations.size()-2;j>=0;j--){
			double nearest=std::abs(line.top()-excitations[j][0]);
			int nearesti=0;
			for(int n = 1; n < excitations[j].size();n++){
				double diff=std::abs(line.top()-excitations[j][n]);
				if(diff<nearest){
					nearest=diff;
					nearesti=n;
				}
			}
			line.push(excitations[j][nearesti]);
		}
		while(!line.empty()){
			eefile << line.top() << " ";
			line.pop();
		}
		eefile<<std::endl;
	}
	eefile.close();
	

	for(int i = 1; i < fcir.eigenvalues.rows();i++){
		std::cout << std::setprecision(6) <<fcir.eigenvalues(i,0)-fcir.eigenvalues(0,0) << std::endl;
	}
	Matrix oEdiff=Matrix::Zero(hfr.operators.size(),hfr.operators.size());
	for(int i = 0; i < hfr.operators.size(); i++){
		int norbs=SlaterDet::codes.norbs;
		oEdiff(i,i)=hfr.E(hfr.operators[i].j%norbs,0)-hfr.E(hfr.operators[i].i%norbs,0);
	}
	std::cout << "SELF ENERGIES CLOSET EIGEN\n";

	/*
	PoleSearch pls(&fcis, & gr, &ic, hfr, fcir, mbptr);
	double mid=1.04751;
	double spread=0.01;
	int steps=500;
	//pls.scan(mid-spread/2.0, spread/100.0,100, "Output/scan");
	pls.scan(5.0/(double)steps,5.0/(double)steps,steps,"Output/scan");
	*/

	fcis.cleanup();
	SlaterDet::cleanStrings();

	libint2::finalize();


	return 0;
}
