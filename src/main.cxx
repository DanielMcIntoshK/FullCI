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
	RPASolver tdhfs(ic);
	RPASolver::RPAResults tdhfr=tdhfs.RPACalc(hfr);
	RPASolver::RPAResults tdhfr2=tdhfs.RPACalcSingTrip(hfr);

	std::cout << "RPA EXCITATION ENERGIES:\n";
	for(int i = 0; i < tdhfr.EE.size();i++){
		std::cout << tdhfr.EE[i]<<" " <<tdhfr2.EE[i]<<std::endl;
	}

	/*
	RPASolver::HRPAResults hrpar=tdhfs.HRPACalc(hfr,0.0001,1000,true);
	std::cout << "HRPA EXCITATION ENERGIES:\n";
	for(int i = 0; i < hrpar.EE.size(); i++){
		std::cout << hrpar.EE[i]<<" " << std::endl;
	}
	*/
	double sumHF=0.0;
	for(int i = 0; i < hfr.nelec; i++){
		sumHF+=2*hfr.E(i,0);
	}
	double test=0.0;
	for(int s1=0; s1 < 2;s1++){
	for(int i = 0; i < hfr.nelec; i++){
		for(int s2=0; s2<2; s2++){
		for(int j = 0; j < hfr.nelec; j++){
			test+=-.5*(ic.mov(i,i,j,j));
			if(s1==s2) test+=.5*(ic.mov(i,j,i,j));
		}}
	}}
	std::cout << "TEST: " << test+sumHF << std::endl;
	
	int ani=hfr.nelec-1,cre=hfr.nelec;
	double sumHF2=sumHF,test2=0.0;
	sumHF2-=hfr.E(ani,0)+hfr.E(cre,0);
	std::vector<int> occ;
	for(int s=0; s <2; s++){
		for(int i = 0; i < hfr.nelec; i++){
			if(s*hfr.norbs+i==ani)continue;
			occ.push_back(s*hfr.norbs+i);
		}
	}
	occ.push_back(cre);
	for(int i = 0; i < occ.size(); i++){
		for(int j = 0; j < occ.size(); j++){
			int ii=occ[i]%hfr.norbs, jj=occ[j]%hfr.norbs;
			int is=occ[i]/hfr.norbs, js=occ[j]%hfr.norbs;

			test2+=-.5*(ic.mov(ii,ii,jj,jj));
			if(is==js) test2+=.5*(ic.mov(ii,jj,ii,jj));
		}
	}
	

	std::cout << "TEST2: " << test2+sumHF2<<std::endl;
	std::cout << "DIFF: " << test2-test << std::endl;
	std::cout << "AC: " << ic.mov(ani,ani,cre,cre)-ic.mov(ani,cre,cre,ani)<<std::endl;
	std::cout << "RATIO: " << (test2-test)/(ic.mov(ani,ani,cre,cre)-ic.mov(ani,cre,cre,ani));

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

	std::ofstream gnum("Output/gnum.vals");
	gnum << "Greens Function Numerically Derived\n";
	for(int i = 0; i <= 5; i++){
		gnum<< i << std::endl;
		for(int p = 2; p <=6; p++){
			for(int q = 2; q <=6; q++){
				gnum<<greensnum[i](p,q)<<" ";
			}
			gnum<<std::endl;
		}
	}
	gnum.close();
	std::ofstream grec("Output/grec.vals");
	gnum << "Greens Function Algebreic Derived\n";
	for(int i = 0; i <= 5; i++){
		grec<< i << std::endl;
		for(int p = 2; p <=6; p++){
			for(int q = 2; q <=6; q++){
				grec<<greens[i](p,q)<<" ";
			}
			grec<<std::endl;
		}
	}
	grec.close();
	std::ofstream mnum("Output/mnum.vals");
	gnum << "Self Energy Numerically Derived\n";
	for(int i = 0; i <= 5; i++){
		mnum<< i << std::endl;
		for(int p = 2; p <=6; p++){
			for(int q = 2; q <=6; q++){
				mnum<<greensnum[i](p,q)<<" ";
			}
			mnum<<std::endl;
		}
	}
	mnum.close();
	std::ofstream mrec("Output/mrec.vals");
	gnum << "Self Energy Algebreic Derived\n";
	for(int i = 0; i <= 5; i++){
		mrec<< i << std::endl;
		for(int p = 2; p <=6; p++){
			for(int q = 2; q <=6; q++){
				mrec<<greens[i](p,q)<<" ";
			}
			mrec<<std::endl;
		}
	}
	mrec.close();
	

	//std::vector<double> excitations;
	std::vector<std::vector<double>> excitations;
	excitations.push_back(tdhfr.EE);

	//excitations[0].resize(tdhfr.EE.rows());
	std::cout << "RPA VALS:\n"; 
	for(int i = 0; i < excitations[0].size(); i++){
		//excitations[0][i]=tdhfr.EE(i,0);
		std::cout << i<< " " << excitations[0][i] << std::endl;
	}
	PoleSearch ps(&fcis,&gr,&ic,hfr,fcir,mbptr);
	for(int i = 1; i <= 4; i++){
		
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

	bool outputEE=true;
	if(outputEE){
		std::ofstream eefile("Output/eevalsN2");
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
	}

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
