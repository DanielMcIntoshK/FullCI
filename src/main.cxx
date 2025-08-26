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
#include "Diagram.h"
#include "SOPPA.h"

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

	FullCISolver fcis;
	FullCISolver::FCIResults fcir =fcis.fci(ic,hfr,1.0);

	std::cout << "FCI E: " << fcir.eigenvalues(0,0) +hfr.enuc << std::endl;
	
	double E=.2;
	std::cout << "ETEST: " << hfr.E(hfr.operators[0].i,0) << " " << hfr.E(hfr.operators[0].j,0)<<std::endl;

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

	int order=20;
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
	
	/*
	PoleSearch pls(&fcis,&gr,&ic,hfr,fcir,mbptr);
	double E_s=0.0;
	double E_f=1.5;
	int steps=10000;
	double dE=(E_f-E_s)/((double)steps);	
	
	pls.mapSelfEnergy(E_s,dE,steps,"Output/new_Mmap");
	*/

	//CALCULATE GREENS FUNCTION RECURSIVELY
	std::cout << "COMPUTING GREENS FUNCTION RECURSIVELY\n";
	int recursiveorder=4;
	std::vector<Matrix> greens=fcis.recursivegreen(recursiveorder,E,hfr,mbptr);
	Matrix rg=greens[0];
	std::cout << "0th order: \n";
	std::cout << greens[0] << std::endl<< std::endl << " " << rg(1,2) << std::endl;
	for(int i = 1; i <=recursiveorder; i++){
		rg+=greens[i];
		std::cout << i << "th order:\n";
		std::cout << greens[i] << std::endl<<std::endl;
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
		//std::cout << i << " " << max << " " << ri << " " << ci << " " << rg(ri,ci)<<" " << greensfull(ri,ci)<<std::endl;
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
	for(int i = 0; i <= recursiveorder; i++){
		val+= greens[i](p,q);
		std::cout << i << ": " << val <<  " " << greens[i](p,q) << " " << greens[i](q,p) <<  std::endl;
	}

	//CALCULATE GREENS FUNCTION NUMERICALLY
	LambdaDeriv ld;
	std::vector<double> grid;
	double gridspace=0.1;
	grid.push_back(0.0);
	int lambdaorder = 2;
	for(int i = 1; i <=lambdaorder/2; i++){
		grid.push_back(gridspace*i);
		grid.push_back(-gridspace*i);
	}
	std::vector<Matrix> greensnum=ld.ComputeGreensNumerical(lambdaorder,0.0,grid,hfr,ic,avocc,E);
	
	std::cout << "GREENS NUMERICAL\n";
	Matrix greensNum=greensnum[0];
	for(int r = 0; r < greensNum.rows(); r++){
	for(int c = 0; c < greensNum.cols(); c++){
		if(std::abs(greensnum[0](r,c))<0.000000000001)greensnum[0](r,c)=0.0;
	}
	}
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
	//return 0;

	//SELF ENERGY CALCULATIONS
	Matrix Mfull=gr.ComputeSelfEnergy(E,hfr,greensfull);

	
	//std::vector<Matrix> senum=gr.ComputeSelfEnergies(E,lambdaorder,hfr,greensnum);
	std::vector<Matrix> senum=ld.ComputeSelfEnergyNumerical(lambdaorder,0.0,grid,hfr,ic,avocc,E);
	std::cout << "SELF ENERGY NUMERICAL\n";
	Matrix seSum=Matrix::Zero(senum[1].rows(),senum[1].cols());;
	for(int i = 1; i < senum.size();i++){
		std::cout << i<<"th:\n"<<senum[i]<<std::endl<<std::endl;
		seSum+=senum[i];
		
		/*
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
		//std::cout << i << " " << max << " "<<ri<<" " << ci<<" "<<Mfull(ri,ci) << " " << seSum(ri,ci)<< std::endl;
		*/
	}
	std::cout << "SUM\n" << seSum<<std::endl<<std::endl;

	std::cout << "OPERATORS:\n";
	for(int i = 0; i < hfr.operators.size(); i++){
		std::cout << i<<" "<<hfr.operators[i].i << " " << hfr.operators[i].j<<std::endl;
	}

	std::cout << "G0:\n" << hfr.getG0(E) << std::endl<<std::endl;

	std::cout << greens[1].rows() << " " << greens[1].cols() << " " << hfr.operators.size()<<std::endl;;
	std::cout << "SELF ENERGY RECURSIVE\n"<<hfr.operators.size()<<std::endl;;
	std::vector<Matrix> serecursive=gr.ComputeSelfEnergies(E,recursiveorder,hfr,greens);

	for(int i = 1; i < senum.size() && i < serecursive.size(); i++){
		Matrix diff=greens[i]-greensnum[i];
		std::cout << "CHECKING SE " << i << std::endl;
		for(int x = 0; x < diff.rows(); x++){
			for(int y = 0; y < diff.cols(); y++){
				if(std::abs(diff(x,y))>0.000000001){
					std::cout << "DIFF ERROR: " << x << ", " << y << ": "<<diff(x,y)<<std::endl;	
				}
				
			}
		}
	}
	
	int o1=3,o2=4;
	std::cout << "G_("<<hfr.operators[o1].i<< " " <<hfr.operators[o1].j
		<<","<<hfr.operators[o2].i<<" " <<hfr.operators[o2].j<<")="
		<< serecursive[2](o1,o2)<<std::endl;
	std::vector<double> gs;
	gs.resize(8);
	int i=hfr.operators[o1].i;
	int in=i%hfr.norbs;
	int j=hfr.operators[o2].i;
	int jn=i%hfr.norbs;
	int a=hfr.operators[o1].j;
	int an=a%hfr.norbs;
	int b=hfr.operators[o2].j;
	int bn=b%hfr.norbs;
	double testv=0.0,testv2=0.0;
	std::cout << "PHLIST:\n";
	for(int i = 0; i < hfr.phlist[0].size(); i++){
		std::cout << hfr.phlist[0][i]<<" ";
	}
	std::cout << "\n";
	for(int i = 0; i < hfr.phlist[1].size(); i++){
		std::cout << hfr.phlist[1][i]<<" ";
	}
	std::cout << "INFO: " << 
		"\ni: " << i/hfr.norbs << " " << i%hfr.norbs<< 
		"\nj: " << j/hfr.norbs << " " << j%hfr.norbs<< 
		"\na: " << a/hfr.norbs << " " << a%hfr.norbs<< 
		"\nb: " << b/hfr.norbs << " " << b%hfr.norbs<< std::endl; 
	for(int k = 0; k < hfr.phlist[0].size(); k++){
		if(hfr.phlist[0][k]==i || hfr.phlist[0][k]==j) continue;
		for(int c = 0; c<hfr.phlist[1].size(); c++){
			if(hfr.phlist[1][c]==a || hfr.phlist[1][c]==b) continue;
			int ks = hfr.phlist[0][k];
			int cs = hfr.phlist[1][c];
			int kn = ks%hfr.norbs,
			    cn = cs%hfr.norbs;
			double eps1=hfr.E(jn,0)+hfr.E(kn,0)-hfr.E(bn,0)-hfr.E(cn,0);
			double delt1=E+hfr.E(an,0)+hfr.E(bn,0)+hfr.E(cn,0)-
				hfr.E(in,0)-hfr.E(jn,0)-hfr.E(kn,0);
			double delt2=E+hfr.E(bn,0)-hfr.E(jn,0);
			double num=ic.movsymphys(j,ks,b,cs,-hfr.norbs)*
				ic.movsymphys(a,cs,i,ks,hfr.norbs);
			testv+=-num/(eps1*delt1*delt2);
			testv2+=-num/(eps1*delt1*delt2);
			std::cout << c<< " " << k << " " <<cn << " " << kn << " " <<cs/hfr.norbs << " " << ks/hfr.norbs << " " <<num << std::endl; 
		}
	}
	std::cout << "TEST Z1G1Z0: " << testv <<" " << testv2<<std::endl;
	for(int k = 0; k < hfr.phlist[0].size(); k++){
		for(int c = 0; c<hfr.phlist[1].size(); c++){
			int kn = hfr.phlist[0][k]%hfr.norbs,
			    cn = hfr.phlist[1][c]%hfr.norbs;
			double eps1=hfr.E(jn,0)+hfr.E(kn,0)-hfr.E(bn,0)-hfr.E(cn,0);
			double delt1=-E+hfr.E(an,0)-hfr.E(in,0);
			double delt2=-E+hfr.E(cn,0)-hfr.E(kn,0);
			double num=ic.movsymphys(j,k,b,c,hfr.norbs)*
				ic.movsymphys(j,c,b,k,hfr.norbs);
			testv+=-num/(eps1*delt1*delt2);
		}
	}

	if(a==b){
	for(int k = 0; k < hfr.phlist[0].size(); k++){
		for(int l = 0; l < hfr.phlist[0].size(); l++){
			for(int c = 0; c<hfr.phlist[1].size(); c++){
				int kn = hfr.phlist[0][k]%hfr.norbs,
				    ln = hfr.phlist[0][l]%hfr.norbs,
				    cn = hfr.phlist[1][c]%hfr.norbs;
				double denom=E+hfr.E(an,0)+hfr.E(cn,0)-hfr.E(kn,0)-hfr.E(ln,0);
				double num=ic.movsymphys(k,l,c,i,hfr.norbs)*
					ic.movsymphys(j,c,k,l,hfr.norbs);
				gs[0]+= -num/(denom*2);
			}
		}
	}
	}
	else{ gs[0]=0.0;}
	gs[1]=0.0;
	for(int k = 0; k < hfr.phlist[0].size(); k++){
		for(int c = 0; c<hfr.phlist[1].size(); c++){
			int kn = hfr.phlist[0][k]%hfr.norbs,
			    cn = hfr.phlist[1][c]%hfr.norbs;
			double denom=E+hfr.E(bn,0)+hfr.E(cn,0)-hfr.E(in,0)-hfr.E(kn,0);
			double num=ic.movsymphys(k,a,b,c,hfr.norbs)*
				ic.movsymphys(j,c,k,i,hfr.norbs);
			gs[2]+=-num/(denom);
		}
	}
	for(int k = 0; k < hfr.phlist[0].size(); k++){
		for(int c = 0; c<hfr.phlist[1].size(); c++){
			int kn = hfr.phlist[0][k]%hfr.norbs,
			    cn = hfr.phlist[1][c]%hfr.norbs;
			double denom=E+hfr.E(an,0)+hfr.E(cn,0)-hfr.E(jn,0)-hfr.E(kn,0);
			double num=ic.movsymphys(j,k,c,i,hfr.norbs)*
				ic.movsymphys(c,a,b,k,hfr.norbs);
			gs[3]+=-num/(denom);
		}
	}
	for(int k = 0; k < hfr.phlist[0].size(); k++){
		for(int l = 0; l<hfr.phlist[0].size(); l++){
			int kn = hfr.phlist[0][k]%hfr.norbs,
			    ln = hfr.phlist[0][l]%hfr.norbs;
			double denom=E+hfr.E(an,0)+hfr.E(bn,0)-hfr.E(kn,0)-hfr.E(ln,0);
			double num=ic.movsymphys(k,l,b,i,hfr.norbs)*
				ic.movsymphys(j,a,k,l,hfr.norbs);
			gs[4]+=+num/(denom*2);
		}
	}
	for(int d = 0; d < hfr.phlist[1].size(); d++){
		for(int c = 0; c<hfr.phlist[1].size(); c++){
			int dn = hfr.phlist[1][d]%hfr.norbs,
			    cn = hfr.phlist[1][c]%hfr.norbs;
			double denom=E+hfr.E(cn,0)+hfr.E(dn,0)-hfr.E(jn,0)-hfr.E(in,0);
			double num=ic.movsymphys(j,a,c,d,hfr.norbs)*
				ic.movsymphys(c,d,b,i,hfr.norbs);
			gs[5]+=+num/(denom*2);
		}
	}
	gs[6]=0.0;
	if(a==b){
	for(int k = 0; k < hfr.phlist[0].size(); k++){
		for(int d = 0; d < hfr.phlist[1].size(); d++){
			for(int c = 0; c<hfr.phlist[1].size(); c++){
				int kn = hfr.phlist[0][k]%hfr.norbs,
				    dn = hfr.phlist[1][d]%hfr.norbs,
				    cn = hfr.phlist[1][c]%hfr.norbs;
				double eps1=hfr.E(in,0)+hfr.E(kn,0)-hfr.E(cn,0)-hfr.E(dn,0);
				double eps2=hfr.E(jn,0)+hfr.E(kn,0)-hfr.E(cn,0)-hfr.E(dn,0);
				double num1=E+hfr.E(an,0)+hfr.E(cn,0)+hfr.E(dn,0)-hfr.E(in,0)
					-hfr.E(jn,0)-hfr.E(kn,0);
				double ints=ic.movsymphys(j,k,c,d,hfr.norbs)*
					ic.movsymphys(c,d,i,k,hfr.norbs);
				
				double num=num1*ints;
				double denom=eps1*eps2;

				gs[7]+= -num/(denom*2);
			}
		}
	}
	}
	else{ gs[7]=0.0;}
	double sum =0.0;
	for(int i = 0; i < gs.size(); i++){
		std::cout << i << " " << gs[i]<<std::endl;
		sum+=gs[i];
	}
	std::cout << "SUM: " << sum<<std::endl;
	
	seSum=Matrix::Zero(seSum.rows(),seSum.cols());
	for(int i = 1; i < serecursive.size();i++){
		std::cout << i<<"th:\n"<<serecursive[i]<<std::endl<<std::endl;
		seSum+=serecursive[i];
		Matrix sepr=seSum;
	}
	std::cout << "SUM\n"<<seSum<<std::endl<<std::endl;

	std::cout << "Full M\n" << Mfull <<std::endl<<std::endl;

	int opk=0, opl=3;
	//opk=3; opl=4;
	opk=2; opl=6;
//	opk=12;opl=16;
//	opk=2; opl=16;
	//opk=2; opl=2;
	//opk=2; opl=3;

	for(int i =0; i < hfr.operators.size()/2; i++){
		for(int j = i+1; j < hfr.operators.size()/2; j++){
			if(hfr.operators[i].i == hfr.operators[j].i &&
					std::abs(greens[2](i,j))>0.00000001){
				std::cout << "POTENTIAL: " << i << " " << j << std::endl;
			}
		}
	}
	
	int opsize=hfr.norbs;
	int ii=hfr.operators[opl].i,ai=hfr.operators[opl].j,
	    ji=hfr.operators[opk].i,bi=hfr.operators[opk].j;

	int ik=ii%opsize,ak=ai%opsize,jk=ji%opsize,bk=bi%opsize;

	std::cout << ii << " " << ai << " : " << ji << " " << bi << std::endl;
	
	std::vector<int> irreduce;
	//Z0G2Z0---------------------------------------
	double ts=0.0;
	std::vector<double> dgv;
	Diagram dgtest(
		std::vector<dnode>{{line(1,bi)},{line(0,ji),line(2)},{line(1),line(3,ai)},{line(2,ii)}},
		std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
	double stest=dgtest.sumfull(&ic,E);
	irreduce.push_back(dgv.size());
	dgv.push_back(stest);
	ts+=stest;
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	dgtest=Diagram(
		std::vector<dnode>{{line(2,bi)},{line(3,ai),line(2)},{line(1),line(0,ji)},{line(1,ii)}},
		std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
	stest=dgtest.sumfull(&ic,E);
	irreduce.push_back(dgv.size());
	dgv.push_back(stest);
	ts+=stest;
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	dgtest=Diagram(
		std::vector<dnode>{{line(2,bi)},{line(0,ji),line(3,ai)},{line(1),line(1)},{line(2,ii)}},
		std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
	stest=dgtest.sumfull(&ic,E);
	dgv.push_back(stest);
	ts+=stest;
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	dgtest=Diagram(
		std::vector<dnode>{{line(1,bi)},{line(2),line(2)},{line(0,ji),line(3,ai)},{line(1,ii)}},
		std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
	stest=dgtest.sumfull(&ic,E);
	dgv.push_back(stest);
	ts+=stest;
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	dgtest=Diagram(
		std::vector<dnode>{{line(1,bi)},{line(2),line(3,ai)},{line(0,ji),line(1)},{line(2,ii)}},
		std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
	stest=dgtest.sumfull(&ic,E);
	dgv.push_back(stest);
	ts+=stest;
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	dgtest=Diagram(
		std::vector<dnode>{{line(2,bi)},{line(2),line(0,ji)},{line(3,ai),line(1)},{line(1,ii)}},
		std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
	stest=dgtest.sumfull(&ic,E);
	dgv.push_back(stest);
	ts+=stest;
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	if(hfr.operators[opk].j==hfr.operators[opl].j){
		dgtest=Diagram(
			std::vector<dnode>{{line(3,ai)},{line(2),line(0,ji)},{line(1),line(1)},{line(2,ii)}},
			std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(stest);
		ts+=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
		dgtest=Diagram(
			std::vector<dnode>{{line(3,ai)},{line(2),line(2)},{line(1),line(0,ji)},{line(1,ii)}},
			std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(stest);
		ts+=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	if(hfr.operators[opk].i==hfr.operators[opl].i){
		dgtest=Diagram(
			std::vector<dnode>{{line(1,ai)},{line(2),line(2)},{line(1),line(3,bi)},{line(0,ii)}},
			std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(stest);
		ts+=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
		dgtest=Diagram(
			std::vector<dnode>{{line(2,ai)},{line(2),line(3,bi)},{line(1),line(1)},{line(0,ii)}},
			std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(stest);
		ts+=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	if((hfr.operators[opk].j==hfr.operators[opl].j)&&
		(hfr.operators[opk].i==hfr.operators[opl].i)){
		dgtest=Diagram(
			std::vector<dnode>{{line(3,ai)},{line(2),line(2)},{line(1),line(1)},{line(0,ii)}},
			std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,true)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(stest);
		ts+=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
		dgtest=Diagram(
			std::vector<dnode>{{line(1),line(1)},{line(0),line(0)}},
			std::vector<resolvent>{resolvent(0,false)}, hfr);
		double ngv=dgtest.sumfull(&ic,E);
		double denomsq=1.0/(E+hfr.E(ai%hfr.norbs,0)-hfr.E(ii%hfr.norbs,0));
		stest=ngv*denomsq*denomsq;
		dgv.push_back(stest);
		ts+=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;

	}
	double z0g2z0sum=ts;
	std::cout << "Z0G2Z0: " << ts << std::endl;
	//Z0G1Z1-----------------------------------------------------
	double z0g1z1sum=0.0;
	dgtest=Diagram(
		std::vector<dnode>{{line(3,ai),line(2)},{line(2,bi)},{line(0),line(1,ji)},{line(0,ii)}},
		std::vector<resolvent>{resolvent(0,false),resolvent(1,true),resolvent(2,true)},hfr);
	stest=dgtest.sumfull(&ic,E);
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	irreduce.push_back(dgv.size());
	dgv.push_back(-stest);
	ts-=stest;
	z0g1z1sum-=stest;
	if(hfr.operators[opk].j==hfr.operators[opl].j){
		dgtest=Diagram(
			std::vector<dnode>{{line(2),line(2)},{line(3,ai)},{line(0),line(1,ji)},{line(0,ii)}},
			std::vector<resolvent>{resolvent(0,false),resolvent(1,true),resolvent(2,true)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(-stest);
		ts-=stest;
		z0g1z1sum-=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	if(hfr.operators[opk].i==hfr.operators[opl].i){
		dgtest=Diagram(
			std::vector<dnode>{{line(2),line(3,ai)},{line(2,bi)},{line(0),line(0)},{line(1,ii)}},
			std::vector<resolvent>{resolvent(0,false),resolvent(1,true),resolvent(2,true)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(-stest);
		ts-=stest;
		z0g1z1sum-=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	if((hfr.operators[opk].j==hfr.operators[opl].j)&&
		(hfr.operators[opk].i==hfr.operators[opl].i)){
		dgtest=Diagram(
			std::vector<dnode>{{line(2),line(2)},{line(3,ai)},{line(0),line(0)},{line(1,ii)}},
			std::vector<resolvent>{resolvent(0,false),resolvent(1,true),resolvent(2,true)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(-stest);
		ts-=stest;
		z0g1z1sum-=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	std::cout << "Z0G1Z1: " << z0g1z1sum << std::endl;
	//Z1G0Z1+----------------------------------------------
	double z1g0z1psum=0.0;
	dgtest=Diagram(
		std::vector<dnode>{{line(2,ai),line(3)},{line(3,bi)},{line(0,ii)},{line(0),line(1,ji)}},
		std::vector<resolvent>{resolvent(0,false),resolvent(1,true),resolvent(2,false)},hfr);
	stest=dgtest.sumfull(&ic,E);
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	irreduce.push_back(dgv.size());
	dgv.push_back(stest);
	ts+=stest;
	z1g0z1psum+=stest;
	if(hfr.operators[opk].j==hfr.operators[opl].j){
		dgtest=Diagram(
			std::vector<dnode>{{line(3),line(3)},{line(2,ai)},{line(0,ji)},{line(0),line(1,ii)}},
			std::vector<resolvent>{resolvent(0,false),resolvent(1,true),resolvent(2,false)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(stest);
		ts+=stest;
		z1g0z1psum+=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	if(hfr.operators[opk].i==hfr.operators[opl].i){
		dgtest=Diagram(
			std::vector<dnode>{{line(2,bi),line(3)},{line(3,ai)},{line(1,ii)},{line(0),line(0)}},
			std::vector<resolvent>{resolvent(0,false),resolvent(1,true),resolvent(2,false)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(stest);
		ts+=stest;
		z1g0z1psum+=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	if((hfr.operators[opk].j==hfr.operators[opl].j)&&
		(hfr.operators[opk].i==hfr.operators[opl].i)){
		dgtest=Diagram(
			std::vector<dnode>{{line(3),line(3)},{line(2,ai)},{line(1,ii)},{line(0),line(0)}},
			std::vector<resolvent>{resolvent(0,false),resolvent(1,true),resolvent(2,false)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(stest);
		ts+=stest;
		z1g0z1psum+=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	std::cout << "Z1G0+Z1: " << z1g0z1psum << std::endl;
	//Z1G0Z1- -----------------------------------------------
	double z1g0z1nsum=0.0;
	dgtest=Diagram(
		std::vector<dnode>{{line(1,ai),line(3)},{line(0,ii)},{line(3,bi)},{line(0),line(2,ji)}},
		std::vector<resolvent>{resolvent(0,false),resolvent(1,true,true),resolvent(2,false)},hfr);
	stest=dgtest.sumfull(&ic,E);
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	irreduce.push_back(dgv.size());
	dgv.push_back(stest);
	ts+=stest;
	z1g0z1nsum+=stest;
	std::cout << "Z1G0-Z1: " << z1g0z1nsum << std::endl;
	//Z1G1Z0------------------------------------------
	double z1g1z0sum=0.0;
	dgtest=Diagram(
		std::vector<dnode>{{line(3,bi)},{line(2,ai),line(3)},{line(1,ii)},{line(0,ji),line(1)}},
		std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,false)},hfr);
	stest=dgtest.sumfull(&ic,E);
	std::cout << "TESTING GRAPH: " << stest << std::endl;
	irreduce.push_back(dgv.size());
	dgv.push_back(-stest);
	ts-=stest;
	z1g1z0sum-=stest;
	if(hfr.operators[opk].j==hfr.operators[opl].j){
		dgtest=Diagram(
			std::vector<dnode>{{line(2,ai)},{line(3),line(3)},{line(1,ii)},{line(1),line(0,ji)}},
			std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,false)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(-stest);
		ts-=stest;
		z1g1z0sum-=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	if(hfr.operators[opk].i==hfr.operators[opl].i){
		dgtest=Diagram(
			std::vector<dnode>{{line(3,bi)},{line(2,ai),line(3)},{line(0,ii)},{line(1),line(1)}},
			std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,false)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(-stest);
		ts-=stest;
		z1g1z0sum-=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	if((hfr.operators[opk].j==hfr.operators[opl].j)&&
		(hfr.operators[opk].i==hfr.operators[opl].i)){
		dgtest=Diagram(
			std::vector<dnode>{{line(2,ai)},{line(3),line(3)},{line(0,ii)},{line(1),line(1)}},
			std::vector<resolvent>{resolvent(0,true),resolvent(1,true),resolvent(2,false)},hfr);
		stest=dgtest.sumfull(&ic,E);
		dgv.push_back(-stest);
		ts-=stest;
		z1g1z0sum-=stest;
		std::cout << "TESTING GRAPH: " << stest << std::endl;
	}
	std::cout << "Z1G1Z0: " << z1g1z0sum<<std::endl;

	double denom1_n = E+hfr.E(hfr.operators[opk].j%hfr.norbs,0)-
			hfr.E(hfr.operators[opk].i%hfr.norbs,0);
	double denom2_n = E+hfr.E(hfr.operators[opl].j%hfr.norbs,0)-
			hfr.E(hfr.operators[opl].i%hfr.norbs,0);
	if((hfr.operators[opk].j==hfr.operators[opl].j)&&
		(hfr.operators[opk].i==hfr.operators[opl].i)){
		dgtest=Diagram(
			std::vector<dnode>{{line(1),line(1)},{line(0),line(0)}},
			std::vector<resolvent>{resolvent(0,false),resolvent(0,false)},hfr);
		double renorm=dgtest.sumfull(&ic,E)/denom1_n;
		ts-=renorm;
		dgv.push_back(-ts);
	}
	std::cout << "GREENS: " << ts << " " << greens[2](opk,opl) <<" " << greensnum[2](opk,opl) << std::endl;
	std::cout << "DIFF: " << greens[2](opk,opl)-greensnum[2](opk,opl) << std::endl;

	dgtest=Diagram(
		std::vector<dnode>{{line(2,bi)},{line(3,ai),line(2)},{line(1),line(0,ji)},{line(1,ii)}},
		std::vector<resolvent>{
			resolvent(0,true),resolvent(1,true,true,false,{ii,ji,ai,bi}),resolvent(2,true)},hfr);
	stest=dgtest.sumfull(&ic,E);
	ts+=stest;
	
	ts=0;
	for(int i = 0; i < dgv.size(); i++){
		ts+=dgv[i];
	}
	std::cout << ts << std::endl;

	std::vector<Diagram> sedgs;
	
	sedgs=std::vector<Diagram>{
		Diagram({{line(2,bi)},{line(0,ji),line(2)},{line(1),line(3,ai)},{line(1,ii)}},
			{resolvent(1,true)},hfr),
		Diagram({{line(1,bi)},{line(3,ai),line(2)},{line(1),line(0,ji)},{line(2,ii)}},
			{resolvent(1,true)},hfr),
		Diagram({{line(2,bi)},{line(3,ai),line(0,ji)},{line(1),line(1)},{line(2,ii)}},
			{resolvent(1,true)},hfr),
		Diagram({{line(1,bi)},{line(2),line(2)},{line(0,ji),line(3,ai)},{line(1,ii)}},
			{resolvent(1,true)},hfr)
		};
	if(hfr.operators[opk].j==hfr.operators[opl].j){
		sedgs.push_back(
			Diagram({{line(3,ai)},{line(0,ji),line(2)},{line(1),line(1)},{line(2,ii)}},
				{resolvent(1,true)},hfr));
		sedgs.push_back(
			Diagram({{line(3,ai)},{line(2),line(2)},{line(0,ji),line(1)},{line(1,ii)}},
				{resolvent(1,true,false,true),
				resolvent(1,false,false,false,{ai,ii}),
				resolvent(1,false,false,false,{ai,ji})},hfr));
		std::cout << sedgs.back().sumfull(&ic,E)<<std::endl;
	}
	if(hfr.operators[opk].i==hfr.operators[opl].i){
		sedgs.push_back(
			Diagram({{line(1,bi)},{line(2),line(2)},{line(1),line(3,ai)},{line(0,ii)}},
				{resolvent(1,true)},hfr));
		sedgs.push_back(
			Diagram({{line(2,bi)},{line(2),line(3,ai)},{line(1),line(1)},{line(0,ii)}},
				{resolvent(1,true,false,true),
				resolvent(1,false,false,false,{ai,ii}),
				resolvent(1,false,false,false,{ai,ji})},hfr));
	}
	double extest=0.0, diff=0.0;
	for(int k = 0; k < hfr.phlist[0].size(); k++){
		int kop=hfr.phlist[0][k];
		for(int c = 0; c < hfr.phlist[1].size(); c++){
		for(int d = 0; d < hfr.phlist[1].size(); d++){
			int cop=hfr.phlist[1][c];
			int dop=hfr.phlist[1][d];

			double ints1=ic.movsymphys(ji,kop,cop,dop,hfr.norbs),
			       ints2=ic.movsymphys(cop,dop,ii,kop,hfr.norbs);

			double denom=-E-hfr.E(kop%hfr.norbs,0)-hfr.E(ai%hfr.norbs,0)
				+hfr.E(cop%hfr.norbs,0)+hfr.E(dop%hfr.norbs,0);

			double eps1=hfr.E(ji%hfr.norbs,0)+hfr.E(kop%hfr.norbs)-
				    hfr.E(cop%hfr.norbs,0)-hfr.E(dop%hfr.norbs,0),
			       eps2=hfr.E(ii%hfr.norbs,0)+hfr.E(kop%hfr.norbs,0)-
				    hfr.E(cop%hfr.norbs,0)-hfr.E(dop%hfr.norbs,0);
			extest+=ints1*ints2/denom;
			diff+=ints1*ints2/(eps1*eps2*denom);
		}}
	}
	double num=(E+hfr.E(ai%hfr.norbs,0)-hfr.E(ii%hfr.norbs,0))*
		  (E+hfr.E(bi%hfr.norbs,0)-hfr.E(ji%hfr.norbs,0));
	std::cout << sedgs[sedgs.size()-1].sumfull(&ic,E) << " " << -extest/2<< " "
	       << sedgs[sedgs.size()-1].sumfull(&ic,E)+extest/2<< " " <<
	       num*diff/2 << std::endl;

	double smtestdiag=0.0;
	for(int i = 0; i < sedgs.size(); i++){
		smtestdiag+=sedgs[i].sumfull(&ic,E);
	}
	double sediff=serecursive[2](opk,opl)-smtestdiag;
	std::cout << "SELFENERGY TEST: " << smtestdiag<< " " << serecursive[2](opk,opl)<<" "<<
		serecursive[2](opk,opl)-smtestdiag << std::endl;

	opk=3; opl=4;
	opk+=hfr.operators.size()/2;
	ii=hfr.operators[opl].i,ai=hfr.operators[opl].j,
	    ji=hfr.operators[opk].j,bi=hfr.operators[opk].i;
	
	sedgs=std::vector<Diagram>{
		Diagram({{line(2,bi),line(1)},{line(0),line(3,ai)},{line(1,ji)},{line(0,ii)}},
			{resolvent(0,false)},hfr),
		Diagram({{line(1),line(1)},{line(3,ai),line(2,bi)},{line(0,ji)},{line(0,ii)}},
			{resolvent(0,false)},hfr),
		Diagram({{line(3,ai),line(2,bi)},{line(0),line(0)},{line(1,ji)},{line(1,ii)}},
			{resolvent(0,false)},hfr),
		Diagram({{line(3,ai),line(1)},{line(0),line(2,bi)},{line(0,ji)},{line(1,ii)}},
			{resolvent(0,false)},hfr)
		};

	smtestdiag=0.0;
	for(int i = 0; i < sedgs.size(); i++){
		smtestdiag-=sedgs[i].sumfull(&ic,E);
	}
	std::cout << "SELFENERGY TEST: " << smtestdiag<< " " << serecursive[2](opk,opl)<<" "<<
		serecursive[2](opk,opl)-smtestdiag << std::endl;
	std::cout << senum[2](opk,opl)<<std::endl;
	
	std::cout << hfr.operators.size() << std::endl;
	for(int i = 0; i < hfr.operators.size(); i++){
		std::cout << i<< " "<<hfr.operators[i].i << " " << hfr.operators[i].j<<std::endl;
	}

	Matrix se2test=Matrix(senum[2].rows(), senum[2].cols());
	Matrix se4test=Matrix(senum[2].rows(), senum[2].cols());
	for(int i = 0; i < se2test.rows(); i++){
	for(int j = 0; j < se2test.cols(); j++){
		se2test(i,j)=0.0;
		se4test(i,j)=0.0;

		int ii=hfr.operators[i].i,ai=hfr.operators[i].j,
	    		ji=hfr.operators[j].i,bi=hfr.operators[j].j;
		bool over=i>=(se2test.rows()/2);
		if(over){
			std::swap(ii,ji);
			std::swap(ai,bi);
			continue;
		}
		if(ai==bi){
			Diagram dg=Diagram(
				{{line(3,ai)},{line(2),line(2)},{line(0,ji),line(1)},{line(1,ii)}},
				{resolvent(1,true,over,true),
				resolvent(1,false,over,false,{ai,ii}),
				resolvent(1,false,over,false,{ai,ji})},hfr);
			se2test(i,j)+=dg.sumfull(&ic,E);
		}
		if(ii==ji){
			Diagram dg=Diagram(
				{{line(2,bi)},{line(2),line(3,ai)},{line(1),line(1)},{line(0,ii)}},
				{resolvent(1,true,over,true),
				resolvent(1,false,over,false,{ai,ii}),
				resolvent(1,false,over,false,{bi,ii})},hfr);
			se2test(i,j)+=dg.sumfull(&ic,E);
		}
		
		if(i==3&&j==6){
			std::cout << ii << " " << ji << " " << ai << " " << bi << std::endl;
		Diagram dg=Diagram(
			{{line(2,bi)},{line(5,ai),line(2)},{line(1),line(1)},{line(4),line(4)},
			{line(0,ji),line(3)},{line(3,ii)}},
			{resolvent(1,false,false,false,{ji,bi}),
			resolvent(2,true,false,false),
			resolvent(3,false,false,false,{ii,ai})},hfr);
		std::cout << "MULT: " << dg.getmultiplicity() << std::endl;
		se4test(i,j)+=dg.sumfull(&ic,E);
		dg=Diagram(
			{{line(4,bi)},{line(2),line(2)},{line(1),line(0,ji)},{line(4),line(5,ai)},
			{line(3),line(3)},{line(1,ii)}},
			{resolvent(1,false,false,false,{ji,bi}),
			resolvent(2,true,false,false),
			resolvent(3,false,false,false,{ii,ai})},hfr);
		se4test(i,j)+=dg.sumfull(&ic,E);
		}
	}}
	std::cout << "G0: " << std::endl << greens[0] << std::endl<<std::endl;
	std::cout << "JUST PSUEDO: " << std::endl << se2test << std::endl<<std::endl;
	std::cout << "REDUCABLE 4th: " << std::endl << se2test*greens[0]*se2test << std::endl<<std::endl;
	//std::cout << "REDUCABLE 4th: " << std::endl << se2test*greens[0] << std::endl;
	std::cout << "4th order S: " <<std::endl<<se4test<<std::endl;

	//std::cout << "SOPPA TEST:\n";
	//SOPPASolver soppasolve(hfr, &ic);
	//Matrix foppa=soppasolve.FOPPA(E);
	//Matrix soppa=soppasolve.SOPPA(E);
	//std::cout << foppa << std::endl <<std::endl;

	//std::cout << serecursive[1] << std::endl;

	return 0;
	/*
	std::ofstream gnum("Output/gnum.vals");
	gnum << "Greens Function Numerically Derived\n";
	for(int i = 0; i <= 5; i++){
		gnum<< i << std::endl;
		for(int p = 2; p <=6; p++){
			for(int q= 2; q <=6; q++){
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
	*/
	

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
	for(int i = 1; i <= order; i++){
		std::cout << "REFINING " << i << "st ORDER\n";
		
		std::vector<double> rexcitations=ps.refinePointsOrder(excitations[i-1],0.000001,i);
	
		for(int i = 0; i < rexcitations.size(); i++){
			std::cout << i << " " << rexcitations[i] << std::endl;
		}
		//excitations.push_back(std::vector<double>());
		
		std::vector<double> nonduplicate;
		nonduplicate.push_back(rexcitations[0]);
		for(int j = 1; j < rexcitations.size();j++){
			if(std::abs(rexcitations[j]-rexcitations[j-1])>0.00001){
				nonduplicate.push_back(rexcitations[j]);
			}
		}
		excitations.push_back(nonduplicate);
		/*
		std::cout << "REFINED "<<i<<"st ORDER\n";
		for(int j = 0; j < rexcitations.size();j++){
			if(j==0 || std::abs(rexcitations[j]-rexcitations[j-1])>0.0001){
				excitations[i].push_back(rexcitations[j]);
				std::cout << j << " " <<  rexcitations[j]<<std::endl;
			}
		}
		*/
		
		//excitations.push_back(rexcitations);
	}
	excitations.push_back(ps.refinePoints(excitations[excitations.size()-1],0.000001));
	std::cout << "REFINED FULL\n";

	for(int i = 0; i < excitations[excitations.size()-1].size();i++){
		std::cout << i << " " <<  excitations[excitations.size()-1][i]<<std::endl;
	}

	bool outputEE=true;
	if(outputEE){
		std::ofstream eefile("Output/eevalsX2");
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
