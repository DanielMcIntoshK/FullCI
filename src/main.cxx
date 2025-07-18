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

		/*
		Matrix G0itemp=hfr.getIndependentG0i();
		sepr=G0itemp-sepr;
		gr.TransformEigen(sepr,hfr);

		//std::cout << i << "th Eigen:\n" << sepr << std::endl << std::endl;

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
		*/
	}
	std::cout << "SUM\n"<<seSum<<std::endl<<std::endl;

	std::cout << "Full M\n" << Mfull <<std::endl<<std::endl;

	std::cout << "GV TEST\n";
	StringMap & sm = SlaterDet::codes;
	for(int i = 0; i < sm.strs.size(); i++){
		std::cout << i << " ";
		for(int j = 0; j < sm.codeblklen; j++){
			for(int k = 0; k < 8; k++){
				if(k+8*j>=sm.norbs) break;

				bool bon=sm.strs[i][j] &(1<<k);
				if(bon) std::cout << "1";
				else std::cout << "0";
			}
		}
		std::cout << std::endl;
	}
	i = 3, j =4, a=5, b=6;
	std::cout << "|0>G1<5|: TARGET: "<<-ic.movsym(i,a,j,b,sm.norbs)/(E*(E+hfr.E(a,0)+hfr.E(b,0)-hfr.E(i,0)-hfr.E(j,0)))<<" ACTUAL: " << fcis.Gm_p[1](0,5)<<std::endl;
	i = 3, j =4, a=5, b=6;
	std::cout << "|3>G1<2|: TARGET: "<<-ic.movsym(a,i,j,b,sm.norbs)/((E+hfr.E(a,0)-hfr.E(i,0))*(E+hfr.E(b,0)-hfr.E(j,0)))<<" ACTUAL: " << fcis.Gm_p[1](3,2)<<std::endl;
	i = 4, j =4, a=5, b=6;
	std::cout << "|1>G1<2|: TARGET: "<<-ic.movsym(a,i,j,b,sm.norbs)/((E+hfr.E(a,0)-hfr.E(j,0))*(E+hfr.E(b,0)-hfr.E(j,0)))<<" ACTUAL: " << fcis.Gm_p[1](1,2)<<std::endl;
	i = 3, j =3, a=5, b=6;
	double sumk=-ic.movsym(a,i,j,b,sm.norbs);
	//for(int k = 0; k < hfr.phlist[0].size(); k++){
	for(int k = 0; k < hfr.nelec; k++){
		int ki = hfr.phlist[0][k];
		if(ki==i)continue;
		sumk+=-ic.movsym(a,ki,ki,b,sm.norbs);
	}
	std::cout <<sumk << std::endl;
	std::cout << "|3>G1<4|: TARGET: "<<-ic.movsym(a,i,j,b,sm.norbs)/((E+hfr.E(a,0)-hfr.E(i,0))*(E+hfr.E(b,0)-hfr.E(j,0)))<<" ACTUAL: " << fcis.Gm_p[1](3,4)<<std::endl;

	int opk=0, opl=3;
	opk=3; opl=4;
	std::cout << "OPP: " << opk << " " << hfr.operators[opk].i << " " << hfr.operators[opk].j << std::endl;
	std::vector<int> valids;
	for(int i = 0; i < fcis.za[1].rows(); i++){
		int strcnt = sm.strs.size();
		int alphaidx=i%strcnt, betaidx=i/strcnt;
		unsigned char * alphacpy=sm.cpystrs[0],*betacpy=sm.cpystrs[1];
		memcpy(alphacpy, sm.strs[alphaidx],sm.codeblklen);
		memcpy(betacpy,sm.strs[betaidx],sm.codeblklen);
		double sign = fcis.opOnSlater(hfr.operators[opk].adjoint(),alphacpy, betacpy);

		if(sign!=0.0)valids.push_back(i);
	}

	int singleexcite=-1;
	std::cout << "OPP: " << opl << " " << hfr.operators[opl].i << " " << hfr.operators[opl].j << std::endl;
	for(int i = 0; i < fcis.za[0].rows(); i++){
		if(std::abs(fcis.za[0](i,opl))>0.0000000001){
			std::cout << sm.printstrphab(i) << " " << sm.printstrab(i)<< " " << fcis.za[0](i,opl)<<std::endl;
			singleexcite=i;
		}
	}

	for(int i = 0; i < valids.size(); i++){
		SlaterDet se=SlaterDet(singleexcite), vi=SlaterDet(valids[i]);
		SlaterCompare scdet=compareSlaterDet(se,vi);

		if(scdet.diff.size()<=4){
			
			std::string prstr=sm.printstrphab(valids[i]);
			std::cout << prstr;

			int loopcount=prstr.size()/4;
			if(loopcount%4!=0)loopcount-=1;
		       
			for(int i = 0; i < 6-loopcount;i++){
				std::cout << "\t";
			}
			std::cout << sm.printstrab(valids[i]) << 
				" " << fcis.za[1](valids[i],opk) << " " << scdet.diff.size() << std::endl;
		}
	}
	std::cout << "12345678901234567890\n";
	std::cout << "\ta\n";
	std::cout << ic.movsymphys(5,12,2,9,sm.norbs)/(hfr.E(2,0)+hfr.E(9%7,0)-hfr.E(5,0)-hfr.E(12%7,0))<< std::endl;
	std::cout << ic.movsymphys(12,5,2,9,sm.norbs)/(hfr.E(2,0)+hfr.E(9%7,0)-hfr.E(5,0)-hfr.E(12%7,0))<< std::endl;
	double sumtest = 0.0;
	//for(int i = 0; i < fcis.Gm_p[0].rows(); i++){
	for(int i = 0; i < valids.size(); i++){
		SlaterDet basese=SlaterDet(0);

		SlaterDet se=SlaterDet(singleexcite), vi=SlaterDet(valids[i]);
		SlaterCompare scdet=compareSlaterDet(se,vi);
		SlaterCompare scdet1=compareSlaterDet(se,basese), scdet2 = compareSlaterDet(vi,basese);
		
		if(scdet.diff.size()<=4){

			sumtest+=fcis.za[0](singleexcite,opl)*fcis.Gm_p[1](singleexcite,valids[i])*fcis.za[1](valids[i],opk);
			//sumtest+=fcis.za[0](j,opk)*fcis.Gm_p[1](j,i)*fcis.za[1](i,opl);
			//std::cout << sm.printstrphab(valids[i]) << " " <<fcis.za[0](singleexcite,opl) << " " << fcis.Gm_p[1](singleexcite,valids[i]) << " " 
			//	<< fcis.za[1](valids[i],opk) << std::endl;
		}
	}
	std::cout << sumtest << std::endl;
	std::cout << (fcis.za[0].transpose()*fcis.Gm_p[1]*fcis.za[1])(opl,opk)<<std::endl;

	std::cout << std::endl;
	double sumdiagram=0.0;
	int opsize=hfr.norbs;
	int ii=hfr.operators[opk].i,ai=hfr.operators[opk].j,
	    ji=hfr.operators[opl].i,bi=hfr.operators[opl].j;
	std::swap(ii,ji); std::swap(ai,bi);

	int ik=ii%opsize,ak=ai%opsize,jk=ji%opsize,bk=bi%opsize;


	
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			//if(ki==ii||ki==ji) continue;
			//if(ci==ai||ci==bi) continue;
			double int1=ic.movsymphys(ci,bi,ji,ki,opsize),
			       int2=ic.movsymphys(ii,ki,ci,ai,opsize);
			double denom1=1.0/(hfr.E(jk,0)+hfr.E(kk,0)-hfr.E(ck,0)-hfr.E(bk,0)),
			       denom2=1.0/(E+hfr.E(ck,0)+hfr.E(ak,0)+hfr.E(bk,0)-hfr.E(jk,0)-hfr.E(ik,0)-hfr.E(kk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));

			sumdiagram-=int1*int2*(denom1*denom2*denom3);
		}
	}
	std::cout << "Z1G1Z0: " << sumdiagram << std::endl;
	sumdiagram=0.0;
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			//if(ki==ii||ki==ji) continue;
			//if(ci==ai||ci==bi) continue;
			double int1=ic.movsymphys(ci,bi,ji,ki,opsize),
			       int2=ic.movsymphys(ii,ki,ci,ai,opsize);
			double denom1=1.0/(E+hfr.E(ak,0)-hfr.E(ik,0)),
			       denom2=1.0/(E+hfr.E(ak,0)+hfr.E(bk,0)+hfr.E(ck,0)-hfr.E(ik,0)-hfr.E(jk,0)-hfr.E(kk,0)),
				denom3=1.0/(hfr.E(ik,0)+hfr.E(kk,0)-hfr.E(ck,0)-hfr.E(ak,0));

			sumdiagram-=int1*int2*(denom1*denom2*denom3);
		}
	}
	std::cout << "Z0G1Z1: " << sumdiagram << std::endl;
	sumdiagram=0.0;
	
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			//if(ki==ii||ki==ji) continue;
			//if(ci==ai||ci==bi) continue;
			double int1=ic.movsymphys(ji,ki,bi,ci,opsize),
			       int2=ic.movsymphys(ai,ci,ii,ki,opsize);
			double denom1=1.0/(hfr.E(ik,0)+hfr.E(kk,0)-hfr.E(ak,0)-hfr.E(ck,0)),
			       denom2=1.0/(hfr.E(jk,0)+hfr.E(kk,0)-hfr.E(bk,0)-hfr.E(ck,0)),
			       denom3=1.0/(E+hfr.E(ak,0)+hfr.E(bk,0)+hfr.E(ck,0)-hfr.E(ik,0)-hfr.E(jk,0)-hfr.E(kk,0));
			sumdiagram+=int1*int2*(denom1*denom2*denom3);
		}
	}
	std::cout << "Z1G+0Z1: " << sumdiagram << std::endl;
	sumdiagram=0.0;
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			//if(ki==ii||ki==ji) continue;
			//if(ci==ai||ci==bi) continue;
			double int1=ic.movsymphys(jk,kk,bk,ck,opsize),
			       int2=ic.movsymphys(ak,ck,ik,kk,opsize);
			double denom1=1.0/(hfr.E(ik,0)+hfr.E(kk,0)-hfr.E(ak,0)-hfr.E(ck,0)),
			       denom2=1.0/(hfr.E(jk,0)+hfr.E(kk,0)-hfr.E(bk,0)-hfr.E(ck,0)),
			       denom3=1.0/(-E+hfr.E(ck,0)-hfr.E(kk,0));
			sumdiagram+=int1*int2*(denom1*denom2*denom3);
		}
	}
	std::cout << "Z1G-0Z1: " << sumdiagram << std::endl;
	sumdiagram=0.0;
	double insize=0.0;
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			double int1=ic.movsymphys(ki,ai,ci,ii,opsize),
			       int2=ic.movsymphys(ji,ci,bi,ki,opsize);
			double denom1=1.0/(E+hfr.E(ak,0)-hfr.E(ik,0)),
			       denom2=1.0/(E+hfr.E(ck,0)-hfr.E(kk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));;
			sumdiagram+=int1*int2*(denom1*denom2*denom3);
			insize+=int1*int2*(denom1*denom2*denom3);
		}
	}
	std::cout << insize << std::endl;
	insize=0.0;
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;
			if(ki==ii||ki==ji) continue;
			if(ci==ai||ci==bi) continue;

			double int1=ic.movsymphys(ji,ki,bi,ci,opsize),
			       int2=ic.movsymphys(ai,ci,ii,ki,opsize);
			double denom1=1.0/(E+hfr.E(ak,0)-hfr.E(ik,0)),
			       denom2=1.0/(E+hfr.E(ck,0)+hfr.E(bk,0)+hfr.E(ak,0)-hfr.E(kk,0)-hfr.E(ik,0)-hfr.E(jk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));;
			sumdiagram+=int1*int2*(denom1*denom2*denom3);
			insize+=int1*int2*(denom1*denom2*denom3);
		}
	}
	std::cout << insize << std::endl;
	insize=0.0;
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			if(ki==ii||ci==bi) continue;

			double int1=ic.movsymphys(ki,ai,bi,ci,opsize),
			       int2=ic.movsymphys(ji,ci,ki,ii,opsize);
			double denom1=1.0/(E+hfr.E(ak,0)-hfr.E(ik,0)),
			       denom2=1.0/(E+hfr.E(ck,0)+hfr.E(bk,0)-hfr.E(kk,0)-hfr.E(ik,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));;
			sumdiagram-=int1*int2*(denom1*denom2*denom3);
			insize-=int1*int2*(denom1*denom2*denom3);
		}
	}
	std::cout << insize << std::endl;
	insize=0.0;
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;
			
			if(ki==ji||ci==ai) continue;

			double int1=ic.movsymphys(ji,ki,ci,ii,opsize),
			       int2=ic.movsymphys(ci,ai,bi,ki,opsize);
			double denom1=1.0/(E+hfr.E(ak,0)-hfr.E(ik,0)),
			       denom2=1.0/(E+hfr.E(ck,0)+hfr.E(ak,0)-hfr.E(kk,0)-hfr.E(jk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));;
			sumdiagram-=int1*int2*(denom1*denom2*denom3);
			insize-=int1*int2*(denom1*denom2*denom3);
		}
	}
	std::cout << insize << std::endl;
	insize=0.0;
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
		for(int d = 0; d < hfr.phlist[1].size();d++){
			int ci=hfr.phlist[1][c],di=hfr.phlist[1][d];
			int ck=ci%opsize,dk=di%opsize;
					
			double int1=ic.movsymphys(ji,ai,ci,di,opsize),
			       int2=ic.movsymphys(ci,di,bi,ii,opsize);
			double denom1=1.0/(E+hfr.E(ak,0)-hfr.E(ik,0)),
			       denom2=1.0/(E+hfr.E(ck,0)+hfr.E(dk,0)-hfr.E(ik,0)-hfr.E(jk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));;
			sumdiagram+=0.5*int1*int2*(denom1*denom2*denom3);
			insize+=0.5*int1*int2*(denom1*denom2*denom3);
		}}
	}
	std::cout << insize << std::endl;
	insize=0.0;
	for(int k = 0; k < hfr.phlist[0].size();k++){
	for(int l = 0; l < hfr.phlist[0].size();l++){
		int ki=hfr.phlist[0][k],li=hfr.phlist[0][l];
		int kk=ki%opsize,lk=li%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			double int1=ic.movsymphys(ki,li,bi,ii,opsize),
			       int2=ic.movsymphys(ji,ai,ki,li,opsize);
			double denom1=1.0/(E+hfr.E(ak,0)-hfr.E(ik,0)),
			       denom2=1.0/(E+hfr.E(ak,0)+hfr.E(bk,0)-hfr.E(kk,0)-hfr.E(lk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));;
			sumdiagram+=0.5*int1*int2*(denom1*denom2*denom3);
			insize+=0.5*int1*int2*(denom1*denom2*denom3);
		}
	}}
	std::cout << insize << std::endl;
	insize=0.0;

	std::cout << "Z0G+2Z0: " << sumdiagram << std::endl;



	/*
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			double int1=ic.movsymphys(ki,ai,ci,ii,opsize),
			       int2=ic.movsymphys(ji,ci,bi,ki,opsize);
			double denom1=1.0/(E-hfr.E(ik,0)+hfr.E(ak,0)),
			       denom2=1.0/(E+hfr.E(ck,0)-hfr.E(kk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));

			sumdiagram+=int1*int2*(denom1*denom2*denom3);
		}
	}
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			double int1=ic.movsymphys(ki,ai,bi,ci,opsize),
			       int2=ic.movsymphys(ji,ci,ki,ii,opsize);
			double denom1=1.0/(E-hfr.E(ik,0)+hfr.E(ak,0)),
			       denom2=1.0/(E+hfr.E(ck,0)+hfr.E(bk,0)-hfr.E(kk,0)-hfr.E(ik,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));

			sumdiagram-=int1*int2*(denom1*denom2*denom3);
		}
	}
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			double int1=ic.movsymphys(ji,ki,ci,ii,opsize),
			       int2=ic.movsymphys(ci,ai,bi,ki,opsize);
			double denom1=1.0/(E-hfr.E(ik,0)+hfr.E(ak,0)),
			       denom2=1.0/(E+hfr.E(ck,0)+hfr.E(ak,0)-hfr.E(kk,0)-hfr.E(jk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));

			sumdiagram-=int1*int2*(denom1*denom2*denom3);
		}
	}
	for(int k = 0; k < hfr.phlist[0].size();k++){
	for(int l = 0; l < hfr.phlist[0].size();l++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;
		int li=hfr.phlist[0][l];
		int lk=li%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			double int1=ic.movsymphys(ki,li,bi,ii,opsize),
			       int2=ic.movsymphys(ji,ai,ki,li,opsize);
			double denom1=1.0/(E-hfr.E(ik,0)+hfr.E(ak,0)),
			       denom2=1.0/(E+hfr.E(ak,0)+hfr.E(bk,0)-hfr.E(kk,0)-hfr.E(lk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));

			sumdiagram+=0.5*int1*int2*(denom1*denom2*denom3);
		}
	}
	}
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;

			double int1=ic.movsymphys(ji,ki,bi,ci,opsize),
			       int2=ic.movsymphys(ai,ci,ii,ki,opsize);
			double denom1=1.0/(E-hfr.E(ik,0)+hfr.E(ak,0)),
			       denom2=1.0/(-E+hfr.E(ck,0)+hfr.E(bk,0)+hfr.E(ak,0)-hfr.E(kk,0)-hfr.E(ik,0)-hfr.E(jk,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));

			sumdiagram+=int1*int2*(denom1*denom2*denom3);
		}
	}
	for(int k = 0; k < hfr.phlist[0].size();k++){
		int ki=hfr.phlist[0][k];
		int kk=ki%opsize;

		for(int c = 0; c < hfr.phlist[1].size();c++){
		for(int d = 0; d < hfr.phlist[1].size();d++){
			int ci=hfr.phlist[1][c];
			int ck=ci%opsize;
			int di=hfr.phlist[1][d];
			int dk=di%opsize;

			double int1=ic.movsymphys(ji,ai,ci,di,opsize),
			       int2=ic.movsymphys(ci,di,bi,ii,opsize);
			double denom1=1.0/(E-hfr.E(ik,0)+hfr.E(ak,0)),
			       denom2=1.0/(E+hfr.E(ck,0)+hfr.E(dk,0)-hfr.E(jk,0)-hfr.E(ik,0)),
			       denom3=1.0/(E+hfr.E(bk,0)-hfr.E(jk,0));

			sumdiagram+=0.5*int1*int2*(denom1*denom2*denom3);
		}
		}
	}
	*/
	
	std::cout << sumdiagram << std::endl;
	std::cout << greens[2](opk,opl)<<std::endl;
	std::cout << sumdiagram-greens[2](opk,opl)<<std::endl;
	
	Matrix V = fcis.CIMat-fcis.H0;
	Matrix E1m=Matrix::Identity(V.rows(),V.cols())*mbptr.energies[1];
	Matrix G2test=fcis.Gm_p[0]*(V-E1m)*fcis.Gm_p[0]*(V-E1m)*fcis.Gm_p[0];
	std::cout << (fcis.za[0].transpose()*G2test*fcis.za[0])(opk,opl)<<std::endl;

	Matrix submat=fcis.Gm_p[0]*(V-E1m)*fcis.Gm_p[0]*fcis.za[0];
	std::vector<int> nonzeroes;
	for(int i = 0; i < fcis.za[0].rows(); i++){
		if(std::abs(submat(i,opl))>0.0000000001){
			//std::cout << i << " " << sm.printstrab(i) << " " << sm.printstrphab(i) << " " << submat(i,opl)<<std::endl;
			nonzeroes.push_back(i);
		}
	}
	
	std::cout << "ATTEMPTING TO SORT:\n";
	for(int i = 0; i < nonzeroes.size(); i++){
		int lowest=i;
		for(int j = i+1; j<nonzeroes.size();j++){
			int lk=nonzeroes[lowest],jk=nonzeroes[j];
			
			std::vector<std::vector<int>> cla=sm.strph(lk%sm.strs.size(),true), clb=sm.strph(lk/sm.strs.size(),false);
			std::vector<std::vector<int>> ja=sm.strph(jk%sm.strs.size(),true), jb=sm.strph(jk/sm.strs.size(),false);

			int clen=cla[0].size()+clb[0].size();
			int jen=ja[0].size()+jb[0].size();

			if(jen<clen){ lowest=j; continue;}
			if(clen<jen){continue;}
				
			int extent=jen;
			int jext=0;
			int cext=0;
			for(int i = 0; i < ja[0].size();i++){
				jext+=std::pow(hfr.nelec,i)*ja[0][i];
			}
			for(int i = 0; i < jb[0].size();i++){
				jext+=std::pow(hfr.nelec,i+ja[0].size())*jb[0][i];
			}
			for(int i = 0; i < ja[1].size();i++){
				jext+=std::pow(hfr.nelec,i+extent)*ja[1][i];
			}
			for(int i = 0; i < jb[1].size();i++){
				jext+=std::pow(hfr.nelec,i+extent+ja[1].size())*jb[1][i];
			}
			for(int i = 0; i < cla[0].size();i++){
				cext+=std::pow(hfr.nelec,i)*cla[0][i];
			}
			for(int i = 0; i < clb[0].size();i++){
				cext+=std::pow(hfr.nelec,i+cla[0].size())*clb[0][i];
			}
			for(int i = 0; i < cla[1].size();i++){
				cext+=std::pow(hfr.nelec,i+extent)*cla[1][i];
			}
			for(int i = 0; i < clb[1].size();i++){
				cext+=std::pow(hfr.nelec,i+extent+cla[1].size())*clb[1][i];
			}

			if(jext<cext){
				lowest=j;
			}
			
		}
		std::swap(nonzeroes[i],nonzeroes[lowest]);
	}

	for(int i = 0; i < nonzeroes.size(); i++){
		std::cout << nonzeroes[i] << " " <<sm.printstrab(nonzeroes[i]) << " " << sm.printstrphab(nonzeroes[i])<<std::endl;
	}

	return 0;
	/*
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
