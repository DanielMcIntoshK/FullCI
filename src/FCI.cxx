#include "FCI.h"
#include <iostream>
#include <cmath>

FullCISolver::FCIResults FullCISolver::fci(ModelParams & mp, IntegralChugger &ic, HartreeFockSolver::HFResults &hf, double lam){
	lambda=lam;
	//lambdaDeriv=(std::abs(lambda-1.0)>.000001);
	lambdaDeriv = true;

	std::cout << "RUNNING FULL CI ENERGY CALCULATION\n";
	
	ints=&ic;
	ints->TransformInts(hf.C);
	ints->TransformFock(hf.C,hf.G);

	std::cout << "BUILDING SLATER DETERMINANT STRINGS\n";
	int N=hf.C.rows();
	int n=mp.nelec/2;
	SlaterDet::buildStrings(N,n);
	
	std::cout << "DECODE TEST:\n";
	for(int i = 0; i < SlaterDet::codes.strs.size(); i++){
		std::cout << i << ": " << SlaterDet::decode(SlaterDet::codes.strs[i]) << std::endl;
	}


	strcnt=SlaterDet::codes.strs.size();
	cisize=strcnt*strcnt;
	bssize=ints->MOH.rows();

	std::cout << "USING DIRECT DIAGONALIZATION OF CI HAMILTONIAN\n";
	computeHamiltonian();	
	
	std::cout << "DIAGONALIZING HAMILTONIAN\n";
	Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(CIMat);

	//std::cout << CIE<<std::endl;

	FCIResults fcir;
	fcir.eigenvalues=eig_solver.eigenvalues();
	fcir.eigenvectors=eig_solver.eigenvectors();
	fcir.fullH=&CIMat;
	return fcir;	
}

void FullCISolver::computeHamiltonian(){
	std::cout << "COMPUTING FCI HAMILTONIAN\n";
	CIMat=Matrix(cisize,cisize);

	for(int r =0; r < CIMat.rows(); r++){
		for(int c=0; c<CIMat.cols(); c++){
			CIMat(r,c)=matrixEl(r,c);
		}
	}
}

double FullCISolver::matrixEl(int r, int c){
	int alpha1idx=r%strcnt, beta1idx=r/strcnt,
	    alpha2idx=c%strcnt, beta2idx=c/strcnt;

	fockTotal=0.0;

	SlaterDet sd1(alpha1idx,beta1idx),sd2(alpha2idx,beta2idx);

	SlaterCompare sc=compareSlaterDet(sd1,sd2);

	double e1 = secondQuantMatel1e(sd1,sd2,sc);
	double e2 = secondQuantMatel2e(sd1,sd2,sc);

	double H0=e1+fockTotal, V=e2-fockTotal;

	double me = H0+lambda*V;

	return me;
}

double FullCISolver::secondQuantMatel1e(SlaterDet & s1, SlaterDet & s2,SlaterCompare & sc){
	double val1e=0.0;
	
	bool output = false;
	
	switch(sc.diff.size()){
	case 0:{
		for(int i = 0; i < sc.share.size(); i++){
			int b = sc.share[i];
			int bi = b%bssize;

			double sqv=secondQuant1e(s1,s2,b,b,bssize);
			val1e+=ints->MOH(bi,bi)*sqv;
			if(lambdaDeriv) fockTotal+=ints->MOG(bi,bi)*sqv;
			
		}
	}break;
	case 2:{
		int p=sc.diff[0],r=sc.diff[1];
		int pi=p%bssize, ri=r%bssize;
		double sqv= secondQuant1e(s1,s2,p,r,bssize);
		val1e=ints->MOH(pi,ri)*sqv;
		if(lambdaDeriv) fockTotal+=ints->MOG(pi,ri)*sqv;
	}break;
	default:break;
	}
	
	return val1e;
}
double FullCISolver::secondQuantMatel2e(SlaterDet & s1, SlaterDet & s2,SlaterCompare & sc){
	double val2e=0.0;

	int n = sc.share.size();

	switch(sc.diff.size()){
	case 0:{
		for(int i = 0; i < sc.share.size(); i++){
			int p =sc.share[i];

			int pi=p%bssize;
			int pspin=p/bssize;
			for(int j = 0; j < sc.share.size(); j++){
				int q =sc.share[j];
				
				int qi=q%bssize;
				int qspin=q/bssize;

				if(qspin!=pspin) continue;
				for(int k = 0; k < sc.share.size(); k++){
					int r =sc.share[k];
				
					int ri = r%bssize;
					int rspin=r/bssize;
					for(int l = 0; l < sc.share.size(); l++){
						int s =sc.share[l];
				
						int si= s%bssize;
						int sspin=s/bssize;

						if(rspin!=sspin)continue;

						val2e+=ints->mov(pi,qi,ri,si)*secondQuant2e(s1,s2,p,r,s,q,bssize);
		}}}}
	}break;
	case 2:{
		int p=sc.diff[0], r=sc.diff[1];
		int pi=p%bssize,ri=r%bssize,
		    pspin=p/bssize,rspin=r/bssize;
		for(int i = 0; i < sc.share.size(); i++){
			int b = sc.share[i];		
			int bi=b%bssize,bspin=b/bssize;

			if(pspin==rspin){
				val2e+=secondQuant2e(s1,s2,p,b,b,r,bssize)*ints->mov(pi,ri,bi,bi);
				val2e+=secondQuant2e(s1,s2,b,p,r,b,bssize)*ints->mov(bi,bi,pi,ri);
			}
			if(pspin==bspin && rspin==bspin){
				val2e+=secondQuant2e(s1,s2,p,b,r,b,bssize)*ints->mov(pi,bi,bi,ri);
				val2e+=secondQuant2e(s1,s2,b,p,b,r,bssize)*ints->mov(bi,ri,pi,bi);
		}}
	}break;

	case 4:{
		int p=sc.diff[0],q=sc.diff[1],r=sc.diff[2],s=sc.diff[3];

		int pi=p%bssize, qi=q%bssize, ri=r%bssize, si=s%bssize,
		    pspin=p/bssize, qspin=q/bssize, rspin=r/bssize, sspin=s/bssize;

		if(pspin==sspin&&qspin==rspin){
			val2e+=secondQuant2e(s1,s2,p,q,r,s,bssize)*ints->mov(pi,si,qi,ri);
			val2e+=secondQuant2e(s1,s2,q,p,s,r,bssize)*ints->mov(qi,ri,pi,si);
		}
		if(qspin==sspin&&pspin==rspin){
			val2e+=secondQuant2e(s1,s2,q,p,r,s,bssize)*ints->mov(qi,si,pi,ri);
			val2e+=secondQuant2e(s1,s2,p,q,s,r,bssize)*ints->mov(pi,ri,qi,si);
		}
	}break;
	default:break;
	}
	
	return 0.5*val2e;
}

void FullCISolver::cleanup(){
	SlaterDet::cleanStrings();
	CIMat=Matrix(0,0);
}

