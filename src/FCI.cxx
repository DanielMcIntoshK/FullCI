#include "FCI.h"
#include <iostream>

FullCISolver::FCIResults FullCISolver::fci(ModelParams & mp, IntegralChugger &ic, HartreeFockSolver::HFResults &hf){
	ints=&ic;
	ints->TransformInts(hf.C);

	int N=hf.C.rows();
	int n=mp.nelec/2;
	SlaterDet::buildStrings(N,n);
	//SlaterDet::printStrings();
	
	strcnt=SlaterDet::codes.strs.size();
	cisize=strcnt*strcnt;

	computeHamiltonian();	
	
	Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(CIMat);

	Matrix CIE=eig_solver.eigenvalues();

	//std::cout << CIE<<std::endl;

	FCIResults fcir;
	fcir.eigenvectors=CIE;
	return fcir;	
}

void FullCISolver::computeHamiltonian(){
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

	SlaterDet sd1(alpha1idx,beta1idx),sd2(alpha2idx,beta2idx);

	SlaterCompare sc=compareSlaterDet(sd1,sd2);

	double e1=matel1e(sc);
	double e2=matel2e(sc);
	double me=e1+e2;

	if((r==0||c==0) && (sc.diff[0].size()+sc.diff[1].size())==2){
		std::cout << "SINGLE EXCITATION GROUND TEST: " << me << std::endl;
		std::cout << e1<< " " << e2 <<std::endl;
	}
	
	return me;
}

double FullCISolver::matel1e(SlaterCompare & sc, bool verbose){
	double val = 0.0;

	int totaldiff=sc.diff[0].size()+sc.diff[1].size();
	switch(totaldiff){
	case 0:
		//Evaluate all doublly occupied orbitals
		for(int b = 0; b < sc.share_do.size(); b++){
			int m = sc.share_do[b];
			val+=2.0*ints->MOH(m,m);
		}
		//Evaluate all singly occupied orbitals
		for(int spin =0; spin < 2; spin++){
			for(int b = 0; b < sc.share_so[spin].size(); b++){
				int m = sc.share_so[spin][b];
				val+=ints->MOH(m,m);
			}
		}
		break;
	case 2:{
		       //Evaluate just the orbitals that differ between slater determinants
			int diffspin=(sc.diff[0].empty())?1:0;
			int i=sc.diff[diffspin][0], j=sc.diff[diffspin][1];

			val=ints->MOH(i,j);
		}break;
	default:
		val=0.0;
	}
	return val;
}

double FullCISolver::matel2e(SlaterCompare & sc, bool verbose){
	double val = 0.0;
	int totaldiff=sc.diff[0].size()+sc.diff[1].size();

	switch(totaldiff){
	case 0:{
		//For all shared double occupied orbitals
		for(int a = 0; a < sc.share_do.size(); a++){
			int i = sc.share_do[a];
			//For all alpha and beta singly occupied orbitals
			for(int d=0; d <2;d++){
			for(int b = 0; b < sc.share_so[d].size(); b++){
				int j = sc.share_so[d][b];
				val+=2.0*ints->mov(i,i,j,j)-ints->mov(i,j,i,j);
			}
			}
			//Itself
			val+=ints->mov(i,i,i,i);
			//All other doubly occupied orbitals
			for(int b = 0; b < sc.share_do.size(); b++){
				int j = sc.share_do[b];

				val+= 4.0*ints->mov(i,i,j,j)-2.0*ints->mov(i,j,j,i);
			}
		}
		//For all the singly occupied alpha orbitals
		for(int a = 0; a < sc.share_so[0].size();a++){
			int i = sc.share_so[0][a];
			//all alpha alpha interactions
			for(int b = 0; b < sc.share_so[0].size();b++){
				int j = sc.share_so[0][b];

				val+=ints->mov(i,i,j,j)-ints->mov(i,j,j,i);
			}
			//all alpha beta interactions
			for(int b = 0; b < sc.share_so[1].size(); b++){
				int j = sc.share_so[1][b];
				val+=ints->mov(i,i,j,j);
			}
		}
		//For all singly occupied beta orbitals
		for(int a = 0; a < sc.share_so[1].size();a++){
			int i = sc.share_so[1][a];
			//all beta beta interactions
			for(int b = 0; b < sc.share_so[1].size();b++){
				int j = sc.share_so[1][b];
				val+=ints->mov(i,i,j,j)-ints->mov(i,j,j,i);
			}
		}
		}break;
	case 2:{
		//Check if diff is in alpha or beta string
		int diffspin=sc.diff[0].empty()?1:0;
		int i = sc.diff[diffspin][0], j = sc.diff[diffspin][1];
		//Interactions with all the doubly occupied orbitals
		for(int b = 0; b < sc.share_do.size(); b++){
			int m = sc.share_do[b];
			val+=2*ints->mov(i,j,m,m)-ints->mov(i,m,m,j);
		}
		//Interactions with singly occupied orbitals of same spin
		for(int b =0; b < sc.share_so[diffspin].size();b++){
			int m = sc.share_so[diffspin][b];
			val+=ints->mov(i,j,m,m)-ints->mov(i,m,m,j);
		}
		//Interactions with singly occupied orbitals of different spin
		int diffspininv=(diffspin+1)%2;
		for(int b = 0; b < sc.share_so[diffspininv].size();b++){
			int m = sc.share_so[diffspininv][b];
			val+=ints->mov(i,j,m,m);
		}
		}break;

	case 4:
	       //If diffs are in alpha and beta strings sperarly
		if(sc.diff[0].size()==2){
			int i=sc.diff[0][0],
			    j=sc.diff[1][0],
			    k=sc.diff[0][1],
			    l=sc.diff[1][1];

			    val=ints->mov(i,k,j,l);
		}
		//Diffs only in alpha or beta strings
		else{
			int diffspin=sc.diff[0].empty()?1:0;

			int i=sc.diff[diffspin][0],
				j=sc.diff[diffspin][1],
				k=sc.diff[diffspin][2],
				l=sc.diff[diffspin][3];
			val=ints->mov(i,k,j,l)-ints->mov(i,l,j,k);
		}
		break;
	}
	return val;
}

double FullCISolver::matel1e(std::vector<int> & diff,std::vector<int> & share){
	double val=0.0;
	switch(diff.size()){
	case 0:
		for(int i = 0; i < share.size();i++) val+=ints->MOH(share[i],share[i]);
		break;
	case 2:
		val=ints->MOH(diff[0],diff[1]);	
		break;
	default:
		val=0.0;
	}
	return val;
}

double FullCISolver::matel2e(std::vector<int> & diff,std::vector<int> & share,bool verbose){

	double val =0.0;
	if(verbose){
		std::cout << "DIFF: ";
		for(auto a: diff) std::cout << a << " ";
		std::cout << std::endl << "SHARE: ";
		for(auto a:share) std::cout << a << " ";
		std::cout << std::endl;
	}
	switch(diff.size()){
	case 0:
		for(int i = 0; i < share.size(); i++){
			for(int j = 0; j < share.size(); j++){
				int m = share[i],n=share[j];
				val+=(ints->mov(m,m,n,n)-ints->mov(m,n,n,m));
				//val+=(ints->mov(i,i,j,j)-ints->mov(i,j,j,i));
			}
		}
		val*=.5;
		break;
	case 2:
		for(int i = 0; i < ints->MOH.rows(); i++){
		//for(int i = 0; i < share.size();i++){
			//int n = share[i];
			int n = i;
			double nv=(ints->mov(diff[0],diff[1],n,n)-ints->mov(diff[0],n,n,diff[1]));
			//nv=(ints->mov(diff[0],n,diff[1],n)-ints->mov(diff[0],n,n,diff[1]));
			if(verbose) std::cout << nv << std::endl;
			val+=nv;
		}
		break;
	case 4:
		//val=ints->mov(diff[0],diff[1],diff[2],diff[3])-ints->mov(diff[0],diff[1],diff[3],diff[2]);
		val=ints->mov(diff[0],diff[2],diff[1],diff[3])-ints->mov(diff[0],diff[3],diff[1],diff[3]);
		break;
	default:
		val=0.0;
	}
	return val;
}

void FullCISolver::cleanup(){
	SlaterDet::cleanStrings();
	CIMat=Matrix(0,0);
}

