#include "FCI.h"
#include "SlaterDet.h"
#include <iostream>

FullCISolver::FCIResults FullCISolver::fci(ModelParams & mp, IntegralChugger &ic, HartreeFockSolver::HFResults &hf){
	ints=&ic;
	ints->TransformInts(hf.C);

	int N=hf.C.rows();
	int n=mp.nelec/2;
	SlaterDet::buildStrings(N,n);
	SlaterDet::printStrings();
	
	strcnt=SlaterDet::codes.strs.size();
	cisize=strcnt*strcnt;

	computeHamiltonian();	
	
	FCIResults fcir;
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

	std::vector<int> diff=sd1|sd2;

	double me= matel1e(diff)+matel2e(diff);

	if(r==0 && diff.size()==2){
		std::cout << diff[0] << " " << diff[1] << std::endl;
		std::cout << "SINGLE EXCITATION GROUND TEST: " << me << std::endl;
	}
	return me;
}

double FullCISolver::matel1e(std::vector<int> & diff){
	double val=0.0;
	switch(diff.size()/2){
	case 0:
		for(int i = 0; i < ints->MOH.rows();i++) val+=ints->MOH(i,i);
		break;
	case 1:
		val=ints->MOH(diff[0],diff[1]);	
		break;
	default:
		val=0.0;
	}
	return val;
}

double FullCISolver::matel2e(std::vector<int> & diff){
	double val =0.0;
	switch(diff.size()/2){
	case 0:
		for(int i = 0; i < ints->MOH.rows(); i++){
			for(int j = 0; j < ints->MOH.rows(); j++){
				val+=(ints->mov(i,j,i,j)-ints->mov(i,j,j,i));
			}
		}
		val*=.5;
		break;
	case 1:
		for(int i = 0; i < ints->MOH.rows();i++){
			val+=(ints->mov(diff[0],diff[1],i,i)-ints->mov(diff[0],i,i,diff[1]));
		}
		break;
	case 2:
		val=ints->mov(diff[0],diff[1],diff[2],diff[3])-ints->mov(diff[0],diff[1],diff[3],diff[2]);
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

