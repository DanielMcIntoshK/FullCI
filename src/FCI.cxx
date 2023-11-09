#include "FCI.h"
#include "SlaterDet.h"
#include <iostream>

FullCISolver::FCIResults FullCISolver::fci(ModelParams & mp, IntegralChugger &ic, HartreeFockSolver::HFResults &hf){
	ints=&ic;
	ints->TransformInts(hf.C);
	for(int i = 0; i < ints->MOH.rows();i++){
	for(int j = 0; j < ints->MOH.rows();j++){
	for(int k = 0; k < ints->MOH.rows();k++){
	for(int l = 0; l < ints->MOH.rows();l++){
		if(ints->mov(i,j,k,l)!=0) std::cout << "NONZERO " << i << " "<<j <<" "<< k <<" " << l << std::endl;
	}}}}

	int N=hf.C.rows();
	int n=mp.nelec/2;
	SlaterDet::buildStrings(N,n);
	SlaterDet::printStrings();
	
	strcnt=SlaterDet::codes.strs.size();
	cisize=strcnt*strcnt;

	computeHamiltonian();	
	
	std::cout << ints->MOH << std::endl;
	std::cout << std::endl;
	std::cout << ints->T+ints->V << std::endl;

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

	std::vector<int> diff=sd1|sd2;

	double me= matel1e(diff)+matel2e(diff);


	if((r==0||c==0) && diff.size()==2){
		std::cout << diff[0] << " " << diff[1] << " " <<matel1e(diff) << " " <<matel2e(diff,true) <<std::endl;
		std::cout << "SINGLE EXCITATION GROUND TEST: " << me << std::endl;
	}
	return me;
}

double FullCISolver::matel1e(std::vector<int> & diff){
	double val=0.0;
	switch(diff.size()){
	case 0:
		for(int i = 0; i < ints->MOH.rows();i++) val+=ints->MOH(i,i);
		break;
	case 2:
		val=ints->MOH(diff[0],diff[1]);	
		break;
	default:
		val=0.0;
	}
	return val;
}

double FullCISolver::matel2e(std::vector<int> & diff,bool verbose){

	double val =0.0;
	switch(diff.size()){
	case 0:
		for(int i = 0; i < ints->MOH.rows(); i++){
			for(int j = 0; j < ints->MOH.rows(); j++){
				val+=(ints->mov(i,j,i,j)-ints->mov(i,j,j,i));
			}
		}
		val*=.5;
		break;
	case 2:
		for(int i = 0; i < ints->MOH.rows();i++){
			double nv=(ints->mov(diff[0],i,diff[1],i)-ints->mov(diff[0],i,i,diff[1]));
			if(verbose) std::cout << nv << std::endl;
			val+=nv;
		}
		break;
	case 4:
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

