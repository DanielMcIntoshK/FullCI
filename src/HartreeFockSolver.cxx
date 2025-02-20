#include "HartreeFockSolver.h"
#include "SlaterDet.h"
#include <stdio.h>
#include <cmath>
#include <array>
#include <algorithm>

void HartreeFockSolver::HFResults::buildOperators(){
	operators.clear();
	for(int s = 0; s < 2; s++){
		for(int i = 0; i < nelec; i++){
			for(int r = nelec; r < norbs; r++){
				operators.push_back(PHOp(i+norbs*s,r+norbs*s));
			}
		}
	}
	for(int s = 0; s < 2; s++){
		for(int i = 0; i < nelec; i++){
			for(int r = nelec; r < norbs; r++){
				operators.push_back(PHOp(r+norbs*s,i+norbs*s));
			}
		}
	}
}

Matrix HartreeFockSolver::HFResults::getG0(double w){
	Matrix G0(operators.size(),operators.size());
	for(int i = 0; i < G0.rows(); i++){
		for(int j = 0; j < G0.cols();j++){
			if(i==j){
				if((operators[i].j%norbs)>(operators[i].i%norbs)){
					G0(i,i)=1.0/(w-(E(operators[i].j%norbs,0)-E(operators[i].i%norbs,0)));
				}
				else{
					G0(i,i)=-1.0/(w-(E(operators[i].j%norbs,0)-E(operators[i].i%norbs,0)));
				}
			}
			else{
				G0(i,j)=0.0;
			}
		}
	}
	return G0;
}

Matrix HartreeFockSolver::HFResults::getG0i(double w){
	Matrix G0i(operators.size(),operators.size());
	for(int i = 0; i < G0i.rows(); i++){
		for(int j = 0; j < G0i.cols();j++){
			if(i==j){
				if((operators[i].j%norbs)>(operators[i].i%norbs)){
					G0i(i,i)=(w-(E(operators[i].j%norbs,0)-E(operators[i].i%norbs,0)));
				}
				else{
					G0i(i,i)=-(w-(E(operators[i].j%norbs,0)-E(operators[i].i%norbs,0)));
				}
			}
			else{
				G0i(i,j)=0.0;
			}
		}
	}
	return G0i;

}

HartreeFockSolver::HartreeFockSolver(){

}

HartreeFockSolver::HFResults HartreeFockSolver::RestrictedHF(ModelParams & params, BasisSet & bs, IntegralChugger & ic){
	HartreeFockSolver::HFResults results;
	
	nelectron = 0;
	for(auto &a: params.atoms) nelectron+=a.atomic_number;
	nelectron -= params.iparams["CHARGE"];

	double nucrepulse = 0.0;
	for(int i = 0; i < params.atoms.size();i++){
		for(int j = i+1; j < params.atoms.size();j++){
			libint2::Atom & a1 = params.atoms[i],
				&a2 = params.atoms[j];
			double dx=a1.x-a2.x, dy=a1.y-a2.y, dz=a1.z-a2.z;
			double rsq=dx*dx+dy*dy+dz*dz;
			double r = std::sqrt(rsq);
			nucrepulse+=a1.atomic_number*a2.atomic_number/r;
		}
	}
	results.enuc=nucrepulse;
	
	/*
	std::cout << "COMPUTING 1 BODY INTEGRALS" << std::endl;
	Matrix S = compute1eints(params,bs,Operator::overlap);
	Matrix T = compute1eints(params,bs,Operator::kinetic);
	Matrix V = compute1eints(params,bs,Operator::nuclear);
	*/
	Matrix S = ic.S;
	Matrix T = ic.T;
	Matrix V = ic.V;

	std::cout << "CALCULATING TRANSFORM MATRIX X" << std::endl;
	Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(S);
	Matrix s = eig_solver.eigenvalues().asDiagonal();
	Matrix U = eig_solver.eigenvectors();
	for(int i = 0; i < S.rows(); i++) s(i,i)=1/std::sqrt(s(i,i));
	Matrix X = U*s*U.transpose();

	Matrix H = T+V;

	T.resize(0,0);
	V.resize(0,0);
	s.resize(0,0);
	U.resize(0,0);

	std::cout << "MAKING INITIAL DENSITY GUESS" << std::endl;
	Matrix D = initialDGuess(params,H,S);

	/*
	std::cout << "COMPUTING 2 BODY INTEGRALS" << std::endl;
	//compute2eints(bs);
	compute2eints(bs);
	*/
	
	std::cout << "BEGINING SCF LOOP" << std::endl;

	int nao = bs.nbf();
	//for(int s = 0; s < bs.size();++s) nao+=bs[s].size();

	int maxiter=params.iparams["SCFCYCLES"];
	double scfconv=params.dparams["SCFCONV"];
	int iter = 0;
	double rmsd = 0.0;
	double ediff = 0.0;
	double ehf = 0.0;
	
	Matrix C;
	Matrix G;
	Matrix F;
	Matrix E;

	do{
		if(++iter > maxiter){
			std::cout << "SCF FAILED TO CONVERGE IN " << iter << " STEPS" << std::endl;
			break;
		}

		double ehf_last = ehf;
		Matrix D_last = D;

		G = computeGMatrix(D,ic);
		F = H+G;

		//F+= computeGMatrix(D);

		Matrix Fx=X.transpose()*F*X;

		//Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(F,S);
		//Eigen::GeneralizedEigenSolver<Matrix> gen_eig_solver(F,S);
		Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(Fx);

		E = eig_solver.eigenvalues();
		Matrix Cx = eig_solver.eigenvectors();
		
		C=X*Cx;

		int ndocc = nelectron/2;
		auto Ct=C.transpose();
		for(int i = 0; i < D.rows(); i++){
			for(int j = 0; j < D.cols(); j++){
				D(i,j)=0;
				for(int a = 0; a < ndocc; a++){
					D(i,j)+=2*C(i,a)*C(j,a);
				}
			}
		}
		//auto C_occ = C.leftCols(ndocc);
		//D=2*C_occ*C_occ.transpose();

		ehf = 0.0;
		for(int i = 0; i < bs.nbf(); i++)
			for(int j = 0; j < bs.nbf(); j++)
				ehf += .5*D(i,j) * (H(i,j)+F(i,j));

		ediff = ehf - ehf_last;
		rmsd = (D-D_last).norm();
		
		std::cout << std::setprecision(10)<<"SCF STEP " << iter << "\tENERGY: " << ehf << "\tEDIFF: " << ediff << std::endl; 
	}while((std::fabs(ediff) > scfconv) || (std::fabs(rmsd) > scfconv));

	results.eelec = ehf;
	results.C = C;
	results.D = D;
	results.E = E;
	results.G = G;
	results.F = F;

	results.nelec=params.nelec/2;
	results.norbs=C.rows();
	results.buildOperators();
	
	return results;
}

Matrix HartreeFockSolver::computeGMatrix(Matrix & D,IntegralChugger &ic){
	int n = ic.tbi.size();

	Matrix G = Matrix::Zero(n,n);
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			for(int s= 0; s < n; s++){
				for(int r = 0; r < n; r++){

					double J = ic(i,j,r,s);
					double K = ic(i,s,r,j);

					G(i,j)+=D(s,r)*(J-.5*K);
				}
			}
		}	
	}
	return G;
}

Matrix HartreeFockSolver::initialDGuess(ModelParams & param, Matrix & H, Matrix & S){
	int nao = 0;
	bool use_soad = true;
	for(const auto & atom: param.atoms){
		int Z = atom.atomic_number;
		if(Z==1||Z==2) nao+=1;
		else if(Z <=10) nao+=5;
		else{use_soad=false;break;}
	}

	Matrix D = Matrix::Zero(S.rows(),S.cols());
	
	if(use_soad){
		int ao_offset = 0;
		for(const auto & atom: param.atoms){
			int Z = atom.atomic_number;
			if(Z==1 || Z==2){
				D(ao_offset,ao_offset)=Z;
				ao_offset+=1;
			}
			else{
				D(ao_offset, ao_offset) = 2;
				D(ao_offset+1, ao_offset+1) = (Z==3)?1:2;
				const double num_electrons_per_2p = (Z>4)?(double)(Z-4)/3:0;
				for(auto xyz=0; xyz<3;++xyz) D(ao_offset+2+xyz, ao_offset+2+xyz)=num_electrons_per_2p;
				ao_offset+=5;
			}
		}
		D *= 0.5;	
	}
	else{
		Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> gen_eig_solver(H,S);
		auto eps = gen_eig_solver.eigenvalues();
		auto C = gen_eig_solver.eigenvectors();
		
		int ndocc=nelectron/2;

		auto C_occ = C.leftCols(ndocc);
		D=C_occ*C_occ.transpose();
	}
	return D;
}


TDHFSolver::TDHFResults TDHFSolver::TDHFCalc(hfresults & hfr, IntegralChugger & ic,bool td){
	StringMap & sm=SlaterDet::codes;

	std::vector<PHOp> operators;
	for(int s = 0; s < 2; s++){
		for(int i = 0; i < sm.nelec; i++){
			for(int r = sm.nelec; r < sm.norbs;r++){
				operators.push_back(PHOp(i+sm.norbs*s,r+sm.norbs*s));
			}	
		}
	}
	std::cout << "OPERATORS: \n";
	for(int i = 0; i < operators.size();i++){
		std::cout << operators[i].i << " " << operators[i].j << std::endl;
	}

	Matrix A(operators.size(), operators.size()),
	       B(operators.size(),operators.size());

	for(int r = 0; r < operators.size(); r++){
		for(int c = 0; c < operators.size(); c++){
			int i=operators[r].i%sm.norbs,
			    a=operators[r].j%sm.norbs,
			    j=operators[c].i%sm.norbs,
			    b=operators[c].j%sm.norbs;

			int ispin=operators[r].i/sm.norbs,
			    jspin=operators[c].i/sm.norbs;
			
			A(r,c)=(r==c)?(hfr.E(a,0)-hfr.E(i,0)):0.0;
			A(r,c)+=ic.mov(a,i,j,b);
			if(ispin==jspin) A(r,c)-=ic.mov(a,b,j,i);
			
			B(r,c)=ic.mov(a,i,b,j);
			if(ispin==jspin) B(r,c)-=ic.mov(a,j,b,i);
		}
	}

	TDHFSolver::TDHFResults tdhfr;
	if(!td){
		Eigen::EigenSolver<Matrix> eigen_solver((A-B)*(A+B));
		std::vector<double> vals;
		for(int i = 0; i < eigen_solver.eigenvalues().rows();i++){
			vals.push_back(std::sqrt(eigen_solver.eigenvalues()(i,0).real()));
		}
		std::sort(vals.begin(),vals.end());

		tdhfr.EE=Matrix(vals.size(),1);
		int count=0;
		for(auto a: vals) tdhfr.EE(count++,0)=a;
	}
	else{
		Eigen::SelfAdjointEigenSolver<Matrix> eigen_solver(A);
		tdhfr.EE=Matrix(eigen_solver.eigenvalues().rows(),1);
		for(int i = 0; i < tdhfr.EE.rows();i++){tdhfr.EE(i,0)=eigen_solver.eigenvalues()(i,0);}
	}
	return tdhfr;

}


