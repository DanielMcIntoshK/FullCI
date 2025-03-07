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
	//Matrix D = initialDGuess(params,H,S);
	Matrix D = Matrix::Zero(H.rows(),H.cols());

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

RPASolver::RPASolver(IntegralChugger & icref):ic{icref}{}

RPASolver::RPAResults RPASolver::RPACalc(hfresults & hr,bool td){
	StringMap & sm=SlaterDet::codes;
	hfr=hr;

	operators.clear();
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

	Matrix A=A0(), B=B0();

	RPASolver::RPAResults tdhfr;
	if(!td){
		Eigen::EigenSolver<Matrix> eigen_solver((A-B)*(A+B));
		std::vector<double> vals;
		for(int i = 0; i < eigen_solver.eigenvalues().rows();i++){
			vals.push_back(std::sqrt(eigen_solver.eigenvalues()(i,0).real()));
		}
		std::sort(vals.begin(),vals.end());
		tdhfr.EE=vals;
		//tdhfr.EE=Matrix(vals.size(),1);
		//int count=0;
		//for(auto a: vals) tdhfr.EE(count++,0)=a;
	}
	else{
		Eigen::SelfAdjointEigenSolver<Matrix> eigen_solver(A);
		Matrix EEtemp=eigen_solver.eigenvalues();
		tdhfr.EE.resize(EEtemp.rows());
		for(int i = 0; i < EEtemp.rows();i++){tdhfr.EE[i]=eigen_solver.eigenvalues()(i,0);}
	}
	return tdhfr;

}

RPASolver::RPAResults RPASolver::RPACalcSingTrip(hfresults & hr, bool TD){
	StringMap & sm=SlaterDet::codes;
	hfr=hr;

	operators.clear();
	for(int i = 0; i < sm.nelec; i++){
		for(int r = sm.nelec; r < sm.norbs;r++){
			operators.push_back(PHOp(i,r));
		}	
	}
	/*
	std::cout << "OPERATORS: \n";
	for(int i = 0; i < operators.size();i++){
		std::cout << i << " "<<operators[i].i << " " << operators[i].j << " " << decode(operators[i].i,operators[i].j)<<std::endl;
	}
	*/

	RPASolver::RPAResults tdhfr;

	Matrix As=A0(0),Bs=B0(0),At=A0(1),Bt=B0(1);

	if(TD){
		Eigen::SelfAdjointEigenSolver<Matrix> eigen_solver(As);
		Matrix EEtemp=eigen_solver.eigenvalues();
		tdhfr.sing.resize(EEtemp.rows());
		for(int i = 0; i < EEtemp.rows();i++){tdhfr.sing[i]=eigen_solver.eigenvalues()(i,0);}
		eigen_solver=Eigen::SelfAdjointEigenSolver<Matrix>(As);
		EEtemp=eigen_solver.eigenvalues();
		tdhfr.trip.resize(EEtemp.rows());
		for(int i = 0; i < EEtemp.rows();i++){tdhfr.trip[i]=eigen_solver.eigenvalues()(i,0);}
	}
	else{
		Eigen::EigenSolver<Matrix> eigen_solver((As-Bs)*(As+Bs));
		Matrix merge=(As-Bs)*(As+Bs);
		tdhfr.sing.resize(eigen_solver.eigenvalues().rows());
		for(int i = 0; i < eigen_solver.eigenvalues().rows();i++){
			tdhfr.sing[i]=std::sqrt(eigen_solver.eigenvalues()(i,0).real());
		}
		eigen_solver=Eigen::EigenSolver<Matrix>((At-Bt)*(At+Bt));
		tdhfr.trip.resize(eigen_solver.eigenvalues().rows());
		for(int i = 0; i < eigen_solver.eigenvalues().rows();i++){
			tdhfr.trip[i]=std::sqrt(eigen_solver.eigenvalues()(i,0).real());
		}

	}
	tdhfr.EE=tdhfr.sing;
	tdhfr.EE.insert(tdhfr.EE.end(),tdhfr.trip.begin(),tdhfr.trip.end());
	std::sort(tdhfr.EE.begin(),tdhfr.EE.end());
	
	return tdhfr;	
}

RPASolver::HRPAResults RPASolver::HRPACalc(hfresults & hfr, double threshold, int terminate,bool USESCF){
	RPASolver::HRPAResults hrpa;
	RPASolver::RPAResults rpa=RPACalcSingTrip(hfr);
	hrpa.EE=rpa.EE;
	
	for(int i = 0; i < 2; i++){
		hrpa.singtrip[i].resize(operators.size());
	}

	std::vector<Matrix> A0s,B0s;
	for(int i = 0; i < 2; i++) {
		A0s.push_back(A0(i));
		B0s.push_back(B0(i));
	}
	hrpa.Cs=initC(B0s);

	threshold=std::pow(threshold,2);

	double sqdiff=0.0;
	int ccount=0;
	do{
		std::vector<double> oldEE=hrpa.EE;

		std::array<Matrix,2> Y,Z;
		for(int i = 0; i < 2; i++){
			Matrix An=A0s[i]+A1(hrpa.Cs,i);
			Matrix Bn=B0s[i]+B1(hrpa.Cs,i);

			Eigen::EigenSolver<Matrix> eigen_solver((An-Bn)*(An+Bn));
			for(int j = 0; j < eigen_solver.eigenvalues().rows(); j++){
				hrpa.singtrip[i][j]=std::sqrt(eigen_solver.eigenvalues()(j,0).real());
			}
			auto YZ_c=eigen_solver.eigenvectors();
			eigen_solver=Eigen::EigenSolver<Matrix>((An+Bn)*(An-Bn));
			auto ZY_c=eigen_solver.eigenvectors();
	
			Matrix YZ(YZ_c.rows(),YZ_c.cols()), ZY(ZY_c.rows(), ZY_c.cols());
			for(int r = 0; r < YZ.rows(); r++){
				for(int c = 0; c <  YZ.cols();c++){
					YZ(r,c)=YZ_c(r,c).real();
					ZY(r,c)=ZY_c(r,c).real();
				}
			}

			Y[i]=0.5*(YZ+ZY);
			Z[i]=0.5*(ZY-YZ);
		}
		for(int i = 0;i < 2; i++){
			hrpa.Cs[i]=constructC(Y[i],Z[i]);
		}

		sqdiff=0.0;
		for(int i = 0; i < hrpa.EE.size(); i++){
			sqdiff+=std::pow(hrpa.EE[i]-oldEE[i],2);
		}
		if(!USESCF)break;
	}while(sqdiff>threshold&&(ccount++<terminate));
	
	if(ccount>=terminate){
		std::cout << "HRPA SCF DID NOT CONVERGE\n";
	}

	return hrpa;
}

Matrix RPASolver::A0(){
	StringMap & sm = SlaterDet::codes;
	Matrix A(operators.size(),operators.size());
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
		}
	}	
	return A;

}

Matrix RPASolver::B0(){
	StringMap & sm = SlaterDet::codes;
	Matrix B(operators.size(),operators.size());
	for(int r = 0; r < operators.size(); r++){
		for(int c = 0; c < operators.size(); c++){
			int i=operators[r].i%sm.norbs,
			    a=operators[r].j%sm.norbs,
			    j=operators[c].i%sm.norbs,
			    b=operators[c].j%sm.norbs;

			int ispin=operators[r].i/sm.norbs,
			    jspin=operators[c].i/sm.norbs;

			B(r,c)=ic.mov(a,i,b,j);
			if(ispin==jspin) B(r,c)-=ic.mov(a,j,b,i);
		}
	}	
	return B;

}

Matrix RPASolver::A1(std::vector<Matrix> C,int S){
	StringMap & sm = SlaterDet::codes;
	Matrix A(operators.size(),operators.size());
	for(int r = 0; r < operators.size(); r++){
		for(int c = 0; c < operators.size(); c++){
			int i=operators[r].i,
			    a=operators[r].j,
			    j=operators[c].i,
			    b=operators[c].j;

			A(r,c)=0.0;
			if(i==j){
			for(int q=sm.nelec; q < sm.norbs; q++){
			for(int u=0; u < hfr.nelec; u++){
			for(int v=0; v < hfr.nelec; v++){
				A(r,c)-=0.5*(ic.mov(a,u,q,v)*C[0](decode(u,b),decode(v,q))+ic.mov(u,b,v,q)*C[0](decode(u,a),decode(v,q)));
			}}}}
			if(i==j){
			for(int p=sm.nelec; p < sm.norbs; p++){
			for(int q=sm.nelec; q < sm.norbs; q++){
			for(int v=0; v < sm.nelec; v++){
				A(r,c)-=0.5*(ic.mov(p,i,q,v)*C[0](decode(j,p),decode(v,q))+ic.mov(j,p,v,q)*C[0](decode(i,p),decode(v,q)));	
			}}}}
		}
	}
	return A;
}

Matrix RPASolver::B1(std::vector<Matrix> C, int S){
	StringMap & sm = SlaterDet::codes;
	Matrix B(operators.size(),operators.size());
	for(int r = 0; r < operators.size(); r++){
		for(int c = 0; c < operators.size(); c++){
			int i=operators[r].i,
			    a=operators[r].j,
			    j=operators[c].i,
			    b=operators[c].j;

			B(r,c)=0.0;
			double sum=0.0;
			double sum2=0.0;
			double sum3=0.0;
			for(int p = sm.nelec; p < sm.norbs; p++){
			for(int u = 0; u < hfr.nelec; u++){
				sum+=ic.mov(a,j,u,p)*C[0](decode(u,p),decode(i,b))+ic.mov(b,i,u,p)*C[0](decode(u,p),decode(j,a));
				sum2+=ic.mov(a,p,u,j)*C[S](decode(i,p),decode(u,b))+ic.mov(b,p,u,i)*C[S](decode(u,a),decode(j,p));
			}
			for(int q = sm.nelec; q<sm.norbs; q++){
				sum3+=ic.mov(a,p,b,q)*C[S](decode(i,p),decode(j,q));
			}
			}
			B(r,c)-=sum*std::pow(-1,S);
			B(r,c)-=sum2;
			B(r,c)+=sum3;
			sum=0.0;
			for(int u = 0; u < sm.nelec; u++){
			for(int v = 0; v < sm.nelec; v++){
				sum+=ic.mov(u,i,v,j)*C[S](decode(u,a),decode(v,b));
			}}
			B(r,c)+=sum;
		}
	}
	return B;
}

Matrix RPASolver::A0(int S){
	Matrix A(operators.size(),operators.size());
	for(int r = 0; r < A.rows(); r++){
		for(int c = 0; c < A.cols(); c++){
			int i=operators[r].i,
			    a=operators[r].j,
			    j=operators[c].i,
			    b=operators[c].j;
			A(r,c)=(a==b&&i==j)?(hfr.E(a,0)-hfr.E(i,0)):0.0;
			double sign=1.0+std::pow(-1.0,(double)S);
			A(r,c)+=sign*ic.mov(a,i,j,b);
			A(r,c)-=ic.mov(a,b,j,i);
		}
	}
	return A;
}

Matrix RPASolver::B0(int S){
	Matrix B(operators.size(),operators.size());
	for(int r = 0; r < B.rows(); r++){
		for(int c = 0; c < B.cols();c++){
			int i=operators[r].i,
			    a=operators[r].j,
			    j=operators[c].i,
			    b=operators[c].j;
			switch(S){
			case 0: B(r,c)=2.0*ic.mov(a,i,b,j)-ic.mov(a,j,b,i);break;
			case 1: B(r,c)=ic.mov(a,j,b,i);break;
			default:B(r,c)=0.0;
			}
		}
	}
	return B;
}

Matrix RPASolver::constructC(Matrix Y, Matrix Z){
	Eigen::FullPivLU<Matrix> lu(Y);
	return Z*lu.inverse();

}

std::vector<Matrix> RPASolver::initC(std::vector<Matrix> Bs){
	for(int k = 0; k < Bs.size(); k++){
		for(int r =0; r < Bs[k].rows(); r++){
		for(int c =0; c < Bs[k].cols(); c++){
			int i = operators[r].i,
			    a = operators[r].j,
			    j = operators[c].i,
			    b = operators[c].j;

			Bs[k](r,c)=-Bs[k](r,c)/(hfr.E(a,0)+hfr.E(b,0)-hfr.E(i,0)-hfr.E(j,0));
		}}

	}
	return Bs; 
}

int RPASolver::decode(int i , int j){
	StringMap & sm = SlaterDet::codes;

	return i*(sm.norbs-sm.nelec)+(j-sm.nelec);
}


