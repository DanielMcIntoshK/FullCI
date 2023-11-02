#include "HartreeFockSolver.h"
#include <stdio.h>
#include <cmath>
#include <array>

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
	MatVec E;

	do{
		if(++iter > maxiter){
			std::cout << "SCF FAILED TO CONVERGE IN " << iter << " STEPS" << std::endl;
			break;
		}

		double ehf_last = ehf;
		Matrix D_last = D;

		G = computeGMatrix(D,ic);
		Matrix F = H+G;

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

		std::cout << "SCF STEP " << iter << "\tENERGY: " << ehf << "\tEDIFF: " << ediff << std::endl; 
	}while((std::fabs(ediff) > scfconv) || (std::fabs(rmsd) > scfconv));

	
	results.eelec = ehf;
	results.C = C;
	results.D = D;
	results.E = E;
	results.G = G;
	
	return results;
}

/*
Matrix HartreeFockSolver::compute1eints(ModelParams & param, BasisSet & bs, libint2::Operator obtype){
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	int n = bs.nbf();
	int nshells=bs.size();

	Matrix intMat(n,n);

	Engine engine(obtype, bs.max_nprim(), bs.max_l(), 0);

	if(obtype == Operator::nuclear){
		std::vector<std::pair<double,std::array<double,3>>> chargepos;
		for(auto & at : param.atoms){
			chargepos.push_back( {(double)at.atomic_number,{{at.x,at.y,at.z}}});
		}
		engine.set_params(chargepos);
	}

	const auto & buf = engine.results();

	for(int i = 0; i < nshells; i++){
		int bf1=bs.shell2bf()[i];
		int n1=bs[i].size();

		for(int j = 0; j <=i; j++){
			int bf2 = bs.shell2bf()[j];
			int n2 = bs[j].size();

			engine.compute(bs[i],bs[j]);

			Eigen::Map<const Matrix> buf_mat(buf[0],n1,n2);
			intMat.block(bf1,bf2,n1,n2)=buf_mat;

			if(i!=j){
				intMat.block(bf2,bf1,n2,n1)=buf_mat.transpose();
			}
		}
	}
	return intMat;
}

void HartreeFockSolver::compute2eints(BasisSet & bs){
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	int n = bs.nbf();
	int nshells=bs.size();

	Engine engine(Operator::coulomb, bs.max_nprim(), bs.max_l());

	const auto & buf = engine.results();

	twobodyints.resize(n);
	for(int a = 0; a < n; a++){
		twobodyints[a].resize(a+1);
		for(int b = 0; b <= a; b++){
			twobodyints[a][b].resize(a+1);
			for(int c = 0; c <= a; c++){
				int dsize=(a==c)?b:c;
				twobodyints[a][b][c].resize(dsize+1);
				for(int d = 0; d<= dsize; d++){
					twobodyints[a][b][c][d]=0.0;
				}
			}
		}
	}

	for(int s1 = 0; s1 < nshells; s1++){
		int bf1_first = bs.shell2bf()[s1];
		int n1 = bs[s1].size();
		for( int s2 = 0; s2 <= s1; s2++){
			int bf2_first = bs.shell2bf()[s2];
			int n2 = bs[s2].size();
			for(int s3 = 0; s3 <= s1; s3++){
				int bf3_first = bs.shell2bf()[s3];
				int n3 = bs[s3].size();

				int s4_max = (s1==s3)?s2:s3;
				for(int s4 = 0; s4 <= s4_max; s4++){
					int bf4_first=bs.shell2bf()[s4];
					int n4 = bs[s4].size();


					engine.compute(bs[s1],bs[s2],bs[s3],bs[s4]);

					const auto* buf_1234=buf[0];
					
					//bool culled=buf_1234==nullptr;
					if(buf_1234==nullptr) continue;

					for(int f1=0,f1234=0; f1<n1;f1++){
					int bf1=f1+bf1_first;
					for(int f2=0; f2<n2;f2++){
					int bf2=f2+bf2_first;
					for(int f3=0; f3<n3;f3++){
					int bf3=f3+bf3_first;
					for(int f4=0; f4<n4;f4++,f1234++){
					int bf4=f4+bf4_first;
						std::array<int,4> co=get2bodyintcord(bf1,bf2,bf3,bf4);
						
						int a=bf1, b=bf2, c=bf3, d=bf4;
						if(b>a)std::swap(a,b);
						if(d>c)std::swap(c,d);
						if(c>a){std::swap(c,a);std::swap(d,b);}
						if(a==c&&d>b){std::swap(a,c);std::swap(b,d);}
						else if(d>c){std::swap(c,d);}
						
						//std::cout << "Computing ("<<bf1<<bf2<<"|"<<bf3<<bf4<<")"<<std::endl;
						double value=buf_1234[f1234];
						twobodyints[co[0]][co[1]][co[2]][co[3]]=value;
					}
					}
					}
					}
				}
			}
		}
	}
}

void HartreeFockSolver::compute2eints_crappy(BasisSet & bs){
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	int n = bs.nbf();
	int nshells=bs.size();

	Engine engine(Operator::coulomb, bs.max_nprim(), bs.max_l());

	const auto & buf = engine.results();

	std::cout << "INITIALIZING LIST\n";
	checklist cl;
	cl.resize(n);
	twobodyints.resize(n);
	for(int a = 0; a < n; a++){
		cl[a].resize(n);
		twobodyints[a].resize(n);
		for(int b = 0; b < n; b++){
			cl[a][b].resize(n);
			twobodyints[a][b].resize(n);
			for(int c = 0; c < n; c++){
				cl[a][b][c].resize(n);
				twobodyints[a][b][c].resize(n);
				for(int d = 0; d< n; d++){
					cl[a][b][c][d]=false;
					twobodyints[a][b][c][d]=0.0;
				}
			}
		}
	}
	std::cout << "COMPUTING INTEGRALS\n";
	for(int s1 = 0; s1 < nshells; s1++){
		int bf1_first = bs.shell2bf()[s1];
		int n1 = bs[s1].size();
		for( int s2 = 0; s2 < nshells; s2++){
			int bf2_first = bs.shell2bf()[s2];
			int n2 = bs[s2].size();
			for(int s3 = 0; s3 < nshells; s3++){
				int bf3_first = bs.shell2bf()[s3];
				int n3 = bs[s3].size();
				for(int s4 = 0; s4 < nshells; s4++){
					int bf4_first=bs.shell2bf()[s4];
					int n4 = bs[s4].size();


					engine.compute(bs[s1],bs[s2],bs[s3],bs[s4]);

					const auto* buf_1234=buf[0];
					
					//bool culled=buf_1234==nullptr;
					if(buf_1234==nullptr) continue;

					for(int f1=0,f1234=0; f1<n1;f1++){
					int bf1=f1+bf1_first;
					for(int f2=0; f2<n2;f2++){
					int bf2=f2+bf2_first;
					for(int f3=0; f3<n3;f3++){
					int bf3=f3+bf3_first;
					for(int f4=0; f4<n4;f4++,f1234++){
					int bf4=f4+bf4_first;
						//std::array<int,4> co=get2bodyintcord(bf1,bf2,bf3,bf4);
						
						int a=bf1, b=bf2, c=bf3, d=bf4;
						if(b>a)std::swap(a,b);
						if(d>c)std::swap(c,d);
						if(c>a){std::swap(c,a);std::swap(d,b);}
						if(a==c&&d>b){std::swap(a,c);std::swap(b,d);}
						else if(d>c){std::swap(c,d);}
						
						//std::cout << "Computing ("<<bf1<<bf2<<"|"<<bf3<<bf4<<")"<<std::endl;
						double value=buf_1234[f1234];
						//twobodyints[co[0]][co[1]][co[2]][co[3]]=value;
						twobodyints[bf1][bf2][bf3][bf4]=value;
						cl[bf1][bf2][bf3][bf4]=true;
					}
					}
					}
					}
				}
			}
		}
	}
	for(int a = 0; a < n; a++){
		cl[a].resize(n);
		for(int b = 0; b < n; b++){
			cl[a][b].resize(n);
			for(int c = 0; c < n; c++){
				cl[a][b][c].resize(n);
				for(int d = 0; d< n; d++){
					if(!cl[a][b][c][d]){

						std::cout << "MISSED: (" << a<<b<<"|"<<c<<d<<")"<<std::endl;
					}
				}
			}
		}
	}
}
*/

Matrix HartreeFockSolver::computeGMatrix(Matrix & D,IntegralChugger &ic){
	int n = ic.tbi.size();

	Matrix G = Matrix::Zero(n,n);
	/*
	for(int a =0; a < twobodyints.size(); a++){
		for(int b = 0; b < twobodyints[a].size(); b++){
			for(int c = 0; c < twobodyints[a][b].size(); c++){
				for(int d =0; d < twobodyints[a][b][c].size(); d++){
					double ival_deg=twobodyints[a][b][c][d];

					G(a,b) += D(c,d) * ival_deg;
					G(c,d) += D(a,b) * ival_deg;
					G(a,c) -= 0.25 * D(b,d) * ival_deg;
					G(b,d) -= 0.25 * D(a,c) * ival_deg;
					G(a,d) -= 0.25 * D(b,c) * ival_deg;
					G(b,c) -= 0.25 * D(a,d) * ival_deg;
				}
			}
		}
	}
	Matrix Gt = G.transpose();
	return 0.5*(G+Gt);
	*/
	for(int i = 0; i < n; i++){
		for(int j = 0; j < n; j++){
			for(int s= 0; s < n; s++){
				for(int r = 0; r < n; r++){
					//auto c1 = get2bodyintcord(i,j,r,s);
					//auto c2 = get2bodyintcord(i,s,r,j);

					//double J=twobodyints[c1[0]][c1[1]][c1[2]][c1[3]];
					//double K=twobodyints[c2[0]][c2[1]][c2[2]][c2[3]];

					double J = ic(i,j,r,s);
					double K = ic(i,s,r,j);

					//double J=twobodyints[i][j][r][s];
					//double K=twobodyints[i][s][r][j];

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
	//return Matrix::Zero(S.rows(),S.cols());
	return D;
}

/*
std::array<int,4> HartreeFockSolver::get2bodyintcord(int a, int b, int c, int d){
	//std::cout << a << " " << b << " " << c <<" " << d << " TO ";
	if(b>a)std::swap(a,b);
	if(d>c)std::swap(c,d);
	if(c>a){std::swap(c,a);std::swap(d,b);}
	if(a==c&&d>b){std::swap(a,c);std::swap(b,d);}
	else if(d>c){std::swap(c,d);}
	//std::cout << a << " " << b << " " << c <<" " << d << std::endl;

	return std::array<int,4>{a,b,c,d};
}
*/
