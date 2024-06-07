#include "Integrals.h"
#include <libint2.hpp>

IntegralChugger::IntegralChugger(BasisSet bs, ModelParams mp):basis{bs}, params{mp}{

}

void IntegralChugger::compute(IntegralType it){
	if(it|OVERLAP){
		S=compute1eints(libint2::Operator::overlap);
	}
	if(it|KINETIC){
		T=compute1eints(libint2::Operator::kinetic);
	}
	if(it|NUCLEAR){
		V=compute1eints(libint2::Operator::nuclear);
	}
	if(it|TWOBODY){
		compute2eints();
	}
}

/*
double IntegralChugger::operator()(int a, int b, int c, int d){
	std::array<int,4> co = get2bodyintcord(a, b, c, d);
	return tbi[co[0]][co[1]][co[2]][co[3]];

}
*/

/*
double IntegralChugger::SlaterInt(SlaterDet s1, SlaterDet s2, Matrix & C, Matrix & E){
	//At the moment I'm only conserned with excitations (For CI)
	if(s1.size()!=s2.size()){return 0;}

	SlaterDet o1, o2;

	for(auto si1=s1.begin(); si1!=s1.end();){
		bool match = false;
		for(auto si2=s2.begin(); si2!=s2.end();){
			if(*si1==*si2){
				match = true;
				si1=s1.erase(si1);
				si2=s2.erase(si2);
				break;
			}
		}
		if(!match){si1++;si2++;}
	}

	if(s1.size() != s2.size()){
		std::cout << "SIZE IMBALANCE THIS SHOULDN'T HAVE HAPPENED\n";
	}

	if(s1.size() == 0) return computeSlater0Diff(o1,o2,C,E);
	else if (s1.size() == 1) return computeSlater1Diff(s1,s2,C,E);
	else if (s1.size() == 2) return computeSlater2Diff(s1,s2,C,E);

	return 0;
}

double IntegralChugger::SlaterInt(SlaterDet s1, SlaterDet s2, Matrix & C, MatVec & E){
	if(s1.norbitals != C.cols() || s2.norbitals != C.cols()|| s1.nelectrons != s2.nelectrons){
		std::cout << "THIS SHOULDN'T HAVE HAPPENED\n";
		return 0;
	}

	int excitediff = std::abs(s1.excitationlevel - s2.excitationlevel);
	
	std::vector<int> diff1;
	std::vector<bool> used;
	used.resize(s2.nelectrons);
	for(int i = 0; i < used.size(); i++){used[i]=false;}
	for(int i =0; i < s1.nelectrons; i++){
		bool match =false;
		for(int j =0; j < s2.nelectrons; j++){
			if(s1[i]==s2[j]&&!used[j]) {
				match = true;
				used[j]=true;
				break;
			}
		}
		if(!match) diff1.push_back(s1[i]);
	}
	std::vector<int> diff2;
	for(int i = 0; i < used.size(); i++){used[i]=false;}
	for(int i =0; i < s2.nelectrons; i++){
		bool match =false;
		for(int j =0; j < s1.nelectrons; j++){
			if(s1[j]==s2[i]&&!used[j]) {
				match = true;
				used[j]=true;
				break;
			}
		}
		if(!match) diff2.push_back(s2[i]);
	}
	if(diff1.size() != diff2.size()) std::cout << "HOW DID YOU DO THIS?\n";

	switch(diff1.size()){
		case 0:
			return computeSlater0Diff(s1,s2,C,E);
			break;
		case 1:
			return computeSlater1Diff(s1,s2,C,E,diff1[0],diff2[0]);
			break;
		case 2:
			return computeSlater2Diff(s1,s2,C,E,diff1[0],diff1[1],diff2[0],diff2[1]);
			break;
		default:
			return 0;
	}
}
*/

Matrix IntegralChugger::compute1eints(libint2::Operator obtype){
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	int n = basis.nbf();
	int nshells=basis.size();

	Matrix intMat(n,n);

	Engine engine(obtype, basis.max_nprim(), basis.max_l(),0);

	if(obtype == Operator::nuclear){
		std::vector<std::pair<double, std::array<double,3>>> chargepos;
		for(auto & at: params.atoms){
			chargepos.push_back({(double)at.atomic_number,{{at.x,at.y,at.z}}});
		}
		engine.set_params(chargepos);
	}

	const auto & buf = engine.results();
	
	for(int i = 0; i < nshells; i++){
		int bf1=basis.shell2bf()[i];
		int n1 = basis[i].size();

		for(int j = 0; j <= i; j++){
			int bf2=basis.shell2bf()[j];
			int n2=basis[j].size();

			engine.compute(basis[i],basis[j]);

			Eigen::Map<const Matrix> buf_mat(buf[0],n1,n2);
			intMat.block(bf1,bf2,n1,n2)=buf_mat;

			if(i!=j){
				intMat.block(bf2,bf1,n2,n1)=buf_mat.transpose();
			}
		}
	}
	return intMat;
}

void IntegralChugger::compute2eints(){
	using libint2::Shell;
	using libint2::Engine;
	using libint2::Operator;

	int n = basis.nbf();
	int nshells=basis.size();

	Engine engine(Operator::coulomb, basis.max_nprim(), basis.max_l());

	const auto & buf = engine.results();

	tbi.resize(n);
	for(int a = 0; a < n; a++){
		tbi[a].resize(a+1);
		for(int b = 0; b <= a; b++){
			tbi[a][b].resize(a+1);
			for(int c = 0; c <= a; c++){
				int dsize=(a==c)?b:c;
				tbi[a][b][c].resize(dsize+1);
				for(int d = 0; d <=dsize; d++){
					tbi[a][b][c][d]=0.0;
				}
			}
		}
	}

	for(int s1= 0; s1 < nshells; s1++){
		int bf1_first = basis.shell2bf()[s1];
		int n1 = basis[s1].size();
		for(int s2=0; s2 <= s1; s2++){
			int bf2_first=basis.shell2bf()[s2];
			int n2 = basis[s2].size();
			for(int s3=0; s3 <= s1; s3++){
				int bf3_first=basis.shell2bf()[s3];
				int n3 = basis[s3].size();

				int s4_max = (s1==s3)?s2:s3;
				for(int s4 = 0; s4 <= s4_max; s4++){
					int bf4_first=basis.shell2bf()[s4];
					int n4 = basis[s4].size();

					engine.compute(basis[s1], basis[s2], basis[s3], basis[s4]);

					const auto * buf_1234=buf[0];

					if(buf_1234==nullptr) continue;

					for(int f1=0,f1234=0; f1<n1;f1++){
					int bf1=f1+bf1_first;
					for(int f2=0; f2<n2;f2++){
					int bf2=f2+bf2_first;
					for(int f3=0; f3<n3; f3++){
					int bf3=f3+bf3_first;
					for(int f4=0; f4<n4; f4++, f1234++){
					int bf4=f4+bf4_first;
						std::array<int,4> co=get2bodyintcord(bf1,bf2,bf3,bf4);
						double value=buf_1234[f1234];
						tbi[co[0]][co[1]][co[2]][co[3]]=value;
					}}}}
				}
			}
		}
	}

}

std::array<int,4> IntegralChugger::get2bodyintcord(int a, int b, int c, int d){
	if(b>a) std::swap(a,b);
	if(d>c) std::swap(c,d);
	if(c>a){std::swap(c,a);std::swap(d,b);}
	if(a==c&&d>b){std::swap(a,c);std::swap(d,b);}
	else if(d>c){std::swap(c,d);}

	return std::array<int,4>{a,b,c,d};
}

double IntegralChugger::tbv(int a, int b, int c, int d){
	auto cord=get2bodyintcord(a,b,c,d);
	return tbi[cord[0]][cord[1]][cord[2]][cord[3]];
}
double IntegralChugger::mov(int a, int b, int c, int d){
	auto cord=get2bodyintcord(a,b,c,d);
	return moi[cord[0]][cord[1]][cord[2]][cord[3]];
}

void IntegralChugger::TransformInts(Matrix & C){
	moi = tbi;

	int n = C.rows();

	std::cout << "TRANSFORMING 1e ints\n";
	MOH=T+V;
	Matrix H=T+V;
	for(int r = 0; r < MOH.rows(); r++){
		for(int c=0; c<MOH.cols();c++){
			double sum =0.0;
			for(int c1=0;c1<n;c1++){
				for(int c2=0; c2<n;c2++){
					sum+=C(c1,r)*C(c2,c)*H(c1,c2);
				}
			}
			MOH(r,c)=sum;
		}
	}
	
	std::cout << "TRANSFORMING 2e ints\n";
	twobodylist t1,t2, test;
	inittwobodylist(t1,n);
	inittwobodylist(t2,n);

	for(int t= 0; t < 4; t++){
		partialtransform(C,t,(t%2==0)?t1:t2,(t%2==0)?t2:(t==3)?moi:t1);
	}
	cleartwobodylist(t1);
	cleartwobodylist(t2);
}

void IntegralChugger::TransformFock(Matrix & C, Matrix & G){
	std::cout << "TRANSFORMING G MATRIX FOR LAMBDA CALCULATION\n";

	int n = C.rows();

	MOG=Matrix(n,n);
	for(int r = 0; r < MOG.rows(); r++){
		for(int c=0; c < MOG.cols();c++){
			double sum = 0.0;
			for(int c1=0; c1<n;c1++){
				for(int c2=0; c2 < n; c2++){
					sum+=C(c1,r)*C(c2,c)*G(c1,c2);
				}
			}
			MOG(r,c)=sum;
		}
	}	
}

void IntegralChugger::partialtransform(Matrix & C, int type,twobodylist & pulllist, twobodylist & addlist){
	int n = C.rows();
	
	for(int i = 0; i < addlist.size();i++){
	for(int j = 0; j < addlist[i].size(); j++){
	for(int k = 0; k < addlist[i][j].size(); k++){
	for(int l = 0; l < addlist[i][j][k].size(); l++){
		addlist[i][j][k][l]=0.0;
		for(int a = 0; a < n; a++){
			double transformval=0.0;
			double lastint=0.0;
			switch(type){
			case 0:
				transformval=C(a,i);
				lastint=tbv(a,j,k,l);
				break;
			case 1:
				transformval=C(a,j);
				lastint=pulllist[i][a][k][l];
				break;
			case 2:
				transformval=C(a,k);
				lastint=pulllist[i][j][a][l];
				break;
			case 3:
				transformval=C(a,l);
				lastint=pulllist[i][j][k][a];
				break;
			default:break;
			}
			addlist[i][j][k][l]+=transformval*lastint;
		}
	}}}}

}

void IntegralChugger::inittwobodylist(twobodylist & tbl, int n){
	tbl.resize(n);
	for(auto & a: tbl){
		a.resize(n);
		for(auto &b: a){
			b.resize(n);
			for(auto & c: b){
				c.resize(n);
			}
		}
	}
}

void IntegralChugger::cleartwobodylist(twobodylist & tbl){
	for(auto & a: tbl){
		for(auto &b: a){
			for(auto & c: b){
				c.clear();
			}
			b.clear();
		}
		a.clear();
	}
	tbl.clear();
}

