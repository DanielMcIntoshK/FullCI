#include "Greens.h"
#include <queue>
#include <cmath>
#include <bitset>
#include <algorithm>
#include <list>
#include <complex>

GreensCalculator::GreensCalculator(){

}


void GreensCalculator::buildOperators(std::vector<double> avocc,HartreeFockSolver::HFResults & hf){
	StringMap & sm=SlaterDet::codes;
	operators.clear();
	for(int s = 0; s < 2; s++){
	for(int i = 0; i < sm.norbs; i++){
		for(int j = 0; j < sm.norbs; j++){
			if((i!=j) && (avocc[i] > avocc[j])){
				operators.push_back(PHOp(i+sm.norbs*s,j+sm.norbs*s));
			}
		}}
	}
	/*
	for(int s =0; s < 2; s++){
		for(int i = 0; i < sm.nelec; i++){
			for(int r = sm.nelec; r < sm.norbs; r++){
				operators.push_back(PHOp(i+sm.norbs*s,r+sm.norbs*s));
			}
		}
	}
	std::cout << "OPERATORS\n";
	for(int i = 0; i < operators.size(); i++){
		std::cout << operators[i].i << " " << operators[i].j << std::endl;
	}
	*/

}

Matrix GreensCalculator::ComputeGreens(double E, IntegralChugger & ic, HartreeFockSolver::HFResults & hf, FullCISolver::FCIResults & fcir){
	ints= &ic;
	Matrix prop=buildProp_Smart(E,fcir);
	return prop;
}

Matrix GreensCalculator::ComputeSelfEnergy(double E, HartreeFockSolver::HFResults & hf, Matrix & G){
	StringMap & sm = SlaterDet::codes;
	Matrix ep(operators.size(),operators.size());
	for(int i = 0; i < ep.rows();i++){
		for(int j = 0; j < ep.cols();j++){
			if(i==j){
				ep(i,i)=hf.E(operators[i].j%sm.norbs,0)-hf.E(operators[i].i%sm.norbs,0);
			}
			else ep(i,j)=0.0;
		}
	}

	Eigen::FullPivLU<Matrix> lu(G);
	Matrix inverse=lu.inverse();
	Matrix G0i=Matrix::Identity(operators.size(),operators.size())*E-ep;
	//std::cout << "G0i for SE\n" <<G0i << std::endl;
	
	return G0i-inverse;
}

Matrix GreensCalculator::buildProp_Smart(double E, FullCISolver::FCIResults & fcir){
	StringMap & sm = SlaterDet::codes;
	
	//std::vector<double> GreensCalculator::computeAvOc_Pure(Matrix state){
	std::vector<double> oc=computeAvOc_Pure(fcir.eigenvectors.col(0));
	/*operators.clear();
	for(int s = 0; s < 2; s++){
	for(int i = 0; i < sm.norbs; i++){
		for(int j = 0; j < sm.norbs; j++){
			//if(oc[i]>oc[j]){
				manifold.push_back(PHOp(i+sm.norbs*s,j+sm.norbs*s));
			//}
		}}
	}
	*/


	Matrix & evecs = fcir.eigenvectors;
	Matrix & evals = fcir.eigenvalues;

	Matrix OM(operators.size(),evecs.rows()), MO(operators.size(),evecs.rows());
	for(int i = 0; i < operators.size();i++){
		for(int j = 0; j < fcir.eigenvectors.rows();j++){
			OM(i,j)=NqM(operators[i],evecs.col(0),evecs.col(j));
			MO(i,j)=NqM(operators[i],evecs.col(j),evecs.col(0));
		}
	}
	std::vector<double> Ep,Em;
	Ep.resize(evecs.rows());Em.resize(evecs.rows());
	for(int i = 0; i < evecs.rows();i++){
		double w0M=evals(i,0)-evals(0,0);
		Ep[i]=1.0/(E+w0M);
		Em[i]=1.0/(E-w0M);
	}

	Matrix prop(operators.size(),operators.size());
	for(int i = 0; i < prop.rows(); i++){
		for(int j = 0; j < prop.cols();j++){
			prop(i,j)=calcPropEl_Smart(i,j,OM,MO,Ep,Em);
		}
	}
	return prop;
}

double GreensCalculator::calcPropEl_Smart(int p, int q,Matrix & OM, Matrix & MO, std::vector<double> & Ep, std::vector<double> & Em){
	double val=0.0;
	for(int n = 1; n < OM.cols();n++){
		val+=(MO(p,n)*MO(q,n))*Em[n]-(OM(q,n)*OM(p,n))*Ep[n];
	}
	return val;
}

std::vector<double> GreensCalculator::computeAvOc_Pure(Matrix state){
	std::vector<double> ocav;

	StringMap & sm=SlaterDet::codes;

	ocav.resize(sm.norbs*2);
	for(int e = 0; e < ocav.size();e++){
		ocav[e]=0;

		int alpha=e/sm.norbs;
		int base=(alpha==0)?e:(e-sm.norbs);

		int block = base/8;
		int scan = 1<<(base%8);
		
		for(int i =0; i < state.size();i++){
			int aid=i%sm.strs.size(),bid=i/sm.strs.size();
			int sd2check=(alpha==0)?aid:bid;
			if(sm.strs[sd2check][block]&scan){
				ocav[e]+=state(i,0)*state(i,0);
			}
		}
	}

	//Might want to change this to restrict operators that excite to a different spin

	return ocav;

}

Matrix GreensCalculator::computeOverlap(Matrix state){
	Matrix S(operators.size(),operators.size());

	for(int r = 0; r < S.rows();r++){
		for(int c = 0; c < S.cols();c++){
			S(r,c)=qq(operators[r],operators[c].adjoint(),state)-qq(operators[c].adjoint(),operators[r],state);	
		}
	}
	return S;
}

Matrix GreensCalculator::computeLMatrix(std::vector<double> & ocavs){
	Matrix lambda(operators.size(),operators.size());
	for(int i = 0; i < operators.size();i++)
		for(int j = 0; j < operators.size();j++)
			lambda(i,j)=0.0;
	for(int i = 0; i < operators.size();i++){
		lambda(i,i)=ocavs[operators[i].j]-ocavs[operators[i].i];
	}
	return lambda;
}

Matrix GreensCalculator::computeAMatrix(Matrix state,Matrix ev,Matrix*fullH){
	Matrix A(operators.size(),operators.size());

	int count = 0;
	for(int r = 0; r < A.rows();r++){
		for(int c = 0; c < A.cols();c++){
			
			PHOp q1=operators[r],q2=operators[c].adjoint();
			A(r,c)=-(qqH(q1,q2,state,ev(0,0))+qqH(q2,q1,state,ev(0,0))-qHq(q1,q2,state,fullH)-qHq(q2,q1,state,fullH));
		}
	}

	return A;
}

Matrix GreensCalculator::computeBMatrix(Matrix state, Matrix ev,Matrix * fullH){
	Matrix B(operators.size(), operators.size());

	for(int r = 0; r < B.rows();r++){
		for(int c = 0; c < B.cols();c++){
			PHOp q1=operators[r],q2=operators[c];
			B(r,c)=qqH(q1,q2,state,ev(0,0))+qqH(q2,q1,state,ev(0,0))-qHq(q1,q2,state,fullH)-qHq(q2,q1,state,fullH);
		}
	}

	return B;
}

double GreensCalculator::qqH(PHOp q1, PHOp q2, Matrix state, double e0){
	return qq(q1,q2,state)*e0;
}

double GreensCalculator::qHq(PHOp q1, PHOp q2, Matrix state, Matrix * fullH){
	return qq(q1,q2,state,fullH);
}

double GreensCalculator::qq(PHOp q1, PHOp q2, Matrix state,Matrix* fullH){
	StringMap & sm = SlaterDet::codes;
	unsigned char * alphacpy=sm.cpystrs[0],*betacpy=sm.cpystrs[1];

	std::priority_queue<chslater>ch1,ch2;

	for(int i = 0; i < state.rows();i++){
		int alphaidx=i%sm.strs.size(),betaidx=i/sm.strs.size();
		memcpy(alphacpy,sm.strs[alphaidx],sm.codeblklen);
		memcpy(betacpy,sm.strs[betaidx],sm.codeblklen);
		
		chslater ch;
		ch.sign = applyPHOp(q1.adjoint(),alphacpy,betacpy);
		if(ch.sign==0.0)continue;
		ch.index = SlaterDet::decode(alphacpy)+sm.strs.size()*SlaterDet::decode(betacpy);
		ch.c=state(i,0);
		ch1.push(ch);
	}
	for(int i = 0; i < state.rows();i++){
		int alphaidx=i%sm.strs.size(),betaidx=i/sm.strs.size();
		memcpy(alphacpy,sm.strs[alphaidx],sm.codeblklen);
		memcpy(betacpy,sm.strs[betaidx],sm.codeblklen);
		
		chslater ch;
		ch.sign = applyPHOp(q2,alphacpy,betacpy);
		if(ch.sign==0.0)continue;
		ch.index = SlaterDet::decode(alphacpy)+sm.strs.size()*SlaterDet::decode(betacpy);
		ch.c=state(i,0);
		ch2.push(ch);
	}

	double value=0.0;
	while(!(ch1.empty()||ch2.empty())){
		chslater  ct1=ch1.top(),ct2=ch2.top();
		if(ct1==ct2){
			double tv=ct1.c*ct2.c*ct1.sign*ct2.sign;
			if(fullH){
				tv*=(*fullH)(ct1.index,ct2.index);
			}
			value+=tv;
			ch1.pop();
			ch2.pop();
		}
		else if(ct1<ct2) ch2.pop();
		else ch1.pop();
	}

	
	return value;
}

double GreensCalculator::NqM(PHOp q, Matrix N, Matrix M){
	StringMap & sm = SlaterDet::codes;
	unsigned char * alphacpy=sm.cpystrs[0],*betacpy=sm.cpystrs[1];

	double val=0;

	for(int i = 0; i < M.rows(); i++){
		int alphaidx=i%sm.strs.size(),betaidx=i/sm.strs.size();
		memcpy(alphacpy,sm.strs[alphaidx],sm.codeblklen);
		memcpy(betacpy,sm.strs[betaidx],sm.codeblklen);

		double sign = applyPHOp(q,alphacpy,betacpy);

		if(sign==0.0) continue;

		int index=SlaterDet::decode(alphacpy)+sm.strs.size()*SlaterDet::decode(betacpy);
		
		val+= sign*M(i,0)*N(index,0);
	}

	return val;
}

//op=a^t_j a_i
double applyPHOp(PHOp op, unsigned char * alpha, unsigned char * beta){
	double sign = 1.0;
	int isalpha=(op.i/SlaterDet::codes.norbs);
	
	int ibase=op.i-SlaterDet::codes.norbs*isalpha,
	    jbase=op.j-SlaterDet::codes.norbs*isalpha;

	unsigned char * sd = (isalpha==0)?alpha:beta;

	int iblock=ibase/8,ioff=ibase%8,
	    jblock=jbase/8,joff=jbase%8;

	if(sd[iblock]&(1<<ioff)){
		int lowercount=(isalpha==0)?0:SlaterDet::codes.nelec;
		for(int b = 0; b < ibase; b++){
			if(sd[b/8]&(1<<(b%8))){
				lowercount++;
			}
		}
		sign*=(lowercount%2==0)?1.0:-1.0;
		sd[iblock]=sd[iblock]^(1<<ioff);
	}
	else return 0.0;

	if(!(sd[jblock]&(1<<joff))){
		int lowercount=(isalpha==0)?0:SlaterDet::codes.nelec;
		for(int b = 0; b < jbase; b++){
			if(sd[b/8]&(1<<(b%8))){
				lowercount++;
			}
		}
		sign*=(lowercount%2==0)?1.0:-1.0;
		sd[jblock]=sd[jblock]|(1<<joff);
	}
	else return 0.0;

	return sign;
}

