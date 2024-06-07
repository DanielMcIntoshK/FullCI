#include "Greens.h"
#include <queue>
#include <cmath>
#include <bitset>

GreensCalculator::GreensCalculator(){

}

void GreensCalculator::computeExcitationSpectra(ModelParams & mp, IntegralChugger & ic, HartreeFockSolver::HFResults & hf, FullCISolver::FCIResults & fcir){
	ints = &ic;

	std::vector<double> ocavs=computeAvOc_Pure(fcir.eigenvectors.col(0));

	StringMap & sm=SlaterDet::codes;

	PHOpp testopp=manifold[3];

	unsigned char * str1=sm.cpystrs[0],*str2=sm.cpystrs[1];
	double value=0.0;
	
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

	//Might want to change this to restrict opperators that excite to a different spin
	for(int s =0; s < 2; s++){
	for(int i = 0; i < sm.norbs;i++){
		PHOpp opp;
		opp.i=i+sm.norbs*s;
		for(int j = 0; j < sm.norbs;j++){
			opp.j=j+sm.norbs*s;
			if(ocav[j]>ocav[i]) manifold.push_back(opp);
		}
	}
	}


	std::cout <<  "MANIFOLD\n";
	for(int i = 0; i < manifold.size(); i++){
		std::cout << manifold[i].i << " " << manifold[i].j << " " << ocav[manifold[i].i] << " " << ocav[manifold[i].j] << std::endl;
	}

	return ocav;

}

Matrix GreensCalculator::computeLMatrix(std::vector<double> & ocavs){
	Matrix lambda(manifold.size(),manifold.size());
	for(int i = 0; i < manifold.size();i++)
		for(int j = 0; j < manifold.size();j++)
			lambda(i,j)=0.0;
	for(int i = 0; i < manifold.size();i++){
		lambda(i,i)=ocavs[manifold[i].j]-ocavs[manifold[i].i];
	}
	return lambda;
}

Matrix GreensCalculator::computeAMatrix(Matrix state,Matrix ev,Matrix*fullH){
	Matrix A(manifold.size(),manifold.size());

	int count = 0;
	for(int r = 0; r < A.rows();r++){
		for(int c = 0; c < A.cols();c++){
			
			PHOpp q1=manifold[r],q2=manifold[c].adjoint();
			A(r,c)=-(qqH(q1,q2,state,ev(0,0))+qqH(q2,q1,state,ev(0,0))-qHq(q1,q2,state,fullH)-qHq(q2,q1,state,fullH));
			if(count++==5){
				std::cout << qqH(q1,q2,state,ev(0,0)) << " : " << qqH(q2,q1,state,ev(0,0))<<std::endl;
				std::cout << qHq(q1,q2,state,fullH) << " : " << qHq(q2,q1,state,fullH)<<std::endl;
			}
		}
	}

	return A;
}

Matrix GreensCalculator::computeBMatrix(Matrix state, Matrix ev,Matrix * fullH){
	Matrix B(manifold.size(), manifold.size());

	for(int r = 0; r < B.rows();r++){
		for(int c = 0; c < B.cols();c++){
			PHOpp q1=manifold[r],q2=manifold[c];
			B(r,c)=qqH(q1,q2,state,ev(0,0))+qqH(q2,q1,state,ev(0,0))-qHq(q1,q2,state,fullH)-qHq(q2,q1,state,fullH);
		}
	}

	return B;
}

double GreensCalculator::qqH(PHOpp q1, PHOpp q2, Matrix state, double e0){
	StringMap & sm = SlaterDet::codes;

	unsigned char * cpy1=sm.cpystrs[0],*cpy2=sm.cpystrs[1];

	bool isalpha=(q1.i/sm.norbs)==0;

	double val=0.0;

	std::priority_queue<chslater> ch1,ch2;

	for(int i = 0; i < state.rows();i++){
		int i_aidx=i%sm.strs.size(), i_bidx=i/sm.strs.size();
		
		int rchange = isalpha?i_aidx:i_bidx;
		int rother  = isalpha?i_bidx:i_aidx;

		memcpy(cpy1,sm.strs[rchange],sm.codeblklen);

		double sign = applyPHOpp(q2,cpy1);

		if(sign==0.0) continue;
		int index = SlaterDet::decode(cpy1);

		chslater ch{index+sm.strs.size()*rother,sign,state(i,0)};
		ch1.push(ch);
	}
	for(int j = 0; j < state.rows();j++){
		int j_aidx=j%sm.strs.size(), j_bidx=j/sm.strs.size();

		int lchange = isalpha?j_aidx:j_bidx;
		int lother  = isalpha?j_bidx:j_aidx;

		memcpy(cpy2,sm.strs[lchange],sm.codeblklen);

		double sign = applyPHOpp(q1.adjoint(),cpy2);

		if(sign==0.0) continue;
		int index = SlaterDet::decode(cpy2);
		
		chslater ch{index+sm.strs.size()*lother,sign,state(j,0)};
		ch2.push(ch);
	}

	while(!(ch1.empty()||ch2.empty())){
		if(ch1.top().index==ch2.top().index){
			val+=ch1.top().sign*ch2.top().sign*ch1.top().c*ch2.top().c;
			ch1.pop();ch2.pop();
		}
		else if(ch1.top()<ch2.top()) ch2.pop();
		else ch1.pop();
	}

	return val*e0;
}

double GreensCalculator::qHq(PHOpp q1, PHOpp q2, Matrix state, Matrix * fullH){
	StringMap & sm = SlaterDet::codes;

	unsigned char * cpy1=sm.cpystrs[0],*cpy2=sm.cpystrs[1];
	bool isalpha=(q1.i/sm.norbs)==0;

	double val=0.0;

	for(int i = 0; i < state.rows();i++){
		int i_aidx=i%sm.strs.size(),i_bidx=i/sm.strs.size();

		int rchange=isalpha?i_aidx:i_bidx;
		int rother=isalpha?i_bidx:i_aidx;

		memcpy(cpy1,sm.strs[rchange],sm.codeblklen);

		double sign1=applyPHOpp(q2,cpy1);
		
		if(sign1==0.0) continue;

		int i_slater_index=SlaterDet::decode(cpy1);

		int i_index=(isalpha)?(i_slater_index+rother*sm.norbs):(rother+i_slater_index*sm.norbs);

		for(int j = 0; j < state.rows();j++){
			int j_aidx=j%sm.strs.size(),j_bidx=j/sm.strs.size();

			int lchange=isalpha?j_aidx:j_bidx;
			int lother=isalpha?j_bidx:j_aidx;

			memcpy(cpy2,sm.strs[lchange],sm.codeblklen);

			double sign2=applyPHOpp(q1.adjoint(),cpy2);

			if(sign2==0.0)continue;

			int j_slater_index=SlaterDet::decode(cpy2);

			int j_index=(isalpha)?(j_slater_index+lother*sm.norbs):(lother+j_slater_index*sm.norbs);

			val+=sign1*sign2*(*fullH)(j_index,i_index);
		}
	}

	return val;
}

//opp=a^t_j a_i
double applyPHOpp(PHOpp opp, unsigned char * sd){
	double sign = 1.0;
	int isalpha=(opp.i/SlaterDet::codes.norbs);
	
	int ibase=opp.i-SlaterDet::codes.norbs*isalpha,
	    jbase=opp.j-SlaterDet::codes.norbs*isalpha;

	int iblock=ibase/8,ioff=ibase%8,
	    jblock=jbase/8,joff=jbase%8;

	if(sd[iblock]&(1<<ioff)){
		int lowercount=(isalpha)?0:SlaterDet::codes.nelec;
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
		int lowercount=(isalpha)?0:SlaterDet::codes.nelec;
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

