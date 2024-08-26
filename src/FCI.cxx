#include "FCI.h"
#include "Greens.h"
#include <iostream>
#include <cmath>
#include <bitset>

FullCISolver::FCIResults FullCISolver::fci(IntegralChugger &ic, HartreeFockSolver::HFResults &hf, double lam){
	lambda=lam;
	//lambdaDeriv=(std::abs(lambda-1.0)>.000001);
	lambdaDeriv = true;

	std::cout << "RUNNING FULL CI ENERGY CALCULATION\n";
	
	ints=&ic;
	//ints->TransformInts(hf.C);
	//ints->TransformFock(hf.C,hf.G);

	//std::cout << "BUILDING SLATER DETERMINANT STRINGS\n";
	//int N=hf.C.rows();
	//int n=mp.nelec/2;
	//SlaterDet::buildStrings(N,n);
	
	strcnt=SlaterDet::codes.strs.size();
	cisize=strcnt*strcnt;
	bssize=ints->MOH.rows();

	std::cout << "USING DIRECT DIAGONALIZATION OF CI HAMILTONIAN\n";
	computeHamiltonian();	
	//std::cout << CIMat.rows() << " " << CIMat.cols() << std::endl;
	
	std::cout << "DIAGONALIZING HAMILTONIAN\n";
	Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(CIMat);

	//std::cout << CIE<<std::endl;

	FCIResults fcir;
	fcir.eigenvalues=eig_solver.eigenvalues();
	fcir.eigenvectors=eig_solver.eigenvectors();
	fcir.fullH=&CIMat;

	fciSuccess=true;

	return fcir;	
}

FullCISolver::MBPTResults FullCISolver::mbpt(HartreeFockSolver::HFResults & hf,int order){
	if(!fciSuccess){
		std::cout << "ERROR: FULLCI HAS NOT BEEN CALCULATED. THIS MUST BE DONE BEFORE MBPT CALCULATION\n";
		return FullCISolver::MBPTResults();
	}
	std::cout << "MBPT HIGH ORDER CALCULATION\n";

	StringMap & sm = SlaterDet::codes;

	std::cout << "CALCULATING H0\n";
	H0=Matrix(CIMat.rows(),CIMat.cols());
	for(int r=0; r < H0.rows(); r++){
		int alpha=r%sm.strs.size(), beta=r/sm.strs.size();
		for(int c = 0; c < H0.cols();c++){
			H0(r,c)=0.0;
			if(r==c){
				for(int i = 0; i < sm.codeblklen;i++){
					for(int b = 0; b < 8; b++){
						if(sm.strs[alpha][i]&(1<<b)) H0(r,c)+=hf.E(i*8+b,0);
						if(sm.strs[beta ][i]&(1<<b)) H0(r,c)+=hf.E(i*8+b,0);
					}
				}
			}
		}
	}
	//std::cout << H0 << std::endl;


	std::cout << "INITIALIZING ZEROTH ORDER WAVEFUNCTION AND ENERGY\n";
	Matrix y_0=Matrix::Zero(CIMat.rows(),1);
	y_0(0,0)=1.0;

	FullCISolver::MBPTResults mbptr;

	double E0=(H0*y_0)(0,0);

	mbptr.wavefunctions.push_back(y_0);
	mbptr.energies.push_back(E0);

	Matrix B=Matrix::Identity(H0.rows(),H0.cols())*E0;

	Matrix H0i=H0-B;

	std::cout << "FORMING INVERSE H0-E0\n";
	for(int i = 1; i <H0i.rows();i++) H0i(i,i)=1.0/H0i(i,i);

	std::cout << "CALCULATING HIGHER ORDER TERMS\n";
	for(int k = 1; k <= order; k++){
		Matrix sigma_k=CIMat*mbptr.wavefunctions[k-1];

		std::cout << sigma_k.transpose()*sigma_k<<std::endl;
		double Ek=sigma_k(0,0);
		if(k==1) Ek-=E0;

		mbptr.energies.push_back(Ek);

		Matrix y_k=Matrix::Zero(sigma_k.rows(),1);

		for(int r = 1; r <= k; r++){
			y_k+=mbptr.energies[r]*mbptr.wavefunctions[k-r];
		}
		y_k+=H0*mbptr.wavefunctions[k-1];
		y_k-=sigma_k;
		y_k=H0i*y_k;

		//y_k/=(y_k.transpose()*y_k)(0,0);

		//std::cout << y_0.transpose()*y_k << std::endl;

		mbptr.wavefunctions.push_back(y_k);
	}
	std::cout << "MBPT COMPLETED\n";

	return mbptr;

}

std::vector<Matrix> FullCISolver::recursivegreen(int order, double E, HartreeFockSolver::HFResults & hfr, FullCISolver::MBPTResults & mbptr){
	if(mbptr.wavefunctions.size()<order){
		std::cout << "ERROR: CANNOT PERFORM GREENS CALCULATION TO " << order << "th ORDER WITH ONLY " << mbptr.wavefunctions.size() << " MBPT ORDER CALCULATION...\nMUST BE GREATER THAN OR EQUAL TO DESIRED GREENS FUNCTION ORDER\n";
		return std::vector<Matrix>();
	}

	//double pm=(pos)?1.0:-1.0;

	StringMap & sm = SlaterDet::codes;
	int nstrings=sm.strs.size()*sm.strs.size();
	
	std::vector<Matrix> Gmiddlep,Gmiddlen;
	Gmiddlep.resize(order+1);
	Gmiddlen.resize(order+1);

	Matrix V=CIMat-H0;

	//First the omega dependent part
	Gmiddlep[0]=Matrix::Zero(nstrings,nstrings);
	Gmiddlen[0]=Matrix::Zero(nstrings,nstrings);
	for(int i = 0; i < nstrings;i++){
		Gmiddlep[0](i,i)=1.0/(E+H0(0,0)-H0(i,i));
		Gmiddlen[0](i,i)=1.0/(E-H0(0,0)+H0(i,i));	
	}

	for(int n = 1; n <=order;n++){
		Gmiddlep[n]=Matrix::Zero(nstrings,nstrings);
		Gmiddlep[n]+=Gmiddlep[n-1]*V*Gmiddlep[0];
		
		Gmiddlen[n]=Matrix::Zero(nstrings,nstrings);
		Gmiddlen[n]+=-Gmiddlen[n-1]*V*Gmiddlen[0];

		for(int k =1; k <= n; k++) {
			Gmiddlep[n]+=-Gmiddlep[n-k]*mbptr.energies[k]*Gmiddlep[0];
			Gmiddlen[n]+= Gmiddlen[n-k]*mbptr.energies[k]*Gmiddlen[0];
		}
	}

	//Now the full Green's function
	std::vector<PHOp> operators;
	for(int s = 0; s<2; s++){
		for(int i = 0; i < sm.norbs; i++){
			for(int j = 0; j < sm.norbs; j++){
				operators.push_back(PHOp(i+sm.norbs*s,j+sm.norbs*s));	
			}
		}
	}

	std::vector<Matrix> G_n;
	G_n.resize(order+1);

	std::vector<Matrix> z,za;
	std::vector<double> D;
	z.resize(order+1);
	za.resize(order+1);
	D.resize(order+1);
	
	for(int n = 0; n <= order; n++){
		std::cout << "CALCULATING G_" << n << std::endl;
		std::cout << "CALCULATING z\n";
		z[n]=Matrix::Zero(nstrings,operators.size());
		za[n]=Matrix::Zero(nstrings,operators.size());
		for(int r = 0; r <z[n].rows(); r++){
			for(int c = 0; c < z[n].cols(); c++){
				double sign=1.0;
				double signa=1.0;

				int index=opOnSlater(operators[c],r);
				int indexa=opOnSlater(operators[c].adjoint(),r);

				if(index>=nstrings){
					index=0;
					sign=0.0;
				}
				if(indexa>=nstrings){
					indexa=0;
					signa=0.0;
				}
				if(index<0){
					index=-index;
					sign=-1.0;
				}
				if(indexa<0){
					indexa=-indexa;
					signa=-1.0;
				}
				z[n](r,c) +=sign*mbptr.wavefunctions[n](index,0);
				za[n](r,c)+=signa*mbptr.wavefunctions[n](indexa,0);
			}
		}
		std::cout << "CALCULATING D\n";
		D[n]=0.0;
		for(int i=0; i <= n; i++){
			D[n]+=(mbptr.wavefunctions[i].transpose()*mbptr.wavefunctions[n-i])(0,0);
		}
		std::cout << "CALCULATING G\n";
		G_n[n]=Matrix::Zero(operators.size(),operators.size());
		Matrix G_f=Matrix::Zero(operators.size(), operators.size()),
		       G_r=Matrix::Zero(operators.size(), operators.size());
		for(int i = 0; i <= n; i++){
		for(int j = 0; j <= n-i; j++){
				G_f+=za[i].transpose()*Gmiddlen[j]*za[n-i-j];
				G_r+=z[i].transpose()*Gmiddlep[j]*z[n-i-j];
				G_n[n]+=G_f-G_r.transpose();
			}
			if(i>0) G_n[n]-=D[i]*G_n[n-i];
		}
	}
	std::cout << std::setprecision(3) <<z[0] << std::endl << std::endl << za[0] << std::endl << std::endl; 

	std::cout << "HARTREE FOCK ENERGIES\n";
	for(int i =0; i < hfr.E.rows(); i++){
		std::cout << i << " " <<  hfr.E(i,0)<<std::endl;

	}
	Matrix G_0=Matrix::Zero(operators.size(),operators.size());
	for(int i = 0; i < G_0.rows();i++){
		int index1=opOnSlater(operators[i],0);
		int index2=opOnSlater(operators[i].adjoint(),0);
		std::cout << "G_0 values: " << operators[i].i << " " << operators[i].j << " ";
		if(index1<nstrings){
			std::cout << "1:"<<E<<":"<<hfr.E(operators[i].i%sm.norbs,0)<<":"<<hfr.E(operators[i].j%sm.norbs)<<":";
			G_0(i,i)+=1.0/(E-hfr.E(operators[i].i%sm.norbs,0)+hfr.E(operators[i].j%sm.norbs,0));
			std::cout <<  G_0(i,i) << " ";
		}
		if(index2<nstrings){
			std::cout << "2:"<<E<<":"<<hfr.E(operators[i].i%sm.norbs,0)<<":"<<hfr.E(operators[i].j%sm.norbs)<<":";
			G_0(i,i)-=1.0/(E+hfr.E(operators[i].i%sm.norbs,0)-hfr.E(operators[i].j%sm.norbs,0));
			std::cout << G_0(i,i);
		}
		std::cout << std::endl;
	}

	std::cout <<"G0 SHOULD BE: " << std::endl;
	for(int i = 0; i < G_0.rows();i++){
		std::cout << operators[i].i << " " << operators[i].j << " "<<G_0(i,i)<<std::endl;
		
	}

	std::cout << "CALCULATED: " << std::endl;
	for(int i = 0; i < G_n[0].rows();i++){
		std::cout << operators[i].i << " " << operators[i].j << " "<<G_n[0](i,i)<<std::endl;
	}
       
	return G_n;
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
	//SlaterDet::cleanStrings();
	CIMat=Matrix(0,0);
}

int FullCISolver::opOnSlater(PHOp op, int det,bool verbose){
	StringMap & sm = SlaterDet::codes;
	int alphaidx=det%sm.strs.size(),betaidx=det/sm.strs.size();
	
	unsigned char * alphacpy = sm.cpystrs[0],*betacpy=sm.cpystrs[1];

	memcpy(alphacpy,sm.strs[alphaidx],sm.codeblklen);
	memcpy(betacpy,sm.strs[betaidx],sm.codeblklen);
	

	double sign = 1.0;

	int isalpha=(op.i/SlaterDet::codes.norbs);

	int ibase = op.i-SlaterDet::codes.norbs*isalpha,
	    jbase = op.j-SlaterDet::codes.norbs*isalpha;

	unsigned char * sd = (isalpha==0)?alphacpy:betacpy;
	
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
	else return sm.strs.size()*sm.strs.size();

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
	else return sm.strs.size()*sm.strs.size();

	return sign*(SlaterDet::decode(alphacpy)+sm.strs.size()*SlaterDet::decode(betacpy));
}

