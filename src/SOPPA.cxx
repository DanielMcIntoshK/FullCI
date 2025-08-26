#include "SOPPA.h"

SOPPASolver::SOPPASolver(hfresults results, IntegralChugger * ints):hfr{results},ic{ints}{
}

Matrix SOPPASolver::FOPPA(double w){
	E=w;
	Matrix A0=ssmat(ABLAB_A0);
	Matrix As=ssmat(ABLAB_Astar);
	Matrix A1=ssmat(ABLAB_A1);
	Matrix B=ssmat(ABLAB_B1);

	std::cout << A0.rows() << " " << A0.cols() << std::endl <<
		As.rows() << " " << As.cols() << std::endl <<
		A1.rows() << " " << A1.cols() << std::endl <<
		B.rows() << " " << B.cols() << std::endl;

	Matrix foppa(A0.rows()+As.rows(),A0.cols()+As.cols());
	for(int r = 0; r < A0.rows(); r++){
	for(int c = 0; c < A0.cols(); c++){
		foppa(r,c)=A0(r,c)-A1(r,c);
		foppa(r+A0.rows(),c+A0.cols())=As(r,c)-A1(r,c);
		foppa(r+A0.rows(),c)=-B(r,c);
		foppa(r,c+A0.cols())=-B(r,c);
	}
	}
	return foppa;
}

Matrix SOPPASolver::SOPPA(double w){
	E=w;
	Matrix A0=ssmat(ABLAB_A0);
	Matrix As=ssmat(ABLAB_A0);
	Matrix A1=ssmat(ABLAB_A0);
	Matrix A2=ssmat(ABLAB_A0);

	Matrix B1=ssmat(ABLAB_B1);
	Matrix B2=ssmat(ABLAB_B2);

	Matrix C1=dsmat(CLAB_C1);

	Matrix D0=ddmat(DLAB_D0);

	Matrix Axt=A0-A1-A2;
	Matrix AxS=As-A1-A2;

	Matrix Bx=B1+B2;

	Matrix tl=Axt-C1.transpose()*(Matrix::Identity(A0.rows(),A0.cols())*E-D0)*C1;
	Matrix br=AxS-C1.transpose()*(-Matrix::Identity(A0.rows(),A0.cols())*E-D0)*C1;


	Matrix soppa(A0.rows()+As.rows(),A0.cols()+As.cols());
	for(int r = 0; r < A0.rows(); r++){
	for(int c = 0; c < A0.cols(); c++){
		soppa(r,c)=tl(r,c);
		soppa(r+A0.rows(),c+A0.cols())=br(r,c);
		soppa(r+A0.rows(),c)=Bx(r,c);
		soppa(r,c+A0.rows())=Bx(r,c);
	}}
	
	return soppa;
}

Matrix SOPPASolver::ssmat(int i){
	Matrix ss(hfr.operators.size()/2, hfr.operators.size()/2);
	
	double (SOPPASolver::*matel)(PHOp, PHOp)=nullptr;
	switch(i){
	case ABLAB_A0:matel=&SOPPASolver::matelA0;break;
	case ABLAB_A1:matel=&SOPPASolver::matelA1;break;
	case ABLAB_A2:matel=&SOPPASolver::matelA2;break;
	case ABLAB_Astar:matel=&SOPPASolver::matelAstar;break;
	case ABLAB_B1:matel=&SOPPASolver::matelB1;break;
	case ABLAB_B2:matel=&SOPPASolver::matelB2;break;
	default:break;
	}

	for(int r = 0; r < ss.rows(); r++){
		for(int c = 0; c < ss.cols(); c++){
			ss(r,c)=(*this.*matel)(hfr.operators[r],hfr.operators[c]);
		}
	}
	return ss;
}

Matrix SOPPASolver::ddmat(int i){
	Matrix dd((hfr.operators.size()*hfr.operators.size())/4, 
			(hfr.operators.size()*hfr.operators.size())/4);
	
	double (SOPPASolver::*matel)(PHOp,PHOp,PHOp, PHOp)=nullptr;
	switch(i){
	case DLAB_D0:matel=&SOPPASolver::matelD0;break;
	case DLAB_D1:matel=&SOPPASolver::matelD1;break;
	default:break;
	}

	for(int r = 0; r < dd.rows(); r++){
		int r1=r%(hfr.operators.size()/2),
		    r2=r/(hfr.operators.size()/2);
		for(int c = 0; c < dd.cols(); c++){
			int c1=c%(hfr.operators.size()/2),
		    	    c2=c/(hfr.operators.size()/2);
			
			dd(r,c)=(*this.*matel)(hfr.operators[r1],hfr.operators[r2],
						hfr.operators[c1],hfr.operators[c2]);
		}
	}
	return dd;
}

Matrix SOPPASolver::dsmat(int i){
	Matrix ds((hfr.operators.size()*hfr.operators.size())/4,
		hfr.operators.size()/2);	
	
	double (SOPPASolver::*matel)(PHOp,PHOp,PHOp)=nullptr;
	switch(i){
	case CLAB_C1:matel=&SOPPASolver::matelC1;break;
	default:break;
	}

	for(int r = 0; r < ds.rows(); r++){
		int r1=r%(hfr.operators.size()/2),
		    r2=r/(hfr.operators.size()/2);
		for(int c = 0; c < ds.cols(); c++){
			ds(r,c)=(*this.*matel)(hfr.operators[r1],hfr.operators[r2],
					hfr.operators[c]);
		}
	}
	return ds;
}

double SOPPASolver::matelA0(PHOp r, PHOp c){
	if((r.i==c.i)&&(r.j==c.j)) return (E+hfr.E(r.j%hfr.norbs,0)-hfr.E(r.i%hfr.norbs,0));
	return 0;
}
double SOPPASolver::matelA1(PHOp r, PHOp c){
	return ic->movsymphys(r.i,c.j,c.i,r.j,hfr.norbs);

}
double SOPPASolver::matelA2(PHOp r, PHOp c){
	double sval=0.0;

	if(r.j==c.j){
		for(int q = 0; q < hfr.phlist[0].size(); q++){
		for(int x = 0; x < hfr.phlist[1].size(); x++){
		for(int y = 0; y < hfr.phlist[1].size(); y++){
			int qop=hfr.phlist[0][q];
			int xop=hfr.phlist[1][x];
			int yop=hfr.phlist[1][y];

			double eps=hfr.E(xop%hfr.norbs,0)+hfr.E(yop%hfr.norbs,0)
				-hfr.E(r.i%hfr.norbs,0)-hfr.E(qop%hfr.norbs,0);
			double num=ic->movsymphys(r.i,qop,xop,yop,hfr.norbs)*
				ic->movsymphys(xop,yop,qop,c.i,hfr.norbs);
			sval+=-num/eps;
		}}}
	}
	if(r.i==c.i){
		for(int q = 0; q < hfr.phlist[0].size(); q++){
		for(int p = 0; p < hfr.phlist[0].size(); p++){
		for(int x = 0; x < hfr.phlist[1].size(); x++){
			int qop=hfr.phlist[0][q];
			int pop=hfr.phlist[0][p];
			int xop=hfr.phlist[1][x];

			double eps=hfr.E(r.j%hfr.norbs,0)+hfr.E(xop%hfr.norbs,0)
				-hfr.E(qop%hfr.norbs,0)-hfr.E(pop%hfr.norbs,0);
			double num=ic->movsymphys(c.j,xop,qop,pop,hfr.norbs)*
				ic->movsymphys(qop,pop,r.j,xop,hfr.norbs);
			sval+=num/eps;
		}}}
	}
	return sval;
}
double SOPPASolver::matelAstar(PHOp r, PHOp c){
	if((r.i==c.i)&&(r.j==c.j)) return (-E+hfr.E(r.j%hfr.norbs,0)-hfr.E(r.i%hfr.norbs,0));
	return 0;
}
double SOPPASolver::matelB1(PHOp r, PHOp c){
	return ic->movsymphys(c.j,r.j,r.i,c.i,hfr.norbs);
}
double SOPPASolver::matelB2(PHOp r, PHOp c){
	double sval=0.0;

	for(int q=0; q < hfr.phlist[0].size(); q++){
	int qop=hfr.phlist[0][q];
	for(int x=0; x < hfr.phlist[1].size(); x++){
	int xop=hfr.phlist[1][x];
		
		double eps1=hfr.E(xop%hfr.norbs,0)+hfr.E(r.j%hfr.norbs,0)
			-hfr.E(c.i%hfr.norbs,0)-hfr.E(qop%hfr.norbs,0);
		double eps2=hfr.E(c.j%hfr.norbs,0)+hfr.E(xop%hfr.norbs,0)
			-hfr.E(r.i%hfr.norbs,0)-hfr.E(qop%hfr.norbs,0);

		double num1=ic->movsymphys(c.j,qop,xop,r.i,hfr.norbs)*
			ic->movsymphys(qop,c.i,xop,r.j,hfr.norbs);
		double num2=ic->movsymphys(r.j,qop,xop,c.i,hfr.norbs)*
			ic->movsymphys(r.i,qop,c.j,xop,hfr.norbs);

		sval-=num1/eps1+num2/eps2;
	}}

	for(int p=0; p < hfr.phlist[0].size(); p++){
	int pop=hfr.phlist[0][p];
	for(int q=0; q < hfr.phlist[0].size(); q++){
	int qop=hfr.phlist[0][q];
		
		double eps=hfr.E(c.j%hfr.norbs,0)+hfr.E(r.j%hfr.norbs,0)
			-hfr.E(pop%hfr.norbs)-hfr.E(qop%hfr.norbs);
		double num=ic->movsymphys(pop,qop,c.i,r.i,hfr.norbs)*
			ic->movsymphys(pop,qop,c.j,r.j,hfr.norbs);

		sval-=0.5*num/eps;
	}}
	
	for(int x=0; x < hfr.phlist[1].size(); x++){
	int xop=hfr.phlist[0][x];
	for(int y=0; y < hfr.phlist[1].size(); y++){
	int yop=hfr.phlist[1][y];
		
		double eps=hfr.E(xop%hfr.norbs,0)+hfr.E(yop%hfr.norbs,0)
			-hfr.E(r.i%hfr.norbs)-hfr.E(c.i%hfr.norbs);
		double num=ic->movsymphys(r.j,c.j,xop,yop,hfr.norbs)*
			ic->movsymphys(r.i,c.i,xop,yop,hfr.norbs);

		sval-=0.5*num/eps;
	}}

	return sval;
}

double SOPPASolver::matelD0(PHOp r1, PHOp r2, PHOp c1, PHOp c2){
	if((r1.i==c1.i)&&(r2.i==c2.i)&&(r1.j==c1.j)&&(r2.j==c2.j)){
		return E+(hfr.E(r2.j%hfr.norbs,0)+hfr.E(r1.j%hfr.norbs,0)
				-hfr.E(r2.i%hfr.norbs,0)-hfr.E(r1.i%hfr.norbs,0));
	}
	return 0.0;
}
double SOPPASolver::matelD1(PHOp r1, PHOp r2, PHOp c1, PHOp c2){
	double sval=0.0;

	sval+=(r1.i==c1.i)? f(r2.i,c2.i,r1,r2,c1,c2):0.0;
	sval+=(r2.i==c2.i)? f(r1.i,c1.i,r1,r2,c1,c2):0.0;
	sval+=(r1.i==c2.i)?-f(r2.i,c1.i,r1,r2,c1,c2):0.0;
	sval+=(r2.i==c1.i)?-f(r1.i,c2.i,r1,r2,c1,c2):0.0;
	sval+=((r1.i==c1.i)&&(r2.i==c2.i))?ic->movsymphys(c1.j,c2.j,r2.j,r1.j,hfr.norbs):0.0;
	sval+=((r2.j==c2.j)&&(r1.j==c1.j))?ic->movsymphys(r1.i,r2.i,c2.i,c1.i,hfr.norbs):0.0;

	return sval;
}

double SOPPASolver::matelC1(PHOp r1, PHOp r2, PHOp c){
	double sval=0.0;

	sval+=(r2.i==c.i)?-ic->movsymphys(c.j,r1.i,r2.j,r1.j,hfr.norbs):0.0;
	sval+=(r1.i==c.i)? ic->movsymphys(c.j,r2.i,r2.j,r1.j,hfr.norbs):0.0;
	sval+=(r2.j==c.j)?-ic->movsymphys(r2.i,r1.i,r1.j,c.i,hfr.norbs):0.0;
	sval+=(r1.j==c.j)? ic->movsymphys(r2.i,r1.i,r2.j,c.i,hfr.norbs):0.0;

	return sval;
}

double SOPPASolver::f(int m, int p, PHOp r1, PHOp r2, PHOp c1, PHOp c2){
	double sval = 0.0;

	sval+=(r1.j==c1.j)? ic->movsymphys(c2.j,m,r2.j,p,hfr.norbs):0.0;
	sval+=(r2.j==c2.j)? ic->movsymphys(c1.j,m,r1.j,p,hfr.norbs):0.0;
	sval+=(c1.j==r2.j)?-ic->movsymphys(c2.j,m,r1.j,p,hfr.norbs):0.0;
	sval+=(r1.j==c2.j)?-ic->movsymphys(c1.j,m,r2.j,p,hfr.norbs):0.0;

	return sval;
}
