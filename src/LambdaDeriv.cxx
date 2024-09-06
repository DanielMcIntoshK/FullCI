#include "LambdaDeriv.h"
#include "FCI.h"

Matrix LambdaDeriv::ComputeSelfEnergyOrder(int M, double x0, std::vector<double> grid, HartreeFockSolver::HFResults hfr, IntegralChugger & ic){
	int N = grid.size()-1;
	std::cout << M  <<" " << N << " " << x0 << " " ;
	for(auto x: grid)std::cout << x << " ";
	std::cout << std::endl;
	std::vector<Matrix> ldet;
	ldet.resize(N+1);
	for(int i = 0; i < ldet.size();i++){
		FullCISolver fcis;
		FullCISolver::FCIResults fcir = fcis.fci(ic,hfr,grid[i]);
		
		GreensCalculator gc;
		ldet[i]=gc.ComputeSelfEnergy(.2,ic,hfr,fcir);
		fcis.cleanup();
	}

	std::vector<Matrix> deltas=generateDeltas(M,N,x0,grid);
	
	
	//for(int n=0; n <=M;n++){
		//std::cout << deltas[n] << std::endl<<std::endl;;
	//}
	//return Matrix(0,0);

	double factorial=1.0;
	for(int n=2;n<=M;n++){
		factorial*=n;
	}

	Matrix selfenergy(ldet[0]);
	for(int i = 0; i < selfenergy.rows();i++){
	for(int j = 0; j < selfenergy.cols();j++){
		selfenergy(i,j)=0.0;
	}
	}
	for(int i =0; i <= N;i++){
		//std::cout << std::endl << ldet[i] << std::endl;

		selfenergy+=ldet[i]*deltas[M](N,i);
		//std::cout << selfenergy << std::endl;
		//int l;
		//std::cin >> l;
		//std::cout << l;
	}
	selfenergy*=1.0/factorial;

	return selfenergy;
}

std::vector<Matrix> LambdaDeriv::ComputeGreensNumerical(int M, double x0, std::vector<double> grid, HartreeFockSolver::HFResults hfr, IntegralChugger & ic, double E){
	int N = grid.size()-1;
	std::vector<Matrix> deltas{generateDeltas(M,N,x0,grid)};
	std::vector<Matrix> ldet;
	ldet.resize(N+1);
	for(int i = 0; i < ldet.size(); i++){
		FullCISolver fcis;
		FullCISolver::FCIResults fcir = fcis.fci(ic,hfr,grid[i]);

		GreensCalculator gc;
		ldet[i]=gc.ComputeGreens(E,ic,hfr,fcir);
		fcis.cleanup();
	}

	std::vector<Matrix> greens;
	greens.resize(M+1);
	greens[0]=ldet[0];
	for(int i = 0; i < ldet[0].rows(); i++){
	for(int j = 0; j < ldet[0].cols();j++){
		if(std::abs(greens[0](i,j))<0.000000001) greens[0](i,j)=0.0;
	}
	}
	double factorial=1.0;
	for(int m = 1; m <= M; m++){
		Matrix gr=Matrix::Zero(ldet[0].rows(),ldet[0].cols());
		factorial*=m;
		for(int n=0; n <=N; n++){
			gr+=ldet[n]*deltas[m](N,n);
				
		}
		gr*=1.0/factorial;
		greens[m]=gr;
	}
	return greens;
}

std::vector<Matrix> LambdaDeriv::generateDeltas(int M,int N, double x0, std::vector<double> & grid){
	std::cout << "COMPUTING DELTAS\n";
	std::vector<Matrix> deltas;
	deltas.resize(M+1);
	std::cout << "INITIALIZING\n";
	for(int i = 0; i <= M; i++){
		deltas[i]=Matrix(N+1,N+1);
		for(int r = 0; r < N+1; r++){
			for(int c = 0; c < N+1;c++){
				deltas[i](r,c)=0.0;
			}
		}
	}
	deltas[0](0,0)=1.0;
	double c1=1.0;
	std::cout << "COMPUTING\n";
	for(int n = 1; n <= N; n++){
		double c2=1.0;
		for(int v = 0; v < n;v++){
			double c3=grid[n]-grid[v];
			c2=c2*c3;
			for(int m = 0; m <= std::min(n,M);m++){
				double pprob=(m==0)?0.0:(m*deltas[m-1](n-1,v));
				deltas[m](n,v)=((grid[n]-x0)*deltas[m](n-1,v)-pprob)/c3;
			}
		}
		for(int m = 0; m <=std::min(n,M);m++){
			double pprob=(m==0)?0.0:(m*deltas[m-1](n-1,n-1));
			deltas[m](n,n)=(c1/c2)*(pprob-(grid[n-1]-x0)*deltas[m](n-1,n-1));
		}
		c1=c2;
	}

	return deltas;

}
