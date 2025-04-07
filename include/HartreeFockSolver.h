#ifndef __HARTREEFOCKSOLVER__H_
#define __HARTREEFOCKSOLVER__H_

#include <vector>
#include <array>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <libint2.hpp>

#include "ModelParams.h"
#include "Integrals.h"

using namespace libint2;

struct PHOp{
	PHOp():i{0},j{0}{}
	PHOp(int a, int b):i{a},j{b}{}
	int i,j;

	PHOp adjoint(){return PHOp(j,i);}
};

//Exactly what it says on the tin
//Class used for solving Hartree Fock
class HartreeFockSolver{
public:
	struct HFResults{
		//The electic energy of the hartree fock ground state
		double eelec;
		//The energy associated with nuclear repulsion
		double enuc;

		//The coeficent matrix, density matrix, the most recently calculated fock matrix - hamiltonian(used for MBPT), and eigenvalues of fock matrix of last iteration
		Matrix C;
		Matrix D;
		Matrix G;
		Matrix F;
		Matrix E;

		int nelec;
		int norbs;

		std::vector<PHOp> operators;

		void buildOperators();
		Matrix getG0(double E);
		Matrix getG0i(double E);
	};
public:
	HartreeFockSolver();

	//Runs A SCF iterative method to minimize 
	HartreeFockSolver::HFResults RestrictedHF(ModelParams & param, BasisSet & bs, IntegralChugger & ic);
private:
	//Takes a density matrix and comptues the G matrix where G=F-H where F is the Fock matrix
	//and H is the core hamiltonian matrix
	Matrix computeGMatrix(Matrix & D, IntegralChugger & ic);
	
	//This makes an initial Guess for the density matrix
	//(At the moment probably more complicated than it nees to be
	Matrix initialDGuess(ModelParams & param,Matrix & H, Matrix &S); 

	int nelectron;
};

typedef HartreeFockSolver::HFResults hfresults;

class RPASolver{
public:
	struct RPAResults{
		std::vector<double> EE;
		std::vector<double> sing;
		std::vector<double> trip;
	};
	struct HRPAResults{
		std::vector<double> EE;
		std::array<std::vector<double>,2> singtrip;
		std::vector<Matrix> Cs;
	};
public:
	RPASolver(IntegralChugger & icref);

	RPASolver::RPAResults RPACalc(hfresults & hfr, bool TammDancoff=false);
	RPASolver::RPAResults RPACalcSingTrip(hfresults & hfr, bool TammDancoff=false);
	RPASolver::HRPAResults HRPACalc(hfresults & hfr, double threshold, int terminate, bool USESCF=true);
private:
	Matrix A0();
	Matrix B0();

	Matrix A0(int S);
	Matrix B0(int S);
	
	Matrix A1(std::vector<Matrix> C, int S);
	Matrix B1(std::vector<Matrix> C, int S);
	Matrix A1new(Matrix Tp, Matrix Th);
	Matrix B1new(Matrix S, Matrix X, int s);

	Matrix S(Matrix C0);
	Matrix X(Matrix CS);
	Matrix Tp(Matrix C0);
	Matrix Th(Matrix C0);


	Matrix constructC(Matrix Y, Matrix Z);
	std::vector<Matrix> initC(std::vector<Matrix> Bs);

	std::vector<PHOp> operators;
	hfresults hfr;
	IntegralChugger & ic;

	int decode(int i, int a);
};

#endif

