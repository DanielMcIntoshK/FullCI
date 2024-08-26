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

#endif

