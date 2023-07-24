#include "HartreeFockSolver.h"
#include <stdio.h>

HartreeFockSolver::HartreeFockSolver(){

}

HartreeFockSolver::HFResults HartreeFockSolver::RestrictedHF(std::vector<Atom> ats, BasisSet bs){
	HartreeFockSolver::HFResults results;
	//results.result="HELLO WORLD!";
	strcpy(results.result, "HELLOWORLD!\0");
	return results;
}
