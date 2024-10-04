#include "PoleSearch.h"
#include <fstream>

PoleSearch::PoleSearch(FullCISolver * fci_, GreensCalculator * gr_,IntegralChugger * ic_):fci{fci_},gr{gr_},ic{ic_}{

}


std::vector<double> PoleSearch::scan(double E_s, double dE, int steps, FullCISolver::FCIResults fcir, 
		FullCISolver::MBPTResults mbptr, HartreeFockSolver::HFResults hfr,  std::string filename){
	
	std::vector<double> avocc=gr->computeAvOc_Pure(fcir.eigenvectors.col(0));
	fci->buildOperators(avocc,hfr);
	gr->buildOperators(avocc,hfr);

	std::ofstream outfile;
	if(filename!="")outfile.open(filename);

	std::vector<double> energies;

	double last = 0.0;

	Matrix orbEdiff=Matrix::Zero(fci->operators.size(),fci->operators.size());
	for(int i = 0; i < fci->operators.size();i++){
		int norbs=SlaterDet::codes.norbs;
		orbEdiff(i,i)=hfr.E(gr->operators[i].j%norbs,0)-hfr.E(gr->operators[i].i%norbs,0);
	}

	for(int i = 0; i < hfr.E.rows();i++){
		std::cout << hfr.E(i,0)<<std::endl;
	}
	std::cout << "HARTREEFOCK DIFF\n";
	for(int i = 0; i < fci->operators.size(); i++){
		int norbs=SlaterDet::codes.norbs;
		std::cout << fci->operators[i].i%norbs << " " << fci->operators[i].j%norbs << " " <<orbEdiff(i,i)<<std::endl;
	}
	//for(int i = 0; i < fcir.eigenvalues.rows(); i++){
	for(int i = 0; i < steps; i++){
		double E = E_s+dE*i;
		Matrix G=gr->ComputeGreens(E,*ic,hfr,fcir);

		Matrix G0=fci->recursivegreen(0,E,hfr,mbptr)[0];
		G-=G0;
		for(int j = 0; j < G0.rows();j++){
			G0(j,j)=1.0/(E-orbEdiff(j,j));
		}
		G+=G0;

		Matrix M = gr->ComputeSelfEnergy(E,hfr,G);

		Eigen::SelfAdjointEigenSolver<Matrix> eig_solver(orbEdiff+M);
		Matrix ems=eig_solver.eigenvalues();
		/*
		double det = 1.0;
		for(int i = 0; i < ems.size(); i++){
			det*=E-ems(i,0);
		}
		double val=det;
		*/
		double lowest = std::abs(E-ems(0,0));
		double val = E-ems(0,0);
		for(int i = 1; i < ems.rows();i++){
			if(std::abs(E-ems(i,0))<lowest){
				lowest=std::abs(E-ems(i,0));
				val=E-ems(i,0);
			}
		}
		
		if(i!=0&&(val<0.0 && last>0.0)){
			energies.push_back(E);
		}
		
		last=val;
		std::cout << E << " " << val << std::endl;

		outfile << E << " " << val << std::endl;
		//}
	}

	std::cout << "POLE SEARCH ENERGIES:\n";
	for(int i = 0; i < energies.size();i++){
		std::cout << energies[i]<<std::endl;
	}

	return energies;

}

