#include "PoleSearch.h"
#include <fstream>

PoleSearch::PoleSearch(FullCISolver * fci_, GreensCalculator * gr_,IntegralChugger * ic_, HartreeFockSolver::HFResults hfr_,
		FullCISolver::FCIResults fcir_,FullCISolver::MBPTResults mbptr_):
	fci{fci_},gr{gr_},ic{ic_},fcir{fcir_},mbptr{mbptr_},hfr{hfr_}{
	
	int norbs=SlaterDet::codes.norbs;
	orbEdiff=Matrix::Zero(hfr.operators.size(),hfr.operators.size());
	for(int i = 0; i < hfr.operators.size();i++) orbEdiff(i,i)=hfr.E(hfr.operators[i].j%norbs,0)-hfr.E(hfr.operators[i].i%norbs,0);
}


std::vector<double> PoleSearch::refinePoints(std::vector<double> EE, double threshold){
	std::vector<double> energies;
	energies.resize(EE.size());
	Matrix G0i=hfr.getIndependentG0i();
	for(int i = 0; i < energies.size();i++){
		energies[i]=EE[i];

		double diff=0.0;
		int n=0;
		do{
			Matrix igreen=gr->ComputeGreens(energies[i], *ic,hfr,fcir);
			Matrix iM=gr->ComputeSelfEnergy(energies[i],hfr,igreen);
			//Matrix md=(Matrix::Identity(G0i.rows(),G0i.cols()))*energies[i]-G0i+iM;
			//Matrix md=iM-G0i;
			Matrix md=iM-G0i;
			gr->TransformEigen(md,hfr);

			Eigen::EigenSolver<Matrix> eig_solver(md);

			double lowest=std::abs(energies[i]-eig_solver.eigenvalues()(0,0).real());
			int lowesti=0;
			for(int k = 1; k < eig_solver.eigenvalues().rows();k++){
				double idiff=std::abs(energies[i]-eig_solver.eigenvalues()(k,0).real());
				if(idiff<lowest){
					lowest=idiff;
					lowesti=k;
				}
			}
			diff=lowest;
			energies[i]=eig_solver.eigenvalues()(lowesti,0).real();
		}while((diff>threshold)&&(threshold>=0.0)&&(n++)<50);
		if(n>50){std::cout << "WARNING THRESHOLD NOT REACHED\n";}
	}
	return energies;
}

std::vector<double> PoleSearch::refinePointsOrder(std::vector<double> EE, double threshold,int order){
	StringMap & sm=SlaterDet::codes;

	std::vector<double> energies;
	energies.resize(EE.size());

	std::vector<PHOp> &operators=hfr.operators;

	Matrix G0i=hfr.getIndependentG0i();
	for(int i = 0; i < energies.size();i++){
		energies[i]=EE[i];

		double diff=0.0;
		int n = 0;
		do{

			std::vector<Matrix> igreens=fci->recursivegreen(order,energies[i],hfr,mbptr);
			std::vector<Matrix> ise=gr->ComputeSelfEnergies(energies[i],order,hfr,igreens);

			Matrix iseSum=ise[1];
			for(int i=2; i <= order;i++){iseSum+=ise[i];}

			//Matrix md=(Matrix::Identity(G0i.rows(),G0i.cols()))*energies[i]-G0i+iseSum;
			//Matrix md=iseSum-G0i;
			Matrix md=iseSum-G0i;
			gr->TransformEigen(md,hfr);

			Eigen::EigenSolver<Matrix> eig_solver(md);

			double lowest=std::abs(energies[i]-eig_solver.eigenvalues()(0,0).real());
			int lowesti=0;
			for(int k = 1; k < eig_solver.eigenvalues().rows();k++){
				double idiff=std::abs(energies[i]-eig_solver.eigenvalues()(k,0).real());
				if(idiff<lowest){
					lowest=idiff;
					lowesti=k;
				}
			}
			
			diff=lowest;
			energies[i]=eig_solver.eigenvalues()(lowesti,0).real();
		}while((diff>threshold)&&(threshold>=0.0)&&(n++)<50);
		if(n>50){std::cout << "WARNING THRESHOLD NOT REACHED ON "<<i<<std::endl;}
	}
	return energies;
	
}


std::vector<double> PoleSearch::scan(double E_s, double dE, int steps,std::string filename){
	
	//std::vector<double> avocc=gr->computeAvOc_Pure(fcir.eigenvectors.col(0));

	std::ofstream outfile;
	if(filename!="")outfile.open(filename);

	std::vector<double> energies;

	double last = 0.0;

	std::cout << "HARTREEFOCK DIFF\n";
	for(int i = 0; i < hfr.operators.size(); i++){
		int norbs=SlaterDet::codes.norbs;
		std::cout << hfr.operators[i].i%norbs << " " << hfr.operators[i].j%norbs << " " <<orbEdiff(i,i)<<std::endl;
	}
	//for(int i = 0; i < fcir.eigenvalues.rows(); i++){
	for(int i = 0; i < steps; i++){
		double E = E_s+dE*i;
		double val=detval(E);

		if(i!=0&&(val<0.0 && last>0.0)){
			energies.push_back(E);
		}
		
		last=val;
		//std::cout << E << " " << val << std::endl;
		

		outfile << E << " " << val << std::endl;
		//}
	}

	std::cout << "POLE SEARCH ENERGIES:\n";
	for(int i = 0; i < energies.size();i++){
		double erefine=refinePoint(energies[i],0.000001);
		std::cout << energies[i]<< " " << erefine<<std::endl;
	}

	return energies;

}

double PoleSearch::refinePoint(double point,double threshold){
	int emergencycount=0;
	double last=0.0,now=point;
	do{
		last=now;

		double p1=detval(now);
		double p2=detval(now+0.00001);

		double slope = (p2-p1)/0.00001;
		double b = p1+-now*slope;

		now=-b/slope;
		//std::cout << std::setprecision(8) << now << " " << std::setprecision(8) <<" " << last<< " " << now-last<<std::endl;

		if(emergencycount++>10)break;
	}while(std::abs(now-last)>threshold);
	if(emergencycount>10){
		std::cout << "DIDNT CONVERGE\n";
	}
	return now;
}

double PoleSearch::detval(double E){
	Matrix G=gr->ComputeGreens(E,*ic,hfr,fcir);
	//Matrix G0=fci->recursivegreen(0,E,hfr,mbptr)[0];
	Matrix G0=hfr.getG0(E);
	G-=G0;
	for(int i = 0; i < G0.rows();i++) G0(i,i)=1.0/(E-orbEdiff(i,i));
	G+=G0;

	Matrix G0i=G0;
	for(int i = 0; i < G0.rows();i++) G0i(i,i)=1.0/G0(i,i);

	//Matrix M = gr->ComputeSelfEnergy(E,hfr,G);
	//Eigen::FullPivLU<Matrix> luinverse(G0i-M);
	
	//return luinverse.determinant();
	return 0;
}

void PoleSearch::mapSelfEnergy(double E_s, double dE, int steps, std::string filename){

	std::vector<std::vector<double>>evs;
	evs.resize(steps);
	int cpercent=-1;
	Matrix G0i=hfr.getIndependentG0i();
	for(int n = 0; n < steps; n++){
		int percent=(int)(100.0*(double)n/(double)steps);
		int fivepercent=5*(percent/5);
		if(fivepercent>cpercent){
			cpercent=fivepercent;
			std::cout << cpercent<<"%"<<std::endl;
		}
		double E_n=E_s+(double)n*dE;
		
		Matrix igreen=gr->ComputeGreens(E_n, *ic,hfr,fcir);
		Matrix iM=gr->ComputeSelfEnergy(E_n,hfr,igreen);
		//Matrix md=(Matrix::Identity(G0i.rows(),G0i.cols()))*E_n-G0i+iM;
		Matrix md=iM-G0i;
		
		gr->TransformEigen(md,hfr);

		evs[n].resize(md.rows());
		for(int i = 0; i < md.rows();i++){
			evs[n][i]=md(i,i);
		}
		continue;
		Eigen::EigenSolver<Matrix> eig_solver(md);

		
		double lowest=std::abs(E_n-eig_solver.eigenvalues()(0,0).real());
		int lowesti=0;
		for(int k = 1; k < eig_solver.eigenvalues().rows();k++){
			double idiff=std::abs(E_n-eig_solver.eigenvalues()(k,0).real());
			if(idiff<lowest){
				lowest=idiff;
				lowesti=k;
			}
		}
		//evs[n]=eig_solver.eigenvalues()(lowesti,0).real();
		
		evs[n].resize(eig_solver.eigenvalues().rows());
		for(int i = 0; i < eig_solver.eigenvalues().rows();i++){
			evs[n][i]=eig_solver.eigenvalues()(i,0).real();
		}
	}
	std::ofstream outfile(filename);
	for(int n = 0; n < steps; n++){
		double E_n = E_s+(double)n*dE;
		outfile<<E_n;
		//outfile << E_n<<" " << evs[n];
		for(int i = 0; i < evs[n].size(); i++){
			outfile<<" "<<evs[n][i];
		}
		outfile << std::endl;
	}
	outfile.close();
}

