#include "MBPT.h"
#include <algorithm>

diagram::diagram(int s, int e, int o):size{s},nelec{e}, norb{o}{
	for(int i = 0; i < size; i++){ structure.push_back(i);}
	nhole=size-1;
	npart=1;
	hpcount.resize(size);
	for(auto &a: hpcount) a= 0;
	chole=0;
	cpart=0;
}

bool diagram::permuteDiagram(){
	bool canonicalyGreater=std::next_permutation(++structure.begin(),structure.end());
	nhole=chole=cpart=0;
	for(auto &a: hpcount) a = 0;
	for(auto a = structure.begin(); a!=structure.end(); ){
		int l = *a;
		a++;
		int r = (a==structure.end())?*structure.begin():*a;
		if(r>l)nhole++;
	}
	npart=size-nhole;
	return canonicalyGreater;
}

bool diagram::incrementhp(){
	int hpc=hpcount.size()-1;
	while(hpc>0){
		hpcount[hpc]++;
		if(hpcount[hpc]>((hpc<nhole)?nelec/2:(norb-nelec)/2)){
			hpcount[hpc]=0;
			hpcount--;
		}
		else break;
	}
	chole=cpart=0;
}

ArbitraryMBPT::ArbitraryMBPT(){

}

ArbitraryMBPT::MBPTResults ArbitraryMBPT::computeMBPT(int degree, int nelec, hfresults & hfr, IntegralChugger & ic){
	MBPTResults mbptr;

	mbptr.degree=degree;
	mbptr.vd.resize(degree);
	mbptr.vd[0]=hfr.eelec+enuc;
	for(int d = 1; d < degree; d++){
		mbptr.vd[d]=0;
		diagram di(degree, nelec);
		do{
			mbptr.vd[d]+=calcDiagram(di, ic, hfr.G);
		}while(di.permuteDiagram());
	}

	return mbptr;
}

double ArbitraryMBPT::calcDiagram(diagram di, IntegralChugger & ic, Matrix G){
	bool notcomplete=true;
	do{
		for(auto a = di.structure.begin(); a !=di.structure.end();){
			auto n = a;
			if(++n==di.structure.end())n=di.structure.begin();
			auto nn = n;
			if(++nn==di.structure.end())nn=di.structure.begin();

			int cv=*a, nv =*n, nnv=*nn;

			bool inhole=cv<nv, outhole = nnv>nv;

			int i = (inhole)?di.hpcount[di.chole]:di.hpcount[di.cpart+di.nhole];
			int j = (outhole)?di.hpcount[di.chole+((inhole)?1:0)]:di.hpcount[di.cpart+((!inhole)?1:0)+di.nhole];

			std::cout << "EVALUATING: " << i << " " << j << std::endl;
		}
	}while(incrementhp());
	return 0.0;
}
