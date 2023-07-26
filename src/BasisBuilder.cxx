#include "BasisBuilder.h"
#include <string>
#include <sstream>
#include <fstream>

BasisBuilder::BasisBuilder(){

}

int BasisBuilder::Build(std::string filename, ModelParams & params){
	std::ifstream in(filename);

	std::vector<libint2::Atom> & atoms = params.atoms;
	if(!in){
		std::cout << "INVALID FILE\n";
		return -1;
	}
	auto elements=libint2::chemistry::get_element_info();

	std::vector<std::vector<libint2::Shell>> elementBasis;
	int maxElement=0;
	while(!in.eof()){
		std::string inputline;
		std::getline(in,inputline);
		if(!inputline.empty()&&inputline[0]!=' '){
			std::stringstream ss(inputline);
			std::string symbol;
			ss >> symbol;
			for(auto & el: elements){
				if(el.symbol==symbol && el.Z>maxElement) {maxElement=el.Z;break;}
			}
		}
	}
	elementBasis.resize(maxElement+1);
	in.close();
	in=std::ifstream(filename);

	int curElement=0;
	while(!in.eof()){
		std::string inputline;
		std::getline(in,inputline);
		if(inputline.empty()) continue;
		
		std::stringstream linestream(inputline);
		
		if(inputline[0]!=' '){
			std::string symbol;
			linestream >> symbol;
			for(auto & el: elements){
				if(el.symbol==symbol) {curElement=el.Z;break;}
			}
			//if(symbol=="Ne") curElement = 0;
		}
		else if(curElement!=0){
			std::string orbitaltype;
			int primcount;
			linestream >> orbitaltype >> primcount;
			
			int l = libint2::Shell::am_symbol_to_l(orbitaltype[0]);
			libint2::svector<double> exponents;
			libint2::svector<double> coef;
			//exponents.resize(primcount);
			//coef.resize(primcount);
			for(int i = 0; i < primcount;i++){
				std::string primline;
				std::getline(in, primline);

				std::stringstream primstream(primline);
				double e, c;
				primstream >> e >> c;

				exponents.push_back(e);
				coef.push_back(c);
				//exponents[i]=e;
				//coef[i]=c;
			}
			//elementBasis[curElement].push_back(libint2::Shell{exponents, {{l,true,coef}}, {0.0,0.0,0.0}})
			elementBasis[curElement].push_back({exponents, {{l,true,coef}},{0.0,0.0,0.0}});
		}
	}
	bs = libint2::BasisSet(atoms,elementBasis,params.cparams["BASISSET"],true);
	in.close();
	/*
	std::vector<libint2::Shell> shells;
	for(auto & at: atoms){
		std::string atSymbol=elements[at.atomic_number].symbol;

		in.seekg(0,in.beg);
		
		bool atomMatch=false;
		while(true){
			std::string inputline;
			std::getline(in,inputline);
			
			if(atomMatch && (inputline[0]!=' '||in.eof())) break;
			if(in.eof() && !atomMatch){
				std::cout << "NO BASIS FOR: " << atSymbol <<std::endl;
				return -1;
			}
			std::stringstream ss(inputline);
			if(!atomMatch){
				std::string symbol;
				ss >> symbol;

				if(inputline[0]!=' ' && symbol == atSymbol){
					atomMatch =true;
				}
			}
			else{
				std::string orbitalType;
				int primcount;
				ss >> orbitalType >> primcount;
				
				int l = libint2::Shell::am_symbol_to_l(orbitalType[0]);
				std::vector<double> exponents;
				std::vector<double> coef;
				for(int i = 0; i < primcount;i++){
					std::string primline;
					std::getline(in, primline);

					std::stringstream primstream(primline);
					double e, c;
					primstream >> e >> c;

					exponents.push_back(e);
					coef.push_back(c);
				}

				shells.push_back({exponents, {{l,true,coef}}, {at.x, at.y, at.z}});
			}
		}
	}
	bs=libint2::BasisSet(atoms,shells);
	*/
	return 0;
}
