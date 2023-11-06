#include "ModelParams.h"
#include <vector>
#include <sstream>
#include <libint2/chemistry/elements.h>

ModelParams::ModelParams(){
	iparams=std::map<std::string,int>{{"CHARGE",0},{"KPOINTS",0},{"JOBTYPE",0},{"PRINT",0},
		{"CUTOFF1",0},{"CUTOFF2",0},{"SCFCYCLES",50},{"DIIS",5},{"MULTIPOLE",0},{"FROZENCORE",0},
		{"CIS_ROOTS",0},{"RPA_ROOTS",0},{"GDIIS",4},{"QPBANDS",0},{"MP2KJOBS",99999},
		{"MP2KJOBE",99999},{"HIGHMP",0},{"HIGHCIS",2},{"HIGHCIE",0},{"HIGHCCS",2},{"HIGHCCE",0},
		{"CCDIIS",3},{"FROZENVIRT",0},{"HIGHCIROOT",1},{"EOMORDERS",0},{"EOMORDERE",0},
		{"EOMROOT",0},{"CCALG",1},{"CCPTROOT",0},{"OEPALG",1},{"POTDUMP",0},{"DENDUMP",0},
		{"RADIAL",0},{"ANGULAR",0},{"SALG",3},{"HIGHGF",0},{"MODULO",1},{"NGRID",3},{"LMAX",-1},
		{"MBGFALG",4},{"CUTOFF3",0},{"POLARALG",1}};
	dparams=std::map<std::string,double>{{"PERIOD",0.0},{"HELIX",0.0},{"INTTOL",1.0e-10},
		{"SCFCONV",1.0e-8},{"SCFRELAX",1.0},{"HFEXCHANGE",1.0},{"SLATER",0.0},{"VWN",0.0},
		{"BECKE88",0.0},{"LYP",0.0},{"MP2",0.0},{"MAXDISK",100.0},{"MAXMEMORY",32.0},
		{"DAVIDSON",1.0e-6},{"OPTCONV",0.4e-3},{"DYSONDAMP",0.0},{"CCCONV",1.0e-6},
		{"CCRELAX",1.0},{"CICONV",1.0e-5},{"CCPTCONV",1.0e-8},{"KERNAL_HF",1.0e99},
		{"KERNAL_S",1.0e99},{"KERNAL_VWN",1.0e99},{"KERNAL_B88",1.0e99},{"KERNAL_LYP",1.0e99},
		{"KERNAL_OEP",1.0e99},{"LINDEP",1.0e-6},{"OMEGA",0.0},{"SINGULAR",1.0e-7},{"TEMP",0.0},
		{"DELTAH",0.3}};
	cparams=std::map<std::string,std::string>{{"JOBNAME","NONE"},{"BASISSET",""},{"UNITS","BOHR"},
		{"AUXILARY","NONE"},{"THEORY","SPECIFIED BELOW"},{"CCTHEORY","NONE"},{"GEMINAL","F12"}};
	lparams=std::map<std::string,bool>{{"DIRECT",true},{"VAPPROX",false},{"SAPPROX",false},
		{"DYSON",false},{"RI_SCF",false},{"RI_MP2",false},{"SINGLET",true},{"TRIPLET",true},
		{"CISQP",false},{"OEP",false},{"IP",false},{"EA",false},{"XCC",false},
		{"SPHERICAL",false},{"SLATER51",false},{"KLI",false},{"AC",false},{"POLAR",false},
		{"POLARX",true},{"POLARY",true},{"POLARZ",true},{"SOS",false},{"RELATIVITY",false},
		{"SINANOGLU",false}};
	//std::cout << iparams.size() + dparams.size() + cparams.size() + lparams.size() << std::endl;
}

int ModelParams::ReadInputFile(std::istream & in){
	std::vector<std::string> cThoeryName{"HF","SNULL","BNULL","SVWN","BLYP",
		"SLYP","BVWN","B3LYP","MP2FC","MP2FULL"};
	
	int nopt = iparams.size()+dparams.size()+cparams.size()+lparams.size();

	std::vector<int> K;
	K.resize(nopt);
	for(auto & a: K) a =0;

	int I,J,L;
	double x,y,z;
	std::string C, D;

	int NSYMM=0;

	std::string inputline;
	do{
		std::getline(in, inputline);
		if(in.eof())break;
		if(inputline.find('!')!=std::string::npos) continue;
		
		L=0;
		if(inputline=="GEOMETRY"){
			//std::cout << "READING GEOMETRY\n";
			NSYMM=0;
			std::string a,x,y,z;

			auto & ellist=libint2::chemistry::get_element_info();

			while(true){
				std::string gline;
				std::getline(in,gline);
				std::stringstream gs(gline);

				gs >> a >> x >> y >> z;

				if(a=="X") break;
				libint2::Atom at;
				
				for(auto el: ellist){
					if(a==el.symbol){
						at.atomic_number=el.Z;
						at.x=std::stod(x);
						at.y=std::stod(y);
						at.z=std::stod(z);
						break;
					}
				}

				atoms.push_back(at);
			}
		}
		else if(inputline=="SYMMETRY"){
			std::cout << "SYMMETRY NOT YET IMPLEMENTED\n";
			return -1;
		}
		else{
			std::string vstr;
			std::getline(in,vstr);
			
			int type = 0;

			if     (iparams.find(inputline)!=iparams.end()){
				iparams[inputline]=std::stoi(vstr);
				type = 1;
			}
			else if(dparams.find(inputline)!=dparams.end()){
				for(auto & c: vstr){
					if(c=='D') c='E';
				}
				dparams[inputline]=std::stod(vstr);
				type=2;
			}
			else if(cparams.find(inputline)!=cparams.end()){
				cparams[inputline]=vstr;
				type=3;
			}
			else if(lparams.find(inputline)!=lparams.end()){
				lparams[inputline]=vstr==".TRUE.";
				type=4;
			}

			//std::cout << type << ": " << inputline << "=" << vstr<<std::endl;
			if(type == 0){
				std::cout << "PARAM NOT FOUND: " << inputline << std::endl;
			}
		}
		
	}while(!in.eof());
	
	nelec=0;
	for(auto &a: atoms) nelec+=a.atomic_number;
	nelec-=iparams["CHARGE"];

	return 0;
}

