#ifndef MODELPARAMS__H_
#define MODELPARAMS__H_
#include <string>
#include <map>
#include <iostream>
#include <libint2/atom.h>

class ModelParams{
public:
	ModelParams();

	std::map<std::string, int>         iparams;
	std::map<std::string, double>      dparams;
	std::map<std::string, std::string> cparams;
	std::map<std::string, bool>        lparams;

	std::vector<libint2::Atom> atoms;

	int ReadInputFile(std::istream & in);
};


#endif 

