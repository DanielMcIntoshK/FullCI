#ifndef BASISBUILDER__H_
#define BASISBUILDER__H_
#include <libint2/basis.h>
#include <iostream>

#include "ModelParams.h"

//extern class libint2::BasisSet;

class BasisBuilder{
public:
	BasisBuilder();

	int Build(std::string filename, ModelParams & params);

	libint2::BasisSet bs;
};

#endif 

