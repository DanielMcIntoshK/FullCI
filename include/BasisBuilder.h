#ifndef BASISBUILDER__H_
#define BASISBUILDER__H_
#include <libint2.hpp>
#include <iostream>

#include "ModelParams.h"

class BasisBuilder{
public:
	BasisBuilder();

	int Build(std::string filename, ModelParams & params);

	libint2::BasisSet bs;
};

#endif 

