#ifndef SLATERDET__H_
#define SLATERDET__H_
#include <vector>

class SlaterDet{
public:
	SlaterDet(int norb, int nelec, int excite=0);

	bool Permute();

	void print();

	int operator[](int i){return occupied[i];}
public:
	std::vector<int> occupied;
	int excitationlevel;

	int norbitals;
	int nelectrons;

	int homo, lumo;
	int max;

	int top;
};

#endif

