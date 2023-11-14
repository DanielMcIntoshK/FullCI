#ifndef SLATERDET__H_
#define SLATERDET__H_
#include <vector>
#include <array>

int rfactorial(int n,int m=1);
int choose(int n, int m);

struct StringMap{
	std::vector<unsigned char*> strs;
	int codeblklen;
};

struct SlaterCompare{
	std::array<std::vector<int>,2> diff, 
		share_so;
	std::vector<int> share_do;
};

class SlaterDet{
public:
	SlaterDet(int alpidx, int betidx);

	std::vector<int> operator^(const SlaterDet & sd) const;
	std::vector<int> operator&(const SlaterDet & sd) const;

	static void buildStrings(int norb,int nelec);
	static void printStrings();
	static void cleanStrings();
public:
	int aidx,bidx;

	static StringMap codes;
};

SlaterCompare compareSlaterDet(SlaterDet & sd1, SlaterDet & sd2);

#endif

