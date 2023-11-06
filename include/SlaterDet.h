#ifndef SLATERDET__H_
#define SLATERDET__H_
#include <vector>

int rfactorial(int n,int m=1);
int choose(int n, int m);

struct StringMap{
	std::vector<unsigned char*> strs;
	int codeblklen;
};

class SlaterDet{
public:
	SlaterDet(int alpidx, int betidx);

	std::vector<int> operator|(const SlaterDet & sd) const;

	static void buildStrings(int norb,int nelec);
	static void printStrings();
	static void cleanStrings();
public:
	int aidx,bidx;

	static StringMap codes;
};

#endif

