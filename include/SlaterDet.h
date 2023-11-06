#ifndef SLATERDET__H_
#define SLATERDET__H_
#include <vector>

struct StringMap{
	std::vector<unsigned char*> strs;
	int codeblklen;
};

class SlaterDet{
public:
	SlaterDet(int alpidx, int betidx);

	std::vector<int> operator|(const SlaterDet & sd) const;

	static void buildStrings(int norb,int nelec, int strcnt);
	
	static void printStrings();
public:
	int aidx,bidx;

	static StringMap codes;
};

#endif

