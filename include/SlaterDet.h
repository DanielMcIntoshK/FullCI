#ifndef SLATERDET__H_
#define SLATERDET__H_
#include <vector>
#include <array>
#include <string>

int rfactorial(int n,int m=1);
int choose(int n, int m);

struct StringMap{
	std::vector<unsigned char*> strs;
	int codeblklen;
	std::array<unsigned char *,4> cpystrs;
	int norbs;
	int nelec;

	std::string strlit(int i);
	std::string printstrab(int i);
	std::vector<std::vector<int>> strph(int i, bool alpha=true);
	std::string printstrphab(int i);
};

struct SlaterCompare{
	std::vector<int> diff, share;
};

class SlaterDet{
public:
	SlaterDet(int strn);
	SlaterDet(int alpidx, int betidx);

	void print();

	static void buildStrings(int norb,int nelec);
	static void buildDecoder();
	static int decode(unsigned char * str);
	static void printStrings();
	static void cleanStrings();

	void copyStrings(int px);
public:
	int aidx,bidx;

	static StringMap codes;
	static std::vector<std::vector<int>> decoder;
};

SlaterCompare compareSlaterDet(SlaterDet & sd1, SlaterDet & sd2);

double secondQuant1e(SlaterDet & sd1, SlaterDet & sd2,int,int,int,bool=false);

double secondQuant2e(SlaterDet & sd1, SlaterDet & sd2, int,int,int,int,int,bool=false);

bool checkOpperator(unsigned char *  str, int opblk, int opidx, bool anihilate);
int countLower(unsigned char *  str, int blk, int idx);

#endif

