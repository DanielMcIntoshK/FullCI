#ifndef SLATERDET__H_
#define SLATERDET__H_
#include <vector>
#include <array>

int rfactorial(int n,int m=1);
int choose(int n, int m);

struct StringMap{
	std::vector<unsigned char*> strs;
	int codeblklen;
	std::array<unsigned char *,4> cpystrs;
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

	void print();

	static void buildStrings(int norb,int nelec);
	static void printStrings();
	static void cleanStrings();

	void copyStrings(int px);
public:
	int aidx,bidx;

	static StringMap codes;
};

SlaterCompare compareSlaterDet(SlaterDet & sd1, SlaterDet & sd2);

double secondQuant1e(SlaterDet & sd1, SlaterDet & sd2,int,int,int);

double secondQuant2e(SlaterDet & sd1, SlaterDet & sd2, int,int,int,int,int);

bool checkOpperator(unsigned char *  str, int opblk, int opidx, bool anihilate);
int countLower(unsigned char *  str, int blk, int idx);

#endif

