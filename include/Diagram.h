#ifndef DIAGRAM__H_
#define DIAGRAM__H_

#include <vector>
#include <list>
#include <map>
#include <set>
#include "HartreeFockSolver.h"

class IntegralChugger;

struct resolvent{
	resolvent(){}
	resolvent(int i_pos, bool b_E=false, bool b_neg=false,bool b_inv=false,
			std::vector<int> ig=std::vector<int>(),
			std::set<int> nl=std::set<int>()):
		pos{i_pos},E{b_E},neg{b_neg},ignore{ig},nloops{nl},invert{b_inv}{}
	int pos=0;
	bool E=true, neg=false;
	std::vector<int> ignore;
	std::set<int> nloops;

	bool invert=false;
};

struct line{
	line():dir{-1},loop{true},hole{false},label{-1}{}
	line(int i_end, int i_label=-1, bool b_hole=false,int i_lp=-1):
		dir{i_end}, loop{i_label==-1}, 
		label{i_label},hole{b_hole},ref{i_lp}{}
	int dir;
	bool loop;
	bool hole;
	int label;
	int ref;
};

typedef std::vector<line> dnode;
typedef std::list<dnode>::iterator dnodeptr;

class Diagram{
public:
	//Diagram(int ne, int no);
	Diagram(std::vector<dnode> outline,std::vector<resolvent> resolve, hfresults res); 
	
	void buildGraph(std::vector<dnode> outline);

	double sumfull(IntegralChugger * ic, double E);
	double sumdiag(std::map<int,int> holelookup,std::map<int,int> partlookup,
			IntegralChugger * ic,double E);

	bool permutemap(std::map<int,int> & pm, bool hole);

	bool disconnected();
	bool reducable();

	int getloops();
	int getholes();
	double getsign();
	double getmultiplicity();

	
	std::list<dnode> diag;
	std::vector<std::vector<resolvent>> rslines;
	std::vector<int> hlooplines,plooplines;

	hfresults hfr;
};

#endif DIAGRAM__H_

