#include "Diagram.h"
#include "Integrals.h"
#include <set>

//Diagram::Diagram(int ne, int no):nelec{ne},norbs{no}{
//
//}

Diagram::Diagram(std::vector<dnode> outline, std::vector<resolvent> resolve,hfresults res):
	hfr{res}{
	buildGraph(outline);
	rslines.resize(diag.size()+1);
	for(int i = 0; i < resolve.size(); i++){
		rslines[resolve[i].pos].push_back(resolve[i]);
	}
}

void Diagram::buildGraph(std::vector<dnode> outline){
	std::vector<dnode> tnodes;
	tnodes.resize(outline.size());
	int runningloop=0;
	for(int i = 0; i < outline.size(); i++){
		for(int j =0; j < outline[i].size(); j++){
			line lref=outline[i][j];
			
			//std::cout << tnodes.size() << " " << i << " " << lref.dir << std::endl;	
			tnodes[i].push_back(
				line(lref.dir,lref.label,lref.dir<i,runningloop));
			tnodes[lref.dir].push_back(
				line(i,lref.label,lref.dir<i,runningloop));
			if(lref.loop) {
				if(lref.dir<i) hlooplines.push_back(runningloop);
				else plooplines.push_back(runningloop);
			}
			runningloop+=1;
		}
	}
	for(int i = 0; i < tnodes.size(); i++){
		diag.push_back(tnodes[i]);
	}
}

double Diagram::sumfull(IntegralChugger * ic, double E){
	
	/*
	for(std::list<dnode>::iterator node = diag.begin(); node!=diag.end(); node++){
		for(int i = 0; i < node->size(); i++){
			line lref=(*node)[i];
			std::cout << lref.dir << " " << lref.label << " ";
			std::cout << lref.ref << " ";
			if(lref.hole) std::cout << "HOLE | ";
			if(!lref.hole) std::cout << "PART | ";
		}
		std::cout << std::endl;
	}
	*/
	std::map<int,int> hlookup, plookup;
	for(int i = 0; i < hlooplines.size(); i++){
		hlookup[hlooplines[i]]=0;
	}	
	for(int i = 0; i < plooplines.size(); i++){
		plookup[plooplines[i]]=0;
	}

	double sum = 0;
	do{
		do{	
			sum+=sumdiag(hlookup, plookup, ic, E);
		}while(permutemap(hlookup,true));
	}while(permutemap(plookup,false));
	return getsign()*sum/getmultiplicity();
}

double Diagram::sumdiag(std::map<int,int> holelookup,std::map<int,int> partlookup, 
		IntegralChugger * ic, double E){
	/*
	std::cout << "GRAPH SETTINGS:\n";
	int posp = 0;
	for(auto itr=diag.begin(); itr!=diag.end(); itr++){
		for(int i = 0; i < itr->size(); i++){
			std::cout << posp;
			if((*itr)[i].hole) std::cout << "<-";
			else std::cout << "->";
			std::cout << (*itr)[i].dir << " | ";
		}
		std::cout << std::endl;
		posp+=1;
	}
	*/
	double numer=1,denom=1;

	int pos=0;
	bool alowed=true;
	std::list<int> chlist, cplist;
	for(auto node=diag.begin(); node!=diag.end(); node++){
		std::vector<int> outs, ins;
		for(int j = 0; j < node->size(); j++){
			line lnref=(*node)[j];
			int label=lnref.label;
			if(lnref.loop){
				if(lnref.hole){
					label=hfr.phlist[0][holelookup[lnref.ref]];
				}
				else{
					label=hfr.phlist[1][partlookup[lnref.ref]];
				}
			}
			if(lnref.dir<pos){
				if(lnref.hole){
					for(auto itr=chlist.begin(); itr!=chlist.end(); itr++){
						if((*itr)==label){
							chlist.erase(itr);
							break;
						}
					}
					outs.push_back(label);
				}
				else{
					for(auto itr=cplist.begin(); itr!=cplist.end(); itr++){
						if((*itr)==label){
							cplist.erase(itr);
							break;
						}
					}
					ins.push_back(label);
				}
			}
			else if(lnref.dir>pos){
				if(lnref.loop){
					for(int rs = 0; rs < rslines[pos].size(); rs++){
						if(rslines[pos][rs].nloops.find(label)!=
								rslines[pos][rs].nloops.end()){
							return 0.0;
						}
					}
				}
				if(lnref.hole){
					chlist.push_back(label);
					ins.push_back(label);
				}
				else{
					cplist.push_back(label);
					outs.push_back(label);
				}
			}
		}
		for(int rs = 0; rs < rslines[pos].size();rs++){
			resolvent & rsps=rslines[pos][rs];
			
			std::list<int> chlist_u=chlist, cplist_u=cplist;
			for(int j = 0; j < rsps.ignore.size(); j++){
				for(auto itr=chlist_u.begin(); itr!=chlist_u.end();){
					if((*itr)==rsps.ignore[j]){
						itr=chlist_u.erase(itr);
						break;
					}
					else{
						itr++;
					}
				}
				for(auto itr=cplist_u.begin(); itr!=cplist_u.end();){
					if((*itr)==rsps.ignore[j]){
						itr=cplist_u.erase(itr);
						break;
					}
					else{
						itr++;
					}
				}
			}

			double val=0.0;
			if(rsps.E){
				val+=(rsps.neg)?-E:E;
			}
			for(auto cp=cplist_u.begin(); cp!=cplist_u.end(); cp++){
				val+=((rsps.E)?1.0:-1.0)*
					hfr.E((*cp)%hfr.norbs,0);
			}
			for(auto ch=chlist_u.begin(); ch!=chlist_u.end(); ch++){
				val+=((rsps.E)?-1.0:1.0)*
					hfr.E((*ch)%hfr.norbs,0);
			}
			if(rsps.invert) val=1.0/val;
			denom*=val;
		}
		if(node->size()==4) {
			numer*=ic->movsymphys(outs[0],outs[1],ins[0],ins[1],hfr.norbs);
		}
		pos+=1;
	}
	return numer/denom;
}

bool Diagram::permutemap(std::map<int,int> & pm, bool hole){
	std::map<int,int>::iterator cpm=pm.begin();

	int loopsize=(hole)?hfr.phlist[0].size():hfr.phlist[1].size();
	while(cpm!=pm.end()){
		cpm->second+=1;
		if(cpm->second>=loopsize){
			cpm->second=0;
			cpm++;
		}
		else{
			break;
		}
	}
	return cpm!=pm.end();

}


bool Diagram::disconnected(){
	return false;
}

bool Diagram::reducable(){
	return false;
}

int Diagram::getloops(){
	std::vector<dnode> outlines, inlines;
	outlines.resize(diag.size());
	inlines.resize(diag.size());
	
	int pos = 0;
	for(auto itr = diag.begin(); itr!=diag.end(); itr++){
		for(int i = 0; i < itr->size(); i++){
			if((*itr)[i].dir > pos){
				if((*itr)[i].hole) inlines[pos].push_back((*itr)[i]);
				else outlines[pos].push_back((*itr)[i]);
			}
			else{
				if(!(*itr)[i].hole) inlines[pos].push_back((*itr)[i]);
				else outlines[pos].push_back((*itr)[i]);
			}
		}
		pos+=1;
	}

	std::set<int> crossedlines;
	//std::list<std::vector<line>> paths;
	int paths=0;
	for(int i = 0; i < outlines.size(); i++){
		for(int j = 0; j < outlines[i].size(); j++){
			if(crossedlines.find(outlines[i][j].ref)!=crossedlines.end()) continue;

			int pos=i;
			int ind=j;
			paths+=1;
			
			line cline=outlines[i][j];
			while(crossedlines.find(cline.ref)==crossedlines.end()){
				crossedlines.insert(cline.ref);
				pos=cline.dir;
				for(int k = 0; k < inlines[pos].size(); k++){
					if(inlines[pos][k].ref==cline.ref){
						ind=k; break;
					}
				}
				cline=outlines[pos][ind];
			}
		}	
	}
	//std::cout << "PATHS: " <<paths << std::endl;
	return paths;
}

int Diagram::getholes(){
	int hcount=0;
	for(auto itr=diag.begin(); itr!=diag.end(); itr++){
		for(int i = 0; i < itr->size(); i++){
			if((*itr)[i].hole) hcount++;
		}
	}
	return hcount/2;
}

double Diagram::getsign(){
	return (((getloops()+getholes())%2)==0)?1.0:-1.0;
}

double Diagram::getmultiplicity(){
	double multiplicity=1.0;
	int pos = 0;
	for(auto itr=diag.begin(); itr!=diag.end(); itr++){
		std::vector<line> looplines;
		for(int i = 0; i < itr->size(); i++){
			if((*itr)[i].dir >pos){
				if((*itr)[i].loop) looplines.push_back((*itr)[i]);
			}
		}
		for(int i = 0; i < looplines.size(); i++){
			int matches=1;
			for(int j = i+1; j < looplines.size(); j++){
				if((looplines[i].hole==looplines[j].hole) &&
					(looplines[i].dir==looplines[j].dir)){
						matches++;
				}
			}
			double factor=1.0;
			for(int mt=matches; mt>1; mt--) factor*=(double)mt;
			multiplicity*=factor;
		}
		pos+=1;
	}
	return multiplicity;
}


