#include "SlaterDet.h"
#include <bitset>
#include <algorithm>
#include <iostream>
#include <cstring>

int rfactorial(int n, int m){
	int v=1;
	for(int i = n; i >m;i--) v*=i;
	return v;
}

int choose(int n, int m){
	return rfactorial(n,n-m)/rfactorial(m);
}

std::string StringMap::strlit(int i){
	std::string strtext="";
	
	unsigned char * str=strs[i];
	for(int block = 0; block < codeblklen; block++){
		unsigned char c = str[block];
		for(int bit = 0; bit < 8; bit++){
			if(bit+block*8>=norbs) break;

			strtext.append((c&1<<bit)?"1":"0");
		}
	}
	return strtext;
}

std::vector<std::vector<int>> StringMap::strph(int i, bool alpha){
	std::vector<int> holes, parts;

	if(alpha) i=i%strs.size();
	else i=i/strs.size();

	unsigned char * str=strs[i];
	for(int block = 0; block < codeblklen; block++){
		unsigned char c = str[block];
		for(int bit = 0; bit < 8; bit++){
			int k = bit+block*8;
			if(k>=norbs)break;

			int kp=k;
			if(!alpha) kp+=norbs;
			
			bool biton=c&1<<bit;
			if(!biton&&k<nelec) holes.push_back(kp);
			else if(biton&&k>=nelec)parts.push_back(kp);
		}
	}
	return std::vector<std::vector<int>>{holes,parts};
}

std::string StringMap::printstrab(int i){
	return strlit(i%strs.size()) + " " + strlit(i/strs.size());
}

std::string StringMap::printstrphab(int i){
	std::vector<std::vector<int>> aph=strph(i, true), bph=strph(i,false);

	std::string strtext="0_";
	for(int i = 0; i < aph[0].size(); i++){
		if(i!=0) strtext+=",";
		strtext+=std::to_string(aph[0][i]);
	}
	for(int i = 0; i < bph[0].size(); i++){
		if(!aph[0].empty()) strtext+=",";
		strtext+=std::to_string(bph[0][i]);
	}
	strtext+="^";
	for(int i = 0; i < aph[1].size(); i++){
		if(i!=0) strtext+=",";
		strtext+=std::to_string(aph[1][i]);
	}
	for(int i = 0; i < bph[1].size(); i++){
		if(!aph[1].empty()) strtext+=",";
		strtext+=std::to_string(bph[1][i]);
	}

	return strtext;
}

SlaterDet::SlaterDet(int strn){
	aidx=strn%SlaterDet::codes.strs.size();
	bidx=strn/SlaterDet::codes.strs.size();
}
SlaterDet::SlaterDet(int alpidx, int betidx){
	aidx=alpidx;
	bidx=betidx;
}

void SlaterDet::print(){
	for(int j = 0; j < SlaterDet::codes.codeblklen;j++){
		//std::bitset<8> pbs(SlaterDet::codes.strs[aidx][SlaterDet::codes.codeblklen-1-j]);
		//std::cout << pbs << " ";
		for(int b=0; b < 8; b++){
			if(SlaterDet::codes.strs[aidx][j]&(1<<b)){
				std::cout << 1;
			}
			else std::cout << 0;
		}
	}
	std::cout << std::endl;
	for(int j = 0; j < SlaterDet::codes.codeblklen;j++){
		//std::bitset<8> pbs(SlaterDet::codes.strs[bidx][SlaterDet::codes.codeblklen-1-j]);
		//std::cout << pbs << " ";
		for(int b=0; b < 8; b++){
			if(SlaterDet::codes.strs[bidx][j]&(1<<b)){
				std::cout << 1;
			}
			else std::cout << 0;
		}
	}
	std::cout << std::endl;
}

void SlaterDet::buildStrings(int norb, int nelec){
	SlaterDet::codes.norbs=norb;
	SlaterDet::codes.nelec=nelec;

	int strcnt=choose(norb,nelec);
	std::vector<bool> permute;
	permute.resize(norb);
	for(int i = 0; i < norb;i++){permute[i]=i<nelec;}

	SlaterDet::codes.codeblklen=norb/8 +1;
	SlaterDet::codes.strs.resize(strcnt);

	std::bitset<8> feedbit;
	int count = 0;
	do{
		SlaterDet::codes.strs[count]=new unsigned char[SlaterDet::codes.codeblklen];
		for(int i = 0; i < SlaterDet::codes.codeblklen;i++){
			feedbit.reset();
			for(int j = 0; j < 8; j++){
				int index = i*8+j;
				feedbit.set(j,(index>=permute.size())?0:permute[index]);
			}
			SlaterDet::codes.strs[count][i]=feedbit.to_ulong();
		}
		count++;
	}while(std::prev_permutation(permute.begin(),permute.end()));
		
	for(int i = 0; i < 4; i++) SlaterDet::codes.cpystrs[i]=new unsigned char[SlaterDet::codes.codeblklen];

	SlaterDet::buildDecoder();
}

void SlaterDet::buildDecoder(){
	SlaterDet::decoder.resize(SlaterDet::codes.norbs-SlaterDet::codes.nelec);
	for(auto & d: SlaterDet::decoder) d.resize(SlaterDet::codes.nelec+1);

	std::vector<std::vector<int>> xvals;
	xvals.resize(SlaterDet::decoder.size()+1);
	for(auto & xv:xvals) xv.resize(SlaterDet::decoder[0].size());

	for(int o = xvals.size()-1;o>=0;o--){
		for(int e=xvals[0].size()-1;e>=0;e--){
			if((e+1)>=xvals[o].size()||(o+1)>=xvals.size()){
				xvals[o][e]=1;
				continue;
			}
			int right = xvals[o][e+1];
			int down = xvals[o+1][e];
			xvals[o][e]=right+down;
		}
	}
	for(int o = 0; o < SlaterDet::decoder.size();o++){
		for(int e = 0; e < SlaterDet::decoder[o].size();e++){
			SlaterDet::decoder[o][e]=((e+1)>=xvals[o].size())?0:xvals[o][e+1];
		}
	}
}

int SlaterDet::decode(unsigned char* str){
	int index = 0; 
	
	int ce=0, co=0;
	for(int b = 0; b < SlaterDet::codes.norbs;b++){
		int codeblock=b/8;
		unsigned char scan = 1<<(b%8);
		if(!(str[codeblock]&scan)){
			index+=SlaterDet::decoder[co][ce];
			co++;
		}
		else{
			ce++;
		}
	}

	return index;
}

void SlaterDet::printStrings(){
	for(int i = 0; i < SlaterDet::codes.strs.size();i++){
		std::cout << i << ": ";
		for(int j = 0; j < SlaterDet::codes.codeblklen;j++){
			std::cout << i << " " << j << " " << SlaterDet::codes.codeblklen-1-j << std::endl; 
			std::bitset<8> pbs(SlaterDet::codes.strs[i][SlaterDet::codes.codeblklen-1-j]);
			std::cout << pbs << " ";
		}
		std::cout << std::endl;
	}
}

void SlaterDet::cleanStrings(){
	std::cout << "CLEANING STRINGS\n";
	for(int i = 0; i < SlaterDet::codes.strs.size();i++){
		delete[] SlaterDet::codes.strs[i];
	}
	SlaterDet::codes.strs.clear();
	for(int i = 0; i < 4; i ++){
		delete[] SlaterDet::codes.cpystrs[i];
	}
	//SlaterDet::codes.cpystrs.clear();
}

StringMap SlaterDet::codes=StringMap();
std::vector<std::vector<int>> SlaterDet::decoder=std::vector<std::vector<int>>();

SlaterCompare compareSlaterDet(SlaterDet & sd1, SlaterDet & sd2){
	auto & strs=SlaterDet::codes.strs;
	int & blklen=SlaterDet::codes.codeblklen;
	int & n=SlaterDet::codes.norbs;
	
	unsigned char * alpha1=strs[sd1.aidx], *alpha2=strs[sd2.aidx],
		*beta1=strs[sd1.bidx], *beta2=strs[sd2.bidx];

	SlaterCompare sc;
	std::vector<int> diff1, diff2;
	for(int i = 0; i < blklen; i++){
		std::bitset<8> a1blk(alpha1[i]),a2blk(alpha2[i]),
			b1blk(beta1[i]),b2blk(beta2[i]);

		std::bitset<8> diffalpha=a1blk^a2blk,diffbeta=b1blk^b2blk;
		std::bitset<8> samealpha=a1blk&a2blk,samebeta=b1blk&b2blk;

		for(int j = 0; j < 8; j++){
			if(diffalpha.test(j)) {((a1blk.test(j))?diff1:diff2).push_back(j+8*i);}
			if(diffbeta.test(j)) {((b1blk.test(j))?diff1:diff2).push_back(j+8*i+n);}
			if(samealpha.test(j)){sc.share.push_back(j+i*8);}
			if(samebeta.test(j)){sc.share.push_back(j+i*8+n);}
		}
	}
	diff1.insert(diff1.end(),diff2.begin(),diff2.end());
	sc.diff=diff1;
	return sc;
}

double secondQuant1e(SlaterDet & sd1, SlaterDet & sd2, int p, int r,int n,bool verbose){
	auto & strs=SlaterDet::codes.strs;
	auto & cpystrs=SlaterDet::codes.cpystrs;
	int & blklen=SlaterDet::codes.codeblklen;

	memcpy(cpystrs[0],strs[sd2.aidx],blklen);
	memcpy(cpystrs[1],strs[sd2.bidx],blklen);
	memcpy(cpystrs[2],strs[sd1.aidx],blklen);
	memcpy(cpystrs[3],strs[sd1.bidx],blklen);

	int pa=(p<n)?0:1, ra=(r<n)?0:1;
	int pblk=(p-n*pa)/8,pidx=(p-n*pa)%8,
	    rblk=(r-n*ra)/8,ridx=(r-n*ra)%8;
	
	int sign=0;

	if(!checkOpperator(cpystrs[ra],rblk,ridx,true)) return 0;
	sign+=countLower(cpystrs[ra],rblk,ridx);
	if(ra==1)sign+=countLower(cpystrs[0],blklen,-1);
	cpystrs[ra][rblk]^=(1<<ridx);
	
	if(!checkOpperator(cpystrs[pa],pblk,pidx,false)) return 0;
	sign+=countLower(cpystrs[pa],pblk,pidx);
	if(pa==1)sign+=countLower(cpystrs[0],blklen,-1);
	cpystrs[pa][pblk]^=(1<<pidx);

	for(int i = 0; i < blklen; i++){
		if(memcmp(cpystrs[0],cpystrs[2],blklen)!=0||
			memcmp(cpystrs[1],cpystrs[3],blklen)!=0) return 0.0;
	}
	return (sign%2==0)?1.0:-1.0;
}

double secondQuant2e(SlaterDet & sd1, SlaterDet & sd2, int p, int q, int r, int s,int n,bool verbose){
	auto & strs=SlaterDet::codes.strs;
	auto & cpystrs=SlaterDet::codes.cpystrs;
	int & blklen=SlaterDet::codes.codeblklen;

	memcpy(cpystrs[0],strs[sd2.aidx],blklen);
	memcpy(cpystrs[1],strs[sd2.bidx],blklen);
	memcpy(cpystrs[2],strs[sd1.aidx],blklen);
	memcpy(cpystrs[3],strs[sd1.bidx],blklen);
	
	int pa=(p<n)?0:1, ra=(r<n)?0:1, qa=(q<n)?0:1, sa=(s<n)?0:1;
	int pblk=(p-n*pa)/8,pidx=(p-n*pa)%8,
	    qblk=(q-n*qa)/8,qidx=(q-n*qa)%8,
	    rblk=(r-n*ra)/8,ridx=(r-n*ra)%8,
	    sblk=(s-n*sa)/8,sidx=(s-n*sa)%8;
	
	int sign = 0;

	if(!checkOpperator(cpystrs[sa],sblk,sidx,true)) return 0;
	sign+=countLower(cpystrs[sa],sblk,sidx);
	sign+=((sa==0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[sa][sblk]^=(1<<sidx);

	if(!checkOpperator(cpystrs[ra],rblk,ridx,true)) return 0;
	sign+=countLower(cpystrs[ra],rblk,ridx);
	sign+=((ra==0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[ra][rblk]^=(1<<ridx);
	
	/*
	if(!checkOpperator(cpystrs[qa],qblk,qidx,false)) return 0;
	sign+=countLower(cpystrs[qa],qblk,qidx);
	sign+=((qa==0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[qa][qblk]^=(1<<qidx);

	if(!checkOpperator(cpystrs[pa],pblk,pidx,false)) return 0;
	sign+=countLower(cpystrs[pa],pblk,pidx);
	sign+=((pa==0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[pa][pblk]^=(1<<pidx);
	*/

	if(!checkOpperator(cpystrs[pa+2],pblk,pidx,true)) return 0;
	sign+=countLower(cpystrs[pa+2],pblk,pidx);
	sign+=((pa==0)?0:countLower(cpystrs[2],blklen,-1));
	cpystrs[pa+2][pblk]^=(1<<pidx);

	if(!checkOpperator(cpystrs[qa+2],qblk,qidx,true)) return 0;
	sign+=countLower(cpystrs[qa+2],qblk,qidx);
	sign+=((qa==0)?0:countLower(cpystrs[2],blklen,-1));
	cpystrs[qa+2][qblk]^=(1<<qidx);

	for(int i = 0; i < blklen; i++){
		if(memcmp(cpystrs[0],cpystrs[2],blklen)!=0||
			memcmp(cpystrs[1],cpystrs[3],blklen)!=0) return 0.0;
	}
	return (sign%2==0)?1.0:-1.0;
}

bool checkOpperator(unsigned char * str, int opblk, int opidx, bool anihilate){
	if(anihilate) return (str[opblk]&(1<<opidx));
	return !(str[opblk]&(1<<opidx));
}
int countLower(unsigned char * str, int blk, int idx){
	int cnt=0;
	for(int i = 0; i < blk; i++){
		cnt+=std::bitset<8>(str[i]).count();
	}
	if(idx>=0){
		for(int i = 0; i < idx; i++){
			if(str[blk]&(1<<i)) cnt++;
		}
	}
	return cnt;
}

double opOnSlater(PHOp op, unsigned char * alpha, unsigned char * beta){
	StringMap & sm = SlaterDet::codes;

	double sign = 1.0;

	int isalpha=(op.i/sm.norbs);

	int ibase=op.i%sm.norbs,
	    jbase=op.j%sm.norbs;

	unsigned char * sd = (isalpha==0)?alpha:beta;

	int iblock=ibase/8,ioff=ibase%8,
	    jblock=jbase/8,joff=jbase%8;

	if(sd[iblock]&(1<<ioff)){
		int lowercount =(isalpha==0)?0:SlaterDet::codes.nelec;
		for(int b = 0; b < ibase; b++){
			if(sd[b/8]&(1<<(b%8))){
				lowercount++;
			}
		}
		sign*=(lowercount%2==0)?1.0:-1.0;
		sd[iblock]=sd[iblock]^(1<<ioff);
	}
	else return 0.0;

	if(!(sd[jblock]&(1<<joff))){
		int lowercount=(isalpha==0)?0:SlaterDet::codes.nelec;
		for(int b = 0; b < jbase; b++){
			if(sd[b/8]&(1<<(b%8))){
				lowercount++;
			}
		}
		sign*=(lowercount%2==0)?1.0:-1.0;
		sd[jblock]=sd[jblock]|(1<<joff);
	}
	else return 0.0;

	return sign;
}

