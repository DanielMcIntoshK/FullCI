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

SlaterDet::SlaterDet(int alpidx, int betidx){
	aidx=alpidx;
	bidx=betidx;
}

std::vector<int> SlaterDet::operator^(const SlaterDet & sd) const{
	std::vector<int> nm1,nm2;

	int nmcnt=0;
	for(int i = 0; i < codes.codeblklen;i++){
		std::bitset<8> c1(codes.strs[aidx][i]),c2(codes.strs[sd.aidx][i]);
		std::bitset<8> check = c1^c2;
		
		for(int j = 0; j < 8;j++){
			if(check.test(j)){
				if(c1.test(j))nm1.push_back(i*8+j);
				else nm2.push_back(i*8+j);
			}
		}
	}
	for(int i = 0; i < codes.codeblklen;i++){
		std::bitset<8> c1(codes.strs[bidx][i]),c2(codes.strs[sd.bidx][i]);
		std::bitset<8> check = c1^c2;
		
		for(int j = 0; j < 8;j++){
			if(check.test(j)){
				if(c1.test(j))nm1.push_back(i*8+j);
				else nm2.push_back(i*8+j);
			}
		}
	}
	nm1.insert(nm1.end(),nm2.begin(),nm2.end());
	return nm1;
}

std::vector<int> SlaterDet::operator&(const SlaterDet & sd) const{
	std::vector<int> share;

	for(int i = 0; i < codes.codeblklen;i++){
		std::bitset<8> c1(codes.strs[aidx][i]),c2(codes.strs[sd.aidx][i]);
		std::bitset<8> check = c1&c2;
		
		for(int j = 0; j < 8;j++){
			if(check.test(j))share.push_back(i*8+j);
		}
	}
	for(int i = 0; i < codes.codeblklen;i++){
		std::bitset<8> c1(codes.strs[bidx][i]),c2(codes.strs[sd.bidx][i]);
		std::bitset<8> check = c1&c2;
		
		for(int j = 0; j < 8;j++){
			if(check.test(j))share.push_back(i*8+j);
		}
	}
	return share;
}

void SlaterDet::print(){
	for(int j = 0; j < SlaterDet::codes.codeblklen;j++){
		std::bitset<8> pbs(SlaterDet::codes.strs[aidx][SlaterDet::codes.codeblklen-1-j]);
		std::cout << pbs << " ";
	}
	std::cout << std::endl;
	for(int j = 0; j < SlaterDet::codes.codeblklen;j++){
		std::bitset<8> pbs(SlaterDet::codes.strs[bidx][SlaterDet::codes.codeblklen-1-j]);
		std::cout << pbs << " ";
	}
	std::cout << std::endl;
}

void SlaterDet::buildStrings(int norb, int nelec){

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
	for(int i = 0; i < SlaterDet::codes.strs.size();i++){
		delete[] SlaterDet::codes.strs[i];
	}
	SlaterDet::codes.strs.clear();
}

StringMap SlaterDet::codes=StringMap();

SlaterCompare compareSlaterDet(SlaterDet & sd1, SlaterDet & sd2){
	SlaterCompare sc;
	std::vector<int> nshared1_a,nshared2_a,
		nshared1_b,nshared2_b;
	for(int i = 0; i < SlaterDet::codes.codeblklen; i++){
		std::bitset<8> c1a(SlaterDet::codes.strs[sd1.aidx][i]),c2a(SlaterDet::codes.strs[sd2.aidx][i]),
			c1b(SlaterDet::codes.strs[sd1.bidx][i]), c2b(SlaterDet::codes.strs[sd2.bidx][i]);
		
		std::bitset<8> checka = c1a^c2a;
		std::bitset<8> checkb = c1b^c2b;

		std::bitset<8> anda = c1a&c2a;
		std::bitset<8> andb = c1b&c2b;
		std::bitset<8> sharedo=anda&andb;
		std::bitset<8> shareso=anda^andb;
		
		for(int j = 0; j < 8;j++){
			if(checka.test(j)){
				if(c1a.test(j))nshared1_a.push_back(i*8+j);
				else nshared2_a.push_back(i*8+j);
			}
			if(checkb.test(j)){
				if(c1b.test(j))nshared1_b.push_back(i*8+j);
				else nshared2_b.push_back(i*8+j);
			}
			if(sharedo.test(j)) sc.share_do.push_back(i*8+j);
			if(shareso.test(j)) sc.share_so[(anda.test(j))?0:1].push_back(i*8+j);
		}
	}
	nshared1_a.insert(nshared1_a.end(),nshared2_a.begin(),nshared2_a.end());
	nshared1_b.insert(nshared1_b.end(),nshared2_b.begin(),nshared2_b.end());
	
	sc.diff[0]=nshared1_a;
	sc.diff[1]=nshared1_b;

	return sc;
}

double secondQuant1e(SlaterDet & sd1, SlaterDet & sd2, int p, int r,int n){
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
	sign+=countLower(cpystrs[ra],rblk,ridx)+((ra!=0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[ra][rblk]^=(1<<ridx);
	
	if(!checkOpperator(cpystrs[pa],pblk,pidx,false)) return 0;	
	sign+=countLower(cpystrs[pa],pblk,pidx)+((pa!=0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[pa][pblk]^=(1<<pidx);

	for(int i = 0; i < blklen; i++){
		if(memcmp(cpystrs[0],cpystrs[2],blklen)!=0||
			memcmp(cpystrs[1],cpystrs[3],blklen)!=0) return 0.0;
	}
	return (sign%2==0)?1.0:-1.0;
}

double secondQuant2e(SlaterDet & sd1, SlaterDet & sd2, int p, int q, int r, int s,int n){
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
	sign+=countLower(cpystrs[sa],sblk,sidx)+((sa!=0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[sa][sblk]^=(1<<sidx);

	if(!checkOpperator(cpystrs[ra],rblk,ridx,true)) return 0;	
	sign+=countLower(cpystrs[ra],rblk,ridx)+((ra!=0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[ra][rblk]^=(1<<ridx);

	if(!checkOpperator(cpystrs[qa],qblk,qidx,false)) return 0;
	sign+=countLower(cpystrs[qa],qblk,qidx)+((qa!=0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[qa][qblk]^=(1<<qidx);
	
	if(!checkOpperator(cpystrs[pa],pblk,pidx,false)) return 0;
	sign+=countLower(cpystrs[pa],pblk,pidx)+((pa!=0)?0:countLower(cpystrs[0],blklen,-1));
	cpystrs[pa][pblk]^=(1<<pidx);

	for(int i = 0; i < blklen; i++){
		if(memcmp(cpystrs[0],cpystrs[2],blklen)!=0||
			memcmp(cpystrs[1],cpystrs[3],blklen)!=0) return 0.0;
	}
	return (sign%2==0)?1.0:-1.0;
}

bool checkOpperator(unsigned char * str, int opblk, int opidx, bool anihilate){
	if(anihilate) return str[opblk]&(1<<opidx);
	return !(str[opblk]&(1<<opidx));
}
int countLower(unsigned char * str, int blk, int idx){
	int cnt=0;
	for(int i = 0; i < blk; i++){
		cnt+=std::bitset<8>(str[blk]).count();
	}
	if(idx>=0){
		for(int i = 0; i < idx; i++){
			if(str[blk]&(1<<i)) cnt++;
		}
	}
	return cnt;
}

