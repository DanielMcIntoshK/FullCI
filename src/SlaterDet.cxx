#include "SlaterDet.h"
#include <bitset>
#include <algorithm>
#include <iostream>

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

std::vector<int> SlaterDet::operator|(const SlaterDet & sd) const{
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

void SlaterDet::buildStrings(int norb, int nelec){

	for(int i = 0; i < 10; i++){
		std::cout << "CHOOSETEST: " << rfactorial(i)<< std::endl;
	}
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

