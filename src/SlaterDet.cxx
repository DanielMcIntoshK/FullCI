#include "SlaterDet.h"
#include <bitset>
#include <algorithm>
#include <iostream>

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

void SlaterDet::buildStrings(int norb, int nelec,int strcnt){
	std::vector<bool> permute;
	permute.resize(norb);
	for(int i = 0; i < norb;i++){permute[i]=i<nelec;}

	SlaterDet::codes.codeblklen=norb/8 +1;
	SlaterDet::codes.strs.resize(strcnt);

	int count = 0;
	do{
		SlaterDet::codes.strs[count]=new unsigned char[SlaterDet::codes.codeblklen];
		for(int i = 0; i < permute.size();i++){
			if(permute[i]) SlaterDet::codes.strs[count][i/8]+=0b1<<(i%8);
		}
		count++;
	}while(std::prev_permutation(permute.begin(),permute.end()));
}

void SlaterDet::printStrings(){
	for(int i = 0; i < SlaterDet::codes.strs.size();i++){
		std::cout << i << ": ";
		for(int j = 0; j < SlaterDet::codes.codeblklen;j++){
			std::bitset<8> pbs(SlaterDet::codes.strs[i][SlaterDet::codes.codeblklen-1-j]);
			std::cout << pbs << " ";
		}
		std::cout << std::endl;
	}
}

StringMap SlaterDet::codes=StringMap();

