#ifndef INTEGRALS__H_
#define INTEGRALS__H_

#include <vector>
#include <array>

#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <libint2.hpp>

#include "ModelParams.h"

using namespace libint2;

//extern class libint2::Operator;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;

typedef std::vector<std::vector<std::vector<std::vector<double>>>> twobodylist;

typedef Eigen::SelfAdjointEigenSolver<Matrix>::RealVectorType MatVec;

//This is a general class used for computing the integrals used in the program
//Originally this was mostly in HartreeFockSolver, but was moved to it's own class
class IntegralChugger{
public:
	enum IntegralType{
		NONE=0b0,
		OVERLAP=0b1,
		KINETIC=0b01,
		NUCLEAR=0b001,
		TWOBODY=0b0001,
		ALL=0b1111
	};
	enum Mode{
		BASIS=0b1,
		MOLECULAR=0b01
	};
public:
	IntegralChugger(BasisSet bs, ModelParams mp);

	//Will compute the integrals defined by IntegralType
	//I don't know why you would do anything but ALL
	void compute(IntegralType it);

	//Nicer way of grabbing two body integrals
	//Makes everything a little prettier
	double operator() (int a, int b, int c, int d){return tbv(a,b,c,d);}
	double tbv(int a, int b, int c, int d);
	double mov(int a, int b, int c, int d);

	void TransformInts(Matrix & C);

	bool isTransformed(){return !moi.empty();}
public:
	//Overlap, kinetic, nuclear matricies
	Matrix S;
	Matrix T;
	Matrix V;

	Matrix MOT;

	//The 4D list of 2body integrals
	twobodylist tbi;
	twobodylist moi;
private:
	BasisSet basis;
	ModelParams params;

	Matrix compute1eints(libint2::Operator obtype);
	void compute2eints();

	std::array<int,4> get2bodyintcord(int a, int b, int c, int d);

	void partialtransform(Matrix & C, int type);
};


#endif

