#ifndef HMM_HPP
#define HMM_HPP

#include <vector>
#include <eigen3/Eigen/Dense>  // eigen: linear algebra library
using namespace std;
using namespace Eigen;


class HMM{
private:
	int n;  // number of states
	int m;  // number of observation in alphabet
	MatrixXd a;  // state transition probability distribution nxn
	MatrixXd b;  // observation symbol probability distribution nxm
	MatrixXd pi; // intitial state distribution

public:
	HMM(int n_, int m_);
	~HMM();

	// Element access
	int N() const;
	int M() const;
	int &N();
	int &M();
	MatrixXd PI() const;
	MatrixXd &PI();
	MatrixXd A() const;
	MatrixXd &A();
	MatrixXd B() const;
	MatrixXd &B();


	// Method
	double Evaluate(const MatrixXi &observation);
	MatrixXd Forward(const MatrixXi &observation);
	MatrixXd Backward(const MatrixXi &observation);

	MatrixXi Decode(const MatrixXi &observation, double & probability);
	double Learn(vector<MatrixXi> &observation, int iterations, double tolerance);

};

#endif
