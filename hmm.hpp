#ifndef HMM_HPP
#define HMM_HPP

#include <vector>
#include <eigen3/Eigen/Dense>  // eigen: linear algebra library

namespace ltk{

using namespace std;
using namespace Eigen;


class HMM{
private:
	int n;  // number of states
	int m;  // number of observation in alphabet
	MatrixXd a;  // state transition probability distribution nxn
	MatrixXd b;  // observation symbol probability distribution nxm
	MatrixXd pi; // initial state distribution

public:
	HMM(int n_, int m_);
	~HMM();

	// Methods
	double Evaluate(const MatrixXi &observation);
	MatrixXd Forward(const MatrixXi &observation, MatrixXd &scale);
	MatrixXd Backward(const MatrixXi &observation, MatrixXd &scale);
	MatrixXi Decode(const MatrixXi &observation, double & probability);
	double Learn(vector<MatrixXi> &observation, int iterations, double tolerance);

	void Reset();

	// private members access
	int N() const;
	int M() const;
	MatrixXd A() const;
	MatrixXd B() const;
	MatrixXd PI() const;
	int &M();
	int &N();
	MatrixXd &A();
	MatrixXd &B();
	MatrixXd &PI();
};

}

#endif
