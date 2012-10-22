/*
 * Hidden Markov Model
 */

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
	MatrixXd pi; // initial state distribution nx1

public:
	HMM(int n_, int m_);
	~HMM();
	// Methods
	double Evaluate(const MatrixXi &observation);
	MatrixXd Forward(const MatrixXi &observation, MatrixXd &scale);
	MatrixXd Backward(const MatrixXi &observation, const MatrixXd &scale);
	MatrixXi Decode(const MatrixXi &observation, double &probability);
	double Learn(vector<MatrixXi> &observations, int iterations, double tolerance);
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

template <class T>
class MyMatrix3x{
private:
	int width, height;
	vector<T> data;
public:
	MyMatrix3x():width(0), height(0){
	}
	MyMatrix3x(int m, int n, int p): width(m), height(n){
		resize(m,n,p);
	}
	void resize(int m, int n, int p){
		data.resize(m*n*p);
	}
	T operator()(int i, int j, int k) const{
		return data[i + j*width + k*width*height];
	}
	T& operator()(int i, int j, int k){
		return data[i + j*width + k*width*height];
	}
};

typedef MyMatrix3x<double> MyMatrix3d;
}
#endif
