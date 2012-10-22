/*
 * hmm.cpp
 *
 *  Created on: Oct 18, 2012
 *      Author: letrungkien7
nnn */

#include "hmm.hpp"
#include <iostream>
#include <cmath>

namespace ltk{
HMM::HMM(int n_, int m_) :
    				n(n_), m(m_) {
	a.resize(n, n);
	b.resize(n, m);
	pi.resize(n, 1);
}

HMM::~HMM(){

}

double HMM::Evaluate(const MatrixXi &observation) {
	MatrixXd c;

	// Forward calculation
	Forward(observation, c);
	// Return likelihood
	return log(c).sum();
}

MatrixXd HMM::Forward(const MatrixXi &observation, MatrixXd &scale) {
	int T = observation.rows();
	MatrixXd fwd(T, n);
	scale.resize(T, 1);

	// 1. Initialization
	fwd.row(0) = (pi.array()*b.col(0).array()).transpose();
	scale(0) = fwd.row(0).sum();
	if(scale(0))
		fwd.row(0)/=scale(0);
	// 2. Induction
	for(int t = 1; t <T; t++){
		fwd.row(t) = fwd.row(t-1)*a;
		for(int i = 0; i<n; ++i)
			fwd(t,i)*=b(i, observation(t));
		scale(t) = fwd.row(t).sum();
		if(scale(t))
			fwd.row(t)/=scale(t);
	}

	// 3. Return
	return fwd;
}

MatrixXd HMM::Backward(const MatrixXi &observation, const MatrixXd &scale){
	int T = observation.rows();
	MatrixXd bwd(T, n);
	// 1. Initialization
	bwd.row(T-1).fill(1.0);
	bwd.row(T-1)/=scale(T-1);

	// 2. Induction
	for(int t = T-2; t>=0; t--){
		for(int i = 0; i< n; ++i){
			double sum = 0;
			for(int j = 0; j<n; ++j){
				sum += a(i,j)*b(j, observation(t))*bwd(t+1, j);
			}
			bwd(t,i)+=sum/scale(t);
		}
	}

	// 3. Return
	return bwd;
}

MatrixXi HMM::Decode(const MatrixXi &observation, double& probability){
	int T = observation.rows();
	MatrixXd delta(T, n);
	MatrixXi phi(T, n);

	double maxWeight;
	double weight;
	int maxState;
	// 1. Initialization
	delta.row(0) = log(pi).transpose() + log(b.col(0)).transpose();

	// 2. Induction
	for(int t = 1; t<T; ++t){
		for(int j = 0; j<n; ++j){
			maxWeight = log((double)delta(t-1, 0))+ log((double)a(0,j));
			for(int i = 1; i < n; ++i){
				weight = log((double)delta(t-1, i)) + log((double)a(i,j));
				if(weight>maxWeight){
					maxWeight = weight;
					maxState = i;
				}
			}

			delta(t, j) = maxWeight + log((double)b(j, observation(t)));
			phi(t, j) = maxState;
		}
	}
	// 3. Find maximum value for time T-1
	maxState = delta.row(T-1).maxCoeff();
	maxWeight = delta(T-1, maxState);

	// 4. Trackback
	MatrixXi path(T,1);
	path(T-1) = maxState;
	for(int t = T-2; t>=0; t++)
		path(t) = phi(path(t+1), t+1);

	probability = exp(maxWeight);
	return path;
}

double HMM::Learn(vector<MatrixXi> &observations, int iterations, double tolerance){
	int K = observations.size();
	int currentIteration = 1;
	bool stop = false;

	// 1. Initialization
	vector<MyMatrix3d> epsilons(K);
	vector<MatrixXd> gamma(K);
	for(int i =0; i<K; i++){
		int T = observations[i].rows();
		epsilons[i].resize(T, n,n);
		gamma[i].resize(T,n);
	}

}
/* Private members access
 *
 */
int HMM::N() const {
	return n;
}
int &HMM::N() {
	return n;
}
int HMM::M() const {
	return m;
}
int &HMM::M() {
	return m;
}
MatrixXd HMM::A() const {
	return a;
}
MatrixXd &HMM::A() {
	return a;
}
MatrixXd HMM::B() const {
	return b;
}
MatrixXd &HMM::B() {
	return b;
}
MatrixXd HMM::PI() const {
	return pi;
}
MatrixXd &HMM::PI() {
	return pi;
}
}
