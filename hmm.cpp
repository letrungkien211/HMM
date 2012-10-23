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

double HMM::Evaluate(const MatrixXi &observation, bool logarithm) {
	MatrixXd scale;
	Forward(observation, scale);
	double prob = 0;
	for(int t=0, T=observation.rows(); t<T; t++)
		prob+=log((double)scale(t));
	return logarithm ? exp(prob) : prob;
}


MatrixXd HMM::Forward(const MatrixXi &observation, MatrixXd &scale) {
	int T = observation.rows();
	MatrixXd fwd(T, n);
	scale.resize(T,1);
	scale.fill(0.0);
	for(int i = 0; i<n; ++i)
		scale(0)+=fwd(0,i) = pi(i)*b(i, observation(0));
	if(scale(0))
		fwd.row(0)/=scale(0);

	for(int t = 1; t<T; ++t){
		for(int j=0; j<n; ++j){
			double sum = 0;
			for(int i = 0; i<n; ++i){
				sum+=fwd(t-1,i)*a(i,j);
			}
			scale(t)+=fwd(t, j) = sum*b(j, observation(t));
		}
		if(scale(t))
			fwd.row(t)/=scale(t);
	}

	return fwd;
}

MatrixXd HMM::Backward(const MatrixXi &observation, const MatrixXd &scale){
}

MatrixXi HMM::Decode(const MatrixXi &observation, double& probability){

}

double HMM::Learn(vector<MatrixXi> &observations, int iterations, double tolerance){
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
