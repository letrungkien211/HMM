/*
 * hmm.cpp
 *
 *  Created on: Oct 18, 2012
 *      Author: letrungkien7
 */

#include "hmm.hpp"
#include <iostream>

HMM::HMM(int n_, int m_) :
n(n_), m(m_) {
	a.resize(n, n);
	b.resize(n, m);
	pi.resize(n, 1);
}

double HMM::Evaluate(const MatrixXi &observation) {
	// 3. Termination
	MatrixXd fwd = Forward(observation);
	return fwd.row(observation.rows()-1).sum();
}

MatrixXd HMM::Forward(const MatrixXi &observation) {
	int T = observation.rows();
	MatrixXd fwd(T, n);
	double scale = 0;

	// 1. Initialization
	fwd.row(0) = (pi.array()*b.col(0).array()).transpose();
	scale = fwd.row(0).sum();
	if(scale)
		fwd.row(0)/=scale;
	// 2. Induction
	for(unsigned int t = 1; t <T; t++){
		fwd.row(t) = fwd.row(t-1)*a;
		for(unsigned int i = 0; i<n; ++i)
			fwd(t,i)*=b(i, observation(t));
		scale = fwd.row(t).sum();
		if(scale)
			fwd.row(t)/=scale;

	}
	return fwd;
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
