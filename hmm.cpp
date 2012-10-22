/*
 * hmm.cpp
 *
 *  Created on: Oct 18, 2012
 *      Author: letrungkien7
 */

#include "hmm.hpp"
#include <iostream>

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
    bwd.row(T-1).fill(1.0/scale(T-1));
  
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
