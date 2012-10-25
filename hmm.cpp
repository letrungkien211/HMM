/*
 * hmm.cpp
 *
 *  Created on: Oct 18, 2012
 *      Author: letrungkien7
 nnn */

#include "hmm.hpp"
#include <iostream>
#include <cmath>
#include <cstdio>
#include <boost/multi_array.hpp>

namespace ltk {

typedef boost::multi_array<double, 3> array_type;

HMM::HMM(int n_, int m_){
	Initialize(n_,m_);
}

HMM::~HMM() {
}

double HMM::Evaluate(const MatrixXi &observation, bool logarithm) {
	MatrixXd scale;
	// Run forward
	Forward(observation, scale);

	// Sum up scaling factors
	double prob = 0;
	for (int t = 0, T = observation.rows(); t < T; t++)
		prob += log((double) scale(t));
	return logarithm ? exp(prob) : prob;
}

MatrixXd HMM::Forward(const MatrixXi &observation, MatrixXd &scale) {
	int T = observation.rows();
	MatrixXd fwd(T, n);
	scale.resize(T, 1);
	scale.fill(0.0);

	// 1. Initialization
	for (int i = 0; i < n; ++i)
		scale(0) += fwd(0, i) = pi(i) * b(i, observation(0));
	if (scale(0))
		fwd.row(0) /= scale(0);

	// 2. Induction
	for (int t = 1; t < T; ++t) {
		for (int j = 0; j < n; ++j) {
			double sum = 0;
			for (int i = 0; i < n; ++i) {
				sum += fwd(t - 1, i) * a(i, j);
			}
			scale(t) += fwd(t, j) = sum * b(j, observation(t));
		}
		if (scale(t))
			fwd.row(t) /= scale(t);
	}

	return fwd;
}

MatrixXd HMM::Backward(const MatrixXi &observation, const MatrixXd &scale) {
	int T = observation.rows();
	MatrixXd bwd(T, n);

	// 1. Initialization
	bwd.row(T - 1).fill(1.0 / (double) scale(T - 1));

	// 2. Induction
	for (int t = T - 2; t >= 0; t--) {
		for (int i = 0; i < n; ++i) {
			double sum = 0;
			for (int j = 0; j < n; ++j)
				sum += a(i, j) * b(j, observation(t + 1)) * bwd(t + 1, j);
			bwd(t, i) = sum / scale(t);
		}
	}

	return bwd;
}

MatrixXi HMM::Decode(const MatrixXi &observation, double& probability) {
	int T = observation.rows();
	MatrixXi path(T,1);
	MatrixXd delta(T,n);
	MatrixXd phi(T,n);

	int maxState;
	double maxWeight;
	double weight;

	// 1. Initialization
	for(int i = 0; i<n; i++){
		delta(0,i) = log((double)pi(i))+log((double)b(i, observation(0)));
		phi(0,i) = 0;
	}
	// 2. Induction
	for(int t = 1; t< T; ++t){
		for(int j = 0; j<n; ++j){
			maxWeight = delta(t-1,0) + log((double)a(0,j));
			maxState = 0;
			for(int i = 0; i< n; ++i){
				weight = delta(t-1, i) + log((double)a(i,j));
				if(weight>maxWeight){
					maxWeight = weight;
					maxState = i;
				}
			}
			delta(t,j) = maxWeight*b(j, observation(t));
			phi(t,j) = maxState;
		}
	}

	// 3. Termination
	maxState = 0;
	maxWeight = delta(0,T-1);
	for(int i = 0; i< n; ++i){
		if(delta(T-1,i)>maxWeight){
			maxState = i;
			maxWeight = delta(T-1,i);
		}
	}

	// 4. Path
	path(T-1) = maxState;
	for(int t = T-1; t >=0; t--)
		path(t) = phi(t+1, path(t+1));

	probability = exp(maxWeight);
	return path;
}

double HMM::Learn(vector<MatrixXi> &observation, int maxIteration, double tolerance) {
	double newLikelihood, oldLikelihood;
	int K = observation.size();
	vector<MatrixXd> gamma(K);
	vector<array_type> epsilon(K);

	//  Initialization
	for (int k = 0; k < K; k++) {
		int T = observation[k].rows();
		gamma[k].resize(T, n);
		boost::array<array_type::index, 3> shape = {{T,n,n}};
		epsilon[k].resize(shape);
	}

	oldLikelihood = -numeric_limits<double>::max();
	newLikelihood = 0.0;
	int currentIteration = 0;

	// Start the loop
	bool stop = false;
	do {
		// Calculate forward, backward, gamma, epsilon
		for (int k = 0; k < K; k++) {
			MatrixXd scale;
			MatrixXi &cobservation = observation[k]; // current observation
			MatrixXd &cgamma = gamma[k]; // current gamma
			array_type &cepsilon = epsilon[k]; // current epsilon
			int T = cobservation.rows();

			// Run forward, backward
			MatrixXd fwd = Forward(cobservation, scale);
			MatrixXd bwd = Backward(cobservation, scale);

			// Calculate gamma
			for (int t = 0; t < T; t++) {
				double s = 0;
				for (int i = 0; i < n; i++) {
					s += cgamma(t, i) = fwd(t, i) * bwd(t, i);
				}
				if (s)
					cgamma.row(t) /= s;
			}

			// Calculate epsilon
			for (int t = 0; t < T - 1; t++) {
				double s = 0;
				for (int i = 0; i < n; i++) {
					for (int j = 0; j < n; j++) {
						s+=cepsilon[t][i][j] = fwd(t,i)*a(i,j)*b(j,cobservation(t+1))*bwd(t+1,j);
					}
				}
				if (s) {
					for (int i = 0; i < n; i++) {
						for (int j = 0; j < n; j++) {
							cepsilon[t][i][j] = cepsilon[t][i][j]/s;
						}
					}
				}
			}
			// Culminate newLikelihood
			for (int t = 0; t < T; ++t)
				newLikelihood += log((double)scale(t));
		}

		newLikelihood /= K; // averaging
		cout << currentIteration << ": " << newLikelihood << endl; // print current likelihood

		// Check convergence
		if (CheckConvergence(oldLikelihood, newLikelihood, currentIteration,
				maxIteration, tolerance)) {
			stop = true;
		} else {
			// Reset variables
			currentIteration++;
			oldLikelihood = newLikelihood;
			newLikelihood = 0.0;

			// Update initial probabilities
			for (int i = 0; i < n; ++i) {
				double sum = 0;
				for (int k = 0; k < K; k++)
					sum += gamma[k](0, i);
				pi(i) = sum / K;
			}

			// Update transition probabilities
			for (int i = 0; i < n; i++) {

				for (int j = 0; j < n; j++) {
					double denominator = 0, numerator = 0;

					for (int k = 0; k < K; k++) {
						int T = observation[k].rows();

						for (int t = 0; t < T - 1; t++) {
							denominator += epsilon[k][t][i][j];
							numerator += gamma[k](t, i);
						}
					}

					a(i,j) = (numerator) ? denominator/numerator : 0.0;
				}
			}

			// Update emission probabilities
			for (int j = 0; j < n; j++) {
				for (int l = 0; l < m; ++l) {
					double denominator = 0, numerator = 0;
					for (int k = 0; k < K; k++) {
						int T = observation[k].rows();
						for (int t = 0; t < T; ++t) {
							if (observation[k](t) == l)
								denominator += gamma[k](t, j);
							numerator += gamma[k](t, j);
						}
					}
					b(j, l) = (numerator) ? denominator / numerator : 0.0;
				}
			}
		}
	}
	while (!stop);

	return newLikelihood;
}

bool HMM::CheckConvergence(double oldLikelihood, double newLikelihood, int currentIteration, int maxIteration, double tolerance) {
	// Update and verify stop criteria
	if (tolerance > 0)
	{
		// Stopping criteria is likelihood convergence
		if (abs(oldLikelihood - newLikelihood) <= tolerance)
			return true;

		if (maxIteration > 0)
		{
			// Maximum iterations should also be respected
			if (currentIteration >= maxIteration)
				return true;
		}
	}
	else
	{
		// Stopping criteria is number of iterations
		if (currentIteration == maxIteration)
			return true;
	}

	// Check if we have reached an invalid state
	if (isnan(newLikelihood) || isinf(newLikelihood))
	{
		return true;
	}

	return false;
}

void HMM::Initialize(int n_, int m_){
	n = n_;
	m = m_;
	a.resize(n, n);
	b.resize(n, m);
	pi.resize(n, 1);

	a.fill(1.0/n);
	b.fill(1.0/m);
	pi.fill(1.0/n);
}

/* Private members access
 *
 */
int HMM::N() const {
	return n;
}

int HMM::M() const {
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
