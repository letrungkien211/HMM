/*
 * Hidden Markov Model
 */

#ifndef HMM_HPP
#define HMM_HPP

#include <vector>
#include <Eigen/Dense>  // eigen: linear algebra library

namespace ltk{

using namespace std;
using namespace Eigen;

class HMM{
private:
	int n;  // Number of states
	int m;  // Number of observation in alphabet
	MatrixXd a;  // State transition probability distribution nxn
	MatrixXd b;  // Observation symbol probability distribution nxm
	MatrixXd pi; // Initial state distribution nx1
	MatrixXd Forward(const MatrixXi &observation, MatrixXd &scale);  // Forward calculation
	MatrixXd Backward(const MatrixXi &observation, const MatrixXd &scale); // Backward calculation
	bool CheckConvergence(double oldLikelihood, double newLikelihood,
			int currentIteration, int maxIteration, double tolerance); // Check convergence in iteration process
public:
	HMM(int n_, int m_); // Constructor
	~HMM();  // Destructor
	// Methods
	double Evaluate(const MatrixXi &observation, bool logarithm = true); // P(O|a,b,pi)
	MatrixXi Decode(const MatrixXi &observation, double &probability); // q= argmax P(O|a,b,pi)
	double Learn(vector<MatrixXi> &observation, int maxIteration, double tolerance); // (a,b,pi) = argmax P(O|a,b,pi)

	void Initialize(int n_, int m_); // Initialization

	// Private members access
	int N() const;
	int M() const;
	MatrixXd A() const;
	MatrixXd B() const;
	MatrixXd PI() const;
	MatrixXd &A();
	MatrixXd &B();
	MatrixXd &PI();
};

}
#endif
