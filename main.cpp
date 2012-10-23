/*
 * main.cpp
 *
 *  Created on: Oct 18, 2012
 *      Author: letrungkien7
 */

#include <iostream>
#include <cstdlib>

#include "hmm.hpp"

using namespace std;
using namespace ltk;

int main(int argc, char* argv[]){
	HMM hmm(2,2);
	hmm.A() << 0.5,0.5,0.0,1.0;
	hmm.B() << 0.5,0.5,0.5,0.5;
	hmm.PI() << 0.75,  0.25;
	MatrixXi ob(2,1);
	MatrixXd scale;

	ob << 1,0;
	cout<< hmm.Forward(ob, scale)<<endl;
	cout<< hmm.Evaluate(ob)<<endl;
	return 0;
}
