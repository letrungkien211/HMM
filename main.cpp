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
	hmm.PI() << 0.5,  0.5, 0.5;
	MatrixXi ob(3,1);
	MatrixXd scale;

	ob << 0,0,1;
	hmm.Forward(ob, scale);
	return 0;
}
