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
	hmm.A() << 0.75,0.25,0.1,0.9;
	hmm.B() << 0.5,0.5,0.5,0.5;
	hmm.PI() << 0.75,  0.25;
	MatrixXi ob(3,1);
	MatrixXd scale;

	vector<MatrixXi> observation(5);
	int len[] = {7,4,10,2,3};
	for(int i=0; i<5; i++){
		observation[i].resize(len[i],1);
	}
	observation[0]<<0,1,1,1,0,1,1;
	observation[1]<< 0,1,1,1;
	observation[2] <<0,1,1,1,0,1,1,1,1,1;
	observation[3]<< 0,1;
	observation[4] <<0,1,1;

	ob << 0,1,1;

	hmm.Learn(observation, 20, 0.001);
	cout <<hmm.A()<<endl<<endl;
	cout <<hmm.B()<<endl<<endl;
	cout<< hmm.Evaluate(ob)<<endl;
	return 0;
}
