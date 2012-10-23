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

int x[] = { 0,1,1,1,1,0,1,1,1,1,
		0,1,1,1,0,1,1,1,1,1,
		0,1,1,1,1,1,1,1,1,1,
		0,1,1,1,1,1       ,
		0,1,1,1,1,1,1       ,
		0,1,1,1,1,1,1,1,1,1 ,
		0,1,1,1,1,1,1,1,1,1 };
int len[] ={10,10,10,6,7,10,10};


int main(int argc, char* argv[]){
	HMM hmm(2,2);
	hmm.A() << 0.75,0.25,0.1,0.9;
	hmm.B() << 0.5,0.5,0.5,0.5;
	hmm.PI() << 0.75,  0.25;
	MatrixXi ob(9,1);
	MatrixXd scale;

	vector<MatrixXi> observation(7);
	for(int i=0; i<7; i++){
		observation[i].resize(len[i],1);
	}
	int index =0;
	for(int i=0; i<7; ++i){
		cout <<index;
		for(int j=0; j<len[i]; j++){
			observation[i](j) = x[index++];
		}
	}

	ob<<0,1,0,1,1,1,1,1,1;

	hmm.Learn(observation, 100, 0.0001);
	cout <<hmm.A()<<endl<<endl;
	cout <<hmm.B()<<endl<<endl;
	cout<< hmm.Evaluate(ob)<<endl;
	return 0;
}
