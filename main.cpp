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
	HMM hmm(3,3);
	hmm.A() << 1,2,3,4,5,6,7,8,9;
	return 0;
}
