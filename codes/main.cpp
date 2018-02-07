// LGA.cpp
// MCM/ICM modeling 
#include "stdafx.h"
#include "cell.h"
#include <iostream>
#include <ctime>
#include <cmath>

using namespace std;

int my_rand(int);

int main() {

	cell solver;

	solver.cmain();
	cout << sizeof(solver) << endl;

	int temp = 2;
	int timeseed = my_rand(2);
	time_t nowtime = time(nullptr);
	//int nowtime = time(nullptr);
	cout << "press any button to exit "<< timeseed << " " << nowtime << endl;
	system("pause");
	return 0;
}

