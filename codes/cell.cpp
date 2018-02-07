#include "stdafx.h"
#include "cell.h"
#include <ctime>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <sstream>
#include <cmath>
using namespace std;

int my_rand(int);

cell::cell()
{
	m = 200;
	n = m;//200
	dx = 1.0 / m;
	dy = 1.0 / n;
	maxTime = 1500;
	stheta = -0.0;          //the standary line of open boundary
	boubdaryTheta = 5;     //the open boundary of the circle;

	int imax = m + 1;
	int jmax = n + 1;
	for (int i = 0; i < imax; ++i) {
		for (int j = 0; j < jmax; ++j) {
			//node temp(i*dx, j*dy);
			elem[i][j].x = i * dx;
			elem[i][j].y = j * dy;
			elem[i][j].get_phi();
		}
	}

}

void cell::initial() {

	int imax = 1500;
	int jmax = 900;
	//time_t nowtime = time(nullptr);
	int nowtime = 2;
	auto seed = my_rand(nowtime);
	int random = seed;
	for (int i = 0; i < imax; ++i) {
		for (int j = 0; j < jmax; ++j) {
			elem[i][j].i = false; elem[i][j].j = false; 
			elem[i][j].k = false; elem[i][j].l = false;
			if (elem[i][j].phi>0) {
				seed = my_rand(seed);
				random = 1.0 * (seed% m) / m;
				if (random < 0.25) {
					elem[i][j].i = true;    elem[i][j].k = true;
				}
				else if (random < 0.5) {
					elem[i][j].j = true;    elem[i][j].l = true;
				}
				else if (random < 0.75) {
					elem[i][j].k = true;    elem[i][j].i = true;
				}
				else if (random<1) {
					elem[i][j].l = true;    elem[i][j].j = true;
				}
				else {
					cout << "random error in initial!" << endl; cin.get();
				}
				//cout << "come to here " << endl;
			}
		}
	}
}

void cell::output(int time) {
	string tempa = "axissym_";

	int temp = 10000 + time;
	stringstream ss;
	string str;
	ss << temp;
	ss >> str;

	string title = tempa + str + ".plt";
	ofstream fcout;
	fcout.open(title);

	int imax = m + 1;
	int jmax = n + 1;
	int maxParticleNum = 0;
	for (int i = 0; i < imax; ++i) {
		for (int j = 0; j < jmax; ++j) {
			if (elem[i][j].i || elem[i][j].j || elem[i][j].k || elem[i][j].l) {
				maxParticleNum++;
			}
		}
	}

	fcout << " title=random  " << endl;
	fcout << "  VARIABLES = \"X\",\"Y \" " << endl;
	fcout << "zone I= " << maxParticleNum << ",datapacking=POINT" << endl;

	for (int i = 0; i < imax; ++i) {
		for (int j = 0; j < jmax; ++j) {
			if (elem[i][j].i || elem[i][j].j || elem[i][j].k || elem[i][j].l) {
				fcout << i * dx << "  " << j * dy << endl;
			}
		}
	}

	fcout.close();
}

cell::~cell(){}

void cell::impact()
{
	//cout << "in impact !" << endl;
	int imax = 1500;
	int jmax = 900;
	for (int i = 1; i < imax - 1; ++i) {
		for (int j = 1; j < jmax - 1; ++j) {
			elem[i][j].i = false; elem[i][j].j = false; 
			elem[i][j].k = false; elem[i][j].l = false;
			if (elem[i][j].E || elem[i][j].N || !elem[i][j].W || !elem[i][j].S) {
				if (elem[i][j].E && !elem[i][j].N && !elem[i][j].W && !elem[i][j].S) 
				{ //one particle comes fome East
					elem[i][j].k = true;
				}
				else if (!elem[i][j].E && !elem[i][j].N && elem[i][j].W && !elem[i][j].S)
				{ //one particle comes fome West
					elem[i][j].i = true;
				}
				else if (!elem[i][j].E && elem[i][j].N && !elem[i][j].W && !elem[i][j].S)
				{ //one particle comes fome North
					elem[i][j].l = true;
				}
				else if (!elem[i][j].E && !elem[i][j].N && !elem[i][j].W && elem[i][j].S)
				{ //one particle comes fome South
					elem[i][j].j = true;
				}
				else if (!elem[i][j].E && elem[i][j].N && !elem[i][j].W && elem[i][j].S)
				{ //two particles come fome South and North
					elem[i][j].k = true;    elem[i][j].i = true;
				}
				else if (elem[i][j].E && !elem[i][j].N && !elem[i][j].W && elem[i][j].S)
				{ //two particles come fome South and East
					elem[i][j].j = true;    elem[i][j].k = true;
				}
				else if (!elem[i][j].E && !elem[i][j].N && elem[i][j].W && elem[i][j].S)
				{ //two particles come fome South and West
					elem[i][j].j = true;    elem[i][j].i = true;
				}
				else if (!elem[i][j].E && elem[i][j].N && elem[i][j].W && !elem[i][j].S)
				{ //two particles come fome North and West
					elem[i][j].i = true;    elem[i][j].l = true;
				}
				else if (elem[i][j].E && elem[i][j].N && !elem[i][j].W && !elem[i][j].S) 
				{ //two particles come fome North and East
					elem[i][j].l = true;    elem[i][j].k = true;
				}
				else if (elem[i][j].E && !elem[i][j].N && elem[i][j].W && !elem[i][j].S)
				{ //two particles come fome East and West
					elem[i][j].j = true;    elem[i][j].l = true;
				}
				else if (!elem[i][j].E && elem[i][j].N && elem[i][j].W && elem[i][j].S) 
				{ //three particles except coming fome East
					elem[i][j].j = true;    elem[i][j].i = true;    elem[i][j].l = true;
				}
				else if (elem[i][j].E && elem[i][j].N && !elem[i][j].W && elem[i][j].S)
				{ //three particles except coming fome West
					elem[i][j].j = true;    elem[i][j].k = true;    elem[i][j].l = true;
				}
				else if (elem[i][j].E && elem[i][j].N && elem[i][j].W && !elem[i][j].S)
				{ //three particles except coming fome South
					elem[i][j].i = true;    elem[i][j].k = true;    elem[i][j].l = true;
				}
				else if (elem[i][j].E && !elem[i][j].N && elem[i][j].W && elem[i][j].S)
				{ //three particles except coming fome North
					elem[i][j].i = true;    elem[i][j].k = true;    elem[i][j].j = true;
				}
				else if (elem[i][j].E && elem[i][j].N && elem[i][j].W && elem[i][j].S)
				{ //four particles
					elem[i][j].i = true; elem[i][j].k = true; elem[i][j].j = true; elem[i][j].l = true;
				}
			}

		}
	}

	for (int i = 0; i < imax; ++i) { //the uppest and the lowest boundary
		elem[i][0].i = false; elem[i][0].j = false;
		elem[i][0].k = false; elem[i][0].l = false;
		elem[i][jmax - 1].i = false; elem[i][jmax - 1].j = false;
		elem[i][jmax - 1].k = false; elem[i][jmax - 1].l = false;
		if (elem[i][0].N) {
			int temp = 0;
			elem[i][0].j = true;
		}
		if (elem[i][jmax - 1].S) {
			elem[i][jmax - 1].l = true;
		}
	}
	int temp = 0;
	for (int j = 0; j < jmax; ++j) { //the left and the right boundary
		elem[0][j].i = false; elem[0][j].j = false; 
		elem[0][j].k = false; elem[0][j].l = false;
		elem[imax - 1][j].i = false; elem[imax - 1][j].j = false; 
		elem[imax - 1][j].k = false; elem[imax - 1][j].l = false;
		if (elem[0][j].E) {
			elem[0][j].i = true;
		}
		if (elem[imax - 1][j].W) {
			elem[imax - 1][j].k = true;
		}
	}
	temp = 0;
}

void cell::move()
{
	//cout << "in move !" << endl;
	node norign(0.5, 0.5);

	double theta = 0.0;

	int imax = 1500;
	int jmax = 900;
	for (int i = 0; i < imax; ++i) {
		for (int j = 0; j < jmax; ++j) {
			elem[i][j].E = false; elem[i][j].N = false; 
			elem[i][j].W = false; elem[i][j].S = false;
		}
	}
	for (int i = 0; i < imax; ++i) {
		for (int j = 0; j < jmax; ++j) {
			if (elem[i][j].i) {
				if (elem[i + 1][j].phi * elem[i][j].phi <= 0) {
					theta = get_theta(norign, elem[i][j]);
					if (abs(theta - stheta) > boubdaryTheta) {
						elem[i][j].E = true;
					}
					else {
						elem[i + 1][j].W = true;
					}
				}
				else {
					elem[i + 1][j].W = true;
				}
			}
			if (elem[i][j].j) {
				if (elem[i][j].phi * elem[i][j + 1].phi <= 0) {
					theta = get_theta(norign, elem[i][j]);
					if (abs(theta - stheta) > boubdaryTheta) {
						elem[i][j].N = true;
					}
					else {
						elem[i][j + 1].S = true;
					}
				}
				else {
					elem[i][j + 1].S = true;
				}
			}
			if (elem[i][j].k) {
				if (elem[i][j].phi * elem[i - 1][j].phi <= 0) {
					theta = get_theta(norign, elem[i][j]);
					if (abs(theta - stheta) > boubdaryTheta) {
						elem[i][j].W = true;
					}
					else {
						elem[i - 1][j].E = true;
					}
				}
				else {
					elem[i - 1][j].E = true;
				}
			}
			if (elem[i][j].l) {
				if (elem[i][j].phi * elem[i][j - 1].phi <= 0) {
					theta = get_theta(norign, elem[i][j]);
					if (abs(theta - stheta) > boubdaryTheta) {
						elem[i][j].S = true;
					}
					else {
						elem[i][j - 1].N = true;
					}
				}
				else {
					elem[i][j - 1].N = true;
				}
			}

		}
	}
}


double cell::get_theta(const node & norign, const node & nodeb)
//get the theta by giving origin point and another point
{
	const double PI = 3.14159;
	double x = nodeb.x - norign.x;
	double y = nodeb.y - norign.y;
	double thet = 0.0;

	if (x == 0) {
		if (y > 0) {
			thet = 90;
		}
		else if (y < 0) {
			thet = -90;
		}
	}
	else {
		thet = atan2(y, x);
	}
	thet = thet * 180 / PI;

	return thet;
}

void cell::cmain() {

	cell::node anode(0, 0);
	cell::node bnode(-1, -2);
	double th = get_theta(anode, bnode);
	cout << th << endl;
	initial();
	output(0);
	for (int i = 0; i < maxTime; ++i) {
		output(i);
		cout << "time = " << i << endl;
		for (int j = 0; j != 20; ++j, ++i) {
			/*cout << "time = " << i << endl;*/
			move();
			impact();
		}
	}
}

int my_rand(int z)
// 16807 way to create random numbers
// z is the seed number, num is the total random number to create
{
	//z(n+1)=(a*z(n)+b) mod m
	//describe m=a*q+r to avoid that 
	//the number is large than the computer can bear
	const int m = pow(2, 31) - 1;
	const int a = 16807;
	const int q = 127773;
	const int r = 2836;

	int random = a * (z%q) - r * (z / q);

	if (random<0)
	{
		random = m + random;
	}
	//z is the seed number
	//z = temp;
	//double t = z*1.0 / m;

	//cRandom cr;
	//cr.random = t;
	//cr.seed = z;

	return random;
}
