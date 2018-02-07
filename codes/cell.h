#pragma once
#include <cmath>
class cell
{
public:
	cell();
	~cell();
	struct node {
		double x, y;//indicate the position of the particle
		bool i, j, k, l, E, N, W, S;
		double phi;
		double get_phi() {
			//the ceter point of the circle is(x0, y0), and ridics is r;
			double x0 = 0.5; double y0 = 0.5;
			double r = 0.2;
			phi = r - sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0));
			return phi;
		}
		node()
		{
			x = 0.0; y = 0.0;
			i = false; j = false; k = false; l = false;
			E = false; N = false; W = false; S = false;
		}

		void operator=(const node temp) {
			x = temp.x; y = temp.y;
			i = temp.i; j = temp.j; k = temp.k; l = temp.l;
			E = temp.E; N = temp.N; W = temp.W; S = temp.S;
		}

		node(double fx, double fy)
		{
			x = fx; y = fy;
			i = false; j = false; k = false; l = false;
			E = false; N = false; W = false; S = false;
		}
		~node() {}
	}elem[1500][900];

	void initial();
	void output(int time);
	void impact();
	void move();
	double get_theta(const node & norign, const node & nodeb);
	void cmain();

private:
	int m;
	int n;
	double dx;
	double dy;
	int maxTime;
	double stheta;//the standary line of open boundary
	int boubdaryTheta;//the open boundary of the circle;
};

