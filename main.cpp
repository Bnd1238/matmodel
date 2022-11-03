#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdio.h>
using namespace std;

double Law(double x)
{
	if (x >= 0)
		return -sqrt(x);
	return 0;
}

double Target(double x0, double t)
{
	if (0 <= t && t < 2 * sqrt(x0))
	{
		return (sqrt(x0) - t / 2) * (sqrt(x0) - t / 2);
	}
	return 0;
}
double RK44(double x, double (*f)(double), double h)
{
	double k1 = f(x);
	double k2 = f(x + h * 1 / 4 * k1);
	double k3 = f(x + h * (3 / 32 * k2 + 9 / 32 * k1));
	double k4 = f(x + h * (1932 / 2197 * k3 - 7200 / 2197 * k2 + 7296 / 2197 * k1));
	double k5 = f(x + h * (439 / 216 * k4 - 8 * k3 + 3680 / 513 * k2 - 845 / 4104 * k1));
	double k6 = f(x + h * (-8 / 27 * k5 + 2 * k4 - 3544 / 2565 * k3 + 1859 / 4104 * k2 - 11 / 40 * k1));
	return h * (k1* 25 / 216 + 1408 / 2565 * k3 + 2197 / 4104* k4-1/5*k5);
}
double RK45(double x, double (*f)(double), double h)
{
	double k1 = f(x);
	double k2 = f(x + h * 1 / 4 * k1);
	double k3 = f(x + h * (3 / 32 * k2 + 9 / 32 * k1));
	double k4 = f(x + h * (1932 / 2197 * k3 -7200 / 2197 * k2 + 7296 / 2197 * k1));
	double k5 = f(x + h * ( 439 / 216 * k4 -8 * k3 + 3680 / 513 * k2 - 845 / 4104 * k1));
	double k6 = f(x + h * (-8 / 27 * k5 + 2 * k4 - 3544 / 2565 * k3 + 1859 / 4104 * k2 - 11 / 40 * k1));
	return h * (k1 * 16 / 135 + 6656 / 12825 * k3 + 28561 / 56430 * k4 - 9 / 50 * k5+2/55*k6);
}

double Trapeze(double x, double (*f)(double), double h)
{
	if (f(x) ==0)
		return 0;
	return h * f(x) + h * h / 4 * -1 / (2 * f(x)) * f(x);
}


int main()
{
	//freopen("a.txt", "w", stdout);
	double h = 1;
	double E = 50;
	double RK = E;
	double T = E;
	double P = E;
	for (double t = 0; Target(50, t) > 0; t = t + h)
	{
		//cout << setprecision(20) << abs(Purpose(50, t) - RK) << endl;;
		//cout  << setprecision(5) << Purpose(50, t) << " " << " " << E <<" "<<T <<endl;
		cout << setprecision(5) <<"Time:"<< t<<" RK:"<<RK <<" Trapeze:"<<T<<" Target:"<<Target(50,t)<<endl;
		RK = RK + RK45(RK, Law, h);
		T = T + Trapeze(T,Law,h);
		if (T <= 0)
			T = 0;
	}
	return 0;
}
