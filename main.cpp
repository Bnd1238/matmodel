#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <vector>
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
	if (f(x) == 0)
		return 0;
	double k1 = f(x);
	double k2 = f(x + h * 1 / 4 * k1);
	double k3 = f(x + h * (3 / 32 * k1 + 9 / 32 * k2));
	double k4 = f(x + h * (1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3));
	double k5 = f(x + h * (439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4));
	double k6 = f(x + h * (-8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5));

	return h * k1 * 25 / 216 + h * 1408 / 2565 * k3 + h * 2197 / 4104 * k4 - h * 1 / 5 * k5;
}
double RK45(double x, double (*f)(double), double h)
{
	if (f(x) == 0)
		return 0;
	double k1 = f(x);
	double k2 = f(x + h * 1 / 4 * k1);
	double k3 = f(x + h * (3 / 32 * k1 + 9 / 32 * k2));
	double k4 = f(x + h * (1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3));
	double k5 = f(x + h * (439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4));
	double k6 = f(x + h * (-8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5));

	return h * k1 * 16 / 135 + h * 6656 / 12825 * k3 + h * 28561 / 56430 * k4 - h * 9 / 50 * k5+ h * 2/55*k6;
}
double Trapeze(double x, double (*f)(double), double h)
{
	if (f(x) ==0)
		return 0;
	return h * f(x) + h * h / 4 * -1 / (2 * f(x)) * f(x);
}

vector<double>start;
vector<double>gamma
{
	1.0,
	0.5 ,
	5.0 / 12.0,
	3.0 / 8.0,
	251.0 / 720.0,
	95.0 / 288.0,
	19087.0 / 60480.0,
	5257.0 / 17280.0,
	1070017.0 / 3628800.0,
	25713.0 / 89600.0,
	26842253.0 / 95800320.0 ,
};
double sup(int i,int n)
{
	switch (i)
	{
		case(0):
			return start[n];
		case(1):
			return start[n]- start[n-1];
		case(2):
			return start[n]-2* start[n-1]+ start[n-2];
		case(3):
			return start[n]-3* start[n-1]+3* start[n-2]- start[n-3];
		case(4):
			return start[n] - 4 * start[n - 1] + 6 * start[n - 2] - 4 * start[n - 3] + start[n - 4];
		case(5):
			return start[n] - 5 * start[n - 1] + 10 * start[n - 2] - 10 * start[n - 3] +5* start[n - 4]- start[n-5];
		case(6):
			return start[n] - 6 * start[n - 1] + 15 * start[n - 2] - 20 * start[n - 3] + 15 * start[n - 4] - 6* start[n - 5]+ start[n-6];
		case(7):
			return start[n] - 7 * start[n - 1] + 21 * start[n - 2] - 35 * start[n - 3] + 35 * start[n - 4] -21* start[n - 5]+7* start[n-6]- start[n-7];
		case(8):
			return start[n] - 8 * start[n - 1] + 28 * start[n - 2] - 56 * start[n - 3] + 70 * start[n - 4] - 56 * start[n - 5] + 28 * start[n - 6] -8* start[n - 7]+ start[n-8];
		case(9):
			return start[n] - 9 * start[n - 1] + 36 * start[n - 2] - 84 * start[n - 3] + 120 * start[n - 4] - 126 * start[n - 5] +84 * start[n - 6] - 36 * start[n - 7] +9* start[n - 8]- start[n-9];
	}
}

double Adams(double h,int n)
{
	double sum = 0;
	for (int i = 0; i < 10; ++i)
	{
		sum = sum + gamma[i] * sup(i,n);
	}
	start.push_back(start[n] + h * sum);
	return h * sum;
}
int main()
{
	freopen("a.txt", "w", stdout);
	double h = 0.01;
	double E = 50;
	double RK = E;
	double T = E;
	double P = E;
	double A = E;
	for (int i = 0; i < 11; ++i)
	{
		start.push_back(Target(50, h * i));
	}
	double n = 10;
	for (double t = 0; t < 5; t = t + h)
	{
		Adams(h, n);
		n++;
	}
	n = 0;
	for (double t = 0; t<5; t = t + h)
	{
		cout <<t<<" "<< setprecision(5) <<T<<" "<< RK<<" "<<start[n]<<" "<<Target(50,t) << " "<<endl;
		RK = RK + RK45(RK, Law, h);
		T = T + Trapeze(T, Law, h);
		n++;
		if (T <= 0)
			T = 0;
		if (RK <= 0)
			RK = 0;
	}
	return 0;
}
