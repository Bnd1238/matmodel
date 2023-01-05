#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <fstream>
using namespace std;

double h;
double E = 54;
double RK = E;
double T = E;


double f(double x)
{
	if (x >= 0)
		return -sqrt(x);
	return 0;
}
double Target(double x0, double t)//Аналитическое решение
{
	if (0 <= t && t < 2 * sqrt(x0))
	{
		return (sqrt(x0) - t / 2) * (sqrt(x0) - t / 2);
	}
	return 0;
}
double RK4(double x)
{
	double k1 = f(x);
	double k2 = f(x + h * 1.0 / 2.0 * k1);
	double k3 = f(x + h * 1.0 / 2.0 * k2);
	double k4 = f(x + h * k3);
	return x + (h * k1 + h * 2 * k2 + h * 2 * k3 + h * k4) / 6;

}
double Trapeze(double x)//Метод неявной тапеции
{
	if (f(x) == 0)
		return 0;
	return x + h * f(x + h / 2.0 * f(x));
}

double RK44(double x, double tol)//Вложенный метод Рунге-Кутты 4(5)
{
	bool ch = 1;
	double y, z;
	double k1;
	double k2;
	double k3;
	double k4;
	double k5;
	double k6;
	
	while(ch)
	{
		ch = 0;
		if (f(x) == 0)
			return 0;
		k1 = f(x);
		k2 = f(x + h * 1.0 / 4.0 * k1);
		k3 = f(x + h * (3.0 / 32.0 * k1 + 9.0/ 32.0 * k2));
		k4 = f(x + h * (1932.0 / 2197.0 * k1 - 7200.0 / 2197.0 * k2 + 7296.0 / 2197.0 * k3));
		k5 = f(x + h * (439.0 / 216.0 * k1 - 8.0 * k2 + 3680.0 / 513.0 * k3 - 845.0 / 4104.0 * k4));
		k6 = f(x + h * (-8.0 / 27.0 * k1 + 2.0 * k2 - 3544.0 / 2565.0 * k3 + 1859.0 / 4104.0 * k4 - 11.0 / 40.0 * k5));
		
		z = x + 16.0 / 13.0  * k1 +  6656.0 / 12825.0  * k3 + 28561.0 / 56430.0  * k4 - 9.0 / 50.0  * k5 + 2.0 / 55.0  * k6;
		y = x + 25.0 / 216.0  * k1 + 1408.0 / 2565.0  * k3 + 2197.0 / 4104.0  * k4 -  k5 / 5.0;
		
		double error = h*abs(z - y);
		z = x + h * k1 * 16.0 / 135.0 + h * 6656.0 / 12825.0 * k3 + h * 28561.0 / 56430.0 * k4 - h * 9.0 / 50.0 * k5 + h * 2.0 / 55.0 * k6;
		if (abs(error  / z) > tol)
		{
			ch = 1;
			h = h * pow(tol*h / error, 1.0 / 4.0);
		}
	}
	return z;
}



double startF[]=
{
	0,			//0
	0,			//1
	0,			//2
	0,			//3
	0,			//4
	0,			//5
	0,			//6
	0,			//7
	0,			//8
	0,			//9
};
long long S = 0,F=0;

double gamma[]=
{
	1.0,						//0
	0.5 ,						//1
	5.0 / 12.0,					//2
	3.0 / 8.0,					//3
	251.0 / 720.0,				//4
	95.0 / 288.0,				//5
	19087.0 / 60480.0,			//6
	5257.0 / 17280.0,			//7
	1070017.0 / 3628800.0,		//8
	25713.0 / 89600.0,			//9
	26842253.0 / 95800320.0 ,	//10
	//4777223.0 / 17418240.0		//11
};
double sup(int i, int n)
{
	if (i == 0)
		return startF[(S+n)%10];
	return sup(i - 1, n) - sup(i - 1, n - 1);
}
void starter()
{
	double t = 0;
	T = E;
	for (int i = 0; i < 10; ++i, t += h)
	{
		cout << setprecision(8) << " " << t << " " << abs(Target(E, t) - T) / Target(E, t) << endl;
		startF[F%10]=f(T);
		F++;
		T = RK4(T);
	}
}
double Adams(double x)
{
	double sum = 0;
	for (int i = 0; i <10; ++i)
	{
		sum = sum + gamma[i] * sup(i, 9);
	}
	x = x + h * sum;
	startF[S%10] = f(x);	
	F++;
	S++;
	return x;
}

int main()
{
	long long calls = 0;
	h = 2;

	bool check = 0;
	double tol = 1;
	double mis = 1e-6;
	cout << "RK4:1\nTrapeze:2\nRK4(5):3\nAB11:4\n";
	int n;
	cin >> n;
	//freopen("a.txt", "w", stdout);

	switch (n)
	{
	case(1):
		h = 0.18;
	
		for (double t = 0; Target(E, t) > 1; t += h)
		{
			cout << setprecision(8) << " " << t << " " << abs(Target(E, t) - T) / Target(E, t) << endl;
			T = RK4(T);
			calls++;
		}
		cout << "calls=" << calls;
		break;
	case(2):
		h = 0.004;
		T = E;
		for (double t = 0; Target(E, t) > 1; t += h)
		{
			cout << setprecision(8) << " " << t << " " << abs(Target(E, t) - T) / Target(E, t) << endl;
			T =Trapeze(T);
			calls++;
		}
		cout << "calls=" << calls;
		break;
	case(3):
		tol = 0.35;
		T = E; 
		h = 1;
		for (double t = 0; Target(E, t) > 1; t += h)
		{
			cout << setprecision(8) << " " << t << " " << abs(Target(E, t) - T) / Target(E, t) << endl;
			//cout << setprecision(8) << " " << t << " " << Target(E, t)<<" "<<T<< endl;
			T = RK44(T, tol);
			calls++;
		}
		cout << "calls=" << calls<<" Rtol="<<tol;
		break;
	case (4):
		h = 0.005;
		T = E;
		starter();
		F = 9;
		for (double t = 10*h; Target(E, t) > 1; t += h)
		{
			cout << setprecision(8) << " " << t << " " << abs(Target(E, t) - T) / Target(E, t) << endl;
			//cout << setprecision(8) << " " << t << " " << Target(E, t) << " " << T << endl;
			T = Adams(T);
			calls++;
		}
		cout << "calls:" << calls+10;
	}
	return 0;
}
