#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <fstream>
using namespace std;

double h;

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
double RK4(double x)// Классический метод Рунге-Кутты 4-ого порядка
{
	double k1 = f(x);
	double k2 = f(x + h * 1.0 / 2.0 * k1);
	double k3 = f(x + h * 1.0 / 2.0 * k2);
	double k4 = f(x + h * k3);
	return x + (h * k1 + h * 2 * k2 + h * 2 * k3 + h * k4) / 6;

}

double err(double x)//Нахождение ошибки с помощью RK4(5)
{
	if (f(x) == 0)
		return 0;
	double k1 = f(x);
	double k2 = f(x + h * 1.0 / 4.0 * k1);
	double k3 = f(x + h * (3.0 / 32.0 * k1 + 9 / 32 * k2));
	double k4 = f(x + h * (1932.0 / 2197.0 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3));
	double k5 = f(x + h * (439.0 / 216.0 * k1 - 8.0 * k2 + 3680.0 / 513.0 * k3 - 845.0 / 4104.0 * k4));
	double k6 = f(x + h * (-8.0 / 27.0 * k1 + 2.0 * k2 - 3544.0 / 2565.0 * k3 + 1859.0 / 4104.0 * k4 - 11.0 / 40.0 * k5));

	return  h * k1 * (16.0 / 135.0-25.0 / 216.0 ) + h * (6656.0 / 12825.0-1408.0 / 2565.0 ) * k3 + h * (28561.0 / 56430.0 -2197.0 / 4104.0 ) * k4 - h * (9.0 / 50.0-1.0 / 5.0 ) * k5 - h * 2.0 / 55.0 * k6;
}
void Hopt(double E,double tol)// поиск оптимального шага
{
	
	double error = err(E);
	if (abs(error) < tol)
		return;
	error = sqrt(pow(error / tol, 2));
	//h = h /2;
	cout << setprecision(20) << error<< endl;
	Hopt(E,tol);
}
double RK45(double x, double tol)//Вложенный метод Рунге-Кутты 4(5)
{
	Hopt(x,tol);
	if (f(x) == 0)
		return 0;
	double k1 = f(x);
	double k2 = f(x + h * 1.0 / 4.0 * k1);
	double k3 = f(x +  (h * 3.0 / 32.0 * k1 + h * 9 / 32 * k2));
	double k4 = f(x +  (h * 1932.0 / 2197.0 * k1 - h * 7200.0 / 2197.0 * k2 + h * 7296.0 / 2197.0 * k3));
	double k5 = f(x +  (h * 439.0 / 216.0 * k1 - h * 8.0 * k2 + h * 3680.0 / 513.0 * k3 - h * 845.0 / 4104.0 * k4));
	double k6 = f(x +  (h * -8.0 / 27.0 * k1 + h * 2.0 * k2 - h * 3544.0 / 2565.0 * k3 + h * 1859.0 / 4104.0 * k4 - h * 11.0 / 40.0 * k5));

	return x + h * k1 * 25.0 / 216.0 + h * 1408.0 / 2565.0 * k3 + h * 2197.0 / 4104.0 * k4 - h * 1.0 / 5.0 * k5;
}
double RK44(double x,double tol)//Вложенный метод Рунге-Кутты 4(5)
{
	Hopt(x, tol);
	if (f(x) == 0)
		return 0;
	double k1 = f(x);
	double k2 = f(x + h * 1.0 / 4.0 * k1);
	double k3 = f(x + h * (3.0 / 32.0 * k1 + 9 / 32 * k2));
	double k4 = f(x + h * (1932.0 / 2197.0 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3));
	double k5 = f(x + h * (439.0 / 216.0 * k1 - 8.0 * k2 + 3680.0 / 513.0 * k3 - 845.0 / 4104.0 * k4));
	double k6 = f(x + h * (-8.0 / 27.0 * k1 + 2.0 * k2 - 3544.0 / 2565.0 * k3 + 1859.0 / 4104.0 * k4 - 11.0 / 40.0 * k5));

	return x + h * k1 * 16.0 / 135.0 + h * 6656.0 / 12825.0 * k3 + h * 28561.0 / 56430.0 * k4 - h * 9.0 / 50.0 * k5+h* 2.0 / 55.0 *k6;
}

double Trapeze(double x)//Метод неявной тапеции
{
	if (f(x) == 0)
		return 0;
	return x+h * f(x) + h * h / 4 * -1 / (2 * f(x)) * f(x);
}

vector<double>startF;
vector<double>gamma
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
};
double sup(int i, int n)//обратные конечные разности
{
	if (i == 0)
		return startF[n];
	return sup(i - 1, n) - sup(i - 1, n - 1);
}
void starter(double E)//разгон методом Рунге-Кутты 4-ого порядка
{
	double time = 0;
	double rk = E;
	for (int i = 0; i < 11; ++i)
	{
		startF.push_back(f(rk));
		rk = RK4(rk);
		time += h;
	}
}
double Adams(double x)//Метод Адамса-Башфорта 11-ого порядка
{
	double sum = 0;
	for (int i = 0; i < 11; ++i)
	{
		sum = sum + gamma[i] * sup(i, 10);
	}
	x = x + h * sum;
	startF.erase(startF.begin(), startF.begin() + 1);
	startF.push_back(f(x));
	return x;
}

int main()
{
	h = 2;
	double E = 50;
	double RK = E;
	double T = E;
	double rk44 = E;
	double A = E;
	for (double t = 0; Target(E, t) > 0; t += h)
	{
		cout << t << setprecision(8) << " " << abs((Target(E, t) - RK)) << " " << abs((Target(E, t) - T)) << endl;
		T = RK44(T,0.0001);
		RK = RK4(RK);
	}
	return 0;
}
