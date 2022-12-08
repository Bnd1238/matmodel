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
double RK4(double x)
{
	double k1 = f(x);
	double k2 = f(x+h*1/2*k1);
	double k3 = f(x+h*1/2*k2);
	double k4 = f(x+h*k3);
	return x + (h*k1 + h*2 * k2 + h*2 * k3 + h*k4) / 6;

}

double err(double x)//Нахождение ошибки с помощью RK4(5)
{
	if (f(x) == 0)
		return 0;
	double k1 = f(x);
	double k2 = f(x + h * 1 / 4 * k1);
	double k3 = f(x + h * (3 / 32 * k1 + 9 / 32 * k2));
	double k4 = f(x + h * (1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3));
	double k5 = f(x + h * (439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4));
	double k6 = f(x + h * (-8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5));

	return  h * k1 * (25 / 216 - 16 / 135) + h * (1408 / 2565 - 6656 / 12825) * k3 + h * (2197 / 4104 - 28561 / 56430) * k4 - h * (1 / 5 - 9 / 50) * k5 - h * 2 / 55 * k6;
}
void Hopt( double E, double tol)// поиск оптимального шага
{
	
	double error = err(E);
	if (abs(error) < tol)
		return ;
	error = sqrt(pow(error / tol, 2));
	h = h * pow(1 / error, 0.2);
	Hopt( E, tol);
}
double RK45(double x,double tol)//Вложенный метод Рунге-Кутты 4(5)
{
	if (f(x) == 0)
		return 0;
	Hopt( x, tol);
	double k1 = f(x);
	double k2 = f(x + h * 1 / 4 * k1);
	double k3 = f(x + h * (3 / 32 * k1 + 9 / 32 * k2));
	double k4 = f(x + h * (1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3));
	double k5 = f(x + h * (439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4));

	return x + h * k1 * 25 / 216 + h * 1408 / 2565 * k3 + h * 2197 / 4104 * k4 - h * 1 / 5 * k5;
}

double Trapeze(double x)//Метод неявной тапеции
{
	if (f(x) ==0)
		return 0;
	return h * f(x) + h * h / 4 * -1 / (2 * f(x)) * f(x);
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
	//4777223.0 / 17418240.0		//11
};
double sup(int i, int n)//обратные конечные разности
{
	if (i == 0)
		return startF[n];
	return sup(i - 1, n) -sup(i - 1, n - 1);
}
void starter(double E)//разгон методом Рунге-Кутты 4
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
double Adams(double x)//Метод Адамса-Башфорта 11 порядка
{
	double sum = 0;
	for (int i = 0; i  <11; ++i)
	{
		sum =sum+ gamma[i] * sup(i,10);
	}
	x = x + h * sum;
	startF.erase(startF.begin(), startF.begin() + 1);
	startF.push_back(f(x));
	return x;
}

int main()
{
	//freopen("a.txt", "w", stdout);
	h = 0.01;
	double E = 50;
	double RK = E;
	double T = E;
	double rk44 = E;
	double A = E;
	
	int i = 0;
	starter(E);
	for (double t = 0; i < 11; t = t + h,i++)
	{
		cout <<t<<" "<< Target(E, t)<<" " <<Target(E,t) <<endl;
		RK = RK4(RK);
	}
	A = RK;
	for (double t = 11*h; Target(E, t)>0; t = t + h)
	{
		A = Adams(A);
		cout <<setprecision(8)<< t << " " << Target(E, t) <<" "<< A << endl;
		RK = RK4(RK);	
	}
/*	ofstream fout;
	fout.open("a.txt");
	for (double t = 0; Target(E, t) > 0 && T > 0; t += Hopt) //вывод со старым шагом
	{
		RK = RK + RK44(RK, Law, h);
		fout << t << " " << setprecision(6) << abs(RK - Target(E, t))/ Target(E, t) << endl;
	}
	fout.close();
	*/
	return 0;
}
