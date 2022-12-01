#include <iostream>
#include <math.h>
#include <iomanip>
#include <stdio.h>
#include <vector>
#include <fstream>
using namespace std;

double Law(double x)
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
double RK44(double x, double (*f)(double), double h)//Вложенный метод Рунге-Кутты 4(5)
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
double err(double x, double (*f)(double), double h,double tol)//Нахождение ошибки с помощью RK4(5)
{
	if (f(x) == 0)
		return 0;
	double k1 = f(x);
	double k2 = f(x + h * 1 / 4 * k1);
	double k3 = f(x + h * (3 / 32 * k1 + 9 / 32 * k2));
	double k4 = f(x + h * (1932 / 2197 * k1 - 7200 / 2197 * k2 + 7296 / 2197 * k3));
	double k5 = f(x + h * (439 / 216 * k1 - 8 * k2 + 3680 / 513 * k3 - 845 / 4104 * k4));
	double k6 = f(x + h * (-8 / 27 * k1 + 2 * k2 - 3544 / 2565 * k3 + 1859 / 4104 * k4 - 11 / 40 * k5));

	return  pow(h * k1 * (25 / 216 - 16 / 135)/tol + h * (1408 / 2565 - 6656 / 12825) * k3/tol + h * (2197 / 4104 - 28561 / 56430) * k4/tol - h * (1 / 5 - 9 / 50) * k5/tol - h * 2 / 55 * k6/tol,2);
}
double Trapeze(double x, double (*f)(double), double h)//Метод неявной тапеции
{
	if (f(x) ==0)
		return 0;
	return h * f(x) + h * h / 4 * -1 / (2 * f(x)) * f(x);
}

vector<double>startF;
vector<double>startX;
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
double sup(int i, int n)//конечные разности
{
	if (i == 0)
		return start[n];
	return sup(i - 1, n) - sup(i - 1, n - 1);
}
void starter(double E, double h)//разгон методом неявной трапеции
{
	double T = E;
	double time = 0;
	for (int i = 0; i < 11; ++i)
	{
		startF.push_back(T);
		startX.push_back(time);
		
		T=T+Trapeze(T,Law,h);
		time += h;
	}
}
void Adams(double h)//Метод Адамса-Башфорта 11 порядка
{
	double sum = 0;
	int siz = startX.size() - 1;
	for (int i = 0; i < 11; ++i)
	{
		sum =sum+ gamma[i] * sup(i,siz);
	}
	startX.push_back(startX[siz] +  h*sum);
	startF.push_back(startF[siz]+Law(startX[siz] + h * sum));
}



double Sopt( double h,double E,double tol)// поиск оптимального шага
{
	double RK = E;
	double error = 0;
	int m = 0;
	double Hopt = 0;
	for (double t = 0; Target(E, t) > 0 && RK > 0; t += h) 
	{
		RK = RK + RK44(RK, Law, h);
		error += err(RK, Law, h, tol);
		++m;
	}
	error = sqrt(error / m);
	return h * pow((1 / error), 1.0 / 5.0);
}

int main()
{
	//freopen("a.txt", "w", stdout);
	double h = 0.1;
	double E = 50;
	double RK = E;
	double T = E;
	double Hopt = h;

	/*ofstream fout;

	fout.open("a.txt");
	for (double t = 0; Target(E, t) > 0 && T > 0; t += Hopt) //вывод со старым шагом
	{
		RK = RK + RK44(RK, Law, h);
		fout << t << " " << setprecision(6) << abs(RK - Target(E, t)) << endl;
	}
	fout.close();*/

	for (int i = 0; i < 100; ++i) //????????С каждой итерациией уточняется оптимальный шаг
	{
		Hopt = Sopt(Hopt, E, 0.001);
	}
	/*RK = E;

	fout.open("b.txt");
	for (double t = 0; Target(E, t) > 0 && T > 0; t += Hopt) //вывод с новым шагом
	{
		RK = RK + RK44(RK, Law, Hopt);
		fout << t << " " << setprecision(6) << abs(RK - Target(E, t))<< endl;
	}
	fout.close();*/
	starter(E, h);
	for(int i=0;i<100;++i)
		Adams(h);
	double t = 0;
	for (int i = 0; i < startF.size(); ++i,t+=h)
		cout <<t<<" "<< startF[i] << endl;
	return 0;
}
