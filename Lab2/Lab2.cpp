#include <iostream>
#include <cmath>
#include <vector>
using namespace std;

const double EPS1 = 10e-9;
const double EPS2 = 10e-10;
const int NIT = 110;

// Система уравнений
vector<double> systemOfEquations(const vector<double>& x) {
	vector<double> result(2);
	result[0] = sin(x[0]) - x[1] - 1.32;
	result[1] = cos(x[1]) - x[0] + 0.85;
	return result;
}

// Якобиан
vector<vector<double>> jacobian(const vector<double>& x) {
	vector<vector<double>> result(2, vector<double>(2));
	result[0][0] = cos(x[0]);
	result[0][1] = -1;
	result[1][0] = -1;
	result[1][1] = -1.0 * sin(x[1]);
	return result;
}

// Реализация метода Ньютона
vector<double> newtonMethod(const vector<double>& vX0, double& delta1, double& delta2) {
	vector<double> vX =vX0;
	int countIter = 0;
	for (int i = 0; i < NIT; i++) {
		vector<double> vf = systemOfEquations(vX);
		vector<vector<double>> vJ = jacobian(vX);

		vector<double> deltaX = { 0.0, 0.0 };
		double determinant = vJ[0][0] * vJ[1][1] - vJ[0][1] * vJ[1][0];

		// Метод крамера для решения системы линейных уравнений
		deltaX[0] = (vJ[1][1] * vf[0] - vJ[0][1] * vf[1]) / determinant;
		deltaX[1] = (-vJ[1][0] * vf[0] + vJ[0][0] * vf[1]) / determinant;

		// Xi+1 = deltaX + Xi
		vX[0] -= deltaX[0];
		vX[1] -= deltaX[1];

		//Критерий выхода
		double delta1, delta2;

		double maxF = numeric_limits<double>::min();
		for (int i = 0; i < vf.size(); i++) {
			if (vf[i] > maxF)
				maxF = vf[i];
		}
		delta1 = maxF;

		double maxValue = numeric_limits<double>::min();
		for (int i = 0; i < vX.size(); i++) {
			if (abs(vX[i]) < 1) {
				if (abs(deltaX[i]) > maxValue) {
					maxValue = abs(deltaX[i]);
				}
			}
			else {
				if (abs(deltaX[i] / vX[i]) > maxValue) {
					maxValue = abs(deltaX[i] / vX[i]);
				}
			}
		}
		delta2 = maxValue;

		countIter++;
		// Проверяем условие сходимости
		if (delta1 <= EPS1 && delta2 <= EPS2 && countIter < NIT) {
			cout << "Решение найдено на итерации : " << countIter << endl;
			return vX;
		}
	}

	cout << "Число итераций больше допустимого" << endl;
	return vX;
}

int main() {
	setlocale(LC_ALL, "Russian");
	double delta1, delta2;
	vector<double> vX0 = { 1.0, 1.0 };

	vector<double> solution = newtonMethod(vX0, delta1, delta2);
	cout << "Дельта 1 = " << delta1 << "  " << "Дельта 2 = " << delta2 << endl;
	cout << "Решение: x1 = " << solution[0] << ", x2 = " << solution[1] << endl;

	return 0;
}
