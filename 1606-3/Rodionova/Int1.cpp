#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

using namespace std;
using namespace std::chrono;

// x^2+y-z
double f(double *X) {
	double x = X[0];
	double y = X[1];
	double z = X[2];

	//return x + y * sin(z);
	return 2 * x + y - 4 * z;
}

double calculate(double *x, double *a, double *b, int N, int i, int n) {
	if (i >= n)
		return f(x);

	double s = 0;
	double h = (b[i] - a[i]) / N;

	for (int j = 0; j < N; j++) {
		x[i] = a[i] + h / 2 + j * h;

		s += calculate(x, a, b, N, i + 1, n);
	}

	return s;
}

// N - разбиения, n - размерность функции
double integrate(double *A, double *B, int n, int N) {
	double s = 0;

	double *x = new double[n];

	s = calculate(x, A, B, N, 0, n);

	double dh = 1;

	for (int i = 0; i < n; i++)
		dh *= (B[i] - A[i]) / N;

	return s * dh;
}

int main() {
	const int n = 3; // размерность интеграла

	double a[n];
	double b[n];
	int N;

	// считываем параметры интегрирования
	for (int i = 0; i < n; i++) {
		cout << "Enter a" << (i + 1) << ", b" << (i + 1) << ": ";
		cin >> a[i] >> b[i];
	}

	cout << "Enter N: ";
	cin >> N; // считываем количество разбиений

	high_resolution_clock::time_point t0 = high_resolution_clock::now(); // начинаем замеры времени
	double I = integrate(a, b, n, N); // вычисляем интеграл
	high_resolution_clock::time_point t1 = high_resolution_clock::now(); // начинаем замеры времени

	auto dt = duration_cast<chrono::microseconds>(t1 - t0); // получаем время выполнения в микросекундах

	cout << "Integrate for " << n << "-space;" << endl;
	cout << "Ranges: " << endl;

	// выводим отрезки интегрирования
	for (int i = 0; i < n; i++)
		cout << "  [" << a[i] << ", " << b[i] << "]" << endl;

	cout << endl;
	cout << "I: " << setprecision(15) << I << endl; // результат
	cout << "time: " << (dt.count() / 1000.0)  << endl; // время
	system("pause");
}