
#include <iostream>
#include <math.h>
#include <omp.h>

#define Eps 0.01



double val_one_dimensional_func(double x);
double val_two_dimensonal_func(double x1, double x2);

void one_dimensional_integral(double a, double b, double h, double *res);
void one_dim_integral_parallel(double a, double b, double h, double *res);

void two_dimensional_integral(double a1, double b1, double a2, double b2, double h, double *res);
void two_dim_integral_parallel(double a1, double b1, double a2, double b2, double h, double *res);
void two_dimensional_integral_fix_point(double a1, double a2, double b2, double h, double *res);

void two_dimensional_integral_multi_step(double a1, double b1, double a2, double b2, double h, double *res);
void two_dim_integral_parallel_multi_step(double a1, double b1, double a2, double b2, double h, double *res);


double val_one_dimensional_func(double x) {
	return sin(x);
}

double val_two_dimensonal_func(double x1, double x2) {
	return sin(x1) * cos(x2);
}

void one_dimensional_integral(double a, double b, double h, double *res) {
	double sum = 0.0;
	for (double i = a; i < b; i += h) {
		sum += val_one_dimensional_func(i) * h;
	}
	*res = sum;
}
void one_dim_integral_parallel(double a, double b, double h, double *res) {
	double total_res = 0.0;
	int i;

#pragma omp parallel for private(i) reduction(+: total_res)
	for (i = 0; i < int((b - a) / h); i++) {
		total_res += val_one_dimensional_func(a + i * h) * h;
	}

	*res = total_res;
}

void two_dim_integral_parallel(double a1, double b1, double a2, double b2, double h, double *res) {
	double sum = 0.0;
	int i, j;
#pragma omp parallel
	{ 
		
#pragma omp for private(i, j) reduction(+: sum)
		for (i = 0; i < int((b1 - a1) / h); i++)
		{
			for (j = 0; j < int((b2 - a2) / h); j++)
			{
				sum += val_two_dimensonal_func(a1 + i * h, a2 + j * h) * h * h;
			}
		}
	}
	*res = sum;
}
void two_dimensional_integral(double a1, double b1, double a2, double b2, double h, double *res) {
	double sum = 0.0;

	for (double i = a1; i < b1; i += h)
	{
		for (double j = a2; j < b2; j += h)
		{
			sum += val_two_dimensonal_func(i, j) * h * h;
		}
	}
	*res = sum;
}

void two_dimensional_integral_fix_point(double a1, double a2, double b2, double h, double *res) {
	double sum = 0.0;

	for (double j = a2; j < b2; j += h)
	{
		sum += val_two_dimensonal_func(a1, j) * h;
	}

	*res = sum;
}

void two_dim_integral_parallel_multi_step(double a1, double b1, double a2, double b2, double h, double *res, int p)
{
	double sum = 0.0;
	double tmpsum = 0.0;
	int i, j;
#pragma omp parallel num_threads(p)
	{
	
			
#pragma omp for private(i, j, tmpsum) reduction(+: sum)
		for (i = 0; i < int((b1 - a1) / h); i++)
		{
			two_dimensional_integral_fix_point(a1 + h * i, a2, b2, h, &tmpsum);
			sum += tmpsum * h;
		}
	}
	*res = sum;
}

void two_dimensional_integral_multi_step(double a1, double b1, double a2, double b2, double h, double *res) {
	double sum = 0.0;
	double tmpsum = 0.0;
	int i;
	for (i = 0; i < int((b1 - a1) / h); i++)
	{
		two_dimensional_integral_fix_point(a1 + h * i, a2, b2, h, &tmpsum);
		sum += tmpsum * h;
	}
	*res = sum;
}

void two_dim_integral_parallel_multi_step(double a1, double b1, double a2, double b2, double h, double * res)
{
}





int main(int argc, char **argv)
{

	double a1 = 0.0, b1 = 5.0, a2 = 0.0, b2 = 5.0;
	double h = 0.001;
	double seqRes = 0.0, parRes = 0.0;
	double t1 = 0.0, t2 = 0.0, parallTimeWork = 0.0, sequentTimeWork = 0.0;



	int p = 4;
	std::cout << "Set a number of threads" << std::endl;
	std:: cin >> p;
	omp_set_num_threads(p);




	for (int i = 1; i < argc; ++i) {
		if ((!strcmp(argv[i], "-a1")) && (i + 1 < argc)) {
			a1 = atof(argv[i + 1]);
		}
		if ((!strcmp(argv[i], "-b1")) && (i + 1 < argc)) {
			b1 = atof(argv[i + 1]);
		}
		if ((!strcmp(argv[i], "-h")) && (i + 1 < argc)) {
			h = atof(argv[i + 1]);
		}
		if ((!strcmp(argv[i], "-a2")) && (i + 1 < argc)) {
			a2 = atof(argv[i + 1]);
		}
		if ((!strcmp(argv[i], "-b2")) && (i + 1 < argc)) {
			b2 = atof(argv[i + 1]);
		}
	}

	if (a1 == INFINITY || b2 == INFINITY) {
		
		t1 = omp_get_wtime();
		one_dimensional_integral(a1, b1, h, &seqRes);
		t2 = omp_get_wtime();
		sequentTimeWork = (t2 - t1) * 1000.;
		std::cout << "Result of sequential algorithm = " << seqRes << std::endl;
		std::cout << "sequentTimeWork = " << sequentTimeWork << " ms" << std::endl;

		t1 = omp_get_wtime();
		one_dim_integral_parallel(a1, b1, h, &parRes);
		t2 = omp_get_wtime();
		parallTimeWork = (t2 - t1) * 1000.;
		std::cout << "Result of parallel algorithm = " << parRes << std::endl;
		std::cout << "parallelTimeWork = " << parallTimeWork << " ms" << std::endl;

	}
	else {

		t1 = omp_get_wtime();
		two_dimensional_integral_multi_step(a1, b1, a2, b2, h, &seqRes);
		t2 = omp_get_wtime();
		sequentTimeWork = (t2 - t1) * 1000.;
		std::cout << "Result of sequential algorithm = " << seqRes << std::endl;
		std::cout << "sequentTimeWork = " << sequentTimeWork << " ms" << std::endl;

		t1 = omp_get_wtime();
		two_dim_integral_parallel_multi_step(a1, b1, a2, b2, h, &parRes);
		t2 = omp_get_wtime();
		parallTimeWork = (t2 - t1) * 1000.;
		std::cout << "Result of parallel algorithm = " << parRes << std::endl;
		std::cout << "parallelTimeWork = " << parallTimeWork << " ms" << std::endl;

		
	}

	system("pause");
}