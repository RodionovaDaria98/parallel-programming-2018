
#include "mpi.h"
#include <iostream>
#include <vector>
#include <ctime>

int main(int argc, char* argv[])
{
	int ProcNum, ProcRank, RecvRank;
	int n;
    short *v1= nullptr;
	short   *v2= nullptr;
	unsigned int res1(0), res2(0), tmp(0);
	double startwtime, endwtime;
	srand(time(0));

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);

	if (ProcRank == 0) {
		//Generate vectors
		std::cin.clear();
		std::cout << "Enter the number of elements" << std::endl;
		while (!(std::cin >> n)) {
			std::cin.clear();
			while (std::cin.get() != '\n');
			std::cout << "Enter the number" << std::endl;
		}
		v1 = new short[n];
		v2 = new short[n];
		for (int i(0); i < n; i++) {
			v1[i] = 1 + rand() % 10;
			v2[i] = 1 + rand() % 10;
		}

		
		//
		startwtime = MPI_Wtime();
	}
	MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if (ProcRank != 0) {
		v1 = new short[n];
		v2 = new short[n];
	}

	MPI_Bcast(v1, n, MPI_SHORT, 0, MPI_COMM_WORLD);
	MPI_Bcast(v2, n, MPI_SHORT, 0, MPI_COMM_WORLD);

	for (int i(ProcRank); i <= n; i += ProcNum) {
		tmp += (v1[i] * v2[i]);
	}

	MPI_Reduce(&tmp, &res2, 1, MPI_UNSIGNED, MPI_SUM, 0, MPI_COMM_WORLD);

	if (ProcRank == 0) {
		endwtime = MPI_Wtime();
		std::cout << "Parallel time is " << endwtime - startwtime << " ms" << std::endl;
		std::cout << "Scalar product = " << res2 << std::endl;
	}
	delete[] v1;
	delete[] v2;
	MPI_Finalize();
	return 0;
}
