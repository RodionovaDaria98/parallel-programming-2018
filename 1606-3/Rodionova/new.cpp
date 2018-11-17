#include <stdlib.h>
#include <cstdio>
#include <mpi.h>
#include <time.h>
#include <iostream>

using namespace std;
const enum Condition { Thinking, Eating, Hungry };

int think() { return Hungry; }
int doneEating() { return Thinking; }
int eat() { return Eating; }

void print_forks(int *forks, int n) {
	for (int i = 0; i < n; i++) {
		cout << "Fork(" << i << ") " << forks[i] << "    ";
	}
	cout << endl << endl;
}

int main(int argc, char* argv[]) {
	int procRank, procSize;
	int condition;
	int *message;
	int n = 8; //number of repeat
	MPI_Status stat;
	double worktime;

	message = new int();

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &procSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
	cout << "Process " << procRank << " reporting start " << endl;
	worktime = MPI_Wtime();

	if (procRank == 0) {
		int *forks;
		cout << "Enter the number of repetitions" << endl;
		cin >> n;
		int repeat = n * (procSize - 1);
		int k = 0;
		forks = new int[procSize - 1];
		for (int i = 0; i < (procSize - 1); ++i) {
			forks[i] = 1;
		}
#ifdef _DEBUG
		print_forks(forks, procSize - 1);
		cout << "Waiter receiving orders " << procRank << endl;
#endif
		while (repeat > 0) {
			MPI_Recv(message, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
			if (stat.MPI_TAG == 0) {
				if ((forks[(stat.MPI_SOURCE - 1)] == 1) && (forks[(stat.MPI_SOURCE) % (procSize - 1)] == 1)) {
					forks[(stat.MPI_SOURCE - 1)] = 0;
					forks[stat.MPI_SOURCE % (procSize - 1)] = 0;
					MPI_Send(message, 1, MPI_INT, stat.MPI_SOURCE, 1, MPI_COMM_WORLD);
#ifdef _DEBUG
					cout << "Waiter gave forks to philosopher " << stat.MPI_SOURCE << endl;
					print_forks(forks, procSize - 1);
#endif
				}
				else {
					MPI_Send(message, 1, MPI_INT, stat.MPI_SOURCE, 0, MPI_COMM_WORLD);
#ifdef _DEBUG
					cout << "Waiter couldn't give forks to philosopher " << stat.MPI_SOURCE << endl;
					print_forks(forks, procSize - 1);
#endif
				}
			}
			else if (stat.MPI_TAG == 1) {
				//spoons recieved
				forks[stat.MPI_SOURCE - 1] = 1;
				forks[(stat.MPI_SOURCE) % (procSize - 1)] = 1;
				MPI_Send(message, 1, MPI_INT, stat.MPI_SOURCE, 2, MPI_COMM_WORLD);
				repeat--;
#ifdef _DEBUG
				cout << "Forks returned from philosopher " << stat.MPI_SOURCE << endl;
				print_forks(forks, procSize - 1);
#endif
			}
		}
#ifdef _DEBUG
		cout << "Philosophers " << " finished" << endl;
#endif
		cout << "Worktime = " << worktime << " sec" << endl;
		MPI_Finalize();
		return 0;
	}
	//oficiant end
	for (int i = 0; i < n; i++) {
#ifdef _DEBUG
		cout << "Philosopher " << procRank << " thinking/hungry" << endl;
#endif
		condition = think();
		//lock
#ifdef _DEBUG
		cout << "Philosopher " << procRank << " trying get forks" << endl;
#endif
		do {
			MPI_Send(message, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
			MPI_Recv(message, 1, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &stat);
			for (int i = 0; i < 10000; i++);
		} while (stat.MPI_TAG == 0); // while not get spoons
#ifdef _DEBUG
		cout << "Philosopher " << procRank << " got forks" << endl;
#endif
		condition = eat();
#ifdef _DEBUG
		cout << "Philosopher " << procRank << " eating, trying to return forks" << endl;
#endif
		MPI_Send(message, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
#ifdef _DEBUG
		cout << "Philosopher " << procRank << " returned forks" << endl;
#endif
		MPI_Recv(message, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &stat);
#ifdef _DEBUG
		cout << "Philosopher " << procRank << " got approvment" << endl;
#endif
		condition = doneEating();
#ifdef _DEBUG
		cout << "Philosopher " << procRank << " done eating" << endl;
#endif
	}
#ifdef _DEBUG
	cout << "Process " << procRank << " finished" << endl;
#endif
	MPI_Finalize();
	return 0;
}