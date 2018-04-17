#include "GRBInterface.h"
#include "KGraph.h"
#include <sstream>
#include <string>
#include <ctime>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>

using namespace std;

int main(int argc, char *argv[])
{
	if (argc < 2)
	{
		cout << "ERROR: Not enough arguments.";
	}
	else if (strcmp(argv[1], "williams") == 0) // Solves minimum spanning tree
	{
		if (argc < 5)
		{
			cout << "error: not enough inputs.";
			return 0;
		}
		KGraph g1(argv[3], argv[3], argv[2]);
		KGraph g2(argv[5], argv[5], argv[4]);
		vector<vector<long>> SP;
		time_t start = clock();
		SP = solveMST(g1, g2);
		cout << "Time in seconds: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}

	else if (strcmp(argv[1], "formulation2") == 0) // Solves formulation2 of districting
	{
		if (argc < 5)
		{
			cout << "error: not enough inputs.";
			return 0;
		}
		KGraph g1(argv[3], argv[3], argv[2]);
		KGraph g2(argv[5], argv[5], argv[4]);
		long V = g1.n;
		KGraph d(argv[7], argv[7], argv[6], V);
		long K = atol (argv[8]);
		long L = atol(argv[9]);
		long U = atol(argv[10]);
		vector<vector<long>> FM2;
		time_t start = clock();
		FM2 = solveFM2(g1, g2, d, K, L, U);
		cout << "Time in seconds: " << (double)(clock() - start) / CLOCKS_PER_SEC << endl;
	}

	else
	{
		cout << "ERROR: Your command is not valid.";
	}
	return EXIT_SUCCESS;
}
