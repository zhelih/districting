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

	else
	{
		cout << "ERROR: Your command is not valid.";
	}
	return EXIT_SUCCESS;
}
