// Copyright Eugene Lykhovyd, 2017-2018.
#include "gurobi_c++.h"
#include "graph.h"
#include "io.h"
#include <algorithm>
#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
#include <cstring>
#include <chrono>
#include <limits>

using namespace std;
//extern const char* gitversion;

int main(int argc, char *argv[])
{
	ios::sync_with_stdio(1);
	//printf("Districting, build %s\n", gitversion);
	if (argc < 9) {
		printf("Not Enough Entries! Please Try Again!", argv[0]);
		return 0;
	}
	// parse command line arguments
	char* dimacs_fname = argv[1];
	char* hash_fname = argv[2];
	char* block_fname = argv[3]; //Get this file from bdistricting.com
	char* distance_fname = argv[4];
	char* population_fname = argv[5];
	int k = std::stoi(argv[6]); //Number of districts
	int b = std::stoi(argv[7]); //Number of blocks based on bdistricting.com
	char* solution_fname = argv[8]; //Final warm-up solution which is in .txt format


	// read inputs
	graph* g = 0;
	vector<vector<int>> dist;
	vector<int> population;
	vector<vector<long long>> assign;
	vector<long long> hash;
	if (read_input_data(dimacs_fname, hash_fname, block_fname, distance_fname, population_fname, k, b, g, dist, population, assign, hash))
		return 1; // failure


	//Clustering Districts
	vector<vector<int>> new_assign(k);

	for (int i = 0; i < g->nr_nodes; i++)
	{
		long long id = hash[i];
		vector<int> count(k,0);
		int max = 0;
		int argmax = -1;
		for (int j = 0; j < k; j++)
		{
			for (int v = 0; v < assign[j].size(); v++)
			{
				if (assign[j][v] == id)
					count[j]++;
			}
		}

		for (long j = 0; j < k; j++)
		{
			if (count[j] > max)
			{
				max = count[j];
				argmax = j;
			}
		}
		new_assign[argmax].push_back(i);
	}


	//Finding Centers
	vector<int> center_node(k, -1);
	for (int i = 0; i < k; i++)
	{	
		int Z_min = std::numeric_limits<int>::max();
		for (int j = 0; j < new_assign[i].size(); j++)
		{
			int center = new_assign[i][j];
			int z = 0;
			for (int u = 0; u < new_assign[i].size(); u++)
			{
				int v = new_assign[i][u];
				z += population[v] * dist[center][v] * dist[center][v];
			}
			if (z < Z_min)
			{
				Z_min = z;
				center_node[i] = center;
			}
		}
	}

	//Find final solution
	vector<int> solution(g->nr_nodes,-1);
	for (int i = 0; i < k; i++)
	{
		int center = center_node[i];
		for (int j = 0; j < new_assign[i].size(); j++)
		{
			int u = new_assign[i][j];
			solution[u] = center;
		}
	}

	printf_solution(solution, solution_fname);
	delete g;
	return 0;
}
