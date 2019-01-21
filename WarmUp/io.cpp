#include "graph.h"
#include "io.h"
#include <vector>
#include <cstdio>
#include "gurobi_c++.h"
#include <iostream>
#include <cstring>
#include <stdio.h>
#include <string.h>

using namespace std;

int read_input_data(const char* dimacs_fname, const char* hash_fname, const char* block_fname, const char* distance_fname, const char* population_fname, int &k, int &b, // INPUTS
	graph* &g, vector<vector<int> >& dist, vector<int>& population, vector<vector<long long>>& assign, vector<long long>& hash) // OUTPUTS
{
	// read dimacs graph
	g = from_dimacs(dimacs_fname);
	if (!g) {
		fprintf(stderr, "Failed to read dimacs graph from %s\n", dimacs_fname);
		return 1;
	}

	// read distances (must be sorted)
	FILE* f = fopen(distance_fname, "r");
	if (!f) {
		fprintf(stderr, "Failed to open %s\n", distance_fname);
		return 1;
	}
	// file contains the first row as a header row, skip it
	// also skip the first element in each row (node id)
	char buf[3000]; //dummy
	fgets(buf, sizeof(buf), f); // skip first line
	dist.resize(g->nr_nodes);
	for (int i = 0; i < g->nr_nodes; ++i)
	{
		dist[i].resize(g->nr_nodes);
		int d; fscanf(f, "%d,", &d); // skip first element
		for (int j = 0; j < g->nr_nodes; ++j) {
			fscanf(f, "%d,", &d);
			dist[i][j] = d;
		}
	}
	fclose(f);

	// read population file
	f = fopen(population_fname, "r");
	if (!f) {
		fprintf(stderr, "Failed to open %s\n", population_fname);
		return 1;
	}
	// skip first line about total population
	fgets(buf, sizeof(buf), f);
	population.resize(g->nr_nodes);
	for (int i = 0; i < g->nr_nodes; ++i) {
		int node, pop;
		fscanf(f, "%d %d ", &node, &pop);
		population[node] = pop;
	}
	fclose(f);

	

	// read block file that is downloaded from bdistricting.com
	f = fopen(block_fname, "r");
	if (!f) {
		fprintf(stderr, "Failed to open %s\n", block_fname);
		return 1;
	}
	assign.resize(k);
	
	for (int i = 0; i < b; ++i)
	{	
		long long block; fscanf(f, "%lld,", &block); // read the block
		int district; fscanf(f, "%d,", &district); // read the district
		string new_block = to_string(block);
		string st = new_block.substr(0, new_block.size() - 4);
		assign[district].push_back(stoll(st));
	}
	fclose(f);

	// read hash file
	f = fopen(hash_fname, "r");
	if (!f) {
		fprintf(stderr, "Failed to open %s\n", hash_fname);
		return 1;
	}

	hash.resize(g->nr_nodes);
	for (int i = 0; i < g->nr_nodes; ++i) {
		int node;
		long long id;
		fscanf(f, "%d %lld ", &node, &id);
		
		string new_id = to_string(id);
		if (new_id[0] == '0')
		{
			new_id.erase(0,1);
		}
		hash[node] = stoll(new_id);
	}
	fclose(f);

	return 0;
}

// construct districts from hess variables
void translate_solution(GRBVar** x, vector<int>& sol, int n)
{
	// translate the solution
	sol.resize(n);

	vector<int> heads(n, 0);
	int cur = 1;
	// firstly assign district number for clusterheads
	for (int i = 0; i < n; ++i)
		if (x[i][i].get(GRB_DoubleAttr_X) > 0.5)
			heads[i] = cur++;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j)
			if (x[i][j].get(GRB_DoubleAttr_X) > 0.5)
				sol[i] = heads[j];

}

// prints the solution <node> <district>
void printf_solution(const vector<int>& sol, const char* fname)
{
	// write solution to file
	if (sol.size() > 0)
	{
		FILE* f;
		if (fname)
			f = fopen(fname, "w");
		else
			f = stderr;

		for (int i = 0; i < sol.size(); ++i)
			fprintf(f, "%d %d\n", i, sol[i]);

		if (fname)
			fclose(f);
	}
}
