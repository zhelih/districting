#ifndef KGRAPH_H
#define KGRAPH_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

extern double copy_time;
using namespace std;

// Much of the KGraph class written by Anurag Verma (student of Sergiy Butenko).
// Modified by Austin Buchanan, and then Hamidreza Validi.

class KGraph
{
	public:
		
		// start new functions by Austin	
			bool IsConnected(vector<long> S);
			bool IsConnected(vector<bool> S);
			bool IsConnected();
		
			long DiameterUnweighted();
			long DiameterWeighted(vector<long> W);
			long DiameterUnweighted(vector<long> S);
			long DiameterUnweighted(vector<bool> S1);

			long LongestShortestPathUnweighted(long origin);
			vector<long> ShortestPathsUnweighted(long origin);
			vector<long> ShortestPathsNodeWeighted(long origin, vector<long> W);
			long LongestShortestPathNodeWeighted(long origin, vector<long> W);
			vector<long> KGraph::ShortestPathsUnweighted(long origin, vector<bool> &S);

			vector< vector<bool> > CreateAdjacenyMatrix();
			bool DeleteNode(long i);  // deletes incident edges (and not the node itself)
			void ComplementGraph(const KGraph &rhs);
		// end new functions by Austin

		vector<long> *adj;  // stores the adj lists. Each list is maintained as a sorted vector.
		vector<long> *pdadj;
		long *degree;       // stores the degree seq
		long n;             // num of nodes
		long m;             // num of edges
		long Delta;         // higest degree. As of now, is instanciated when graph is read or copied, not maintained afterwards
		string name;        // name of the graph. Could be anything.

		/* Constructors */
		KGraph();
		KGraph(string nm);
		KGraph(long n);		
		KGraph(string nm, string file, string type);
		KGraph(const KGraph &rhs);

		/* File IO utility functions */
		void ReadDIMACSGraph(string file);
		void ReadDIMACSColorGraph(string file);
		void ReadDIMACSGraphParallel(string file);
		void ReadDATGraph(string file);
		void ReadSNAPGraph(string file);
		void ReadPrimalGraph(string file);
		void ReadDualGraph(string file);
		bool WriteGVizGraph(vector <bool> F, string file);
		bool WriteGVizGraph(string file);
		void WriteSNAPGraph(string outfile);
		void WriteDIMACSGraph(string outfile);

		/* General purpose utility functions */
		bool CheckValid();
		void Duplicate(const KGraph &rhs);
		void DuplicateConnected(const KGraph &rhs,  map<long, long> &node_map);
		bool AddEdge(long i, long j, bool reverseToo, bool safe);
		bool DeleteEdge(long i, long j, bool reverseToo, bool safe);
		long ConnectedVertices(); // number of non-isolated nodes

		/* common neighbor functions */
		long CommonNeighbors(long u, vector<long> &v);
		long CommonNeighbors(long u, long v);
		bool HasKCommonNeighbors(long u, long v, long k);
		vector<long> CommonNeighborsList(long u, long v);
		vector<long> CommonNeighborsList(long u, vector<long> &v);		
		
		/* functions for induced subgraphs */
		void FindInducedGraph(vector<bool> &S);
		KGraph CreateInducedGraph(vector<long> &S);
		KGraph CreateInducedGraph(vector<long> &S, vector<long> &ReverseMap);
		void FindConnectedComponents(vector< vector< long> > &components, vector<long> &degreeZero);

		/* Destructors */
		void clear();
		~KGraph();		
		void ReadIEEEGraph(string file);
};

#endif
