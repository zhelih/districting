#include "KGraph.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <queue>
#include <set>
#include <map>
#include <omp.h>
# include <string.h>
using namespace std;


bool KGraph::IsConnected(vector<bool> S)
{
	/* Is the subgraph induced by S connected? */
	long u, v;
	vector<bool> reached(n, false);
	vector<long> children, parents;

	for (long i = 0; i < n; i++) //find a root node
		if (S[i]) {
			children.push_back(i);
			reached[i] = true;
			break;
		}

	for (; !children.empty(); ) { //do BFS
		parents = children;
		children.clear();
		for (long i = 0; i < parents.size(); i++) { //for each parent, examine the children
			u = parents[i];
			for (long j = 0; j < degree[u]; j++) {
				v = adj[u][j];
				if (S[v] && !reached[v]) { //can only use vertices in S
					reached[v] = true;
					children.push_back(v);
				}
			}
		}
	}
	for (long i = 0; i < n; i++) //if a vertex in S hasn't been reached, return false
		if (S[i] && !reached[i])
			return false;

	return true;
}
bool KGraph::IsConnected(vector<long> S1)
{
	vector<bool> S(n, false); //convert S1 to bool
	for (long i = 0; i < S1.size(); i++)
		S[S1[i]] = true;
	return IsConnected(S);
}

bool KGraph::IsConnected()
{
	//Is G connected?
	vector<bool> S(n, true);
	return IsConnected(S);
}

long KGraph::DiameterUnweighted()
{
	/* Returns the (unweighted) diameter of a connected graph.
	If the graph is not connected it returns n.
	Solves a series of shortest paths problems. */
	vector<long> ShortestPaths = ShortestPathsUnweighted(0); //check to see if the graph is connected
	for (long i = 0; i < n; i++)
		if (ShortestPaths[i] == n)
			return n;

	long diameter = 0, temp_longest;
	for (long i = 0; i < n; i++) { //solve shortest paths problem, originating from node i
		temp_longest = LongestShortestPathUnweighted(i);
		if (temp_longest > diameter)
			diameter = temp_longest;
	}
	return diameter;
}

long KGraph::DiameterWeighted(vector<long> W) {
	vector<long> ShortestPaths = ShortestPathsUnweighted(0); //check to see if the graph is connected
	if (*max_element(ShortestPaths.begin(), ShortestPaths.end()) == n) return (long)LONG_MAX;

	long diameter = -1, temp_longest;
	for (long i = 0; i < n; i++) { //solve shortest paths problem, originating from node i
		temp_longest = LongestShortestPathNodeWeighted(i, W);
		if (temp_longest > diameter)
			diameter = temp_longest;
	}
	return diameter;
}

long KGraph::DiameterUnweighted(vector<long> S)
{
	/* returns the diameter of the graph induced by S */
	KGraph g = CreateInducedGraph(S);
	return g.DiameterUnweighted();
}

long KGraph::DiameterUnweighted(vector<bool> S1)
{
	vector<long> S;
	for (long i = 0; i < n; i++)
		if (S1[i])
			S.push_back(i);
	return DiameterUnweighted(S);
}

vector<long> KGraph::ShortestPathsUnweighted(long origin)
{
	vector<bool> S(n, true);
	return ShortestPathsUnweighted(origin, S);
}


vector<long> KGraph::ShortestPathsUnweighted(long origin, vector<bool> &S)
{
	/*Finds the shortest paths from node v to all other nodes in graph G[S].
	Assumes the graph is connected.
	Performs BFS.*/
	long u, v;
	vector<long> dist(n, n); //shortest distance from origin node to each other node. dist[i] = n means i not reachable
	if (!S[origin]) return dist;  // if origin not in S, return infinities.
	vector<bool> reached(n, false);
	vector<long> children, parents;

	children.push_back(origin);
	dist[origin] = 0; //the origin node is distance 0 from itself
	reached[origin] = true;

	for (long d = 1; !children.empty(); d++) { //for each distance
		parents = children;
		children.clear();
		for (long i = 0; i < parents.size(); i++) { //for each parent, examine the children
			u = parents[i];
			for (long j = 0; j < degree[u]; j++) {
				v = adj[u][j];
				if (!reached[v] && S[v]) {
					reached[v] = true;
					dist[v] = d;
					children.push_back(v);
				}
			}
		}
	}
	return dist;
}

//vector<long> KGraph::ShortestPathsUnweighted(long origin)
//{ 
//	/*Finds the shortest paths from node v to all other nodes.
//	Assumes the graph is connected. 
//	Performs BFS.*/
//	long u, v;
//	vector<long> dist(n,n); //shortest distance from origin node to each other node. dist[i] = n means i not reachable
//	vector<bool> reached(n,false);
//	vector<long> children,parents;
//	
//	children.push_back(origin);
//	dist[origin] = 0; //the origin node is distance 0 from itself
//	reached[origin] = true;
//
//	for(long d = 1; !children.empty(); d++) { //for each distance
//		parents = children;
//		children.clear();
//		for(long i=0; i<parents.size(); i++) { //for each parent, examine the children
//			u = parents[i];
//			for(long j=0; j<degree[u]; j++)	{
//				v = adj[u][j];
//				if(!reached[v])	{
//					reached[v] = true;
//					dist[v] = d;
//					children.push_back(v);
//				}
//			}
//		}
//	}
//	return dist;
//}
long KGraph::LongestShortestPathUnweighted(long origin)
{
	vector<long> SP = ShortestPathsUnweighted(origin);
	return (long)*max_element(SP.begin(), SP.end());
}

vector<long> KGraph::ShortestPathsNodeWeighted(long origin, vector<long> W) {

	long v;
	long infty = long(floor((double)0.4*LONG_MAX));
	vector<long> dist(n, infty); //shortest distance from origin node to each other node. dist[i] = \infty means i not reachable
	vector<bool> Q(n, true); // the set of vertices whose distance labels are not yet permanent
	//W[origin] = 0;   // for convenience, since we are not counting the weight of the first vertex in the path as part of its length
	dist[origin] = 0; //the origin node is distance 0 from itself
	long minDistance;
	long minVertex;

	// do node-weighted version of Dijkstra's shortest path algorithm.
	for (long i = 0; i < n; i++)
	{
		// find a vertex u from Q of minimum distance label dist[u]
		minDistance = (long)LONG_MAX;
		minVertex = -1;
		for (long v = 0; v < n; v++)
		{
			if (!Q[v]) continue;
			if (dist[v] < minDistance)
			{
				minDistance = dist[v];
				minVertex = v;
			}
		}

		// remove minVertex from Q
		Q[minVertex] = false;

		// update distance labels of neighbors
		for (long j = 0; j < degree[minVertex]; j++)
		{
			v = adj[minVertex][j];
			if (Q[v] && dist[minVertex] + W[minVertex] < dist[v])
			{
				dist[v] = dist[minVertex] + W[minVertex];
			}
		}
	}
	return dist;
}

long KGraph::LongestShortestPathNodeWeighted(long origin, vector<long> W)
{
	vector<long> SP = ShortestPathsNodeWeighted(origin, W);
	return (long)*max_element(SP.begin(), SP.end());
}

vector< vector<bool> > KGraph::CreateAdjacenyMatrix()
{
	/* Creates the associated adjacency matrix. Not recommended for large graphs. */
	long v;
	vector<bool> tempAdjZeros(n, false);
	vector< vector<bool> > adjMatrix(n, tempAdjZeros);
	for (long i = 0; i < n; i++)
	{
		for (long j = 0; j < degree[i]; j++)
		{
			v = adj[i][j];
			adjMatrix[i][v] = true;
		}
	}
	return adjMatrix;
}

bool KGraph::DeleteNode(long i)
{
	/* Deletes all edges which are incident to node i.
	This does NOT actually remove the node per se.
	The value g.n remains the same. */
	for (long j = 0; j < degree[i]; j++)
	{
		long v = adj[i][j];
		DeleteEdge(v, i, false, false);
	}
	m -= degree[i];
	adj[i].clear();
	degree[i] = 0;
	return true;
}


void KGraph::ComplementGraph(const KGraph &rhs)
{
	/* Makes the graph the complement of rhs graph.*/
	if (n > 0)
		clear();
	n = rhs.n;
	m = n*(n - 1) / 2 - rhs.m;
	name = rhs.name;
	degree = new long[n];
	adj = new vector<long>[n];
	vector<long>::iterator it; //position within adjacency list
	for (long i = 0; i < n; i++)
	{
		degree[i] = n - rhs.degree[i] - 1;	//update degrees
		if (degree[i] == 0) continue;
		it = rhs.adj[i].begin();
		for (long j = 0; j < n; j++) //add edges.
		{
			if (j == i) continue; //no self-loops
			if (it == rhs.adj[i].end())
				adj[i].push_back(j);
			else if (*it == j)
				it++;
			else adj[i].push_back(j);
		}
	}
}

/* Default constructor. Does nothing but initializing n and m. */
KGraph::KGraph()
{
	n = 0;
	m = 0;
	Delta = 0;
}

/* Default constructor, but names the graph in addition. */
KGraph::KGraph(string nm)
{
	name = nm;
	n = 0;
	m = 0;
	Delta = 0;
}

/* Useful constructor. Names the graph, and reads it from a files. */
KGraph::KGraph(string nm, string file, string type)
{
	n = 0;
	m = 0;
	Delta = 0;
	name = nm;
	if (type == "dimacs")
		ReadDIMACSGraph(file);
	else if (type == "snap_d")
		ReadSNAPGraph(file);
	else if (type == "primal")
		ReadPrimalGraph(file);
	else if (type == "dual")
		ReadDualGraph(file);
	else if (type == "dimacs_color")
		ReadDIMACSColorGraph(file);
	else if (type == "DAT")
		ReadDATGraph(file);
	else if (type == "ieee")
		ReadIEEEGraph(file);
	else
		cerr << "Format " << type << " not found.\n";
}

/* Default constructor, but initializes the number of nodes and the vectors in addition */
KGraph::KGraph(long nodes)
{
	n = nodes;
	m = 0;
	Delta = 0;
	degree = new long[n];
	memset(degree, 0, sizeof(long)*n);
	adj = new vector<long>[n];
}

/* Copy constructor */
KGraph::KGraph(const KGraph &rhs)
{
	Duplicate(rhs);
}

/* Copying function. Makes the calling graph same as the passed graph. */
void KGraph::Duplicate(const KGraph &rhs)
{
	if (n > 0)
		clear();
	n = rhs.n;
	m = rhs.m;
	name = rhs.name;
	Delta = rhs.Delta;
	degree = new long[n];
	adj = new vector<long>[n];
	long i = 0;
	//#pragma omp parallel for	
	for (i = 0; i < n; i++)
	{
		degree[i] = rhs.degree[i];
		adj[i] = rhs.adj[i];
	}
}

/* Copying function. Makes the calling graph same as the passed graph, but removes isolated vertices
   and renumbers the remaining vertices, making n smaller. */
void KGraph::DuplicateConnected(const KGraph &rhs, map<long, long> &node_map)
{
	if (n > 0)
		clear();

	map<long, long> node;
	for (long i = 0; i < rhs.n; i++)
		if (rhs.degree[i] > 0)
		{
			node[i] = n;
			node_map[n] = i;
			n++;
		}
	m = rhs.m;
	name = rhs.name;
	degree = new long[n];
	adj = new vector<long>[n];


	for (long i = 0; i < rhs.n; i++)
	{
		long t = 0;
		if (rhs.degree[i] <= 0)
			continue;
		t = node.find(i)->second;
		degree[t] = rhs.degree[i];
		adj[t].resize(rhs.degree[i]);
		for (long s = 0; s < rhs.degree[i]; s++)
		{
			adj[t][s] = node.find(rhs.adj[i][s])->second;
		}
	}
}

long KGraph::ConnectedVertices()
{
	/* Returns the number of nodes which have positive degree. */
	long connectedVertices = 0;
	for (long i = 0; i < n; i++)
		if (degree[i] > 0)
			connectedVertices++;
	return connectedVertices;
}

/* reverseToo: if false, j will be added to i's adj, but not the other way round.
 * safe: if true, a check will be performed to make sure the edge does not already exists. */
bool KGraph::AddEdge(long i, long j, bool reverseToo, bool safe)
{
	vector<long>::iterator it;
	if (degree[i] == 0 || j > adj[i][degree[i] - 1])
		adj[i].push_back(j);
	else
	{
		it = lower_bound(adj[i].begin(), adj[i].end(), j);
		if (!safe)
			adj[i].insert(it, j);
		else
		{
			if (*it != j)
				adj[i].insert(it, j);
			else
				return false;
		}
	}
	degree[i]++;

	if (reverseToo)
	{
		if (degree[j] == 0)
			adj[j].push_back(i);
		else
		{
			it = lower_bound(adj[j].begin(), adj[j].end(), i);
			if (!safe)
				adj[j].insert(it, i);
			else
			{
				if (*it != i)
					adj[j].insert(it, i);
				else return false;
			}
		}
		degree[j]++;
	}
	return true;
}

/* reverseToo: if false, j will be removed from i's adj, but not the other way round.
 * safe: if true, a check will be performed to make sure the edge exists. */
bool KGraph::DeleteEdge(long i, long j, bool reverseToo, bool safe)
{
	vector<long>::iterator it = lower_bound(adj[i].begin(), adj[i].end(), j);
	if (!safe)
		adj[i].erase(it);
	else
	{
		if (it != adj[i].end() && *it == j)
			adj[i].erase(it);
		else return false;
	}
	degree[i]--;
	if (reverseToo)
	{
		it = lower_bound(adj[j].begin(), adj[j].end(), i);
		if (!safe)
			adj[j].erase(it);
		else
		{
			if (it != adj[j].end() && *it == i)
				adj[j].erase(it);
			else return false;
		}
		degree[j]--;
	}
	return true;
}

bool KGraph::CheckValid()
{
	long m1 = 0;
	for (long i = 0; i < n; i++)
		m1 += degree[i];
	if (m1 != 2 * m)
	{
		cerr << "ERROR: " << m1 << " != " << 2 * m << endl;
		return false;
	}
	for (long i = 0; i < n; i++)
	{
		for (long j = 0; j < degree[i] - 1; j++)
		{
			if (adj[i][j] >= adj[i][j + 1])
			{
				cerr << "ERROR: " << "adj[i][j]=" << adj[i][j] << " >= " << adj[i][j + 1] << "adj[i][j+1]" << endl;
				return false;
			}
		}
	}
	return true;
}

/* Returns the number of common neighbors node u has with v */
long KGraph::CommonNeighbors(long u, long v)
{
	long t = 0;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = degree[v];
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &adj[v].front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1 < q1end && q < qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t++;
			q++;
			q1++;
		}
		else if (t1 > t2)
			q1++;
		else
			q++;
	}
	return t;
}

/* Returns the number of common neighbors node u has with v */
bool KGraph::HasKCommonNeighbors(long u, long v, long k)
{
	long t = 0;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = degree[v];
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &adj[v].front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1 < q1end && q < qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t++;
			if (t >= k)
				return true;
			q++;
			q1++;
		}
		else if (t1 > t2)
			q1++;
		else
			q++;
	}
	if (t >= k)
		return true;
	return false;
}

/* Returns the sorted list of common neighbors node u has with v */
long KGraph::CommonNeighbors(long u, vector<long> &v)
{
	long t = 0;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = v.size();
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &v.front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1 < q1end && q < qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t++;
			q++;
			q1++;
		}
		else if (t1 > t2)
			q1++;
		else
			q++;
	}
	return t;
}

/* Returns the sorted list of common neighbors node u has with v */
vector<long> KGraph::CommonNeighborsList(long u, long v)
{
	vector<long> t;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = degree[v];
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &adj[v].front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1 < q1end && q < qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t.push_back(t1);
			q++;
			q1++;
		}
		else if (t1 > t2)
			q1++;
		else
			q++;
	}
	return t;
}

/* Returns the sorted list of common neighbors node u has with v */
vector<long> KGraph::CommonNeighborsList(long u, vector<long> &v)
{
	vector<long> t;
	long q = 0, q1 = 0;
	long qend = degree[u], q1end = v.size();
	long t1, t2;
	const long* ptr1 = (q1end != 0) ? &v.front() : NULL;
	const long* ptr = (qend != 0) ? &adj[u].front() : NULL;

	while (q1 < q1end && q < qend)
	{
		t1 = ptr[q];
		t2 = ptr1[q1];
		if (t1 == t2)
		{
			t.push_back(t1);
			q++;
			q1++;
		}
		else if (t1 > t2)
			q1++;
		else
			q++;
	}
	return t;
}

/* Finds the subgraph induced by all the nodes in S. S is sorted */
void KGraph::FindInducedGraph(vector<bool> &S)
{
	for (long i = 0; i < n; i++)
		if (!S[i])
			DeleteNode(i);
}

KGraph KGraph::CreateInducedGraph(vector<long> &S, vector<long> &ReverseMap)
{
	/* Finds the subgraph induced by all the nodes in S. S is sorted */
	unsigned long S_size = S.size();
	KGraph g(S_size);
	ReverseMap.resize(n, -1);
	for (long i = 0; i < S_size; i++)
		ReverseMap[S[i]] = i;
	for (unsigned long i = 0; i < S_size; i++)
	{
		g.adj[i] = CommonNeighborsList(S[i], S);
		g.degree[i] = g.adj[i].size();
		for (long j = 0; j < g.degree[i]; j++) //relabel the vertices for the new, smaller graph
			g.adj[i][j] = ReverseMap[g.adj[i][j]];
		g.m += g.degree[i];
	}
	g.m /= 2;
	return g;
}


KGraph KGraph::CreateInducedGraph(vector<long> &S)
{
	/* Finds the subgraph induced by all the nodes in S. S is sorted */
	unsigned long S_size = S.size();
	KGraph g(S_size);
	vector<long> R(n, -1); //inverse sorting
	for (long i = 0; i < S_size; i++)
		R[S[i]] = i;
	for (unsigned long i = 0; i < S_size; i++)
	{
		g.adj[i] = CommonNeighborsList(S[i], S);
		g.degree[i] = g.adj[i].size();
		for (long j = 0; j < g.degree[i]; j++) //relabel the vertices for the new, smaller graph
			g.adj[i][j] = R[g.adj[i][j]];
		g.m += g.degree[i];
	}
	g.m /= 2;
	return g;
}

/* This function finds all the connected components in a graph and places them in separate clusters.
 * Nodes with degree 0 are placed in a vector called degreeZero.
 * Arguments:
		* clusters: the vector in which the components are stored. Could be non-empty, signifying
					previously found components (For consistency, vertices in old clusters should have degree 0).
		* degreeZero: array where singleton nodes are placed. Cleared before being filled. */
void KGraph::FindConnectedComponents(vector< vector< long> > &clusters, vector<long> &degreeZero)
{
	//cerr<<m<<" edges to start off in FindConnectedComponents."<<endl;
	long v;
	degreeZero.clear();
	bool* label = new bool[n];
	for (long i = 0; i < n; i++)
		label[i] = false;
	for (long i = 0; i < n; i++)
	{
		if (degree[i] != 0 && label[i] == false)
		{
			vector<long> cluster;
			cluster.push_back(i);
			label[i] = true;
			long c = 0;
			while (c != cluster.size())
			{
				long j = cluster[c];
				if (label[j] == false)
				{
					cluster.push_back(j);
					label[j] = true;
				}
				for (long t = 0; t < degree[j]; t++)
				{
					v = adj[j][t];
					if (label[v] == false)
					{
						cluster.push_back(v);
						label[v] = true;
					}
				}
				c++;
			}
			sort(cluster.begin(), cluster.end());
			clusters.push_back(cluster);
		}
	}
	for (long i = 0; i < n; i++)
		if (label[i] == false)
			degreeZero.push_back(i);
	delete[] label;
}

/* Reads a graph from a file in the DIMACS format.
 * The format is described at the bottom of http://www.cc.gatech.edu/dimacs10/downloads.shtml */
void KGraph::ReadDIMACSGraph(string file)
{
	if (n > 0)
		clear();
	cerr << "ReadDIMACSGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	long m1;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File " << file << " not found\n";
		exit(-1);
	}

	bool lineread = false;
	string line;
	while (!lineread)
	{
		getline(input, line);
		if (line[0] != '%')
			lineread = true;
	}
	istringstream nLine = istringstream(line);
	nLine >> n >> m1;

	cerr << n << " nodes, " << m1 << " edges.\n";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << n / 100000 << " dots: ";
	for (long i = 0; i < n; i++)
	{
		lineread = false;
		line.clear();
		while (!lineread)
		{
			getline(input, line);
			if (line[0] != '%')
				lineread = true;
		}
		if ((i + 1) % 100000 == 0)
			cerr << ".";
		if (line == "") continue;
		istringstream iLine = istringstream(line);
		v = -1;
		while (!iLine.eof())
		{
			iLine >> u;
			if (u != v)
			{
				adj[i].push_back(u - 1);
				degree[i]++;
				m++;
				//if(i>(u-1) && !binary_search(adj[u-1].begin(), adj[u-1].end(), i))
				//	cerr<<i<<"-"<<u<<" found, but "<<u<<"-"<<i<<"wasn't\n";
			}
			v = u;
		}
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		if (degree[i] != adj[i].size()) { cerr << " Error on line number " << i << endl; exit(0); }
		Delta = max(degree[i], Delta);
	}
	cerr << endl;
	m = m / 2;
	if (m1 != m)
	{
		cerr << "WRONG DATA!!!!!!!!!!!!!!!!!!!!! " << m << " != " << m1 << endl;
		exit(0);
	}
}

/* Reads a graph from a file in the DIMACS-2 format (clique/coloring challenges).
 * The format is described at  */
void KGraph::ReadDIMACSColorGraph(string file)
{
	if (n > 0)
		clear();
	cerr << "ReadDIMACSColoringGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	long m1;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File " << file << " not found\n";
		exit(-1);
	}

	bool lineread = false;
	string line;
	while (!lineread)
	{
		getline(input, line);
		if (line != "" && line[0] != 'c')
			lineread = true;
	}
	istringstream nLine = istringstream(line);
	string p;
	nLine >> p >> temp >> n >> m1;

	cerr << n << " nodes, " << m1 << " edges.\n";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m1 / 1000000 << " dots: ";
	for (long i = 0; i < m1; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		getline(input, line);
		istringstream nLine = istringstream(line);
		nLine >> temp;
		if (temp == "e")
		{
			nLine >> u >> v;
			if (u == v)
				continue;
			adj[u - 1].push_back(v - 1);
			adj[v - 1].push_back(u - 1);
		}
		else i--;
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i < n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	if (m1 != m)
	{
		cerr << "Possible error in ReadDirectedGraphFromFile: " << m1 << "!=" << m << "\n";
		//	exit(0);
	}
}

/* Reads a graph from a file in the DIMACS format.
 * The format is described at the bottom of http://www.cc.gatech.edu/dimacs10/downloads.shtml */
void KGraph::ReadDIMACSGraphParallel(string file)
{
	if (n > 0)
		clear();
	cerr << "ReadDIMACSGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	long m1;
	string temp;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File " << file << " not found\n";
		exit(-1);
	}

	bool lineread = false;
	string line;
	while (!lineread)
	{
		getline(input, line);
		if (line[0] != '%')
			lineread = true;
	}
	istringstream nLine = istringstream(line);
	nLine >> n >> m1;

	cerr << n << " nodes, " << m1 << " edges.\n";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << n / 100000 << " dots: ";
	int nthreads = 0;
	int tid = 0;
	long lineNum = 0;

#pragma omp parallel shared(input,lineNum) private(line,tid)
	{
		long u, v;
		long share = 0;
		tid = nthreads++;
		cout << "Number of threads = " << tid << endl;
		string line;
		bool lineread;
		long i = 0;

		while (i < n && lineNum < n)
		{
			lineread = false;
			line.clear();
#pragma omp critical
			if (lineNum < n)
			{
#pragma omp flush (lineNum)
				while (!lineread && lineNum < n)
				{
					getline(input, line);
					if (line[0] != '%')
						lineread = true;
				}
				i = lineNum;
				lineNum++;
			}

			if ((i + 1) % 100000 == 0)
				cerr << ".";

			share++;
			if (i < n)
			{
				//cerr<<tid<<"------Using \t"<<i<<"\t"<<lineNum<<endl;
				if (line == "") continue;
				istringstream iLine = istringstream(line);
				v = -1;

				while (!iLine.eof())
				{
					iLine >> u;
					if (u != v)
					{
						adj[i].push_back(u - 1);
						degree[i]++;
#pragma omp atomic
						m++;
						//if(i>(u-1) && !binary_search(adj[u-1].begin(), adj[u-1].end(), i))
						//	cerr<<i<<"-"<<u<<" found, but "<<u<<"-"<<i<<"wasn't\n";
					}
					v = u;
				}
				sort(adj[i].begin(), adj[i].end());
				adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
				if (degree[i] != adj[i].size()) { cerr << " Error on line number " << i << " " << degree[i] << " " << adj[i].size() << endl; exit(0); }
#pragma omp critical(delta)
				Delta = max(degree[i], Delta);
			}
		}
		cerr << tid << " share = " << (double)share / (double)n << endl;
	}
	cerr << endl;
	m = m / 2;
	if (m1 != m)
	{
		cerr << "WRONG DATA!!!!!!!!!!!!!!!!!!!!! " << m << " != " << m1 << endl;
		exit(0);
	}
}
// reads the .dat graphs from Simonetti et al (2011) The Minimum Connected Dominating Set Problem: Formulation, Valid Inequalities and a Branch-and-Cut Algorithm
void KGraph::ReadDATGraph(string file)
{
	cerr << "ReadDATGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	input >> n >> m;
	cerr << n << " nodes, " << m << " edges suggested. ";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m / 1000000 << " dots: ";
	for (long i = 0; i < m; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		input >> u >> v;
		if (u == v)
			continue;
		v--; u--;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i < n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	cerr << m << " edges read\n";
}

void KGraph::ReadDualGraph(string file)
{
	cerr << "ReadDualGraph " << endl;
	n = 0;
	m = 0;
	string temp;

	long Lnum1, Lnum2;
	long z;
	ifstream input;
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	// READ IN NUMBER OF NODES
	input >> temp >> n;
	cerr << "# of Nodes in Dual Graph: " << n << ".\n";
	input >> temp >> m;
	cerr << "# of Edges in Dual Graph: " << m << ".\n";
	pdadj = new vector<long>[m];
	adj = new vector<long>[n];
	edge = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(long));
	// READ IN THE EDGES
	for (long k = 0; k < m; k++)
	{
		input >> temp;
		for (long i = 0; i < temp.length(); i++)
		{
			if (temp[i] == ',')
				z = i;
		}
		string num1, num2;
		//cerr << z << endl;
		for (long i = 0; i < temp.length(); i++)
		{
			if (i == z)
				continue;
			if (i < z && temp[i] != '(')
			{
				num1.push_back(temp[i]);
			}
			if (i > z && temp[i] != ')')
			{
				num2.push_back(temp[i]);
			}
		}
		Lnum1 = stol(num1);
		Lnum2 = stol(num2);

		pdadj[k].push_back(Lnum1);
		pdadj[k].push_back(Lnum2);
		if (Lnum1 == Lnum2)
		{
			adj[Lnum1].push_back(Lnum2);
			edge[Lnum1].push_back(k);
		}
		else
		{
			adj[Lnum1].push_back(Lnum2);
			adj[Lnum2].push_back(Lnum1);
			edge[Lnum1].push_back(k);
			edge[Lnum2].push_back(k);
		}
		degree[Lnum1]++;
		degree[Lnum2]++;
	}
}


void KGraph::ReadPrimalGraph(string file)
{
	cerr << "ReadPrimalGraph " << endl;
	n = 0;
	long n_f = 0;
	string temp;
	long Lnum1, Lnum2;
	long z;
	ifstream input;
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	// READ IN NUMBER OF NODES
	input >> temp >> n;
	//n = n_f - 1;
	cerr << "# of Nodes in Primal Graph: " << n << ".\n";
	input >> temp >> m;
	cerr << "# of Edges in Primal Graph: " << m << ".\n";
	pdadj = new vector<long>[m];
	adj = new vector<long>[n];
	edge = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(long));
	// READ IN THE EDGES
	for (long k = 0; k < m; k++)
	{
		input >> temp;
		for (long i = 0; i < temp.length(); i++)
		{
			if (temp[i] == ',')
				z = i;
		}
		string num1, num2;

		for (long i = 0; i < temp.length(); i++)
		{
			if (i == z)
				continue;
			if (i < z && temp[i] != '(')
			{
				num1.push_back(temp[i]);
			}
			if (i > z && temp[i] != ')')
			{
				num2.push_back(temp[i]);
			}
		}

		Lnum1 = stol(num1);
		//cerr << Lnum1 << endl;
		Lnum2 = stol(num2);
		//cerr << Lnum2 << endl;
		pdadj[k].push_back(Lnum1);
		pdadj[k].push_back(Lnum2);

		if (Lnum1 == Lnum2)
		{
			adj[Lnum1].push_back(Lnum2);
			edge[Lnum1].push_back(k);
		}
		else
		{
			adj[Lnum1].push_back(Lnum2);
			adj[Lnum2].push_back(Lnum1);
			edge[Lnum1].push_back(k);
			edge[Lnum2].push_back(k);
		}
		degree[Lnum1]++;
		degree[Lnum2]++;
	}
}

/* Reads a graph from a file in the SNAP format.
 * The format is described at http://snap.stanford.edu/data/index.html
 * Note: if (i,j) is read from file, both (i,j) and (j,i) are added to the graph to make sure its undirected*/
void KGraph::ReadSNAPGraph(string file)
{
	cerr << "ReadSNAPGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	input >> n >> temp >> m >> temp;

	cerr << n << " nodes, " << m << " edges suggested. ";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m / 1000000 << " dots: ";
	for (long i = 0; i < m; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		input >> u >> v;
		if (u == v)
			continue;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i < n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	cerr << m << " edges read\n";
}

/* Writes the graph to a file in the SNAP format.
 * The format is described at http://snap.stanford.edu/data/index.html
 * Note: if (i,j) is read from file, both (i,j) and (j,i) are added to the graph to make sure its undirected*/
void KGraph::WriteSNAPGraph(string file)
{
	cerr << "WriteSNAPGraph ";
	ofstream output;
	output.open(file.c_str(), ios::out);
	if (!output.is_open())
	{
		cout << "File " << file << " could not be opened!!!\n";
		return;
	}

	output << n << " nodes, " << m << " edges.\n";

	cerr << "Wait till " << n / 100000 << " dots: ";
	for (long i = 0; i < n; i++)
	{
		if ((i + 1) % 100000 == 0)
			cerr << ".";
		for (long j = 0; j < degree[i]; j++)
			output << i << "\t" << adj[i][j] << endl;
	}
}

/* Writes the graph to a file in the DIMACS-10 format.
 * The format is described at the bottom of http://www.cc.gatech.edu/dimacs10/downloads.shtml */
void KGraph::WriteDIMACSGraph(string file)
{
	cerr << "WriteDIMACSGraph: " << file << endl;
	ofstream output;
	output.open(file.c_str(), ios::out);
	if (!output.is_open())
	{
		cout << "File " << file << " could not be opened!!!\n";
		return;
	}

	output << n << " " << m << endl;

	cerr << "Wait till " << n / 100000 << " dots: ";
	for (long i = 0; i < n; i++)
	{
		if ((i + 1) % 100000 == 0)
			cerr << ".";
		for (long j = 0; j < degree[i]; j++)
			output << adj[i][j] + 1 << " ";
		output << endl;
	}
}

/* Writes the graph in the format such that GraphViz can plot it.
 * No position is specified, so graphviz determines whats best */
bool KGraph::WriteGVizGraph(string file)
{
	vector<bool> Final(n, false);
	return WriteGVizGraph(Final, file);
}
bool KGraph::WriteGVizGraph(vector<bool> Final, string file)
{
	string outfile = file + ".viz";
	ofstream gviz;
	gviz.open(outfile.c_str(), ios::out);
	gviz << "graph test\n{\nnode [shape=circle, pin=true, fontsize=1];\n";
	double scale = 10;
	for (long i = 0; i < n; i++)
		gviz << i + 1 << ";\n";
	long temp;
	for (long i = 0; i < n; i++)
	{
		for (long j = 0; j < degree[i]; j++)
		{
			temp = adj[i][j];
			if (i > temp) continue;
			if (Final[i] && Final[temp])
			{
				gviz << i + 1 << " -- " << temp + 1
					<< "[" << "color = black" << "]"
					<< ";\n";
			}
			else if (Final[i] || Final[temp])
			{
				gviz << i + 1 << " -- " << temp + 1
					<< "[" << "color = gray55" << "]"
					<< ";\n";
			}
			else
			{
				gviz << i + 1 << " -- " << temp + 1
					<< "[" << "color = gray92" << "]"
					<< ";\n";
			}
		}
	}
	gviz << "}";
	gviz.close();
	return true;
}

void KGraph::clear()
{
	//cout<<"===Clearing DG "<<name<<"\n";
	if (n > 0)
	{
		for (long i = 0; i < n; i++)
			adj[i].clear();
		delete[]adj;
		delete[]degree;
		n = 0;
		m = 0;
	}
}

KGraph::~KGraph()
{
	//cout<<"===Destructing DG "<<name<<"\n";
	if (n > 0)
	{
		for (long i = 0; i < n; i++)
			adj[i].clear();
		delete[]adj;
		delete[]degree;
		n = 0;
		m = 0;
	}
}
void KGraph::ReadIEEEGraph(string file)
{
	cerr << "ReadIEEEGraph ";
	m = 0;
	n = 0;
	Delta = 0;
	string temp;
	long u, v;
	ifstream input;
	char* t = "";
	input.open(file.c_str(), ios::in);
	if (!input.is_open())
	{
		cout << "File not found\n";
		exit(-1);
	}
	input >> n >> m;

	cerr << n << " nodes, " << m << " edges suggested. ";

	adj = new vector<long>[n];
	degree = new long[n];
	memset(degree, 0, n * sizeof(int));

	cerr << "Wait till " << m / 1000000 << " dots: ";
	for (long i = 0; i < m; i++)
	{
		if ((i + 1) % 1000000 == 0)
			cerr << ".";
		input >> u >> v;
		if (u == v)
			continue;
		v--; u--;
		adj[u].push_back(v);
		adj[v].push_back(u);
	}
	cerr << endl;
	m = 0;
	for (long i = 0; i < n; i++)
	{
		sort(adj[i].begin(), adj[i].end());
		adj[i].erase(unique(adj[i].begin(), adj[i].end()), adj[i].end());
		degree[i] = adj[i].size();
		m += degree[i];
		Delta = max(degree[i], Delta);
		if (adj[i].size() != degree[i])
		{
			cerr << "Error in ReadDirectedGraphFromFile\n";
			exit(0);
		}
	}
	m = m / 2;
	cerr << m << " edges read\n";
}
