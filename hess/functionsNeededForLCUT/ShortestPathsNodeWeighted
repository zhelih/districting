// Much of this code is written by Anurag Verma (student of Sergiy Butenko).
// Modified by Hamidreza Validi.

vector<int> graph::ShortestPathsNodeWeighted(int origin, const vector<int>& population, vector<bool> &deletedNodes, int total_pop)
{
    int v, u;
    int infty = total_pop + 1;
    vector<int> dist(nr_nodes, infty); //shortest distance from origin node to each other node. dist[i] = \infty means i not reachable
    vector<bool> Q(nr_nodes, true); // the set of vertices whose distance labels are not yet permanent
    dist[origin] = population[origin]; //the distance from the origin node is equal to its population
    int minDistance;
    int minVertex;
    
    //make Q[v]= false for deleted nodes
    for (int i = 0; i < nr_nodes; ++i)
    {
        if (!deletedNodes[i]) continue;
        Q[i] = false;
    }
    
    // do node-weighted version of Dijkstra's shortest path algorithm.
    for (int i = 0; i < nr_nodes; ++i)
    {
        // find a vertex u from Q of minimum distance label dist[u]
        minDistance = infty;
        minVertex = -1;
        for (int v = 0; v < nr_nodes; ++v)
        {
            if (!Q[v]) continue;
            if (dist[v] < minDistance)
            {
                minDistance = dist[v];
                minVertex = v;
            }
        }
        
        if (minVertex == -1) continue;
        
        // remove minVertex from Q
        Q[minVertex] = false;

        // update distance labels of neighbors
        for (int v : nb(minVertex))
        {
            if (Q[v] && dist[minVertex] + population[v] < dist[v])
            {
                dist[v] = dist[minVertex] + population[v];
            }
        }
    }
    return dist;
}
