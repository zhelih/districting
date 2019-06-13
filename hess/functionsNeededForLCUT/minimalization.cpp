vector<int> minimalization(graph* g, vector<int>& C_0, const vector<int>& population, int a, int b, int U, int total_pop)
{
    vector<int> minimalSet;

    vector<int> distanceFrom_a(g->nr_nodes);

    vector<bool> deletedNodes(g->nr_nodes, false);

    for (int i = 0; i < C_0.size(); ++i)
    {
        int cur = C_0[i];
        deletedNodes[cur] = true;
    }

    for (int i = 0; i < C_0.size(); ++i)
    {
        int cur = C_0[i];
        deletedNodes[cur] = false;
        distanceFrom_a = g->ShortestPathsNodeWeighted(a, population, deletedNodes, total_pop);
        if (distanceFrom_a[b] <= U)
        {
            minimalSet.push_back(cur);
            deletedNodes[cur] = true;
        }
    }
    return minimalSet;
}
