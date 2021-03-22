// https://www.techiedelight.com/least-cost-path-weighted-digraph-using-bfs/
#include <iostream>
#include <vector>
#include <queue>
using namespace std;

struct Edge {
	int source, dest, weight;
};

class Graph
{
public:
	vector<vector<int>> adjList;
	Graph(vector<Edge> const &edges, int N, int maxWeight)
	{
		adjList.resize(maxWeight*N);

		for (auto &edge: edges)
		{
			int v = edge.source;
			int u = edge.dest;
			int weight = edge.weight;
    	for(int intermediate = 0; intermediate < weight - 1; intermediate++)
			{
				adjList[v + intermediate*N].push_back(v + (intermediate+1)*N);
				adjList[v + (intermediate+1)*N].push_back(v + intermediate*N);
			}
			adjList[v + (weight-1)*N].push_back(u);
			adjList[u].push_back(v + (weight-1)*N);
			
		}
	}
};

void printPath(vector<int> const &predecessor, int v, int &cost, int N)
{
	if (v < 0) {
		return;
	}
	printPath(predecessor, predecessor[v], cost, N);
	cost++;
	if (v < N) {
		cout << v << " ";
	}
}

void BFS(Graph const &graph, int source, int dest, int N)
{
	vector<bool> discovered(3*N, false);
	discovered[source] = true;
	vector<int> predecessor(3*N, -1);
	queue<int> q;
	q.push(source);
  int foundPath = 0;
	while (!q.empty())
	{
		int curr = q.front();
		q.pop();
		if (curr == dest)
		{
			int cost = -1;
			cout << "The least-cost path between " << source << " and " << dest << " is ";
      printPath(predecessor, dest, cost, N);
      cout << "having cost " << cost;
      foundPath = 1;
		}
    for (int v: graph.adjList[curr])
		{
			if (!discovered[v])
			{
				discovered[v] = true;
				q.push(v);
    		predecessor[v] = curr;
			}
		}
	}
  if (foundPath == 0)
    cout << "Found no path";    
}

int main()
{
	int N = 5;
	vector<Edge> edges =
	{
		{0, 1, 3}, {0, 4, 1}, {1, 2, 2}, {1, 3, 3},
		{1, 4, 8}, {4, 2, 2}, {4, 3, 8}
	};
	Graph graph(edges, N, 20);
	int source = 0, dest = 2;
	BFS(graph, source, dest, 5);
	return 0;
}