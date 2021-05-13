#include <omp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <bits/stdc++.h>

#define INFTY   1e8

using namespace std;

class Graph
{
  int V;
public:
  list< pair<int, int> > *adj;
  vector<int> path;

  Graph(int V);
  int vertNum();
  void addEdge(int u, int v, int w);
  void modifyEdge(int u, int v, int w);
  void Dijkstras_Algo(int s, int dest);
  void BellmanFord_Algo(int s);
  void FloydWarshall_Algo(int s);
};

int Graph::vertNum()
{
  return this->V;
}
Graph::Graph(int V)
{
  this->V = V;
  adj = new list< pair<int, int> >[V];
}
void Graph::addEdge(int u, int v, int w)
{
  adj[u].push_back(make_pair(v, w));
  adj[v].push_back(make_pair(u, w));
}
void Graph::modifyEdge(int u, int v, int w)
{
  list< pair<int, int> >::iterator it;
  for (it = this->adj[u].begin(); it != this->adj[u].end(); it++)
  {
    if ((*it).first == v)
      (this->adj)[u].erase(it);
    break;
  }
  for (it = this->adj[v].begin(); it != this->adj[v].end(); it++)
  {
    if ((*it).first == u)
      (this->adj)[v].erase(it);
    break;
  }
  
  adj[u].push_back(make_pair(v, w));
  adj[v].push_back(make_pair(u, w));
}

int main(int argv, char **argc)
{
  const int V    = strtol(argc[1], (char **)NULL, 10), E = 3*V;
  int numThreads = strtol(argc[2], (char **)NULL, 10);
  int *graph = new int[V*V];
  int *dist  = new int[V*V];
  Graph g(V);
  
  g.addEdge(0,2,1);
  for (int j = 1; j < E; j++)
    g.addEdge(rand()%(V) + 0, rand()%(V) + 0, rand()%(50) + 1);
  
  omp_set_num_threads(numThreads);
  
  double t1 = omp_get_wtime();
  
  for(int i=0; i < V; i++)
  {
    for(int j=0; j<V; j++)
      graph[i*V + j] = 1e9;
    
    list< pair<int, int> >::iterator it;
    for(it = g.adj[i].begin(); it != g.adj[i].end(); it++)
      graph[i*V + (*it).first] = (*it).second;
  }
  
  #pragma omp parallel
  for(int mid = 0; mid < V; mid++)
  {
    #pragma omp for
    for (int i = 0; i < V*V; i++)
      dist[(i/V)*V + (i % V)] = min( graph[(i / V)*V + (i % V)], graph[(i / V)*V + mid] + graph[mid*V + (i % V)] );
    
    #pragma omp for
    for(int i=0; i < V*V; i++)
      graph[(i / V)*V + (i % V)] = dist[(i / V)*V + (i % V)];
  }
  double t2 = omp_get_wtime() - t1;
  cout<<"Time taken for OpenMP implementation with (V = \t"<<V<<") and (numThreads = \t"<<numThreads<<") = "<<t2<<endl;
  
  delete graph;
  delete dist;
  return 0;
}
