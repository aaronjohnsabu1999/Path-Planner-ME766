#include <omp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <bits/stdc++.h>

#define INFTY   1e8

#ifdef __CUDACC__
#define CUDA_CALLABLE_MEMBER __host__ __device__
#else
#define CUDA_CALLABLE_MEMBER
#endif 

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

__global__ void matrixMul(Graph g)
{
  for(int mid = 0; mid < n_v; mid++)
	{
		for (int i = 0; i < n_v*n_v; i++)
		{			
				list< pair<int, int> >::iterator it;
				int w1=1e8;
				int w2=1e8;
				int w3=1e8;
				
				int start = i / n_v;
				int end = i % n_v;
				for (it = g.adj[start].begin(); it != g.adj[start].end(); it++)
				{
					if ((*it).first == end)
					(*it).second = w1;
					if((*it).first == mid)
					(*it).second = w2;
				}
				
				for (it = g.adj[mid].begin(); it != g.adj[mid].end(); it++)
				{
					if ((*it).first == end)
					(*it).second = w3;

				}
				
				int w = min(w1, w2+w3);
				h.modifyEdge(start,end,w);

		}

		#pragma omp parallel for shared(g)
		for(int i=0; i < n_v; i++)
		{
			list< pair<int, int> >::iterator it;
			for(it = h.adj[i].begin(); it != dis.adj[i].end(); it++)
			{
				g.modifyEdge(i, (*it).first, (*it).second);
			}
			
		}
	}
	
}

int main(int argv, char **argc)
{
	int V = strtol(argc[1], (char **)NULL, 10), E = 3*V;
	
	Graph g(V);
	Graph h(V);
	g.addEdge(0,2,1);
	for (int j = 1; j < E; j++) {
		g.addEdge(rand()%(V) + 0, rand()%(V) + 0, rand()%(50 - 1 + 1) + 1);
	}


	int numThreads = strtol(argc[2], (char **)NULL, 10), threadNum;
	omp_set_num_threads(numThreads);

	int n_v = g.vertNum();
	
  int threads = 10 * 10;
  int blocks  = (N + threads - 1) / threads;
  double t1 = omp_get_wtime();
  dim3 THREADS (threads, threads);
  dim3 BLOCKS  ( blocks,  blocks);
	matrixMul<<<BLOCKS, THREADS>>>(g);
	cudaDeviceSynchronize();
  double t2 = omp_get_wtime() - t1;
	cout<<"Total Time taken for \t"<<V<<" vertices = "<<t2<<" seconds."<<endl;

	return 0;
}
