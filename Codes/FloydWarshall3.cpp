#include <omp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <bits/stdc++.h>

#define	INFTY	1e8

using namespace std;

int main(int argv, char **argc)
{
	int V = strtol(argc[1], (char **)NULL, 10), E = 0;
	int numThreads = strtol(argc[2], (char **)NULL, 10), threadNum;
	int *g = new int(V*V);
	int *h = new int(V*V);
	
	omp_set_num_threads(numThreads);
	
	for (int start = 0; start < V; start++)
		for (int end = 0; end < V; end++)
			if (rand()%((V*V)/(3*V)))
			{
				g[start*V + end] = rand()%500 + 1;
				E++;
			}
			else
				g[start*V + end] = INFTY;
	
	// #pragma omp parallel
	for (int mid = 0; mid < V; mid++)
	{
		// #pragma omp for
		for (int startEnd = 0; startEnd < V*V; startEnd++)
		{
			int start = startEnd / V;
			int end   = startEnd % V;
			h[start*V + end] = min(g[start*V + end], g[start*V + mid] + g[mid*V + end]);
		}
		
		// #pragma omp for
		for (int startEnd = 0; startEnd < V*V; startEnd++)
			g[startEnd] = h[startEnd];
	}
	
	return 0;
}