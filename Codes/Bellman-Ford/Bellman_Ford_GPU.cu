#include <bits/stdc++.h>
#include <stdio.h>
#include <climits>
#include <random>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;
  
struct Edge 
{
    int src, dest, weight;
};
    
__global__ void SetupDist(int V, int *dist, int src) 
{	
	// Initialize Distances as INT_Max
    int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    for (int i = index; i < V; i += stride)
        dist[i] = INT_MAX;
    dist[src] = 0;
}
  
__global__ void BellmanFord(int V, int E, struct Edge* edges, int *dist) {
  // Kernel with Grid Stride Loop
  int index = blockIdx.x * blockDim.x + threadIdx.x;
  int stride = blockDim.x * gridDim.x;
  
 // Step 2: Relax all edges |V| - 1 times
  for (int j = index; j < E; j += stride) {
    int u      = edges[j].src;
    int v      = edges[j].dest;
    int weight = edges[j].weight;

 // Atomic to avoid two threads writing simultaneously
	if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
	  atomicMin(&dist[v], dist[u] + weight);
  }
  return;
}
  
__global__ void Check_Neg_Cycle(int E, struct Edge* edges, int *dist) 
{

  int index = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = blockDim.x * gridDim.x;

    int flag = 1;
    for (int i = index; i < E; i += stride) {
        int u = edges[i].src;
        int v = edges[i].dest;
        int weight = edges[i].weight;
        if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
            flag = -1;
        }
    }
    if(flag == -1)
        printf("Graph contains negative weight cycle");

return;
}

// Driver program to test above functions
int main()
{
  vector<vector<int>> vect;
    fstream disFile;
    disFile.open("../USA-road-d_NY.txt");

    if (disFile.is_open()){ 
        string tp;

        while(getline(disFile, tp)){ 
            vector<int> pair;
            if(tp[0] == 'a'){
                int a;
                string str(tp.begin()+1, tp.begin()+tp.size() );
                stringstream ss (str);
                while ((ss >> a))
                    pair.push_back (a);
                vect.push_back(pair);
            }
        }

        disFile.close(); 
   }

    int V = 264346; // Number of vertices in graph
    int E = vect.size(); // Number of edges in graph

    int *host_dist;
    Edge *host_edges;
    host_edges    = (Edge*)malloc(E*sizeof(Edge));
    host_dist     = (int*)malloc(V*sizeof(int));

	Edge *gpu_edges;
	int *gpu_dist;

	// allocate memory
	cudaMalloc(&gpu_edges, E * sizeof(Edge));
	cudaMalloc(&gpu_dist, V * sizeof(int));

    for (int j = 1; j < E; j++) 
    {
    host_edges[j].src    = vect[j][0]-1;
    host_edges[j].dest   = vect[j][1]-1;
    host_edges[j].weight = vect[j][2];
  	}
  
    cout<<host_edges[10].weight;
	
	// copy over to device/GPU
	cudaMemcpy(gpu_edges, host_edges, E*sizeof(Edge), cudaMemcpyHostToDevice);
	cudaMemcpy(gpu_dist, host_dist, V*sizeof(int), cudaMemcpyHostToDevice);

	// Setup threads
    int blockSize = 256;
    int numBlocksE = (E + blockSize - 1) / blockSize;

    int numBlocksV = (V + blockSize - 1) / blockSize;

    SetupDist<<<numBlocksV, blockSize>>>(V, gpu_dist, 0);
    cudaDeviceSynchronize();

    for (int i = 1; i <= V - 1; i++) 
        {
            BellmanFord<<<numBlocksE, blockSize>>>(V, E, gpu_edges, gpu_dist);
            cudaDeviceSynchronize();
        }

Check_Neg_Cycle<<<numBlocksE, blockSize>>>(E, gpu_edges, gpu_dist);
cudaDeviceSynchronize();

  // Copy results to host
  cudaMemcpy(host_dist, gpu_dist, V*sizeof(int), cudaMemcpyDeviceToHost);

  // cleanup memory
  cudaFree(gpu_edges);
  cudaFree(gpu_dist);
  free(host_dist);
  free(host_edges);
  
  return 0;
}
