#include <bits/stdc++.h>
  
struct Edge {
    int src, dest, weight;
};
  
struct Graph {
// V-> Number of vertices, E-> Number of edges
    int V, E;
  
    struct Edge* edge;
};
  
struct Graph* createGraph(int V, int E)
{
    struct Graph* graph = new Graph;
    graph->V = V;
    graph->E = E;
    graph->edge = new Edge[E];
    return graph;
}
  
void printArr(int dist[], int n)
{
    printf("Vertex   Distance from Source\n");
    for (int i = 0; i < n; ++i)
        printf("%d \t\t %d\n", i, dist[i]);
}
  
void BellmanFord(struct Graph* graph, int src)
{
    int V = graph->V;
    int E = graph->E;
    int dist[V];
  
    // Step 1: Initialize distances from src to all other vertices as INFINITE
    #pragma omp parallel for
    for (int i = 0; i < V; i++)
        dist[i] = INT_MAX;

    dist[src] = 0;
  
    // Step 2: Relax all edges |V| - 1 times
    
    for (int i = 1; i <= V - 1; i++) 
    {   
        #pragma omp parallel for // Order in which they are done doesn't matter
        for (int j = 0; j < E; j++) 
        {
            int u = graph->edge[j].src;
            int v = graph->edge[j].dest;
            int weight = graph->edge[j].weight;
            if (dist[u] != INT_MAX && dist[u] + weight < dist[v])
               {
               #pragma omp atomic write
                dist[v] = dist[u] + weight;
            }
        }
    }
    #pragma omp barrier
    // Step 3: check for negative-weight cycles
    int flag = 1;
    #pragma omp parallel for
    for (int i = 0; i < E; i++) 
    {
        int u = graph->edge[i].src;
        int v = graph->edge[i].dest;
        int weight = graph->edge[i].weight;
        if (dist[u] != INT_MAX && dist[u] + weight < dist[v]) {
            flag = -1;
        }
    }
    #pragma omp barrier
    if(flag == -1)
    printf("Graph contains negative weight cycle");
    //printArr(dist, V);
    
    return;
}
  
// Driver program to test above functions
int main()
{
    int V = 5000; // Number of vertices in graph
    int E = 800; // Number of edges in graph
    struct Graph* graph = createGraph(V, E);
  
    for (int j = 1; j < E; j++) 
    {
    graph->edge[j].src    = rand()%(V) + 0;
    graph->edge[j].dest   = rand()%(V) + 0;
    graph->edge[j].weight = rand()%(50 - 1 + 1) + 1;
  }
  
    BellmanFord(graph, 0);
  
    return 0;
}
