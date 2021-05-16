#include<bits/stdc++.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

// This class represents a directed graph
// using adjacency list representation
class Graph
{
  int V;  // No. of vertices
public:
  list< pair<int, int> > *adj;  // List of all edge pairs
  vector<int> path;

  Graph(int V);
  int vertNum();
  void addEdge(int u, int v, int w);
  void modifyEdge(int u, int v, int w);
  void Dijkstras_Algo(int s, int dest);
  void BellmanFord_Algo(int s);
  void FloydWarshall_Algo(int s);
};
// Return number of vertices in the graph
int Graph::vertNum()
{
  return this->V;
}
// Constructor: Allocate memory for adjacency list
Graph::Graph(int V)
{
  this->V = V;
  adj = new list< pair<int, int> >[V];
}
// Add an edge to the graph
void Graph::addEdge(int u, int v, int w)
{
  adj[u].push_back(make_pair(v, w));
  adj[v].push_back(make_pair(u, w));
}
// Modify an edge in the graph
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
// Shortest paths from src to dest or to all other vertices
void Graph::Dijkstras_Algo(int src, int dest)
{
  
}

void Graph::BellmanFord_Algo(int src)
{
  
}

void Graph::FloydWarshall_Algo(int src)
{
  
}

int main()
{
  /*
  int V = 9;
  Graph g(V);

  // making a graph
  g.addEdge(0, 1, 4);
  g.addEdge(0, 7, 8);
  g.addEdge(1, 2, 8);
  g.addEdge(1, 7, 11);
  g.addEdge(2, 3, 7);
  g.addEdge(2, 8, 2);
  g.addEdge(2, 5, 4);
  g.addEdge(3, 4, 9);
  g.addEdge(3, 5, 14);
  g.addEdge(4, 5, 10);
  g.addEdge(5, 6, 2);
  g.addEdge(6, 7, 1);
  g.addEdge(6, 8, 6);
  g.addEdge(7, 8, 7);
  */
  
  fstream disFile;
  disFile.open("./USA-road-d_NY.txt");
  vector<vector<int>> vect;
  
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
  
  Graph g(264346);
  cout<<"Number of edges: "<<vect.size()<<'\n';
    for (int i = 0; i < vect.size(); i++) {
      g.addEdge(vect[i][0]-1,vect[i][1]-1,vect[i][2]);
    }
  
  // Code for algorithm  
}
