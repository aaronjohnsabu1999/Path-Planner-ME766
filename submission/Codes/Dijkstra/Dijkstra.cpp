#include<bits/stdc++.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <omp.h>

# define INF 0x3f3f3f3f
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
  void Dijkstras_Algo(int s, int dest);
  void Dijkstras_Algo_omp(int s, int dest);
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
// Shortest paths from src to dest or to all other vertices
void Graph::Dijkstras_Algo(int src, int dest)
{
	// Create a set to store vertices that are being preprocessed
	set< pair<int, int> > setds;

	// Create a vector for distances and initialize all distances as infinite (INF)
	vector<int> dist(V, INF);
	vector<int> parent(V, 0); 
	parent[src] = -1;
	// Insert source itself in Set and initialize its distance as 0.
	setds.insert(make_pair(0, src));
	dist[src] = 0;

	// Looping till all shortest distance are finalized
	// then setds will become empty
	double start, end;
	start = omp_get_wtime();	
	while (!setds.empty())
	{
		// The first vertex in Set is the minimum distance
		// vertex, extract it from set.
		pair<int, int> tmp = *(setds.begin());
		setds.erase(setds.begin());

		// vertex label is stored in second of pair (it
		// has to be done this way to keep the vertices
		// sorted distance (distance must be first item
		// in pair)
		int u = tmp.second;

		// 'i' is used to get all adjacent vertices of a vertex
		list< pair<int, int> >::iterator i;
		for (i = adj[u].begin(); i != adj[u].end(); ++i)
		{
			// Get vertex label and weight of current adjacent of u.
			int v = (*i).first;
			int weight = (*i).second;

			// If there is shorter path to v through u.
			if (dist[v] > dist[u] + weight)
			{
				parent[v] = u;
				/* If distance of v is not INF then it must be in
					our set, so removing it and inserting again
					with updated less distance.
					Note : We extract only those vertices from Set
					for which distance is finalized. So for them,
					we would never reach here. */
				if (dist[v] != INF)
					setds.erase(setds.find(make_pair(dist[v], v)));
				// Updating distance of v
				dist[v] = dist[u] + weight;
				setds.insert(make_pair(dist[v], v));
			}
		}
	}
	end = omp_get_wtime();
	printf("Serial Running time: %f ms\n", (end - start)*1000);
	int j;
	for (int i = 0; i < V; ++i){
		j = i+ 1;
		ofstream myfile ("USA_NY_out.txt", ios::app);
		myfile << j << " \t\t " << dist[i] << endl;
		myfile.close();		
	}
}

int search(int *a, int item, int low, int high)
{
    if (high <= low)
        return (item > a[low]) ?
                (low + 1) : low;
 
    int mid = (low + high) / 2;
 
    if (item == a[mid])
        return mid + 1;
 
    if (item > a[mid])
        return search(a, item,
                            mid + 1, high);
    return search(a, item, low,
                        mid - 1);
}

void sort(int *a1, int *a2, int index)
{
	int i, loc, j, k, selected2, selected1;

	for(i = 1; i < index; i++)
	{
		j = i - 1;
		selected2 = a2[i];
		selected1 = a1[i];

		loc = search(a2, selected2, 0, j);

		while( j >= loc)
		{
			a2[j+1] = a2[j];
			a1[j+1] = a1[j];
			j--;
		}
		a2[j+1] = selected2;
		a1[j+1] = selected1;
	}

}

void erase(int start, int *a1, int *a2, int index)
{
	for(int i=start; i < index + 2; i++){
		a1[i] = a1[i+1];
		a2[i] = a2[i+1];
	}
}

void Graph::Dijkstras_Algo_omp(int src, int dest)
{
	int vertices = this->V;
	int i,v,weight,d;
	bool * connected = new bool[vertices];


	for(i = 0; i<vertices;i++)
	{
		connected[i] = false;
	}
	connected[src] = 0;

	int *dist = new int[vertices];
	for(i=0; i<vertices; i++){
		dist[i] = INF;
	}
	dist[src] = 0;

	vector<list< pair<int, int> >::iterator> partial_element;

	vector< vector<list< pair<int, int> >::iterator> > element;

	list< pair<int, int> >::iterator k;
	
	for(i=0; i<vertices; i++)
	{
		partial_element.clear();
		for (k = adj[i].begin(); k != adj[i].end(); ++k)
		{
			partial_element.push_back(k);
		}
		element.push_back(partial_element);		
	}

	for (k = adj[src].begin(); k != adj[src].end(); ++k)
	{
		v = (*k).first;
		weight = (*k).second;
		dist[v] = weight;
	}

	int first, thread, n_thread, last, min_d, min_node, step;
	bool check = true;

	double start, end;
	start = omp_get_wtime();

	for(step=0;step<vertices;step++)
	{
		#pragma omp parallel private(i, weight, k, first, last, thread, step, min_d, min_node) shared(connected,d,dist,v,n_thread,element,check)
		{
			thread = omp_get_thread_num();
			n_thread = omp_get_num_threads();
			first = (thread*vertices)/n_thread;
			last = ((thread+1)*vertices)/n_thread;

			#pragma omp single
			{
				d = INF;
				v = -1;
			}

			min_d = INF;
			min_node = -1;
			for (int i = first; i <= last; i++ )
			{
				if ( !connected[i] && dist[i] < min_d)
				{
				  min_d = dist[i];
				  min_node = i;
				}
			}

			#pragma omp critical
			{
				if(min_d < d)
				{
					d = min_d;
					v = min_node;
				}
			}

			#pragma omp barrier

			#pragma omp single
			{
				if(v != -1)
				{
					connected[v] = true;
				}
			}
			#pragma omp barrier					
		}

		if(v != -1)
		{
			int j;
			#pragma omp parallel for schedule(runtime) private(i, j, weight, k) shared(dist,v,element)
			for(j =0; j < element[v].size(); j++)
			{
				k = element[v].at(j);
				i = (*k).first;
				if((!connected[i]))
				{
					weight = (*k).second;
					if(dist[i] > dist[v] + weight)
					{
						dist[i] = dist[v] + weight;
					}
				}
			}				
		}
	}
    
	end = omp_get_wtime();
	printf("Parallel Running time: %f ms\n", (end - start)*1000);	

	delete[] connected;

	int j;
	for (int i = 0; i < V; ++i){
		j = i+ 1;
		ofstream myfile ("USA_NY_parallel_out.txt", ios::app);
		myfile << j << " \t\t " << dist[i] << endl;
		myfile.close();			
	}
	delete[] dist;
}
void Graph::BellmanFord_Algo(int src)
{
  
}
void Graph::FloydWarshall_Algo(int src)
{
  
}

int main()
{
	fstream disFile;
	disFile.open("USA-road-d_NY.txt");
	
	vector< vector<int> > vect;
	
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
	
	g.Dijkstras_Algo(0, 100);
	g.Dijkstras_Algo_omp(0, 100);
	return 0;
}