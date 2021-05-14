#include<bits/stdc++.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <bits/stdc++.h>
#include <omp.h>

# define INF 0x3f3f3f3f
using namespace std;

// This class represents a directed graph using
// adjacency list representation
class Graph
{
	int V; // No. of vertices

	// In a weighted graph, we need to store vertex
	// and weight pair for every edge
	//list< pair<int, int> > *adj;

public:
	list< pair<int, int> > *adj;

	vector<int> path;

	Graph(int V); // Constructor

	// function to add an edge to graph
	void addEdge(int u, int v, int w);

	// prints shortest path from s
	void Dijkstras_Algo(int s, int dest);
	void Dijkstras_Algo_omp(int s, int dest);
	void BellmanFord_Algo(int s);
	void FloydWarshall_Algo(int s);
};

// Allocates memory for adjacency list
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

// Prints shortest paths from src to all other vertices.. dest is the maximum weigth on a edge in the graph

void Graph::Dijkstras_Algo(int src, int dest)
{
	// Create a set to store vertices that are being
	// prerocessed
	set< pair<int, int> > setds;

	// Create a vector for distances and initialize all
	// distances as infinite (INF)
	vector<int> dist(V, INF);
	vector<int> parent(V, 0); 
	parent[src] = -1;
	// Insert source itself in Set and initialize its
	// distance as 0.
	setds.insert(make_pair(0, src));
	dist[src] = 0;

	/* Looping till all shortest distance are finalized
	then setds will become empty */
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
			// Get vertex label and weight of current adjacent
			// of u.
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
	//vector<vector<int>> vec(V);
	// Print shortest distances stored in dist[]
	//printf("Vertex Distance from Source\n");
	int j;
	for (int i = 0; i < V; ++i){
		//printf("%d \t\t %d\n", i, dist[i]);
		j = i+ 1;
		ofstream myfile ("USA_NY_out.txt", ios::app);
		//ofstream myfile ("tmp_out.txt", ios::app);
		myfile << j << " \t\t " << dist[i] << endl;
		myfile.close();		
		//j = i;
		//printf("%d  ", j);
		//while(parent[j]!=-1){
		//	vec[i].push_back(j);
		//	j = parent[j];
			//printf("%d  ", j);
		//}
		//printf("\n");
		//reverse(vec[i].begin(), vec[i].end());
		//for(int k=0; k< vec[i].size(); k++){
		//	printf("%d  ", vec[i][k]);
		//}
		//printf("\n");		
	}
	//path = vec[dest];
	//return_path(vec[dest]);
	/*vector<int> z = vec[dest];
	for(int it=0; it<path.size();it++){
		printf("%d  ", path[it]);
	}
	return z;*/
}


// void Graph::Dijkstras_Algo_omp_old(int src, int dest)
// {
// 	// Create a set to store vertices that are being
// 	// prerocessed
// 	set< pair<int, int> > setds;

// 	// Create a vector for distances and initialize all
// 	// distances as infinite (INF)
// 	vector<int> dist(V, INF);
// 	vector<int> parent(V, 0); 
// 	parent[src] = -1;
// 	// Insert source itself in Set and initialize its
// 	// distance as 0.
// 	setds.insert(make_pair(0, src));
// 	dist[src] = 0;

// 	/* Looping till all shortest distance are finalized
// 	then setds will become empty */
// 	std::vector<list< pair<int, int> >::iterator> element;
// 	int count = 0;
// 	int count1 = 0;
// 	double start, end;
// 	start = omp_get_wtime();
// 	while (!setds.empty())
// 	{
// 		cout<< count++ <<endl;
// 		// The first vertex in Set is the minimum distance
// 		// vertex, extract it from set.
// 		pair<int, int> tmp = *(setds.begin());
// 		setds.erase(setds.begin());

// 		// vertex label is stored in second of pair (it
// 		// has to be done this way to keep the vertices
// 		// sorted distance (distance must be first item
// 		// in pair)
// 		int u = tmp.second;

// 		// 'i' is used to get all adjacent vertices of a vertex
// 		list< pair<int, int> >::iterator i;
		
// 		int v;
// 		int weight;

// 		element.clear();
// 		for (i = adj[u].begin(); i != adj[u].end(); ++i)
// 		{
// 			element.push_back(i);
// 		}

// 		i = element.at(0);
// 		int j;
// 		// std::vector<int> a1;
// 		// std::vector<int> a2;
// 		int *b1 = new int[element.size()];
// 		int *b2 = new int[element.size()];
// 		count1 = 0;

// 		#pragma omp parallel for schedule(runtime) private(i, j, v, weight) shared(setds,dist,parent,u,b1,b2,count1)
// 		//#pragma omp parallel for shared(setds, a1, a2)
// 		//#pragma omp parallel for schedule(runtime) private(j)
// 		for (j = 0; j < element.size(); ++j)
// 		{
// 			//cout<< j << endl;
// 			// Get vertex label and weight of current adjacent
// 			// of u.
			
// 			//cout<< "1" << endl;
// 			i = element.at(j);
// 			//cout<< "2" << endl;
// 			v = (*i).first;
// 			//cout<< "3" << endl;
// 			weight = (*i).second;
// 			//cout<< "4" << endl;

// 			// If there is shorter path to v through u.
// 			if (dist[v] > dist[u] + weight)
// 			{
// 				//parent[v] = u;
// 				/* If distance of v is not INF then it must be in
// 					our set, so removing it and inserting again
// 					with updated less distance.
// 					Note : We extract only those vertices from Set
// 					for which distance is finalized. So for them,
// 					we would never reach here. */
// 				if (dist[v] != INF)
// 					setds.erase(setds.find(make_pair(dist[v], v)));

// 				// Updating distance of v
// 				dist[v] = dist[u] + weight;
// 				#pragma omp critical
// 				{
// 					//setds.insert(make_pair(dist[v], v));
// 					// a1.push_back(dist[v]);
// 					// a2.push_back(v);
// 					b1[count1] = dist[v];
// 					b2[count1] = v;
// 					count1++;
// 				}
// 				//setds.insert(make_pair(3, 4));
// 			}
			
// 		}
// 		#pragma omp barrier
		
// 		// for(j=0;j<a1.size();j++){
// 		// 	setds.insert(make_pair(a1.at(j), a2.at(j)));
// 		// }
		
// 		for(j=0;j<count1;j++){
// 			setds.insert(make_pair(b1[j], b2[j]));
// 		}
// 		delete[] b1;
// 		delete[] b2;	

// 	}
// 	end = omp_get_wtime();
// 	printf("Parallel Running time: %f ms\n", (end - start)*1000);
// 	// vector<vector<int>> vec(V);
// 	// Print shortest distances stored in dist[]
// 	// printf("Vertex Distance from Source\n");
// 	int j;
// 	for (int i = 0; i < V; ++i){
// 		//printf("%d \t\t %d\n", i, dist[i]);
// 		j = i+ 1;
// 		ofstream myfile ("USA_NY_parallel_out.txt", ios::app);
// 		//ofstream myfile ("tmp_parallel_out.txt", ios::app);
// 		myfile << j << " \t\t " << dist[i] << endl;
// 		myfile.close();		
// 		//j = i;
// 		//printf("%d  ", j);
// 		//while(parent[j]!=-1){
// 		//	vec[i].push_back(j);
// 		//	j = parent[j];
// 			//printf("%d  ", j);
// 		//}
// 		//printf("\n");
// 		//reverse(vec[i].begin(), vec[i].end());
// 		//for(int k=0; k< vec[i].size(); k++){
// 		//	printf("%d  ", vec[i][k]);
// 		//}
// 		//printf("\n");		
// 	}
// 	//path = vec[dest];
// 	//return_path(vec[dest]);
// 	/*vector<int> z = vec[dest];
// 	for(int it=0; it<path.size();it++){
// 		printf("%d  ", path[it]);
// 	}
// 	return z;*/
// }
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

// void erase(int start, int *a1, int *a2, int *ptr)
// {
// 	for(int i=start; i < *ptr; i++){
// 		a1[i] = a1[i+1];
// 		a2[i] = a2[i+1];
// 	}
// 	*ptr = *ptr - 1;
// }


void erase(int start, int *a1, int *a2, int index)
{
	for(int i=start; i < index + 2; i++){
		a1[i] = a1[i+1];
		a2[i] = a2[i+1];
	}
}


// void Graph::Dijkstras_Algo_omp(int src, int dest)
// {
// 	int vertices = this->V;
// 	int *dist = new int[vertices];
// 	for(int i=0; i<vertices; i++){
// 		dist[i] = INF;
// 	}
// 	dist[src] = 0;
// 	int n_thread,my_thread;

// 	int *a1 = new int[vertices];//conatins vertices with unverified min
// 	int *a2 = new int[vertices];//distance for above vertex
// 	int index = -1;//last valid index for above two
// 	index++;
// 	a1[index] = src;
// 	a2[index] = dist[src];
// 	//sort(a1,a2,index);
// 	std::vector<list< pair<int, int> >::iterator> element;
// 	int tmp;
// 	double start, end;
// 	start = omp_get_wtime();

// 	while(index>=0)
// 	{
// 		sort(a1, a2, index);
// 		int u = a1[0];
// 		//erase(0,a1,a2,&index);
// 		erase(0,a1,a2,index);
// 		index--;

// 		int j, v, weight;
// 		list< pair<int, int> >::iterator i;

// 		element.clear();
// 		for (i = adj[u].begin(); i != adj[u].end(); ++i)
// 		{
// 			element.push_back(i);
// 		}

// 		//i = element.at(0);

// 		//#pragma omp parallel for schedule(runtime) private(i, j, v, weight, tmp) shared(dist,u,a1,a2,index)
// 		for(j=0; j < element.size(); j++)
// 		{
// 			i = element.at(j);

// 			v = (*i).first;

// 			weight = (*i).second;

// 			if(dist[v] > dist[u] + weight)
// 			{
// 				if(dist[v] != INF)
// 				{
// 					//#pragma omp critical
// 					{
// 						tmp = search(a1, v, 0, index);
// 						//erase(tmp, a1, a2, &index);
// 						erase(tmp, a1, a2, index);
// 						index--;
// 					}
// 				}

// 				dist[v] = dist[u] + weight;
// 				//#pragma omp barrier

// 				//#pragma omp critical
// 				{
// 					index++;
// 					a1[index] = v;
// 					a2[index] = dist[v];
// 				}
// 				//#pragma omp barrier
				
// 				//#pragma omp critical
// 				{
// 					sort(a1, a2, index);
// 				}
// 			}

// 		}

// 	} 

// 	end = omp_get_wtime();
// 	printf("Parallel Running time: %f ms\n", (end - start)*1000);
	
// 	int j;
// 	for (int i = 0; i < V; ++i){
// 		j = i+ 1;
// 		ofstream myfile ("USA_NY_parallel_out.txt", ios::app);
// 		//ofstream myfile ("tmp_parallel_out.txt", ios::app);
// 		myfile << j << " \t\t " << dist[i] << endl;
// 		myfile.close();			
// 	}

// }


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

	#pragma omp parallel private(weight, k, first, last, thread, step, min_d, min_node) shared(connected,d,dist,v,n_thread,adj,element,check)
	{
		thread = omp_get_thread_num();
		n_thread = omp_get_num_threads();
		first = (thread*vertices)/n_thread;
		last = ((thread+1)*vertices)/n_thread;

		for(step=1;step<vertices;step++)
		//while(check)
		//for(step=first;step<last;step++)
		{
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

			if(v != -1)
			{
				for (int i = first; i <= last; i++ )
				{
					if ( !connected[i] )
					{
						int j;
						//for(k = adj[v].begin(); k != adj[v].end(); ++k)
						for(j = 0; j < element[v].size(); ++j)
						{
							k = element[v].at(j);
							if(i == (*k).first)
							{
								weight = (*k).second;
								if(weight < INF)
								{
									if(dist[i] > dist[v] + weight)
									{
										dist[i] = dist[v] + weight;
									}
								}
							}
						}
					  // if ( ohd[mv][i] < i4_huge )
					  // {
					  //   if ( mind[mv] + ohd[mv][i] < mind[i] )  
					  //   {
					  //     mind[i] = mind[mv] + ohd[mv][i];
					  //   }
					  // }
					}
				}				
			}

			#pragma omp barrier

			// check = false;
			// #pragma omp single
			// {
			// 	for(int i=0; i< vertices;i++)
			// 	{
			// 		if(connected[i] = false)
			// 		{
			// 			check = true;
			// 		}
			// 	}
			// }

			// #pragma omp barrier

		}
	}

	end = omp_get_wtime();
	printf("Parallel Running time: %f ms\n", (end - start)*1000);	

	delete[] connected;

	int j;
	for (int i = 0; i < V; ++i){
		j = i+ 1;
		ofstream myfile ("USA_NY_parallel_out.txt", ios::app);
		//ofstream myfile ("tmp_parallel_out.txt", ios::app);
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

	return 0;
	*/
// 	fstream disFile;
// 	disFile.open("../USA-road-d_NY.txt");
	
// 	vector<vector<int>> vect;
	
// 	if (disFile.is_open()){ 
// 		string tp;
// //		int count =15;	
// 		while(getline(disFile, tp)){ 
// //			getline(disFile, tp);
// 			vector<int> pair;
// 			if(tp[0] == 'a'){
// 				int a;
// 				string str(tp.begin()+1, tp.begin()+tp.size() );
// 				stringstream ss (str);
// 				while ((ss >> a))
// 					pair.push_back (a);
// 				vect.push_back(pair);
// 				//cout<<vect.size();
// 			}
// //			count --;
// 		}
// 		disFile.close(); 
//    }
   
// //    for (int i = 0; i < vect.size(); i++) {
// //        for (int j = 0; j < vect[i].size(); j++)
// //            cout << vect[i][j] << " ";
// //        cout << endl;
// //    }

// 	Graph g(264346);

// 	cout<<"Number of edges: "<<vect.size()<<'\n';
	
//     for (int i = 0; i < vect.size(); i++) {
//     	g.addEdge(vect[i][0]-1,vect[i][1]-1,vect[i][2]);
//     }
	
// 	g.Dijkstras_Algo(0, 100);
// 	g.Dijkstras_Algo_omp(0, 100);
	
	// int V = 25;
	// Graph g(V);
	// int i,j;
	// for(int x=0; x<5; x++){
	// 	for(int y=0; y<5;y++){
	// 		if(g.adj[(5*y+x)].empty()==true){
	// 		if(y==0){
	// 			if(x==0){
	// 				g.addEdge(0,1,1);
	// 				g.addEdge(0,5,1);
	// 			}
	// 			else if(x==4){
	// 				g.addEdge(4,3,1);
	// 				g.addEdge(4,9,1);					
	// 			}
	// 			else{
	// 				i = 5*y+x;
	// 				j = 5*(y+1)+x;
	// 				g.addEdge(i,j,1);
	// 				j = i+1;
	// 				g.addEdge(i,j,1);
	// 				j = i-1;
	// 				g.addEdge(i,j,1);
	// 			}
	// 		}
	// 		else if(y==4){
	// 			if(x==0){
	// 				g.addEdge(20,21,1);
	// 				g.addEdge(20,15,1);
	// 			}
	// 			else if(x==4){
	// 				g.addEdge(24,23,1);
	// 				g.addEdge(24,19,1);					
	// 			}
	// 			else{
	// 				i = 5*y+x;
	// 				j = 5*(y-1)+x;
	// 				g.addEdge(i,j,1);
	// 				j = i+1;
	// 				g.addEdge(i,j,1);
	// 				j = i-1;
	// 				g.addEdge(i,j,1);
	// 			}				
	// 		}
	// 		else if(x==0){
	// 			i = 5*y+x;
	// 			j = 5*(y-1)+x;
	// 			g.addEdge(i,j,1);
	// 			j = 5*(y+1)+x;
	// 			g.addEdge(i,j,1);
	// 			j = i+1;
	// 			g.addEdge(i,j,1);				
	// 		}
	// 		else if(x==4){
	// 			i = 5*y+x;
	// 			j = 5*(y-1)+x;
	// 			g.addEdge(i,j,1);
	// 			j = 5*(y+1)+x;
	// 			g.addEdge(i,j,1);
	// 			j = i-1;
	// 			g.addEdge(i,j,1);				
	// 		}
	// 		else{
	// 			i = 5*y+x;
	// 			//if(g.adj[i].empty()==true){
	// 			j = 5*(y-1)+x;
	// 			g.addEdge(i,j,1);
	// 			j = 5*(y+1)+x;
	// 			g.addEdge(i,j,1);
	// 			j = i-1;
	// 			g.addEdge(i,j,1);
	// 			j = i+1;
	// 			g.addEdge(i,j,1);
	// 			//}	
	// 		}
	// 		}			
	// 	}
	// }

	// //g.Dijkstras_Algo(0, 5);
	// g.Dijkstras_Algo_omp(0, 5);

	int V = 10000;
	int E = 50000;
	Graph g(V);
	g.addEdge(0,2,1);
	for (int j = 1; j < E; j++) {
		g.addEdge(rand()%(V) + 0, rand()%(V) + 0, rand()%(50 - 1 + 1) + 1);
	}
	g.Dijkstras_Algo(0, 5);	
	g.Dijkstras_Algo_omp(0, 5);

	return 0;
}
