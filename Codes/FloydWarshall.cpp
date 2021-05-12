#include <omp.h>
#include <vector>
#include <fstream>
#include <iostream>
#include <bits/stdc++.h>

#define INFTY   1e8

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
void Graph::Dijkstras_Algo(int src, int dest)
{


}
void Graph::BellmanFord_Algo(int src)
{


}
void Graph::FloydWarshall_Algo(int src)
{
  
}

int main(int argv, char **argc)
{
	/*
	fstream disFile;
	disFile.open("./USA-road-d_NY.txt");
	
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
	*/
    
	int V = strtol(argc[1], (char **)NULL, 10), E = 3*V;
	Graph g(V);
	g.addEdge(0,2,1);
	for (int j = 1; j < E; j++) {
		g.addEdge(rand()%(V) + 0, rand()%(V) + 0, rand()%(50 - 1 + 1) + 1);
	}
	/*
	// Code for generating predefined edges in a 25-vertex graph -- Too short for comparison
	int i,j;
	for(int x=0; x<5; x++){
		for(int y=0; y<5;y++){
			if(g.adj[(5*y+x)].empty()==true){
			if(y==0){
				if(x==0){
					g.addEdge(0,1,1);
					g.addEdge(0,5,1);
				}
				else if(x==4){
					g.addEdge(4,3,1);
					g.addEdge(4,9,1);
				}
				else{
					i = 5*y+x;
					j = 5*(y+1)+x;
					g.addEdge(i,j,1);
					j = i+1;
					g.addEdge(i,j,1);
					j = i-1;
					g.addEdge(i,j,1);
				}
			}
			else if(y==4){
				if(x==0){
					g.addEdge(20,21,1);
					g.addEdge(20,15,1);
				}
				else if(x==4){
					g.addEdge(24,23,1);
					g.addEdge(24,19,1);
				}
				else{
					i = 5*y+x;
					j = 5*(y-1)+x;
					g.addEdge(i,j,1);
					j = i+1;
					g.addEdge(i,j,1);
					j = i-1;
					g.addEdge(i,j,1);
				}
			}
			else if(x==0){
				i = 5*y+x;
				j = 5*(y-1)+x;
				g.addEdge(i,j,1);
				j = 5*(y+1)+x;
				g.addEdge(i,j,1);
				j = i+1;
				g.addEdge(i,j,1);
			}
			else if(x==4){
				i = 5*y+x;
				j = 5*(y-1)+x;
				g.addEdge(i,j,1);
				j = 5*(y+1)+x;
				g.addEdge(i,j,1);
				j = i-1;
				g.addEdge(i,j,1);
			}
			else{
				i = 5*y+x;
				//if(g.adj[i].empty()==true){
				j = 5*(y-1)+x;
				g.addEdge(i,j,1);
				j = 5*(y+1)+x;
				g.addEdge(i,j,1);
				j = i-1;
				g.addEdge(i,j,1);
				j = i+1;
				g.addEdge(i,j,1);
				//}
			}
			}
		}
	}
	*/
	/*
	// Code for importing NYC Road Map -- Impossible to solve in a lifetime
	Graph g(264346);
	for (int i = 0; i < vect.size(); i++) {
		g.addEdge(vect[i][0]-1,vect[i][1]-1,vect[i][2]);
	}
	cout<<"Number of edges: "<<vect.size()<<'\n';
	*/
	
	// Floyd-Warshall Algorithm
	int numThreads = strtol(argc[2], (char **)NULL, 10), threadNum;
	int start, end, mid;
	// vector< pair<int,int> > xStart, xMid, xEnd;
	omp_set_num_threads(numThreads);
	
	double t1 = omp_get_wtime();
	// #pragma omp parallel for private(start, end) shared(g)
	for (mid = 0; mid < g.vertNum(); mid++)
    {
		double tM1 = omp_get_wtime();
		threadNum = omp_get_thread_num();
		
		vector< pair<int,int> > xMid;
		xMid.clear();
		xMid.reserve(g.adj[mid].size());
		copy(g.adj[mid].begin(), g.adj[mid].end(), xMid.begin());
		
		#pragma omp parallel for private(start) shared(g)
		for (start = 0; start < g.vertNum(); start++)
		{
			double tS1 = omp_get_wtime();
			
			vector< pair<int,int> > xStart;
			xStart.clear();
			xStart.reserve(g.adj[start].size());
			copy(g.adj[start].begin(), g.adj[start].end(), xStart.begin());
			
			#pragma omp parallel for private(end)
			for (end = 0; end < g.vertNum(); end++)
			{
				double tE1 = omp_get_wtime();
				
				vector< pair<int,int> > xEnd;
				xEnd.clear();
				xEnd.reserve(g.adj[end].size());
				copy(g.adj[end].begin(), g.adj[end].end(), xEnd.begin());
				
				int  posStartEnd = -1, posStartMid = -1, posMidEnd = -1;
				int  valStartEnd,      valStartMid,      valMidEnd;
				int doneStartEnd,     doneStartMid,     doneMidEnd;
				
				for (auto itStart = xStart.begin(); itStart != xStart.end(); ++itStart)
				{
					if (doneStartEnd == 0)
					{
						posStartEnd += 1;
						valStartEnd = (*itStart).second;
					}
					if (doneStartMid == 0)
					{
						posStartMid += 1;
						valStartMid = (*itStart).second;
					}
					if ((*itStart).first == end)
						doneStartEnd = 1;
					if ((*itStart).first == mid)
						doneStartMid = 1;
					if (doneStartEnd == 1 && doneStartMid == 1)
						break;
				}
				
				for (auto itMid = xMid.begin(); itMid != xMid.end(); ++itMid)
				{
					if (doneMidEnd == 0)
					{
						posMidEnd += 1;
						valMidEnd = (*itMid).second;
					}
					if ((*itMid).first == end)
						doneMidEnd = 1;
					if (doneMidEnd == 1)
						break;
				}
				pair<int,int> replaceNode = make_pair(posStartEnd, valStartEnd);
				(g.adj)[start].remove(replaceNode);
				(g.adj)[start].push_back(make_pair(posStartEnd, min(valStartEnd, valStartMid + valMidEnd)));
				
				double tE2 = omp_get_wtime() - tE1;
				// cout<<"      "<<end<<": "<<tE2<<endl;
			}
			
			double tS2 = omp_get_wtime() - tS1;
			// cout<<"    "<<start<<": "<<tS2<<endl;
		}
		double tM2 = omp_get_wtime() - tM1;
		// cout<<"  "<<mid<<": "<<tM2<<endl;
    }
	double t2 = omp_get_wtime() - t1;
	cout<<"Total Time taken for \t"<<V<<" vertices = "<<t2<<" seconds."<<endl;
	
	return 0;
}