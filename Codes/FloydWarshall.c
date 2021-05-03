#include <omp.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define R       6371000.0
#define PI      3.1415926535897932
#define INFTY   1e20
#define deg2rad PI/180.0

#ifndef min
#define min(a,b)  (((a) < (b)) ? (a) : (b))
#endif

int     numStates     =  35;
int     numDistricts  = 594;
float  *Lattitudes    = {0};
float  *Longitudes    = {0};
float **Distances     = {0};

// State capitals
int **capitals = {  {1, 0},   {1, 8},  {1, 32},  {1, 51},  {1, 87}, {1, 100}, {1, 114}, {1, 117}, {1, 118},
                  {1, 120}, {1, 121}, {1, 130}, {1, 100}, {2, 170,      175}, {2, 183,      191}, {1, 211},
                  {1, 217}, {1, 253}, {1, 256}, {1, 263}, {1, 316}, {2, 342,      347}, {1, 349}, {1, 355},
                  {1, 364}, {1, 389}, {4, 401,      402,      403,      404}, {1, 100}, {1, 438}, {1, 454},
                  {1, 459}, {1, 491}, {1, 537}, {1, 566}, {1, 585}};

// Starting district number of each state
int **districtStarts = {  0,   2,  25,  40,  63, 100, 101,
                        117, 118, 120, 121, 123, 148, 167,
                        179, 193, 215, 242, 256, 257, 305,
                        339, 348, 355, 363, 371, 401, 405,
                        422, 454, 458, 488, 492, 562, 575}

// Check if a value is present in an array of length array[0] starting from position 1
bool valInArr(int value, int *array)
{
  for(int i = 1; i <= array[0]; i++) 
    if(array[i] == value) 
      return 1; 
  return 0;
}

// Calculate distance between two locations given their lattitudes and longitudes - Haversian formula
float calcDist(float latt1, float long1, float latt2, float long2)
{
  float a = pow(sin((latt2-latt1)*deg2rad/2), 2) + cos(latt1*deg2rad)*cos(latt2*deg2rad)*pow(sin((long2-long1)*deg2rad/2), 2);
  return 2*R*atan2(sqrt(1-a), sqrt(a));
}

// Find state number for a particular district
int state(int district)
{
  for (int state = 1; state < numStates; state++)
    if (district < districtStarts[state])
      return state-1;
  return numStates-1;
}

// Check if two districts are connected
int connected(int start, int end)
{
  if (valInArr(start, capitals[state(start)]))
    if (valInArr(state(end), neighbors[state(start)]))
      return 1;
  if (valInArr(end, capitals[state(end)-1]))
    if (valInArr(state(end), neighbors[state(start)]))
      return 1;
  return 0;
}

// Import to Lattitudes and Longitudes from a CSV file
void dataImport()
{
  FILE *fp = fopen("IndiaDistricts.csv", "r");
  
  if (!fp)
    printf("Can't open file\n");
  
  else
  {
    char buffer[1024];
    
    int district = -1;
    int datatype = 0;
  
    while (fgets(buffer, 1024, fp))
    {
      datatype = 0;
      district++;
      char* value = strtok(buffer, ",");
      
      while (value)
      {
        if (datatype == 2)
          Lattitudes[district] = value;
        if (datatype == 3)
          Longitudes[district] = value;
        datatype++;
      }
      fclose(fp);
    }
}

int main(int argc, char **argv)
{
  int numThreads = strtol(argv[1], (char **)NULL, 10), threadNum;
  int start, end, mid;
  
  // Compute Distances
  dataImport();
  for (start = 0; start < N; start++)
    for (end = 0; end < N; end++)
      if (connected(start, end) == 0)
        Distances[start][end] = INFTY;
      else if(start != end)
        Distances[start][end] = calcDist(Lattitudes[start], Longitudes[start], Lattitudes[end], Longitudes[end]);
      else
        Distances[start][end] = 0.0;
  
  // Start Timing
  double t1 = omp_get_wtime();
  
  // Set number of threads for OpenMP
  omp_set_num_threads(numThreads);
  
  // Perform Floyd-Warshall Algorithm
  for(threadNum = 1; threadNum <= 10; threadNum++)
  {
    #pragma omp parallel shared(Distances)
    for (mid = 0; mid < N; mid++)
    {
      float *distMid = Distances[mid];
      #pragma omp parallel for private(start, end) schedule(dynamic)
      for (start = 0; start < N; start++)
      {
        int *distStart = Distances[start];
        for (end = 0; end < N; end++)
          distStart[end] = min(distStart[end], distStart[mid] + distMid[end]);
      }
    }
  }
  // Stop Timing
  double t2 = omp_get_wtime() - t1;
  
  // Exit
  return 0;
}