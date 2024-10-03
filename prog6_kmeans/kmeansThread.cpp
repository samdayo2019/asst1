#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <thread>
#include <fstream> //included to allow for creation of file to input profiled data over multiple iterations.

#include "CycleTimer.h"

using namespace std;

#define THREADS 8 //hardcode max threads to 8: quad-core + hyper-threading

typedef struct {
  // Control work assignments
  int start, end;

  // Shared by all functions
  double *data;
  double *clusterCentroids;
  int *clusterAssignments;
  double *currCost;
  int M, N, K;

  int threadID; // add new member for threadID.
  int starting_dp;
} WorkerArgs;


/**
 * Checks if the algorithm has converged.
 * 
 * @param prevCost Pointer to the K dimensional array containing cluster costs 
 *    from the previous iteration.
 * @param currCost Pointer to the K dimensional array containing cluster costs 
 *    from the current iteration.
 * @param epsilon Predefined hyperparameter which is used to determine when
 *    the algorithm has converged.
 * @param K The number of clusters.
 * 
 * NOTE: DO NOT MODIFY THIS FUNCTION!!!
 */
static bool stoppingConditionMet(double *prevCost, double *currCost,
                                 double epsilon, int K) {
  for (int k = 0; k < K; k++) {
    if (abs(prevCost[k] - currCost[k]) > epsilon)
      return false;
  }
  return true;
}

/**
 * Computes L2 distance between two points of dimension nDim.
 * 
 * @param x Pointer to the beginning of the array representing the first
 *     data point.
 * @param y Poitner to the beginning of the array representing the second
 *     data point.
 * @param nDim The dimensionality (number of elements) in each data point
 *     (must be the same for x and y).
 */
double dist(double *x, double *y, int nDim) {
  double accum = 0.0;
  for (int i = 0; i < nDim; i++) {
    accum += pow((x[i] - y[i]), 2);
  }
  return sqrt(accum);
}

/**
 * Assigns each data point to its "closest" cluster centroid.
 * N = number of dimensions (features) for each data point. 
 * M = total number of data points (samples) 
 * K = total number of clusters (cluster centres)
 */
void computeAssignments(WorkerArgs *const args) {
  double *minDist = new double[args->M]; // minDist array should match number of datapoints processed which is m

  // Initialize arrays
  /*
    These array initializations are very simple operations that can be completed in constant-time. 
    As a result, their overall time complexity is linear O(M). Creating threads to parallelize these operations, and the
      associated overhead in synchronizing them (in joining), spawning them, and starting them would likely take longer 
      than computing these operations serially. This also applies to having two seperate threads initialize the seperate 
      arrays in parallel.
  */
  for (int m = 0; m < args->M; m++) {
    minDist[m] = 1e30;
    args->clusterAssignments[args->starting_dp + m] = -1;
  }
  // Assign datapoints to closest centroids

  for (int k = args->start; k < args->end; k++) {
    for (int m = 0; m < args->M; m++) {
      double d = dist(&args->data[(args->starting_dp + m)* args->N],
                      &args->clusterCentroids[k * args->N], args->N);
      if (d < minDist[m]) {
        minDist[m] = d;
        args->clusterAssignments[args->starting_dp + m] = k;
      }
    }
  }

  free(minDist);
}

// function hard-coded for 8 threads.
void computeAssignmentsWithThreads(WorkerArgs *const args){
  std::thread worker_threads[THREADS]; // init 8 threads
  WorkerArgs work_args[THREADS];  // init 8 args structs for each thread.

  int work_split = args->M / THREADS;

  // Update per-thread arguments. Simply copying over args and adding thread id member.
  // all threads point to same clusterAssignment, clusterCentroids, currCost arrays, and data arrays.
  for(int i = 0; i < THREADS; i++){
    work_args[i].threadID = i;
    work_args[i].clusterAssignments = args->clusterAssignments; 
    work_args[i].clusterCentroids = args->clusterCentroids;
    work_args[i].currCost = args->currCost;
    work_args[i].data = args->data;
    work_args[i].start = args->start;
    work_args[i].end = args->end;
    work_args[i].K = args->K; // each thread processes all k clusters
    work_args[i].N = args->N; // number of features remains unchanged as dimensionality required for euclidean distance.

    work_args[i].starting_dp = work_args[i].threadID*work_split;

    work_args[i].M = (work_args[i].threadID == THREADS - 1) ? args->M - work_args[i].starting_dp : work_split; // we need to split the M among clusters. We begin with simple static decomp

  }

  for(int i = 1; i < THREADS; i++){
    worker_threads[i] = std::thread(computeAssignments, &work_args[i]);
  }


  computeAssignments(&work_args[0]);

  for(int i = 1; i < THREADS; i++){
    worker_threads[i].join();
  }

}

/**
 * Given the cluster assignments, computes the new centroid locations for
 * each cluster.
 */
void computeCentroids(WorkerArgs *const args) {
  int *counts = new int[args->K];

  // Zero things out
  for (int k = 0; k < args->K; k++) {
    counts[k] = 0;
    for (int n = 0; n < args->N; n++) {
      args->clusterCentroids[k * args->N + n] = 0.0;
    }
  }


  // Sum up contributions from assigned examples
  for (int m = 0; m < args->M; m++) {
    int k = args->clusterAssignments[m];
    for (int n = 0; n < args->N; n++) {
      args->clusterCentroids[k * args->N + n] +=
          args->data[m * args->N + n];
    }
    counts[k]++;
  }

  // Compute means
  for (int k = 0; k < args->K; k++) {
    counts[k] = max(counts[k], 1); // prevent divide by 0
    for (int n = 0; n < args->N; n++) {
      args->clusterCentroids[k * args->N + n] /= counts[k];
    }
  }

  free(counts);
}

/**
 * Computes the per-cluster cost. Used to check if the algorithm has converged.
 */
void computeCost(WorkerArgs *const args) {
  double *accum = new double[args->K];

  // Zero things out
  for (int k = 0; k < args->K; k++) {
    accum[k] = 0.0;
  }

  // Sum cost for all data points assigned to centroid
  for (int m = 0; m < args->M; m++) {
    int k = args->clusterAssignments[m];
    accum[k] += dist(&args->data[m * args->N],
                     &args->clusterCentroids[k * args->N], args->N);
  }

  // Update costs
  for (int k = args->start; k < args->end; k++) {
    args->currCost[k] = accum[k];
  }

  free(accum);
}

/**
 * Computes the K-Means algorithm, using std::thread to parallelize the work.
 *
 * @param data Pointer to an array of length M*N representing the M different N 
 *     dimensional data points clustered. The data is layed out in a "data point
 *     major" format, so that data[i*N] is the start of the i'th data point in 
 *     the array. The N values of the i'th datapoint are the N values in the 
 *     range data[i*N] to data[(i+1) * N].
 * @param clusterCentroids Pointer to an array of length K*N representing the K 
 *     different N dimensional cluster centroids. The data is laid out in
 *     the same way as explained above for data.
 * @param clusterAssignments Pointer to an array of length M representing the
 *     cluster assignments of each data point, where clusterAssignments[i] = j
 *     indicates that data point i is closest to cluster centroid j.
 * @param M The number of data points to cluster.
 * @param N The dimensionality of the data points.
 * @param K The number of cluster centroids.
 * @param epsilon The algorithm is said to have converged when
 *     |currCost[i] - prevCost[i]| < epsilon for all i where i = 0, 1, ..., K-1
 */
void kMeansThread(double *data, double *clusterCentroids, int *clusterAssignments,
               int M, int N, int K, double epsilon) {

  // Used to track convergence
  double *prevCost = new double[K];
  double *currCost = new double[K];

  // The WorkerArgs array is used to pass inputs to and return output from
  // functions.
  WorkerArgs args;
  args.data = data;
  args.clusterCentroids = clusterCentroids;
  args.clusterAssignments = clusterAssignments;
  args.currCost = currCost;
  args.M = M;
  args.N = N;
  args.K = K;
  args.starting_dp = 0;
  

  // Initialize arrays to track cost
  for (int k = 0; k < K; k++) {
    prevCost[k] = 1e30;
    currCost[k] = 0.0;
  }

  float start_assignment = 0;
  float stop_assignment = 0;
  float start_centroids = 0;
  float stop_centroids  = 0;
  float start_cost = 0;
  float stop_cost = 0;

  float diff_assignment = 0; 
  float diff_centroids = 0; 
  float diff_cost = 0; 
  
  ofstream outputFile("Average_Profiled_Elapsed_Times.csv");

  outputFile << "Assignment" << ","  << "Centroids" << "," << "Cost" << endl;

  /* Main K-Means Algorithm Loop */
  int iter = 0;
  while (!stoppingConditionMet(prevCost, currCost, epsilon, K)) {
    // Update cost arrays (for checking convergence criteria)
    for (int k = 0; k < K; k++) {
      prevCost[k] = currCost[k];
    }

    // Setup args struct
    args.start = 0;
    args.end = K;


    // embedd profiling code to time execution time of each of the three functions
    start_assignment = CycleTimer::currentSeconds();
    computeAssignmentsWithThreads(&args);
    stop_assignment = CycleTimer::currentSeconds();

    start_centroids = CycleTimer::currentSeconds();
    computeCentroids(&args);
    stop_centroids = CycleTimer::currentSeconds();

    start_cost = CycleTimer::currentSeconds();
    computeCost(&args);
    stop_cost = CycleTimer::currentSeconds();

    // printf("%d ", iter);

    iter++;
    diff_assignment += (stop_assignment - start_assignment);
    diff_centroids += (stop_centroids - start_centroids);
    diff_cost += (stop_cost - start_cost);

  }

  diff_assignment /= iter; 
  diff_centroids /= iter; 
  diff_cost /= iter; 

  outputFile << diff_assignment << "," << diff_centroids << "," << diff_cost << endl;
  outputFile.close(); // finish writing into the output file. Will be opened and processed by analysis.py

  printf("Elapse time for Assignment function: %f.\nElapsed time for Centroids function: %f.\nElapsed time for Cost function: %f.\nIter was: %d\n", diff_assignment, diff_centroids, diff_cost, iter);

  free(currCost);
  free(prevCost);
}
