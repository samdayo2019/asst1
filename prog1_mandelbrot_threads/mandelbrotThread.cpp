#include <stdio.h>
#include <thread>

#include "CycleTimer.h"

typedef struct {
    float x0, x1;
    float y0, y1;
    unsigned int width;
    unsigned int height;
    int maxIterations;
    int* output;
    int threadId;
    int numThreads;
} WorkerArgs;


extern void mandelbrotSerial(
    float x0, float y0, float x1, float y1,
    int width, int height,
    int startRow, int numRows,
    int maxIterations,
    int output[]);


//
// workerThreadStart --
//
// Thread entrypoint.
void workerThreadStart(WorkerArgs * const args) {

    // TODO FOR CS149 STUDENTS: Implement the body of the worker
    // thread here. Each thread should make a call to mandelbrotSerial()
    // to compute a part of the output image.  For example, in a
    // program that uses two threads, thread 0 could compute the top
    // half of the image and thread 1 could compute the bottom half.

    double start = CycleTimer::currentSeconds();

    //-------PART 1.1------------//
    

    // /*
    //     The height field controls the number of rows we process in the image
    //     The startrows field controls where in the image we begin
    //     Those two variables need to be adjusted in order to select a subset of the rows (ie., top vs bottom half) 
    //         and to assign the different portions to the different threads. 
    //     The starting row should be tied to the thread id number.
        
    //     1. Calculate number of rows we want to split into based on the number of threads. 
    //     2. Determine the starting row we for each thread to begin at.
    //     3. Ensure that we do not exceed the number of threads for the last row.
    // */

    // // insert timing code to measure execution time for each thread

    // int num_rows_split = args->height/args->numThreads; // splitting image up
    // int starting_row = args->threadId*num_rows_split; // offset image start

    // // account for non-equal splits (when height not divisible by num threads) for last thread.
    // int number_rows_process = 0;
    // if(args->threadId == args->numThreads-1){
    //     number_rows_process = args->height - starting_row; 
    // }
    // else{
    //     number_rows_process = num_rows_split;
    // }

    // // call MandelbrotSerial using the updated row and heights
    // mandelbrotSerial(
    // args->x0, args->y0, args->x1, args->y1,
    // args->width, args->height,
    // starting_row, number_rows_process,
    // args->maxIterations,
    // args->output);

    // replace with printf at the end 
    //printf("Hello world from thread %d. Number of rows computed %d. Starting row %d. Execution took %f\n", args->threadId, number_rows_process, starting_row, end-start);


    //----------------------------------//

    //-------PART 1.4-------------//

    /*
    To improve performance speed up without added synchronization, we should employ a work decomp stragey that distributes the work more evenly.
        We can assign each thread rows that are not contiguous, such that on average the required number of iterations to solve per thread is more uniform.

    Here we implement a round-robin decomposition policy based on the rows. 
    */

    int starting_row = args->threadId; 
    int rows_computed = 0;
    while(unsigned(starting_row) < args->height){
        // continuously call mandelbrot one row at a time so we can do non-contiguous rows
        mandelbrotSerial(
        args->x0, args->y0, args->x1, args->y1,
        args->width, args->height,
        starting_row, 1,
        args->maxIterations,
        args->output);

        rows_computed++;
        starting_row += args->numThreads;
    }


    /*
    Is it possible to further improve performance by separating the image into blocks instead to be processed instead of rows?
    */

    double end = CycleTimer::currentSeconds();
    
    printf("Hello world from thread %d. Number of rows computed %d. Execution took %f\n", args->threadId, rows_computed, end-start);

}

//
// MandelbrotThread --
//
// Multi-threaded implementation of mandelbrot set image generation.
// Threads of execution are created by spawning std::threads.
void mandelbrotThread(
    int numThreads,
    float x0, float y0, float x1, float y1,
    int width, int height,
    int maxIterations, int output[])
{
    static constexpr int MAX_THREADS = 32;

    if (numThreads > MAX_THREADS)
    {
        fprintf(stderr, "Error: Max allowed threads is %d\n", MAX_THREADS);
        exit(1);
    }

    // Creates thread objects that do not yet represent a thread.
    std::thread workers[MAX_THREADS];
    WorkerArgs args[MAX_THREADS];

    for (int i=0; i<numThreads; i++) {
      
        // TODO FOR CS149 STUDENTS: You may or may not wish to modify
        // the per-thread arguments here.  The code below copies the
        // same arguments for each thread
        args[i].x0 = x0;
        args[i].y0 = y0;
        args[i].x1 = x1;
        args[i].y1 = y1;
        args[i].width = width;
        args[i].height = height;
        args[i].maxIterations = maxIterations;
        args[i].numThreads = numThreads;
        args[i].output = output;
      
        args[i].threadId = i;
    }

    // Spawn the worker threads.  Note that only numThreads-1 std::threads
    // are created and the main application thread is used as a worker
    // as well.
    for (int i=1; i<numThreads; i++) {
        workers[i] = std::thread(workerThreadStart, &args[i]);
    }
    
    workerThreadStart(&args[0]);

    // join worker threads
    for (int i=1; i<numThreads; i++) {
        workers[i].join();
    }
}

