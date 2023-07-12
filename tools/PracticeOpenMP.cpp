#include <stdio.h>
#include <iostream>
#include <omp.h>

using namespace std;

static long num_steps = 1000000000;
double step;

#define NUM_THREADS 4

int main()
{
    int i, nthreads;
    double x, pi = 0.0, sum = 0.0;

    step = 1.0/(double)num_steps;
    double start_time = omp_get_wtime();

    omp_set_num_threads(NUM_THREADS);

    
    #pragma omp parallel
    {   
        #pragma omp for reduction(+:sum) schedule(auto)
        for (int i = 0 ; i < num_steps; i++ )
        {
            x = (i + 0.5) * step;
            sum += 4.0/(1.0 + x*x);
        }
        
    }
    pi+= sum * step; 

    double end_time = omp_get_wtime();
    double running_time = end_time - start_time;

    cout << "Approximate value of pi: " << pi << endl;
    cout << "Running time: " << running_time << " seconds" << endl;

    return 0;
}
