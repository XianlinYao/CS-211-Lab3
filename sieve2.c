/*
 *   Sieve of Eratosthenes
 *
 *   Programmed by Michael J. Quinn
 *
 *   Last modification: 7 September 2001
 */

#include "mpi.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define MIN(a,b)  ((a)<(b)?(a):(b))

int main(int argc, char **argv)
{
    double elapsedTime;
    int pid;
    int psize;
    unsigned long long int n;
    unsigned long long int low_value, high_value;
    unsigned long long int size, proc0_size;
    char *marked;
    unsigned long long int i, index, first;
    unsigned long long int prime;
    unsigned long long int count, global_count;
    unsigned long long int local_first, local_prime, local_size;
    char *local_marked;

    /**
     * Initialize MPI
     */
    MPI_Init(&argc, &argv);
    MPI_Barrier(MPI_COMM_WORLD);
    elapsedTime = -MPI_Wtime();
    MPI_Comm_rank(MPI_COMM_WORLD, &pid);
    MPI_Comm_size(MPI_COMM_WORLD, &psize);
    if (argc != 2)
    {
        if (pid == 0)
            printf("Command line syntax error.\n");
        MPI_Finalize();
        exit(1);
    }

    /**
     * Parameters Initialization
     */
    n = atoll(argv[1]);
    low_value = 2 + BLOCK_LOW(pid, psize, n-1);
    high_value = 2 + BLOCK_HIGH(pid, psize, n-1);
    // size = BLOCK_SIZE(pid, psize, n-1);
    low_value = low_value + (low_value + 1) % 2;
    high_value = high_value - (high_value + 1) % 2;
    size = (high_value - low_value) / 2 + 1;
    local_size  = (int)sqrt((double)(n)) - 1;
    
    /**
     * process 0 must holds all primes used
     */
    proc0_size = (n/2 - 1) / psize;
    if ((2 + proc0_size) < (int) sqrt((double) n/2))
    {
        if (pid == 0)
            printf("Too many processes.\n");
        MPI_Finalize();
        exit(1);
    }

    /**
     * Allocation
     */
    marked = (char*) malloc(size);
    local_marked = (char *) malloc (local_size);
    if (marked == NULL || local_marked == NULL)
    {
        printf("PID: %d - Cannot allocate enough memory.\n", pid);
        MPI_Finalize();
        exit(1);
    }

    /**
     * Core Function
     */
    local_prime = 2;
    for (i = 0; i < local_size; i++)
        local_marked[i] = 0;
    index = 0;
    do
    {
        local_first = local_prime * local_prime - 2;
        for (i = local_first; i < local_size; i += local_prime)
            local_marked[i] = 1;
        while (local_marked[++index] == 1);
        local_prime = 2 + index;
    } while (local_prime * local_prime <= n);
    

    for (i = 0; i < size; i++)
        marked[i] = 0;
    index = 0;
    prime = 3;
    do
    {
        if (prime * prime > low_value)
            first = (prime * prime - low_value) / 2;
        else
        {
            if ((low_value % prime) == 0)
                first = 0;
            else
                first = (prime - (low_value % prime) + low_value / prime % 2 * prime) / 2;
        }
        for (i = first; i < size; i += prime)
            marked[i] = 1;
        while(local_marked[++index] == 1);
        prime = index + 2;
    } while (prime * prime <= n);
    count = 0;
    for (i = 0; i < size; i++)
        if (marked[i] == 0)
            count++;
    if (pid == 0)
        count++;    // 2
    if (psize > 1)
        MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    
    elapsedTime += MPI_Wtime();

    /**
     * Print the results
     */
    if (pid == 0)
        printf("The total number of primes: %lld, time: %.6lf s in %d processes.\n", global_count, elapsedTime, psize);
    MPI_Finalize();
    return 0;
}
