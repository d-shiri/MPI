// NOTE1: compile the file with -lm -ldl flags
// NOTE2: run the file as mpirun -np 4 a.out n.txt

#include <stdio.h>
#include <math.h>
#include "mpi.h"

int main(int argc, char* argv[]){
    int p, rank, n, tag;
    int dest, source;
    double h, start, finish, integral, integral_sum;
    int section;
    MPI_Status status;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0){
        FILE *fptr = fopen(argv[1], "r");
        if (fptr == NULL){printf("can't open the file");}
        else if (fptr){fscanf(fptr, "%d", &n);} 
        else {n = 16;}
        tag = 0;
        for (dest = 1; dest < p; ++dest){ // distribute to rest
            MPI_Send(&n, 1, MPI_INT, dest, tag, MPI_COMM_WORLD);
        }

    }
    else{ // get n from process 0
        tag = 0; 
        source = 0;
        MPI_Recv(&n, 1, MPI_INT, source, tag, MPI_COMM_WORLD, &status);

    }
    // at this point all processes have n
    int a, b;//integral bounds
    a = 0; b = 1;
    h = 1.0/n;
    section = n/p;
    start = rank * section;
    finish = (rank + 1) * section;
    double Trap(int start, int finish, double h);
    
    integral = Trap(start, finish, h);
    if (rank != 0) { 
        tag = 0;
        MPI_Send(&integral, 1, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD); 
    } else {
        integral_sum = integral;

        for (source = 1; source < p; source++) {
            MPI_Recv(&integral, 1, MPI_DOUBLE, source, 0, MPI_COMM_WORLD, &status);
            integral_sum += integral;
        }
        integral_sum *= 4;
    }
    if (rank == 0){
        printf("number of trapz:%d\nfrom a:%d to b:%d\nintegral: %f\n", n, a, b, integral_sum);
    }
    
    MPI_Finalize();
    return 0;
}
double Trap(int start, int finish, double h){
    // Trapezoidal rule
    double integral;
    int i;
    integral = ((1.0/(1 + pow(start*h, 2))) + (1.0/(1 + pow(finish*h, 2)))) * h / 2.0;
    for (i = start + 1; i < finish; ++i) {
            integral += h * (1.0/(1 + pow(i*h, 2)));
    }
    return integral;    
}

