#include <matrix/matrix.hpp>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <unistd.h>
#include <mpi.h>

using namespace std;
    
int main(int argc, char** argv){
    int size = 10000;
    double * result = (double *)malloc(sizeof(double) * size);
    double * result1 = NULL;
    double * result2 = NULL;
    int my_rank, p;
    for(int i = 0; i < size; i++){
        result[i] = i/100.0;
    }

    cout << "COO" << endl; 
    COO myMatrix3;
    if(my_rank == 0){
        myMatrix3.initialize("data/smg2s.mtx");
    }
    

    MPI_Init(&argc, &argv);

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    
    /* spmv_mpi(&myMatrix3, result, size, &result2);

    cout << "result_MPI :" << endl;
    for(int i = 0; i < size; i++){
        cout << result2[i] << ' ';
    }
    cout << endl;

    free(result2);
    result2 = NULL;

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "____CSR____" <<endl;
    CSR myMatrix;
    myMatrix.initialize("data/smg2s.mtx");

    spmv_mpi(&myMatrix, result, size, &result2);
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "result_MPI CSR " << my_rank << " :" << endl;
    for(int i = 0; i < size; i++){
        cout << result2[i] << ' ';
    }
    cout << endl;

    free(result2);
    result2 = NULL; */

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "____ELL____" <<endl;
    ELL myMatrix2;
    myMatrix2.initialize("data/smg2s10000.mtx");

    free(result2);
    result2 = NULL;

    spmv_mpi(&myMatrix2, result, size, &result2);

    cout << "result_MPI ELL " << my_rank << " :" << endl;
    for(int i = 0; i < 8; i++){
        cout << result2[i] << ' ';
    }
    cout << endl;

    
    free(result2);
    result2 = NULL;
    cout << "____SGP____" <<endl;
    SGP myMatrix1;
    if(my_rank == 0){
        myMatrix1.initialize("data/smg2s10000.mtx");
        //myMatrix1.data();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "CC" <<endl;

    spmv_mpi(&myMatrix1, result, size, &result2);

    cout << "result_MPI SGP " << my_rank << " :" << endl;
    for(int i = 0; i < 8; i++){
        cout << result2[i] << ' ';
    }
    cout << endl; 


    MPI_Finalize();
    return 0;
}
