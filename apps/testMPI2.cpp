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
#include <ctime>
#include <sys/time.h>

#include <chrono>

using namespace std;
    
double my_gettimeofday(){
    struct timeval tmp_time;
    gettimeofday(&tmp_time, NULL);
    return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

int main(int argc, char** argv){
    int size = 217918;
    string file = "data/pwtk.mtx";
    double start;
    int elapsed_seconds;
    double * result = NULL;
    double * result1 = NULL;
    double * result2 = NULL;
    int my_rank, p;
    MPI_Init(&argc, &argv);

    cout << file << endl;

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if(my_rank == 0){
        result = (double *)malloc(sizeof(double) * size);
        for(int i = 0; i < size; i++){
            result[i] = i/(1000/3.0);
        }
    }
    cout << "COO" << endl; 
    COO myMatrix3;
    myMatrix3.initialize(file);
    

    
    
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    spmv_mpi(&myMatrix3, result, size, &result2);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS COO : " <<  my_gettimeofday() - start << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    cout << "result_MPI :" << endl;
    for(int i = 0; i < 10; i++){
        cout << result2[i] << ' ';
    }
    cout << endl;

    free(result2);
    result2 = NULL;

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "____CSR____" <<endl;
    CSR myMatrix;
    myMatrix.initialize(file);

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    spmv_mpi(&myMatrix, result, size, &result2);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS CSR : " << my_gettimeofday() - start << endl;
    }

    cout << "result_MPI CSR " << my_rank << " :" << endl;
    for(int i = 0; i < 10; i++){
        cout << result2[i] << ' ';
    }
    cout << endl;

    free(result2);
    result2 = NULL;

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "____ELL____" <<endl;
    ELL myMatrix2;
    myMatrix2.initialize(file);

    free(result2);
    result2 = NULL;

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    spmv_mpi(&myMatrix2, result, size, &result2);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS ELL : " << my_gettimeofday() - start << endl;
    }

    cout << "result_MPI ELL " << my_rank << " :" << endl;
    for(int i = 0; i < 10; i++){
        cout << result2[i] << ' ';
    }
    cout << endl;

    
    free(result2);
    result2 = NULL;
    cout << "____SGP____" <<endl;
    SGP myMatrix1;
    if(my_rank == 0){
        myMatrix1.initialize(file);
        //myMatrix1.data();
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "CC" <<endl;

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    spmv_mpi(&myMatrix1, result, size, &result2);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS SGP : " << my_gettimeofday() - start << endl;
    }
    cout << "result_MPI SGP " << my_rank << " :" << endl;
    for(int i = 0; i < 10; i++){
        cout << result2[i] << ' ';
    }
    cout << endl; 

    MPI_Finalize();
    return 0;
}
