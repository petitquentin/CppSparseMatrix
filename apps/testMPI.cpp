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

    int my_rank, p;
    int size = 5;
    MPI_Init(&argc, &argv);

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    double * result = NULL;
    if(my_rank == 0){
        cout << "JE PASSE LA " << my_rank << endl;
        result = (double *)malloc(sizeof(double) * size);
        result[0] = 1;
        result[1] = 2;
        result[2] = 3;
        result[3] = 4;
        result[4] = 5;
    }
    double * result1 = NULL;
    double * result2 = NULL;
    cout << sizeof(p) << endl;
    

    double * vec = (double *)malloc(sizeof(double) * 8);
    for(int i = 0; i < 8; i++){
        vec[i] = i + 1;
    }

    cout << "COO" << endl; 
    COO myMatrix3;
    myMatrix3.initialize("data/test2.mtx");
    //myMatrix3.print();
    //myMatrix3.spmv(result, size, &result1);
    MPI_Barrier(MPI_COMM_WORLD);
    /* cout << "result :" << endl;
    for(int i = 0; i < size; i++){
        cout << result1[i] << ' ';
    }
    cout << endl; */

    
    
    spmv_mpi(&myMatrix3, result, size, &result2);

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
    myMatrix.initialize("data/test2.mtx");

    spmv_mpi(&myMatrix, result, size, &result2);
    MPI_Barrier(MPI_COMM_WORLD);
    cout << "result_MPI CSR " << my_rank << " :" << endl;
    for(int i = 0; i < size; i++){
        cout << result2[i] << ' ';
    }
    cout << endl;

    free(result2);
    result2 = NULL;

    MPI_Barrier(MPI_COMM_WORLD);
    cout << "____ELL____" <<endl;
    ELL myMatrix2;
    myMatrix2.initialize("data/test3.mtx");
    //myMatrix2.spmv(result, size, &result2);
    /* cout << "result :" << endl;
    for(int i = 0; i < size; i++){
        cout << result2[i] << ' ';
    }
    cout << endl; */

    free(result2);
    result2 = NULL;

    spmv_mpi(&myMatrix2, vec, 8, &result2);

    cout << "result_MPI ELL " << my_rank << " :" << endl;
    for(int i = 0; i < 8; i++){
        cout << result2[i] << ' ';
    }
    cout << endl;

    
    free(result2);
    result2 = NULL;
    cout << "____SGP____" <<endl;
    SGP myMatrix1;
    myMatrix1.initialize("data/test3.mtx");
    MPI_Barrier(MPI_COMM_WORLD);

    spmv_mpi(&myMatrix1, vec, 8, &result2);

    cout << "result_MPI SGP " << my_rank << " :" << endl;
    for(int i = 0; i < 8; i++){
        cout << result2[i] << ' ';
    }
    cout << endl;


    MPI_Finalize();
    return 0;
}
