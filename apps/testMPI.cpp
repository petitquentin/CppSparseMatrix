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
    vector<double> result(5);
    vector<double> result1;
    vector<double> result2;
    int my_rank, p;
    cout << sizeof(p) << endl;
    result[0] = 1;
    result[1] = 2;
    result[2] = 3;
    result[3] = 4;
    result[4] = 5;

    cout << "COO" << endl; 
    COO myMatrix3;
    myMatrix3.initialize("data/test2.mtx");
    //myMatrix3.print();
    result1 = myMatrix3.spmv(result);
    cout << "result :" << endl;
    for(int i = 0; i < result1.size(); i++){
        cout << result1[i] << ' ';
    }
    cout << endl;

    MPI_Init(&argc, &argv);

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    result2 = spmv_mpi(&myMatrix3, result);

    cout << "result_MPI :" << endl;
    for(int i = 0; i < result2.size(); i++){
        cout << result2[i] << ' ';
    }
    cout << endl;

    MPI_Barrier(MPI_COMM_WORLD);
    CSR myMatrix;
    myMatrix.initialize("data/test2.mtx");

    result = spmv_mpi(&myMatrix, result);

    MPI_Finalize();
    cout << "result_MPI CSR " << my_rank << " :" << endl;
    for(int i = 0; i < result.size(); i++){
        cout << result[i] << ' ';
    }
    cout << endl;
    return 0;
}
