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
    result[0] = 1;
    result[1] = 2;
    result[2] = 3;
    result[3] = 4;
    result[4] = 5;

    cout << "COO" << endl; 
    COO myMatrix3;
    myMatrix3.initialize("data/test2.mtx");
    myMatrix3.print();
    result1 = myMatrix3.spmv(result);
    cout << "result :" << endl;
    for(int i = 0; i < result1.size(); i++){
        cout << result1[i] << ' ';
    }
    cout << endl;

    MPI_Init(&argc, &argv);

    result1 = spmv_mpi(myMatrix3, result);

    MPI_Finalize();
    cout << "result_MPI :" << endl;
    for(int i = 0; i < result1.size(); i++){
        cout << result1[i] << ' ';
    }
    cout << endl;

    return 0;
}
