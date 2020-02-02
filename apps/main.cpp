#include <matrix/matrix.hpp>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <unistd.h>

using namespace std;
    
int main(int argc, char** argv){
    printf("COUCOU\n");
    int size = 5;
    double * result = (double*)malloc(sizeof(double) * size);
    double * result1 = NULL;
    result[0] = 1;
    result[1] = 2;
    result[2] = 3;
    result[3] = 4;
    result[4] = 5;

    cout << "COO" << endl; 
    COO myMatrix3;
    myMatrix3.initialize("data/test2.mtx");
    myMatrix3.print();
    myMatrix3.spmv(result, size, &result1);
    cout << "result :" << endl;
    for(int i = 0; i < size; i++){
        cout << result1[i] << ' ';
    }
    cout << endl;

    cout << "CSR" << endl; 
    CSR myMatrix;
    myMatrix.initialize("data/test2.mtx");
    myMatrix.print();myMatrix.spmv(result, size, &result1);
    cout << "result :" << endl;
    for(int i = 0; i < size; i++){
        cout << result1[i] << ' ';
    }
    cout << endl; 
    
    cout << "ELLPACK" << endl; 
    ELL myMatrix2;
    myMatrix2.initialize("data/test2.mtx");
    myMatrix2.print();
     myMatrix2.spmv(result, size, &result1);
    cout << "result :" << endl;
    for(int i = 0; i < size; i++){
        cout << result1[i] << ' ';
    }
    cout << endl;

    

    cout << endl << "SGP" << endl;
    SGP myMatrix4;
    myMatrix4.initialize("data/test2.mtx");
    myMatrix4.print();
    myMatrix4.data();
    myMatrix3.spmv(result, size, &result);
    cout << "result :" << endl;
    for(int i = 0; i < size; i++){
        cout << result1[i] << ' ';
    }
    cout << endl;

    return 0;
}
