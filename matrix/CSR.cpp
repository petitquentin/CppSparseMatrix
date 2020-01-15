#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <matrix/matrix.hpp>
include "mpi.h"

using namespace std;


void CSR::initialize(string path){
    vector<double> values;
    vector<long int> row;
    vector<long int> col;

    read_mtx_file(path, values, row, col);

    //Initialize ptr
    ptr.clear();
    val.clear();
    ind.clear();

    ptr.push_back(0);
    for(long int i = 0; i < MNL[0]; i++){
        ptr.push_back(0);
    }
    for(long int n = 0; n < MNL[2]; n++){
        ptr[row[n]]++;
        val.push_back(0);
        ind.push_back(0);
    }
    for(long int i = 0, cumsum = 0; i< MNL[0]; i++){
        long int temp = ptr[i];
        ptr[i] = cumsum;
        cumsum += temp;
    }
    ptr[ptr.size()-1] = MNL[2];

    for(long int n = 0; n < MNL[2]; n++){
        long int r = row[n];
        long int dest = ptr[r];
        //cout << dest <<endl;
        ind[dest] = col[n];
        val[dest] = values[n];

        ptr[r]++;
    }

    for(long int i = 0, last = 0; i <= MNL[0]; i++){
        long int temp = ptr[i];
        ptr[i] = last;
        last = temp;
    }

};

void CSR::print(){
    long int c = 0;
    for(long int i = 1; i < ptr.size(); i++){
        long int start = ptr[i-1];
        long int end = ptr[i];
        vector<long int>::const_iterator first = ind.begin() + start;
        vector<long int>::const_iterator last = ind.begin() + end;
        vector<long int> row(first,last);
        for(int j = 0; j < ptr.size(); j++){
            if(count(row.begin(), row.end(), j) == 0)
                cout << '0' << ' ';
            else{
                cout << val[c] << '(' << c << ')' << ' ';
                c++;
            }
        }
        std::cout << std::endl;
        

    }
    //print ind
    cout << "print ind" << endl;
    for(int i = 0; i < ind.size(); i++){
        cout << ind[i] << ' ';
    }
    cout << endl;
    //print val
    cout << "print val" << endl;
    for(int i = 0; i < val.size(); i++){
        cout << val[i] << ' ';
    }
    cout << endl;
    //print ptr
    cout << "print ptr" << endl;
    for(int i = 0; i < ptr.size(); i++){
        cout << ptr[i] << ' ';
    }
    cout << endl;
    cout << endl;
    cout << val.size() << ' ' << ind.size() << ' ' << ptr.size() << ' ' << ptr[ptr.size()-1] << ' ' << endl;
};


vector<double> CSR::spmv(vector<double> denseVector){
    typedef vector<double> RowDouble;
    RowDouble result(MNL[1]);
    if(MNL[1] != denseVector.size()){
        for(int i = 0; i < result.size(); i++){
            result[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
        return result;
    }
    for(int i = 0; i < result.size(); i++){
        result[i] = 0;
    }
    for(long int i = 0; i < ptr.size()-1; i++){
        long int start = ptr[i];
        long int end = ptr[i+1];
        vector<long int>::const_iterator firstInd= ind.begin() + start;
        vector<long int>::const_iterator lastInd = ind.begin() + end;
        vector<double>::const_iterator firstVal = val.begin() + start;
        vector<double>::const_iterator lastVal = val.begin() + end;
        vector<long int> rowInd(firstInd,lastInd);
        vector<long int> rowVal(firstVal, lastVal);
        for(int j = 0; j < rowInd.size(); j++){
            result[i] += denseVector[rowInd[j]] * rowVal[j]; 
        }
    }
    
    /* for(long int i = 0; i < ptr.size(); i++){
        result[i] += denseVector[ind[i]] * val[i];
    } */
    return result;
}

vector<double> CSR::spmv_mpi(vector<double> denseVector){
    int my_rank;
    int p;
    int tag=0;

    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    typedef vector<double> RowDouble;
    RowDouble result(MNL[1]);
    if(MNL[1] != denseVector.size()){
        for(int i = 0; i < result.size(); i++){
            result[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
        return result;
    }
    for(int i = 0; i < result.size(); i++){
        result[i] = 0;
    }
    for(long int i = my_rank * ceil((ptr.size() - 1)/p); i < min((ptr.size() - 1), (my_rank + 1) * ceil((ptr.size() - 1)/p); i++){
        long int start = ptr[i];
        long int end = ptr[i+1];
        vector<long int>::const_iterator firstInd= ind.begin() + start;
        vector<long int>::const_iterator lastInd = ind.begin() + end;
        vector<double>::const_iterator firstVal = val.begin() + start;
        vector<double>::const_iterator lastVal = val.begin() + end;
        vector<long int> rowInd(firstInd,lastInd);
        vector<long int> rowVal(firstVal, lastVal);
        for(int j = 0; j < rowInd.size(); j++){
            result[i] += denseVector[rowInd[j]] * rowVal[j]; 
        }
    }
    MPI_Reduce(&result, &result, sizeof(result), MPI_DOUBLE, MPI_SUM, tag, MPI_COMM_WORLD);

    MPI_Finalize();
    if(my_rank ==0){
        return result;
    }
    
}