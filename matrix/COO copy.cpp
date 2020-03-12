#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <mpi.h>
#include <matrix/matrix.hpp>

using namespace std;

void COO::initialize(string path)
{
    val.clear();
    row.clear();
    col.clear();
    read_mtx_file(path, val, row, col);
};

void COO::print(){
    typedef vector<double> RowDouble;
    RowDouble actualRow(MNL[1]);
    for(int i = 0; i < MNL[0]; i ++){
        for(int j = 0; j < MNL[1]; j++){
            actualRow[j] = 0;
        }
        for(int j = 0; j < row.size(); j++){
            if(row[j] == i){
                actualRow[col[j]] = val[j];
            }
            
        }
        cout << endl;
        for(int j = 0; j < actualRow.size(); j++){
            cout << actualRow[j] << ' ';
        }
    }
    cout << endl;
};

vector<double> COO::spmv(vector<double> denseVector){
    
    typedef vector<double> RowDouble;
    RowDouble result(MNL[1]);
    if(MNL[1] != denseVector.size()){
        for(long int i = 0; i < result.size(); i++){
            result[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
        return result;
    }
    for(long int i = 0; i < result.size(); i++){
        result[i] = 0;
    }
    for(long int i = 0; i < row.size(); i++){
        result[row[i]] += val[i]*denseVector[col[i]];
    }
    return result;
}

/* vector<double> COO::spmv(vector<double> denseVector){
    
    
    RowDouble result(MNL[1]);
    if(MNL[1] != denseVector.size()){
        for(long int i = 0; i < result.size(); i++){
            result[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
        return result;
    }
    for(long int i = 0; i < result.size(); i++){
        result[i] = 0;
    }
    for(long int i = 0; i < row.size(); i++){
        result[row[i]] += val[i]*denseVector[col[i]];
    }
    return result;
} */

int COO::getRowSize(){
    return row.size();
}

long int * COO::getRow(){
    return row.data();
}
long int * COO::getCol(){
    return col.data();
}

double * COO::getVal(){
    return val.data();
}



/* vector<double> COO::spmv_mpi(vector<double> denseVector){
    int my_rank;
    int p;
    int tag=0;

    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    typedef vector<double> RowDouble;result.size()
    RowDouble result[MNL[1]];
    if(MNL[1] != denseVector.size()){
        for(long int i = 0; i < result.size(); i++){
            result[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
        return result;
    }
    for(long int i = 0; i < result.size(); i++){
        result[i] = 0;
    }
    for(long int i = my_rank * ceil(row.size()/p); i < min(row.size(), (my_rank + 1) * ceil(row.size()/p); i++){
        result[row[i]] += val[i]*denseVector[col[i]];
    }
    MPI_Reduce(&result, &result, sizeof(result), MPI_DOUBLE, MPI_SUM, tag, MPI_COMM_WORLD);

    MPI_Finalize();
    if(my_rank ==0){
        return result;
    }
    
} */