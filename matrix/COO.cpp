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
    if(val != NULL){
        free(val);
        val = NULL;
    }
    if(row != NULL){
        free(row);
        row = NULL;
    }
    if(col != NULL){
        free(col);
        col = NULL;
    }


    vector<double> values;
    vector<long int> rowInit;
    vector<long int> colInit;
    read_mtx_file(path, values, rowInit, colInit);

    
    val = (double *)malloc(sizeof(double) * MNL[2]);
    col = (long int *)malloc(sizeof(long int) * MNL[2]);
    row = (long int *)malloc(sizeof(long int) * MNL[2]);

    //val = values.data();
    //col = colInit.data();
    //row = rowInit.data();
    copy(values.data(), values.data() + values.size(), val);
    copy(colInit.data(), colInit.data() + colInit.size(), col);
    copy(rowInit.data(), rowInit.data() + rowInit.size(), row);
    

    /* cout << "Col" << endl;
    for(int i = 0; i < MNL[2]; i++){
        cout << col[i] << ' ';
    }
    cout << endl;

    cout << "row" << endl;
    for(int i = 0; i < MNL[2]; i++){
        cout << row[i] << ' ';
    }
    cout << endl;

    cout << "Val" << endl;
    for(int i = 0; i < MNL[2]; i++){
        cout << val[i] << ' ';
    }
    cout << endl; */
};

void COO::print(){
    typedef vector<double> RowDouble;
    RowDouble actualRow(MNL[1]);
    for(int i = 0; i < MNL[0]; i ++){
        for(int j = 0; j < MNL[1]; j++){
            actualRow[j] = 0;
        }
        for(int j = 0; j < MNL[2]; j++){
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

void COO::spmv(double * denseVector, int sizeDenseVector, double ** result){
    if(*result != NULL){
        *result = (double *)realloc(*result, sizeDenseVector * sizeof(double));
    }else{
        *result = (double *)malloc(sizeof(double) * sizeDenseVector);
    }
    if(MNL[1] != sizeDenseVector){
        for(long int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
    }else{
        for(int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = 0;
        }
        for(int i = 0; i < MNL[2]; i++){
            (*result)[row[i]] += val[i]*denseVector[col[i]];
        }
    }
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
    return MNL[2];
}

long int * COO::getRow(){
    return row;
}
long int * COO::getCol(){
    return col;
}

double * COO::getVal(){
    return val;
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