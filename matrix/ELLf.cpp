#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <matrix/matrix.hpp>
#include "mpi.h"

using namespace std;


void ELLf::initialize(string path){
    vector<double> values;
    vector<long int> row;
    vector<long int> col;

    read_mtx_file(path, values, row, col);
    
    

    int maxInd = 0;
    int nbValuesInLine;
    int maxIndTab[MNL[0]];
    for(int j = 0; j <MNL[0]; j++){
        maxIndTab[j] = 0;
    }
    for(int j = 0; j <MNL[2]; j++){
        maxIndTab[row[j] - 1]++;
    }
    for(int j = 0; j <MNL[0]; j++){
        if(maxIndTab[j] > maxInd){
            maxInd = maxIndTab[j];
        }
    }
    maxColInd = maxInd;

    //Initialization of colInd 
    vector<long int> initialization;
    if(colInd != NULL){
        free(colInd);
    }
    if(val != NULL){
        free(val);
    }
    colInd = (long int *)malloc(sizeof(long int) * maxInd * MNL[0]);
    val = (float *)malloc(sizeof(float) * maxInd * MNL[0]);
    for(int i = 0; i < MNL[0]; i++){
        for(int j = 0; j < maxInd; j++){
            colInd[i*maxInd+j] = -1;
            val[i*maxInd+j] = nan("0");
        }
        initialization.push_back(0);
    }

    for(int i = 0; i< MNL[2]; i ++){
        colInd[row[i]*maxInd+initialization[row[i]]] = col[i];
        val[row[i]*maxInd+initialization[row[i]]] = ((float)(values[i]));
        initialization[row[i]] = initialization[row[i]] + 1; 
    }
    /* cout << "ColInd : " << endl;
    for(int i = 0; i < MNL[0]; i++){
        for(int j = 0; j < maxInd; j++){
            cout << colInd[i*maxInd + j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
    cout << "Val : " << endl;
    for(int i = 0; i < MNL[0]; i++){
        for(int j = 0; j < maxInd; j++){
            cout << val[i*maxInd+j] << ' ';
        }
        cout << endl;
    } */
};

void ELLf::print(){
    typedef vector<float> RowDouble;
    RowDouble row(MNL[1]);
    for(int i = 0; i < MNL[0]; i ++){
        for(int j = 0; j < MNL[1]; j++){
            row[j] = 0;
        }
        for(int j = 0; j < maxColInd; j++){
            if(colInd[i*maxColInd + j] != -1){
                row[colInd[i*maxColInd + j]] = val[i*maxColInd + j];
            }
            
        }
        cout << endl;
        for(int j = 0; j < row.size(); j++){
            cout << row[j] << ' ';
        }
    }
    cout << endl;
};

void ELLf::spmv(double * denseVector, int sizeDenseVector, double ** result){
    if(*result != NULL){
        *result = (double *)realloc(*result, sizeDenseVector * sizeof(double));
    }else{
        *result = (double *)malloc(sizeof(double) * sizeDenseVector);
    }
    if(MNL[1] != sizeDenseVector){
        for(int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
    }else{
        for(long int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = 0;
        }
        for(long int i = 0; i < MNL[0]; i++){
            for(long int j = 0; j < maxColInd; j++){
                if(colInd[i*maxColInd + j] != -1){
                    (*result)[i] +=  ((double)(val[i*maxColInd + j])) * denseVector[(int) (colInd[i* maxColInd + j])];
                    //result[i] +=  val[i][j] * denseVector[colInd[i][j]];
                }
            }
        }
    }
}

void ELLf::spmv(float * denseVector, int sizeDenseVector, double ** result){
    if(*result != NULL){
        *result = (double *)realloc(*result, sizeDenseVector * sizeof(double));
    }else{
        *result = (double *)malloc(sizeof(double) * sizeDenseVector);
    }
    if(MNL[1] != sizeDenseVector){
        for(int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
    }else{
        for(long int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = 0;
        }
        for(long int i = 0; i < MNL[0]; i++){
            for(long int j = 0; j < maxColInd; j++){
                if(colInd[i*maxColInd + j] != -1){
                    (*result)[i] +=  ((double)(val[i*maxColInd + j] * (denseVector[(int) (colInd[i* maxColInd + j])])));
                    //result[i] +=  val[i][j] * denseVector[colInd[i][j]];
                }
            }
        }
    }
}

void ELLf::spmv(float * denseVector, int sizeDenseVector, float ** result){
    if(*result != NULL){
        *result = (float *)realloc(*result, sizeDenseVector * sizeof(float));
    }else{
        *result = (float *)malloc(sizeof(float) * sizeDenseVector);
    }
    if(MNL[1] != sizeDenseVector){
        for(int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
    }else{
        for(long int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = 0;
        }
        for(long int i = 0; i < MNL[0]; i++){
            for(long int j = 0; j < maxColInd; j++){
                if(colInd[i*maxColInd + j] != -1){
                    (*result)[i] +=  val[i*maxColInd + j] * (denseVector[(int) (colInd[i* maxColInd + j])]);
                    //result[i] +=  val[i][j] * denseVector[colInd[i][j]];
                }
            }
        }
    }
}

long int ELLf::sizeColInd(int num){
    if(num == 0){
        return MNL[0];
    }else {
        if(num == 1){
            return maxColInd;
        }else{
            return -1;
        }
    }
}

float * ELLf::getVal(){
    return val;
}

long int * ELLf::getColInd(){
    return colInd;
}

float ELLf::norm(){
    float norm = 0;

    for(long int i = 0; i < MNL[0]; i++){
        for(long int j = 0; j < maxColInd; j++){
            if(colInd[i*maxColInd + j] != -1){
                norm += val[i*maxColInd + j];
            }

        }
    }
    return norm;
}


void ELLf::freeData(){
    if(colInd != NULL){
        free(colInd);
        colInd = NULL;
    }
    if(colInd != NULL){
        free(val);
        val = NULL;
    }
}
/* vector<double> ELL::spmv_mpi(vector<double> denseVector){
    int my_rank;
    int p;
    int tag=0;

    MPI_Status status;
    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    typedef vector<double> RowDouble;
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