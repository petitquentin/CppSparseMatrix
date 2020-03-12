#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <matrix/matrix.hpp>
#include <mpi.h>

using namespace std;


void CSRf::initialize(string path){
    vector<double> values;
    vector<long int> row;
    vector<long int> col;

    read_mtx_file(path, values, row, col);

    //Initialize ptr
    if(ptr != NULL){
        free(ptr);
    }
    if(val != NULL){
        free(val);
    }
    if(ind != NULL){
        free(ind);
    }

    ptr = (long int *)malloc(sizeof(long int) * MNL[0]+1);
    ind = (long int *)malloc(sizeof(long int) * MNL[2]);
    val = (float *)malloc(sizeof(float) * MNL[2]);
    
    for(long int i = 0; i < MNL[0]+1; i++){
        ptr[i] = 0;
    }
    for(long int n = 0; n < MNL[2]; n++){
        ptr[row[n]]++;
        val[n] = 0;
        ind[n] = 0;
    }
    for(long int i = 0, cumsum = 0; i< MNL[0]; i++){
        long int temp = ptr[i];
        ptr[i] = cumsum;
        cumsum += temp;
    }
    ptr[MNL[0]] = MNL[2];

    for(long int n = 0; n < MNL[2]; n++){
        long int r = row[n];
        long int dest = ptr[r];
        //cout << dest <<endl;
        ind[dest] = col[n];
        val[dest] = ((float)(values[n]));

        ptr[r]++;
    }

    for(long int i = 0, last = 0; i <= MNL[0]; i++){
        long int temp = ptr[i];
        ptr[i] = last;
        last = temp;
    }

};

void CSRf::print(){
    long int c = 0;
    for(long int i = 1; i < MNL[0]+1; i++){
        long int start = ptr[i-1];
        long int end = ptr[i];
        for(int j = 0; j < MNL[1]; j++){
            int it = 0;
            for(int k = start; k < end; k++){
                if(ind[k] == j){
                    it++;
                }
            }
            if(it == 0){
                cout << '0' << ' ';
            }else{
                cout << val[c] << '(' << c << ')' << ' ';
                c++;
            }
        }
        std::cout << std::endl;
        

    }
    //print ind
    cout << "print ind" << endl;
    for(int i = 0; i < MNL[2]; i++){
        cout << ind[i] << ' ';
    }
    cout << endl;
    //print val
    cout << "print val" << endl;
    for(int i = 0; i < MNL[2]; i++){
        cout << val[i] << ' ';
    }
    cout << endl;
    //print ptr
    cout << "print ptr" << endl;
    for(int i = 0; i < MNL[0] + 1; i++){
        cout << ptr[i] << ' ';
    }
    cout << endl;
};


void CSRf::spmv(double * denseVector, int sizeDenseVector, double ** result){
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
        for(int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = 0;
        }
        for(long int i = 0; i < MNL[0]; i++){
            long int start = ptr[i];
            long int end = ptr[i+1];
            for(int j = start; j < end; j++){
                (*result)[i] += denseVector[ind[j]] * ((double)(val[j])); 
            }
        }
    }
    
    /* for(long int i = 0; i < ptr.size(); i++){
        result[i] += denseVector[ind[i]] * val[i];
    } */
}

void CSRf::spmv(float * denseVector, int sizeDenseVector, double ** result){
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
        for(int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = 0;
        }
        for(long int i = 0; i < MNL[0]; i++){
            long int start = ptr[i];
            long int end = ptr[i+1];
            for(int j = start; j < end; j++){
                (*result)[i] += ((double)(denseVector[ind[j]] * val[j])); 
            }
        }
    }
    
    /* for(long int i = 0; i < ptr.size(); i++){
        result[i] += denseVector[ind[i]] * val[i];
    } */
}

void CSRf::spmv(float * denseVector, int sizeDenseVector, float ** result){
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
        for(int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = 0;
        }
        for(long int i = 0; i < MNL[0]; i++){
            long int start = ptr[i];
            long int end = ptr[i+1];
            for(int j = start; j < end; j++){
                (*result)[i] += denseVector[ind[j]] * val[j]; 
            }
        }
    }
    
    /* for(long int i = 0; i < ptr.size(); i++){
        result[i] += denseVector[ind[i]] * val[i];
    } */
}

int CSRf::getValSize(){
    return MNL[2];
}
long int * CSRf::getPtr(){
    return ptr;
}
long int * CSRf::getInd(){
    return ind;
}

float * CSRf::getVal(){
    return val;
}

float CSRf::norm(){
    float norm = 0;

    for(long int i = 0; i < MNL[2]; i++){
        norm += val[i];    
    }
    return norm;
}

void CSRf::freeData(){
    if(val != NULL){
        free(val);
        val = NULL;
    }
    if(ind != NULL){
        free(ind);
        ind = NULL;
    }
    if(ptr != NULL){
        free(ptr);
        ptr = NULL;
    }
}

/* vector<double> CSR::spmv_mpi(vector<double> denseVector){
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
    
} */