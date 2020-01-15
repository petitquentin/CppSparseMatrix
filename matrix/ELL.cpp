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


void ELL::initialize(string path){
    vector<double> values;
    vector<long int> row;
    vector<long int> col;

    read_mtx_file(path, values, row, col);
    
    int maxInd = 0;
    int nbValuesInLine;
    for (int i = 1; i <= MNL[0]; i++){
        nbValuesInLine = count(row.begin(), row.end(), i);
        if(nbValuesInLine > maxInd){
            maxInd = nbValuesInLine;
        }
    }

    //Initialization of colInd 
    typedef vector<long int> RowLongInt;
    typedef vector<double> RowDouble;
    vector<long int> initialization;
    colInd.clear();
    val.clear();
    for(int i = 0; i < MNL[0]; i++){
        RowLongInt rowLongInt(maxInd);
        RowDouble rowDouble(maxInd);
        for(int j = 0; j < maxInd; j++){
            rowLongInt[j] = -1;
            rowDouble[j] = nan("0");
        }
        colInd.push_back(rowLongInt);
        val.push_back(rowDouble);
        initialization.push_back(0);
    }

    for(int i = 0; i< row.size(); i ++){
        colInd[row[i]][initialization[row[i]]] = col[i];
        val[row[i]][initialization[row[i]]] = values[i];
        initialization[row[i]] = initialization[row[i]] + 1; 
    }
    cout << "ColInd : " << endl;
    for(int i = 0; i < colInd.size(); i++){
        for(int j = 0; j < colInd[i].size(); j++){
            cout << colInd[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
    cout << "Val : " << endl;
    for(int i = 0; i < val.size(); i++){
        for(int j = 0; j < val[i].size(); j++){
            cout << val[i][j] << ' ';
        }
        cout << endl;
    }
};

void ELL::print(){
    typedef vector<double> RowDouble;
    RowDouble row(MNL[1]);
    for(int i = 0; i < MNL[0]; i ++){
        for(int j = 0; j < MNL[1]; j++){
            row[j] = 0;
        }
        for(int j = 0; j < colInd[i].size(); j++){
            if(colInd[i][j] != -1){
                row[colInd[i][j]] = val[i][j];
            }
            
        }
        cout << endl;
        for(int j = 0; j < row.size(); j++){
            cout << row[j] << ' ';
        }
    }
    cout << endl;
};

vector<double> ELL::spmv(vector<double> denseVector){
    typedef vector<double> RowDouble;
    RowDouble result(MNL[1]);
    if(MNL[1] != denseVector.size()){
        for(int i = 0; i < result.size(); i++){
            result[i] = NULL;
        }
        cout << "the vector is not of the right size" << endl;
        return result;
    }
    for(long int i = 0; i < result.size(); i++){
        result[i] = 0;
    }
    for(long int i = 0; i < colInd.size(); i++){
        for(long int j = 0; j < colInd[i].size(); j++){
            if(colInd[i][j] != -1){
                result[i] +=  val[i][j] * denseVector[colInd[i][j]];
            }
        }
    }
    return result;
}

vector<double> ELL::spmv_mpi(vector<double> denseVector){
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
    
}