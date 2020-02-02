#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <matrix/matrix.hpp>

using namespace std;


void SGP::initialize(string path){
    vector<double> values;
    vector<long int> row;
    vector<long int> col;

    read_mtx_file(path, values, row, col);
    
    int nbValuesRow;
    int nbValuesCol;
    for (int i = 1; i <= MNL[0]; i++){
        nbValuesRow = 0;
        nbValuesCol = 0;
        for(int j = 0; j <= MNL[2]; j++){
            if(row[j] == i){
                nbValuesRow++;
            }
            if(col[j] == i){
                nbValuesCol++;
            }
        }
        if(nbValuesRow > maxInd){
            maxInd = nbValuesRow;
        }
        if(nbValuesCol > maxInd){
            maxInd = nbValuesCol;
        }
    }
    

    typedef vector<long int> RowLongInt;
    typedef vector<double> RowDouble;
    vector<long int> initializationAcc;
    vector<long int> initializationAcr;
    if(acc != NULL){
        free(acc);
        acc = NULL;
    }
    if(acr != NULL){
        free(acr);
        acr = NULL;
    }
    if(ai != NULL){
        free(ai);
        ai = NULL;
    }
    if(aj != NULL){
        free(aj);
        aj = NULL;
    }
    if(ic != NULL){
        free(ic);
        ic = NULL;
    }
    if(jc != NULL){
        free(jc);
        jc = NULL;
    }
    acr = (double *)malloc(sizeof(double) * MNL[0] * maxInd);
    acc = (double *)malloc(sizeof(double) * MNL[0] * maxInd);
    ic = (long int*)malloc(sizeof(long int) * MNL[0] * maxInd);
    jc = (long int*)malloc(sizeof(long int) * MNL[0] * maxInd);
    ai = (long int*)malloc(sizeof(long int) * MNL[0] * maxInd);
    aj = (long int*)malloc(sizeof(long int) * MNL[0] * maxInd);

    for(int i = 0; i < MNL[0]; i++){
        initializationAcc.push_back(0);
        initializationAcr.push_back(0);
    }
    for(int i = 0; i < maxInd * MNL[0]; i++){
        ai[i] = -1;
        aj[i] = -1;
        ic[i] = -1;
        jc[i] = -1;
        acr[i] = nan("0");
        acc[i] = nan("0");
    }

    for(int i = 0; i< MNL[2]; i ++){
        aj[row[i] * maxInd + initializationAcc[row[i]]] = col[i];
        acc[row[i] * maxInd + initializationAcc[row[i]]] = values[i];
        initializationAcc[row[i]] = initializationAcc[row[i]] + 1; 
    }
    for(int i = 0; i< MNL[2]; i ++){
        ai[initializationAcr[col[i]] * MNL[0] + col[i]] = row[i];
        acr[initializationAcr[col[i]] * MNL[0] + col[i]] = values[i];
        initializationAcr[col[i]] = initializationAcr[col[i]] + 1;
    }

    //Initialization of Ic and Jc
    for(int i = 0; i < maxInd; i++){
        for(int j = 0; j < MNL[0]; j++){
            if(ai[i * MNL[0] + j] != -1){
                long int result = -1;
                for(long int k  = 0; k < maxInd; k++){
                    if(aj[ai[i*MNL[0] + j] + k] == j){
                        result = k;
                    }
                    if(k != -1){
                        jc[i*MNL[0] + j] = result;
                    }else{
                        cout << "We'd a problem" << endl;
                    }
                }
                //auto result = find(aj[ai[i][j]].begin(), aj[ai[i][j]].end(), j);
                //if(result != aj[ai[i* MNL[0] + j]].end()){
            }
            if(aj[j * maxInd + i] != -1){
                long int k = 0;
                while(ai[k * MNL[0] + (aj[j* maxInd + i])] != j){
                    k++;
                }
                ic[j*maxInd + i] = k;
            }
            
        }
    }
};

void SGP::data(){
    cout << "---" << endl;
    cout << "Acc" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < MNL[0]; i++){
        for(long int j = 0; j < maxInd; j++){
            cout << acc[i * maxInd + j] << ' ';
        }
        cout << endl;
    }
    cout << endl;


    cout << "---" << endl;
    cout << "Acr" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < maxInd; i++){
        for(long int j = 0; j < MNL[0]; j++){
            cout << acr[i * MNL[0] + j] << ' ';
        }
        cout << endl;
    }
    cout << endl;

    cout << "---" << endl;
    cout << " Ai" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < maxInd; i++){
        for(long int j = 0; j < MNL[0]; j++){
            cout << ai[i * MNL[0] + j] << ' ';
        }
        cout << endl;
    }
    cout << endl;

    cout << "---" << endl;
    cout << " Aj" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < MNL[0]; i++){
        for(long int j = 0; j < maxInd; j++){
            cout << aj[i * maxInd + j] << ' ';
        }
        cout << endl;
    }
    cout << endl;

    cout << "---" << endl;
    cout << " Jc" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < maxInd; i++){
        for(long int j = 0; j < MNL[0]; j++){
            cout << jc[i * MNL[0] + j] << ' ';
        }
        cout << endl;
    }
    cout << endl;

    cout << "---" << endl;
    cout << " Ic" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < MNL[0]; i++){
        for(long int j = 0; j < maxInd; j++){
            cout << ic[i * maxInd + j] << ' ';
        }
        cout << endl;
    }
    cout << endl;

}

void SGP::print(){
    typedef vector<double> RowDouble;
    RowDouble row(MNL[1]);
    for(int i = 0; i < MNL[0]; i ++){
        for(long int j = 0; j < MNL[1]; j++){
            row[j] = 0;
        }
        for(int j = 0; j < maxInd; j ++){
            if(aj[i * maxInd + j] != -1){
                row[aj[i * maxInd + j]] = acc[i * maxInd + j];
            }
        }
        /* for(int j = 0; j < MNL[0]; j++){
            if(aj[i * maxInd + j] != -1){
                row[aj[i * maxInd + j]] = acc[i * maxInd + j];
            }
            
        } */
        cout << endl;
        for(int j = 0; j < row.size(); j++){
            cout << row[j] << ' ';
        }
    }
    cout << endl;
};

void SGP::spmv(double * denseVector, int sizeDenseVector, double ** result){
    if(*result != NULL){
        *result = (double *)realloc(*result, sizeDenseVector * sizeof(double));
    }else{
        *result = (double *)malloc(sizeof(double) * sizeDenseVector);
    }
    if(MNL[1] != sizeDenseVector){
        for(int i = 0; i < sizeDenseVector; i++){
            (*result)[i] = 0;
        }
        cout << "the vector is not of the right size" << endl;
    }else{
        for(long int i = 0; i < sizeDenseVector; i ++){
            (*result)[i] = 0;
        }

        for(long int i = 0; i < MNL[0]; i++){
            for(long int j = 0; j < maxInd; j++){
                if(aj[i * maxInd + j] != -1){
                    (*result)[i] += acc[i * maxInd + j] * denseVector[aj[i * maxInd + j]];
                }
            }
        }
    }
    /* for(long int i = 0; i < result.size(); i++){
        result[i] = 0;
    }
    for(long int i = 0; i < colInd.size(); i++){
        for(long int j = 0; j < colInd[i].size(); j++){
            if(colInd[i][j] != -1){
                result[i] +=  val[i][j] * denseVector[colInd[i][j]];
            }
        }
    } */
}