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
    
    int maxInd = 0;
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
    acc.clear();
    acr.clear();
    ai.clear();
    aj.clear();
    ic.clear();
    jc.clear();

    for(int i = 0; i < MNL[0]; i++){
        RowLongInt rowLongInt(maxInd);
        RowDouble rowDouble(maxInd);
        for(int j = 0; j < maxInd; j++){
            rowLongInt[j] = -1;
            rowDouble[j] = nan("0");
        }
        aj.push_back(rowLongInt);
        ic.push_back(rowLongInt);
        acc.push_back(rowDouble);
        initializationAcc.push_back(0);
        initializationAcr.push_back(0);
    }
    for(int i = 0; i < maxInd; i++){
        RowLongInt rowLongInt(MNL[0]);
        RowDouble rowDouble(MNL[0]);
        for(int j = 0; j < MNL[0]; j++){
            rowLongInt[j] = -1;
            rowDouble[j] = nan("0");
        }
        ai.push_back(rowLongInt);
        jc.push_back(rowLongInt);
        acr.push_back(rowDouble);
    }

    for(int i = 0; i< MNL[2]; i ++){
        aj[row[i]][initializationAcc[row[i]]] = col[i];
        acc[row[i]][initializationAcc[row[i]]] = values[i];
        initializationAcc[row[i]] = initializationAcc[row[i]] + 1; 
    }
    for(int i = 0; i< MNL[2]; i ++){
        ai[initializationAcr[col[i]]][col[i]] = row[i];
        acr[initializationAcr[col[i]]][col[i]] = values[i];
        initializationAcr[col[i]] = initializationAcr[col[i]] + 1;
    }

    //Initialization of Ic and Jc
    for(int i = 0; i < ai.size(); i++){
        for(int j = 0; j < aj.size(); j++){
            if(ai[i][j] != -1){
                auto result = find(aj[ai[i][j]].begin(), aj[ai[i][j]].end(), j);
                if(result != aj[ai[i][j]].end()){
                    jc[i][j] = distance(aj[ai[i][j]].begin(), result);
                }else{
                    cout << "We'd a problem" << endl;
                }
            }
            if(aj[j][i] != -1){
                long int k = 0;
                while(ai[k][aj[j][i]] != j){
                    k++;
                }
                ic[j][i] = k;
            }
            
        }
    }
};

void SGP::data(){
    cout << "---" << endl;
    cout << "Acc" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < acc.size(); i++){
        for(long int j = 0; j < acc[i].size(); j++){
            cout << acc[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;


    cout << "---" << endl;
    cout << "Acr" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < acr.size(); i++){
        for(long int j = 0; j < acr[i].size(); j++){
            cout << acr[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;

    cout << "---" << endl;
    cout << " Ai" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < ai.size(); i++){
        for(long int j = 0; j < ai[i].size(); j++){
            cout << ai[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;

    cout << "---" << endl;
    cout << " Aj" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < aj.size(); i++){
        for(long int j = 0; j < aj[i].size(); j++){
            cout << aj[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;

    cout << "---" << endl;
    cout << " Jc" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < jc.size(); i++){
        for(long int j = 0; j < jc[i].size(); j++){
            cout << jc[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;

    cout << "---" << endl;
    cout << " Ic" << endl;
    cout << "---" << endl;
    for(long int i = 0; i < ic.size(); i++){
        for(long int j = 0; j < ic[i].size(); j++){
            cout << ic[i][j] << ' ';
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
        for(int j = 0; j < aj[i].size(); j++){
            if(aj[i][j] != -1){
                row[aj[i][j]] = acc[i][j];
            }
            
        }
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

        for(long int i = 0; i < acc.size(); i++){
            for(long int j = 0; j < acc[0].size(); j++){
                if(aj[i][j] != -1){
                    (*result)[i] += acc[i][j] * denseVector[aj[i][j]];
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