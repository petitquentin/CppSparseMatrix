#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <matrix/matrix.hpp>
#include <mpi.h>

#define TAGROW   1
#define TAGCOL   2
#define TAGIND   2
#define TAGPTR   1
#define TAGVAL   3
#define TAGSIZE  4
#define TAGINDEX 5
#define TAGAI  200
#define TAGACR 201
#define TAGJC  202
#define TAGAJ  204
#define TAGIC  205
#define TAGACC 203
#define TAGREDUCE   500
#define TAGSPREAD   501
#define TAGGET      502
#define TAGEXCHANGE 503
#define TAGVECTOR   808

using namespace std;

//Matrix values : float AND denseVector : double
void spmv_mpi(COOf *matrix, double * denseVector, int sizeDenseVector, double ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    long int MNL1;
    int rowSize;
    if(my_rank == 0){
        MNL1 = matrix->getMNL(1);
        rowSize = matrix->getRowSize();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&MNL1, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    

    typedef vector<double> RowDouble;
    //RowDouble result(MNL1);
    if(*result != NULL){
        free(*result);
        *result = NULL;
    }
    *result = (double *)malloc(sizeof(double)*MNL1);
    if(my_rank == 0){
        /* vector<long int> row = matrix.getRow();
        vector<long int> col = matrix.getCol();
        vector<double> val = matrix.getVal(); */
        long int * row;
        long int * col;
        float * val;
        row = (long int *)malloc(sizeof(long int)*rowSize);
        col = (long int *)malloc(sizeof(long int)*rowSize);
        val = (float *)malloc(sizeof(float)* rowSize);
        row = matrix->getRow();
        col = matrix->getCol();
        val = matrix->getVal();

        
        //Check if denseVector is in the same size
        if(MNL1 == sizeDenseVector){
            for(long int i = 0; i < MNL1; i++){
                (*result)[i] = 0;
            }
            //Send data to processus
            for(int i = 1; i < p; i++){
                int sizeData = max(0, (int)(min(rowSize, (int)((i + 1) * ceil(rowSize*1.0/p))) - i * ceil(rowSize*1.0/p)));
                MPI_Send(&row[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_LONG, i, TAGROW, MPI_COMM_WORLD);
                MPI_Send(&val[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_FLOAT, i, TAGVAL, MPI_COMM_WORLD);
                MPI_Send(&col[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_LONG, i, TAGCOL, MPI_COMM_WORLD);
            }

            for(int i = 1; i < p; i++){
                MPI_Send(denseVector, sizeDenseVector, MPI_DOUBLE, i, TAGVECTOR, MPI_COMM_WORLD);
            }
            for(long int i = 0; i < min(rowSize, (int)(ceil(rowSize*1.0/p))); i++){
                (*result)[row[i]] += ((double)(val[i]))*denseVector[col[i]];
            }
            free(row);
            free(col);
            free(val);
        }else{
        //If not in the same size, stop the function
            cout << "the vector is not of the right size" << endl;
            free(*result);
            *result = NULL;
        }
    }else{
        if(MNL1 == sizeDenseVector){
            for(long int i = 0; i < MNL1; i++){
                (*result)[i] = 0;
            }
            int sizeData = max(0, (int)(min(rowSize, (int)((my_rank + 1) * ceil(rowSize*1.0/p))) - my_rank * ceil(rowSize*1.0/p)));
            /* vector<long int> row(sizeData);
            vector<double> val(sizeData);
            vector<long int> col(sizeData); */
            long int * row = (long int *)malloc(sizeof(long int)*sizeData);
            long int * col = (long int *)malloc(sizeof(long int)*sizeData);
            float * val = (float *)malloc(sizeof(float) * sizeData);
            MPI_Recv(row,  sizeData, MPI_LONG, 0, TAGROW, MPI_COMM_WORLD, &status);
            MPI_Recv(val, sizeData, MPI_FLOAT, 0, TAGVAL, MPI_COMM_WORLD, &status);
            MPI_Recv(col, sizeData, MPI_LONG, 0, TAGCOL, MPI_COMM_WORLD, &status);
            double * vec = NULL;
            vec = (double *)malloc(sizeof(double) * sizeDenseVector);
            MPI_Recv(vec, sizeDenseVector, MPI_DOUBLE, 0, TAGVECTOR, MPI_COMM_WORLD, &status);
            for (long int  i = 0; i < sizeData; i++){
                (*result)[row[i]] += ((double)(val[i]))*vec[col[i]];
            }

            free(row);
            free(col);
            free(val);
            free(vec);
        }else{
            free(*result);
            *result = NULL;
        }
    }
    
    if(MNL1 == sizeDenseVector){
        //double * resultReduce;
        //resultReduce = (double *)malloc(sizeof(double)*MNL1);
        //MPI_Reduce(result, resultReduce, MNL1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        double * result1 = (double*)malloc(sizeof(double) * MNL1);
        MPI_Allreduce(*result, result1, MNL1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        free(*result);
        *result = result1;
        //free(result);
        //free(resultReduce);
    }
    
    

}

void spmv_mpi(CSRf * matrix, double * denseVector, int sizeDenseVector, double ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    long int MNL1;
    int valSize;
    if(my_rank == 0){
        MNL1 = matrix->getMNL(1);
        valSize = matrix->getValSize();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&MNL1, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&valSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(result != NULL){
        free(*result);
        *result = NULL;
    }
    

    typedef vector<double> RowDouble;
    //RowDouble result(MNL1);
    
    if(my_rank == 0){
        (*result) = (double *)malloc(sizeof(double)*MNL1);
        long int * ptr;
        long int * ind;
        float * val;
        ptr = (long int *)malloc(sizeof(long int)*(MNL1+1));
        ind = (long int *)malloc(sizeof(long int)*valSize);
        val = (float *)malloc(sizeof(float)* valSize);
        ptr = matrix->getPtr();
        ind = matrix->getInd();
        val = matrix->getVal();

        
        //Check if denseVector is in the same size
        if(MNL1 == sizeDenseVector){
            for(long int i = 0; i < MNL1; i++){
                (*result)[i] = 0;
            }
            //Send data to processus
            long int * indexProcessus;
            indexProcessus = (long int *)malloc(sizeof(long int)*p);
            long int actualProcessus = 0;
            long int lastInd = 0;
            for(long int i = 1; i < MNL1+1; i++){
                actualProcessus++;
                if(actualProcessus >= p){
                    actualProcessus = p;
                    //cout << "Break it" << endl;
                    break;
                }
                while(ptr[i]-ptr[lastInd] < ceil(valSize/p) and i<MNL1){
                    i++;
                }
                long int sizeData = i-lastInd;
                MPI_Send(&sizeData, 1, MPI_LONG, actualProcessus, TAGSIZE, MPI_COMM_WORLD);
                MPI_Send(&ptr[lastInd], sizeData + 1, MPI_LONG, actualProcessus, TAGPTR, MPI_COMM_WORLD);
                indexProcessus[actualProcessus] = lastInd;
                MPI_Send(&ind[ptr[lastInd]], ptr[i] - ptr[lastInd], MPI_LONG, actualProcessus, TAGIND, MPI_COMM_WORLD);
                MPI_Send(&val[ptr[lastInd]], ptr[i] - ptr[lastInd], MPI_FLOAT, actualProcessus, TAGVAL, MPI_COMM_WORLD);
                lastInd = i;
            }
            for(int i = 1; i < actualProcessus; i++){
                MPI_Send(denseVector, sizeDenseVector, MPI_DOUBLE, i, TAGVECTOR, MPI_COMM_WORLD);
            }
            if(actualProcessus >= p){
                for(int i = lastInd-1; i < MNL1; i++){
                    for(int j = 0; j < ptr[i+1]-ptr[i]; j++){
                        (*result)[i] += denseVector[ind[j+ptr[i]]] * ((double)(val[j+ptr[i]])); 
                        //result[i] +=  denseVector[ind[j+ptr[i]-1]] * val[j+ptr[i]-1];
                    }
                }
                actualProcessus = actualProcessus - 1;
            }else{
                long int sizeData = -1;
                for(int i = actualProcessus+1; i<p; i++){
                    MPI_Send(&sizeData, 1, MPI_LONG, i, TAGSIZE, MPI_COMM_WORLD);
                }
            }
            int length = 0;
            for(long int i = 0; i < actualProcessus; i++){
                MPI_Recv(&length, 1, MPI_LONG, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
                MPI_Recv(&((*result)[indexProcessus[status.MPI_SOURCE]]), length, MPI_DOUBLE, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            free(ptr);
            free(ind);
            free(val);
            free(indexProcessus);
        }else{
        //If not in the same size, stop the function
            cout << "the vector is not of the right size" << endl;

            free(ptr);
            free(ind);
            free(val);
        }
    }else{
        if(MNL1 == sizeDenseVector){
            
            
            long int sizePtr;
            MPI_Recv(&sizePtr, 1, MPI_LONG, 0, TAGSIZE, MPI_COMM_WORLD, &status);
            if(sizePtr != -1){
                long int * ptr = (long int *)malloc(sizeof(long int)*(sizePtr+1));
                MPI_Recv(ptr, sizePtr+1, MPI_LONG, 0, TAGPTR, MPI_COMM_WORLD, &status);
                long int sizeData = ptr[sizePtr]-ptr[0];
                long int * ind = (long int *)malloc(sizeof(long int)*sizeData);
                float * val = (float *)malloc(sizeof(float) * sizeData);
                MPI_Recv(ind, sizeData, MPI_LONG, 0, TAGIND, MPI_COMM_WORLD, &status);
                MPI_Recv(val, sizeData, MPI_FLOAT, 0, TAGVAL, MPI_COMM_WORLD, &status);
                *result = (double *)malloc(sizeof(double)*sizePtr);
                double * vec = NULL;
                vec = (double *)malloc(sizeof(double) * sizeDenseVector);
                MPI_Recv(vec, sizeDenseVector, MPI_DOUBLE, 0, TAGVECTOR, MPI_COMM_WORLD, &status);
                for(long int i = 0; i < sizePtr; i++){
                    (*result)[i] = 0;
                }
                for (long int  i = 0; i < sizePtr; i++){
                    for(long int j = 0; j < ptr[i+1]-ptr[i]; j++){
                        (*result)[i] +=  vec[ind[j+ptr[i]-ptr[0]]] * ((double)(val[j+ptr[i]-ptr[0]]));
                    }
                }
                MPI_Send(&sizePtr, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                MPI_Send(*result, sizePtr, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                free(*result);
                free(ptr);
                free(ind);
                free(val);
                free(vec);
                
            }else{
                cout << my_rank << " : Ok, I don't work" << endl;
            }
            *result = (double *)malloc(MNL1 * sizeof(double));
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
}

void spmv_mpi(ELLf *matrix, double * denseVector, int sizeDenseVector, double ** result){
    int my_rank;
    int p;
    int tag=0;
    long int MNL1;
    int rowSize;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(my_rank == 0){
        MNL1 = matrix->getMNL(1);
        rowSize = matrix->sizeColInd(1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&MNL1, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    

    //typedef vector<double> RowDouble;
    //RowDouble result(MNL1);
    if(*result != NULL){
        free(*result);
        *result = NULL;
    }
    if(my_rank == 0){
        /* vector<long int> row = matrix.getRow();
        vector<long int> col = matrix.getCol();
        vector<double> val = matrix.getVal(); */
        float * val;
        long int * indCol;
        indCol = (long int *)malloc(sizeof(long int) * matrix->sizeColInd(0) * matrix->sizeColInd(1));
        val = (float *)malloc(sizeof(float) * matrix->sizeColInd(0) * matrix->sizeColInd(1));
        
        *result = (double *)malloc(sizeof(double)*MNL1);
        indCol = matrix->getColInd();
        val = matrix->getVal();
        //cout << "VALEUR 0 : " << val[0] << endl;
    
        //Check if denseVector is in the same size
        if(MNL1 == sizeDenseVector){
            for(long int i = 0; i < MNL1; i++){
                (*result)[i] = 0;
            }
            //Send data to processus
            long int * indexProcessus;
            indexProcessus = (long int *)malloc(sizeof(long int)*p);
            long int actualProcessus = 0;
            long int lastInd = 0;
            for(long int i = 1; i < MNL1+1; i++){
                actualProcessus++;
                if(actualProcessus >= p){
                    actualProcessus = p;
                    break;
                }
                i = min(MNL1 - 1, (long)(lastInd + ceil(1.0*MNL1/p)));
                long int sizeData = i-lastInd;
                MPI_Send(&sizeData, 1, MPI_LONG, actualProcessus, TAGSIZE, MPI_COMM_WORLD);
                indexProcessus[actualProcessus] = lastInd;
                MPI_Send(&indCol[lastInd*rowSize], (i-lastInd)*rowSize, MPI_LONG, actualProcessus, TAGIND, MPI_COMM_WORLD);
                MPI_Send(&val[lastInd*rowSize], (i-lastInd)*rowSize, MPI_FLOAT, actualProcessus, TAGVAL, MPI_COMM_WORLD);
                //cout << "Envoie de " << sizeData << " à " << actualProcessus << " (indexProcessus = " << indexProcessus[actualProcessus] << ") length = "<< ptr[i]- ptr[lastInd] << endl;
                lastInd = i;
            }
            for(int i = 1; i < actualProcessus; i++){
                MPI_Send(denseVector, sizeDenseVector, MPI_DOUBLE, i, TAGVECTOR, MPI_COMM_WORLD);
            }
            if(actualProcessus >= p){
                for(int i = lastInd-1; i < MNL1; i++){
                    for(int j = 0; j < rowSize; j++){
                        if(indCol[j+i*rowSize] != -1){
                            (*result)[i] += denseVector[indCol[j+i*rowSize]] * ((double)(val[j+i*rowSize]));
                        }
                         
                        //result[i] +=  denseVector[ind[j+ptr[i]-1]] * val[j+ptr[i]-1];
                        //result[i] +=  val[i*maxColInd + j] * denseVector[colInd[i* maxColInd + j]];
                    }
                }
                actualProcessus = actualProcessus - 1;
            }else{
                long int sizeData = -1;
                for(int i = actualProcessus+1; i<p; i++){
                    MPI_Send(&sizeData, 1, MPI_LONG, i, TAGSIZE, MPI_COMM_WORLD);
                }
            }
            int length = 0;
            for(long int i = 0; i < actualProcessus; i++){
                MPI_Recv(&length, 1, MPI_LONG, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
                MPI_Recv(&((*result)[indexProcessus[status.MPI_SOURCE]]), length, MPI_DOUBLE, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }else{
        //If not in the same size, stop the function
            cout << "the vector is not of the right size (size = "<< sizeDenseVector << " MNL = "<< MNL1<< ")" << endl;
        }
    }else{
        if(MNL1 == sizeDenseVector){
            
            
            long int nbRowProc;
            MPI_Recv(&nbRowProc, 1, MPI_LONG, 0, TAGSIZE, MPI_COMM_WORLD, &status);
            if(nbRowProc != -1){
                //MPI_Recv(ptr, sizePtr+1, MPI_LONG, 0, TAGPTR, MPI_COMM_WORLD, &status);
                long int sizeData = nbRowProc * rowSize;
                long int * ind = (long int *)malloc(sizeof(long int)*sizeData);
                float * val = (float *)malloc(sizeof(float) * sizeData);
                MPI_Recv(ind, sizeData, MPI_LONG, 0, TAGIND, MPI_COMM_WORLD, &status);
                MPI_Recv(val, sizeData, MPI_FLOAT, 0, TAGVAL, MPI_COMM_WORLD, &status);
                *result = (double *)malloc(sizeof(double)*nbRowProc);
                for(long int i = 0; i < nbRowProc; i++){
                    (*result)[i] = 0;
                }
                double * vec = NULL;
                vec = (double *)malloc(sizeof(double) * sizeDenseVector);
                MPI_Recv(vec, sizeDenseVector, MPI_DOUBLE, 0, TAGVECTOR, MPI_COMM_WORLD, &status);
                for (long int  i = 0; i < nbRowProc; i++){
                    for(long int j = 0; j < rowSize; j++){
                        if(ind[i*rowSize+j] != -1){
                            (*result)[i] +=  vec[ind[i*rowSize+j]] * ((double)(val[i*rowSize+j]));
                        }
                    }
                }
                MPI_Send(&nbRowProc, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                MPI_Send(*result, nbRowProc, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                free(*result);
                free(ind);
                free(val);
                free(vec);
                
            }else{
                cout << my_rank << " : Oh ok, I don't work" << endl;
            }
            *result = (double *)malloc(MNL1 * sizeof(double));
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
    

}


//Matrix values : double AND denseVector : float

void spmv_mpi(COOf *matrix, float * denseVector, int sizeDenseVector, float ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    long int MNL1;
    int rowSize;
    if(my_rank == 0){
        MNL1 = matrix->getMNL(1);
        rowSize = matrix->getRowSize();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&MNL1, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    

    typedef vector<double> RowDouble;
    //RowDouble result(MNL1);
    if(*result != NULL){
        free(*result);
        *result = NULL;
    }
    *result = (float *)malloc(sizeof(float)*MNL1);
    if(my_rank == 0){
        /* vector<long int> row = matrix.getRow();
        vector<long int> col = matrix.getCol();
        vector<double> val = matrix.getVal(); */
        long int * row;
        long int * col;
        float * val;
        row = matrix->getRow();
        col = matrix->getCol();
        val = matrix->getVal();

        
        //Check if denseVector is in the same size
        if(MNL1 == sizeDenseVector){
            for(long int i = 0; i < MNL1; i++){
                (*result)[i] = 0;
            }
            //Send data to processus
            for(int i = 1; i < p; i++){
                int sizeData = max(0, (int)(min(rowSize, (int)((i + 1) * ceil(rowSize*1.0/p))) - i * ceil(rowSize*1.0/p)));
                MPI_Send(&row[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_LONG, i, TAGROW, MPI_COMM_WORLD);
                MPI_Send(&val[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_FLOAT, i, TAGVAL, MPI_COMM_WORLD);
                MPI_Send(&col[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_LONG, i, TAGCOL, MPI_COMM_WORLD);
            }

            for(int i = 1; i < p; i++){
                MPI_Send(denseVector, sizeDenseVector, MPI_FLOAT, i, TAGVECTOR, MPI_COMM_WORLD);
            }
            for(long int i = 0; i < min(rowSize, (int)(ceil(rowSize*1.0/p))); i++){
                (*result)[row[i]] += val[i]*denseVector[col[i]];
            }
            /* free(row);
            free(col);
            free(val); */
        }else{
        //If not in the same size, stop the function
            cout << "the vector is not of the right size" << endl;
            free(*result);
            *result = NULL;
        }
    }else{
        if(MNL1 == sizeDenseVector){
            for(long int i = 0; i < MNL1; i++){
                (*result)[i] = 0;
            }
            int sizeData = max(0, (int)(min(rowSize, (int)((my_rank + 1) * ceil(rowSize*1.0/p))) - my_rank * ceil(rowSize*1.0/p)));
            /* vector<long int> row(sizeData);
            vector<double> val(sizeData);
            vector<long int> col(sizeData); */
            long int * row = (long int *)malloc(sizeof(long int)*sizeData);
            long int * col = (long int *)malloc(sizeof(long int)*sizeData);
            float * val = (float *)malloc(sizeof(float) * sizeData);
            MPI_Recv(row,  sizeData, MPI_LONG, 0, TAGROW, MPI_COMM_WORLD, &status);
            MPI_Recv(val, sizeData, MPI_FLOAT, 0, TAGVAL, MPI_COMM_WORLD, &status);
            MPI_Recv(col, sizeData, MPI_LONG, 0, TAGCOL, MPI_COMM_WORLD, &status);
            float * vec = NULL;
            vec = (float *)malloc(sizeof(float) * sizeDenseVector);
            MPI_Recv(vec, sizeDenseVector, MPI_FLOAT, 0, TAGVECTOR, MPI_COMM_WORLD, &status);
            for (long int  i = 0; i < sizeData; i++){
                (*result)[row[i]] += val[i]*vec[col[i]];
            }

            free(row);
            free(col);
            free(val);
            free(vec);
        }else{
            free(*result);
            *result = NULL;
        }
    }
    
    if(MNL1 == sizeDenseVector){
        //double * resultReduce;
        //resultReduce = (double *)malloc(sizeof(double)*MNL1);
        //MPI_Reduce(result, resultReduce, MNL1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        float * result1 = (float*)malloc(sizeof(float) * MNL1);
        MPI_Allreduce(*result, result1, MNL1, MPI_FLOAT, MPI_SUM, MPI_COMM_WORLD);
        free(*result);
        *result = result1;
        //free(result);
        //free(resultReduce);
    }
    
    

}

void spmv_mpi(CSRf * matrix, float * denseVector, int sizeDenseVector, float ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    long int MNL1;
    int valSize;
    if(my_rank == 0){
        MNL1 = matrix->getMNL(1);
        valSize = matrix->getValSize();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&MNL1, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&valSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if(result != NULL){
        free(*result);
        *result = NULL;
    }
    

    typedef vector<double> RowDouble;
    //RowDouble result(MNL1);
    
    if(my_rank == 0){
        (*result) = (float *)malloc(sizeof(float)*MNL1);
        long int * ptr;
        long int * ind;
        float * val;
        ptr = matrix->getPtr();
        ind = matrix->getInd();
        val = matrix->getVal();

        
        //Check if denseVector is in the same size
        if(MNL1 == sizeDenseVector){
            for(long int i = 0; i < MNL1; i++){
                (*result)[i] = 0;
            }
            //Send data to processus
            long int * indexProcessus;
            indexProcessus = (long int *)malloc(sizeof(long int)*p);
            long int actualProcessus = 0;
            long int lastInd = 0;
            for(long int i = 1; i < MNL1+1; i++){
                actualProcessus++;
                if(actualProcessus >= p){
                    actualProcessus = p;
                    //cout << "Break it" << endl;
                    break;
                }
                while(ptr[i]-ptr[lastInd] < ceil(valSize/p) and i<MNL1){
                    i++;
                }
                long int sizeData = i-lastInd;
                MPI_Send(&sizeData, 1, MPI_LONG, actualProcessus, TAGSIZE, MPI_COMM_WORLD);
                MPI_Send(&ptr[lastInd], sizeData + 1, MPI_LONG, actualProcessus, TAGPTR, MPI_COMM_WORLD);
                indexProcessus[actualProcessus] = lastInd;
                MPI_Send(&ind[ptr[lastInd]], ptr[i] - ptr[lastInd], MPI_LONG, actualProcessus, TAGIND, MPI_COMM_WORLD);
                MPI_Send(&val[ptr[lastInd]], ptr[i] - ptr[lastInd], MPI_FLOAT, actualProcessus, TAGVAL, MPI_COMM_WORLD);
                lastInd = i;
            }
            for(int i = 1; i < actualProcessus; i++){
                MPI_Send(denseVector, sizeDenseVector, MPI_FLOAT, i, TAGVECTOR, MPI_COMM_WORLD);
            }
            if(actualProcessus >= p){
                for(int i = lastInd-1; i < MNL1; i++){
                    for(int j = 0; j < ptr[i+1]-ptr[i]; j++){
                        (*result)[i] += denseVector[ind[j+ptr[i]]] * val[j+ptr[i]]; 
                        //result[i] +=  denseVector[ind[j+ptr[i]-1]] * val[j+ptr[i]-1];
                    }
                }
                actualProcessus = actualProcessus - 1;
            }else{
                long int sizeData = -1;
                for(int i = actualProcessus+1; i<p; i++){
                    MPI_Send(&sizeData, 1, MPI_LONG, i, TAGSIZE, MPI_COMM_WORLD);
                }
            }
            int length = 0;
            for(long int i = 0; i < actualProcessus; i++){
                MPI_Recv(&length, 1, MPI_LONG, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
                MPI_Recv(&((*result)[indexProcessus[status.MPI_SOURCE]]), length, MPI_FLOAT, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_FLOAT, 0, MPI_COMM_WORLD);
            free(indexProcessus);
        }else{
        //If not in the same size, stop the function
            cout << "the vector is not of the right size" << endl;

            free(ptr);
            free(ind);
            free(val);
        }
    }else{
        if(MNL1 == sizeDenseVector){
            
            
            long int sizePtr;
            MPI_Recv(&sizePtr, 1, MPI_LONG, 0, TAGSIZE, MPI_COMM_WORLD, &status);
            if(sizePtr != -1){
                long int * ptr = (long int *)malloc(sizeof(long int)*(sizePtr+1));
                MPI_Recv(ptr, sizePtr+1, MPI_LONG, 0, TAGPTR, MPI_COMM_WORLD, &status);
                long int sizeData = ptr[sizePtr]-ptr[0];
                long int * ind = (long int *)malloc(sizeof(long int)*sizeData);
                float * val = (float *)malloc(sizeof(float) * sizeData);
                MPI_Recv(ind, sizeData, MPI_LONG, 0, TAGIND, MPI_COMM_WORLD, &status);
                MPI_Recv(val, sizeData, MPI_FLOAT, 0, TAGVAL, MPI_COMM_WORLD, &status);
                *result = (float *)malloc(sizeof(float)*sizePtr);
                float * vec = NULL;
                vec = (float *)malloc(sizeof(float) * sizeDenseVector);
                MPI_Recv(vec, sizeDenseVector, MPI_FLOAT, 0, TAGVECTOR, MPI_COMM_WORLD, &status);
                for(long int i = 0; i < sizePtr; i++){
                    (*result)[i] = 0;
                }
                for (long int  i = 0; i < sizePtr; i++){
                    for(long int j = 0; j < ptr[i+1]-ptr[i]; j++){
                        (*result)[i] +=  vec[ind[j+ptr[i]-ptr[0]]] * val[j+ptr[i]-ptr[0]];
                    }
                }
                MPI_Send(&sizePtr, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                MPI_Send(*result, sizePtr, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
                free(*result);
                free(ptr);
                free(ind);
                free(val);
                free(vec);
                
            }else{
                cout << my_rank << " : Ok, I don't work" << endl;
            }
            *result = (float *)malloc(MNL1 * sizeof(float));
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
    }
}

void spmv_mpi(ELLf *matrix, float * denseVector, int sizeDenseVector, float ** result){
    int my_rank;
    int p;
    int tag=0;
    long int MNL1;
    int rowSize;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(my_rank == 0){
        MNL1 = matrix->getMNL(1);
        rowSize = matrix->sizeColInd(1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&MNL1, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rowSize, 1, MPI_INT, 0, MPI_COMM_WORLD);

    

    //typedef vector<double> RowDouble;
    //RowDouble result(MNL1);
    if(*result != NULL){
        free(*result);
        *result = NULL;
    }
    if(my_rank == 0){
        /* vector<long int> row = matrix.getRow();
        vector<long int> col = matrix.getCol();
        vector<double> val = matrix.getVal(); */
        float * val;
        long int * indCol;
        
        *result = (float *)malloc(sizeof(float)*MNL1);
        indCol = matrix->getColInd();
        val = matrix->getVal();
        //cout << "VALEUR 0 : " << val[0] << endl;
    
        //Check if denseVector is in the same size
        if(MNL1 == sizeDenseVector){
            for(long int i = 0; i < MNL1; i++){
                (*result)[i] = 0;
            }
            //Send data to processus
            long int * indexProcessus;
            indexProcessus = (long int *)malloc(sizeof(long int)*p);
            long int actualProcessus = 0;
            long int lastInd = 0;
            for(long int i = 1; i < MNL1+1; i++){
                actualProcessus++;
                if(actualProcessus >= p){
                    actualProcessus = p;
                    break;
                }
                i = min(MNL1 - 1, (long)(lastInd + ceil(1.0*MNL1/p)));
                long int sizeData = i-lastInd;
                MPI_Send(&sizeData, 1, MPI_LONG, actualProcessus, TAGSIZE, MPI_COMM_WORLD);
                indexProcessus[actualProcessus] = lastInd;
                MPI_Send(&indCol[lastInd*rowSize], (i-lastInd)*rowSize, MPI_LONG, actualProcessus, TAGIND, MPI_COMM_WORLD);
                MPI_Send(&val[lastInd*rowSize], (i-lastInd)*rowSize, MPI_FLOAT, actualProcessus, TAGVAL, MPI_COMM_WORLD);
                //cout << "Envoie de " << sizeData << " à " << actualProcessus << " (indexProcessus = " << indexProcessus[actualProcessus] << ") length = "<< ptr[i]- ptr[lastInd] << endl;
                lastInd = i;
            }
            for(int i = 1; i < actualProcessus; i++){
                MPI_Send(denseVector, sizeDenseVector, MPI_FLOAT, i, TAGVECTOR, MPI_COMM_WORLD);
            }
            if(actualProcessus >= p){
                for(int i = lastInd-1; i < MNL1; i++){
                    for(int j = 0; j < rowSize; j++){
                        if(indCol[j+i*rowSize] != -1){
                            (*result)[i] += denseVector[indCol[j+i*rowSize]] * val[j+i*rowSize];
                        }
                         
                        //result[i] +=  denseVector[ind[j+ptr[i]-1]] * val[j+ptr[i]-1];
                        //result[i] +=  val[i*maxColInd + j] * denseVector[colInd[i* maxColInd + j]];
                    }
                }
                actualProcessus = actualProcessus - 1;
            }else{
                long int sizeData = -1;
                for(int i = actualProcessus+1; i<p; i++){
                    MPI_Send(&sizeData, 1, MPI_LONG, i, TAGSIZE, MPI_COMM_WORLD);
                }
            }
            int length = 0;
            for(long int i = 0; i < actualProcessus; i++){
                MPI_Recv(&length, 1, MPI_LONG, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
                MPI_Recv(&((*result)[indexProcessus[status.MPI_SOURCE]]), length, MPI_FLOAT, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_FLOAT, 0, MPI_COMM_WORLD);
            free(indexProcessus);
        }else{
        //If not in the same size, stop the function
            cout << "the vector is not of the right size (size = "<< sizeDenseVector << " MNL = "<< MNL1<< ")" << endl;
        }
    }else{
        if(MNL1 == sizeDenseVector){
            
            
            long int nbRowProc;
            MPI_Recv(&nbRowProc, 1, MPI_LONG, 0, TAGSIZE, MPI_COMM_WORLD, &status);
            if(nbRowProc != -1){
                //MPI_Recv(ptr, sizePtr+1, MPI_LONG, 0, TAGPTR, MPI_COMM_WORLD, &status);
                long int sizeData = nbRowProc * rowSize;
                long int * ind = (long int *)malloc(sizeof(long int)*sizeData);
                float * val = (float *)malloc(sizeof(float) * sizeData);
                MPI_Recv(ind, sizeData, MPI_LONG, 0, TAGIND, MPI_COMM_WORLD, &status);
                MPI_Recv(val, sizeData, MPI_FLOAT, 0, TAGVAL, MPI_COMM_WORLD, &status);
                *result = (float *)malloc(sizeof(float)*nbRowProc);
                for(long int i = 0; i < nbRowProc; i++){
                    (*result)[i] = 0;
                }
                float * vec = NULL;
                vec = (float *)malloc(sizeof(float) * sizeDenseVector);
                MPI_Recv(vec, sizeDenseVector, MPI_FLOAT, 0, TAGVECTOR, MPI_COMM_WORLD, &status);
                for (long int  i = 0; i < nbRowProc; i++){
                    for(long int j = 0; j < rowSize; j++){
                        if(ind[i*rowSize+j] != -1){
                            (*result)[i] +=  vec[ind[i*rowSize+j]] * val[i*rowSize+j];
                        }
                    }
                }
                MPI_Send(&nbRowProc, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                MPI_Send(*result, nbRowProc, MPI_FLOAT, 0, tag, MPI_COMM_WORLD);
                free(*result);
                free(ind);
                free(val);
                free(vec);
                
            }else{
                cout << my_rank << " : Oh ok, I don't work" << endl;
            }
            *result = (float *)malloc(MNL1 * sizeof(float));
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_FLOAT, 0, MPI_COMM_WORLD);
        }
    }
    

}

    
