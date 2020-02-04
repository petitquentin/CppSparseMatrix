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

using namespace std;

void spmv_mpi(COO *matrix, double * denseVector, int sizeDenseVector, double ** result){
    int my_rank;
    int p;
    int tag=0;
    long int MNL1 = matrix->getMNL(1);
    int rowSize = matrix->getRowSize();

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

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
        double * val;
        row = (long int *)malloc(sizeof(long int)*rowSize);
        col = (long int *)malloc(sizeof(long int)*rowSize);
        val = (double *)malloc(sizeof(double)* rowSize);
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
                MPI_Send(&val[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_DOUBLE, i, TAGVAL, MPI_COMM_WORLD);
                MPI_Send(&col[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_LONG, i, TAGCOL, MPI_COMM_WORLD);
            }
            for(long int i = 0; i < min(rowSize, (int)(ceil(rowSize*1.0/p))); i++){
                (*result)[row[i]] += val[i]*denseVector[col[i]];
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
            double * val = (double *)malloc(sizeof(double) * sizeData);
            MPI_Recv(row,  sizeData, MPI_LONG, 0, TAGROW, MPI_COMM_WORLD, &status);
            MPI_Recv(val, sizeData, MPI_DOUBLE, 0, TAGVAL, MPI_COMM_WORLD, &status);
            MPI_Recv(col, sizeData, MPI_LONG, 0, TAGCOL, MPI_COMM_WORLD, &status);
            for (long int  i = 0; i < sizeData; i++){
                (*result)[row[i]] += val[i]*denseVector[col[i]];
            }

            free(row);
            free(col);
            free(val);
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

void spmv_mpi(CSR * matrix, double * denseVector, int sizeDenseVector, double ** result){
    int my_rank;
    int p;
    int tag=0;
    long int MNL1 = matrix->getMNL(1);
    int valSize = matrix->getValSize();
    if(result != NULL){
        free(*result);
        *result = NULL;
    }
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    typedef vector<double> RowDouble;
    //RowDouble result(MNL1);
    
    if(my_rank == 0){
        (*result) = (double *)malloc(sizeof(double)*MNL1);
        long int * ptr;
        long int * ind;
        double * val;
        ptr = (long int *)malloc(sizeof(long int)*(MNL1+1));
        ind = (long int *)malloc(sizeof(long int)*valSize);
        val = (double *)malloc(sizeof(double)* valSize);
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
                    cout << "Break it" << endl;
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
                MPI_Send(&val[ptr[lastInd]], ptr[i] - ptr[lastInd], MPI_DOUBLE, actualProcessus, TAGVAL, MPI_COMM_WORLD);
                lastInd = i;
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
                double * val = (double *)malloc(sizeof(double) * sizeData);
                MPI_Recv(ind, sizeData, MPI_LONG, 0, TAGIND, MPI_COMM_WORLD, &status);
                MPI_Recv(val, sizeData, MPI_DOUBLE, 0, TAGVAL, MPI_COMM_WORLD, &status);
                *result = (double *)malloc(sizeof(double)*sizePtr);
                for(long int i = 0; i < sizePtr; i++){
                    (*result)[i] = 0;
                }
                for (long int  i = 0; i < sizePtr; i++){
                    for(long int j = 0; j < ptr[i+1]-ptr[i]; j++){
                        (*result)[i] +=  denseVector[ind[j+ptr[i]-ptr[0]]] * val[j+ptr[i]-ptr[0]];
                    }
                }
                MPI_Send(&sizePtr, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                MPI_Send(*result, sizePtr, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                free(*result);
                free(ptr);
                free(ind);
                free(val);
                
            }else{
                cout << my_rank << " : Ok, I don't work" << endl;
            }
            *result = (double *)malloc(MNL1 * sizeof(double));
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
}

void spmv_mpi(ELL *matrix, double * denseVector, int sizeDenseVector, double ** result){
    int my_rank;
    int p;
    int tag=0;
    long int MNL1 = matrix->getMNL(1);
    int rowSize = matrix->sizeColInd(1);

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

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
        double * val;
        long int * indCol;
        indCol = (long int *)malloc(sizeof(long int) * matrix->sizeColInd(0) * matrix->sizeColInd(1));
        val = (double *)malloc(sizeof(double) * matrix->sizeColInd(0) * matrix->sizeColInd(1));
        
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
                MPI_Send(&val[lastInd*rowSize], (i-lastInd)*rowSize, MPI_DOUBLE, actualProcessus, TAGVAL, MPI_COMM_WORLD);
                //cout << "Envoie de " << sizeData << " à " << actualProcessus << " (indexProcessus = " << indexProcessus[actualProcessus] << ") length = "<< ptr[i]- ptr[lastInd] << endl;
                lastInd = i;
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
                double * val = (double *)malloc(sizeof(double) * sizeData);
                MPI_Recv(ind, sizeData, MPI_LONG, 0, TAGIND, MPI_COMM_WORLD, &status);
                MPI_Recv(val, sizeData, MPI_DOUBLE, 0, TAGVAL, MPI_COMM_WORLD, &status);
                *result = (double *)malloc(sizeof(double)*nbRowProc);
                for(long int i = 0; i < nbRowProc; i++){
                    (*result)[i] = 0;
                }
                for (long int  i = 0; i < nbRowProc; i++){
                    for(long int j = 0; j < rowSize; j++){
                        if(ind[i*rowSize+j] != -1){
                            (*result)[i] +=  denseVector[ind[i*rowSize+j]] * val[i*rowSize+j];
                        }
                    }
                }
                MPI_Send(&nbRowProc, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                MPI_Send(*result, nbRowProc, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                free(*result);
                free(ind);
                free(val);
                
            }else{
                cout << my_rank << " : Oh ok, I don't work" << endl;
            }
            *result = (double *)malloc(MNL1 * sizeof(double));
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(*result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }
    }
    

}

void spmv_mpi(SGP *matrix, double * denseVector, int sizeDenseVector, double ** result){
    int my_rank;
    int p;
    int my_x, my_y, sizeX, sizeY;
    int tag=0;
    long int MNL1;
    int maxInd;
    MPI_Status status;
    MPI_Request request;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    if(my_rank == 0){
        MNL1 = matrix->getMNL(1);
        maxInd = matrix->getMaxInd();
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&MNL1, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&maxInd, 1, MPI_INT, 0, MPI_COMM_WORLD);

    double * acr = NULL;
    long int * ai = NULL;
    long int * jc = NULL;
    double * acc = NULL;
    long int * aj = NULL;
    long int * ic = NULL;
    double * vecX = NULL;


    typedef vector<double> RowDouble;
    //RowDouble result(MNL1);
    if(*result != NULL){
        free(*result);
        *result = NULL;
    }
    int xGrid;
    int yGrid;
    if(p <= MNL1 * 2){
        xGrid = ceil(1.0*MNL1/p);
        yGrid = maxInd;
    }else{
        xGrid = 1;
        yGrid = ceil(maxInd/(1.0*p/MNL1));
    }
    //cout << "xGrid = " << xGrid << " yGrid = " << yGrid << " MNL1 = " << MNL1 << endl;
    //cout << "COUCOU " << my_rank << endl;
    if(sizeDenseVector == MNL1){
        if(my_rank == 0){
            /* vector<long int> row = matrix.getRow();
            vector<long int> col = matrix.getCol();
            vector<double> val = matrix.getVal(); */
            
            long int * aic = (long int *)malloc(sizeof(long int) * MNL1 * maxInd);
            long int * ajc = (long int *)malloc(sizeof(long int) * MNL1 * maxInd);
            long int * icc = (long int *)malloc(sizeof(long int) * MNL1 * maxInd);
            long int * jcc = (long int *)malloc(sizeof(long int) * MNL1 * maxInd);
            double * acrc = (double *)malloc(sizeof(double) * MNL1 * maxInd);
            double * accc = (double *)malloc(sizeof(double) * MNL1 * maxInd);

            aic = matrix->getAi();
            ajc = matrix->getAj();
            icc = matrix->getIc();
            jcc = matrix->getJc();
            acrc = matrix->getAcr();
            accc = matrix->getAcc();

            for(int i = 0; i < maxInd; i ++){
                for(int j = 0; j < MNL1; j ++){
                    int k = i*MNL1+j;
                    int actualProcessus = floor(i/yGrid) * (MNL1/xGrid) + floor(j/xGrid);
                    if(actualProcessus != 0){
                        //cout << actualProcessus << endl;
                        MPI_Send(&k, 1, MPI_INT, actualProcessus, TAGSIZE, MPI_COMM_WORLD);
                        MPI_Send(&aic[k], 1, MPI_LONG, actualProcessus, TAGAI, MPI_COMM_WORLD);
                        MPI_Send(&jcc[k], 1, MPI_LONG, actualProcessus, TAGJC, MPI_COMM_WORLD);
                        MPI_Send(&acrc[k], 1, MPI_DOUBLE, actualProcessus, TAGACR, MPI_COMM_WORLD);
                        MPI_Send(&ajc[j*maxInd+i], 1, MPI_LONG, actualProcessus, TAGAJ, MPI_COMM_WORLD);
                        MPI_Send(&icc[j*maxInd+i], 1, MPI_LONG, actualProcessus, TAGIC, MPI_COMM_WORLD);
                        MPI_Send(&accc[j*maxInd+i], 1, MPI_DOUBLE, actualProcessus, TAGACC, MPI_COMM_WORLD);
                    }
                }
            }

            for(int i = 0; i < floor(maxInd/yGrid); i++){
                for(int j = 0; j < MNL1; j++){
                    if(floor(i * MNL1/xGrid) + floor(j/xGrid) != 0){
                        MPI_Send(&denseVector[j], 1, MPI_DOUBLE, floor(i * MNL1/xGrid) + floor(j/xGrid), j%xGrid, MPI_COMM_WORLD);
                    }
                }
            }
            sizeX = xGrid;
            sizeY = yGrid;
            my_x = 0;
            my_y = 0;
            acr = (double *)malloc(sizeof(double) * sizeX * sizeY);
            ai = (long int *)malloc(sizeof(long int) * sizeX * sizeY);
            jc = (long int *)malloc(sizeof(long int) * sizeY * sizeX);
            acc = (double *)malloc(sizeof(double) * sizeX * sizeY);
            aj = (long int *)malloc(sizeof(long int) * sizeX * sizeY);
            ic = (long int *)malloc(sizeof(long int) * sizeY * sizeX);
            vecX = (double *)malloc(sizeof(double) * sizeX);
            for(int i = 0; i < sizeY; i++){
                for(int j = 0; j < sizeX; j++){
                    acr[i*sizeX + j] = acrc[i*MNL1 + j];
                    ai[i*sizeX + j] = aic[i*MNL1 + j];
                    jc[i*sizeX + j] = jcc[i*MNL1 + j];
                    acc[i*sizeX + j] = accc[j*maxInd+i];
                    aj[i*sizeX + j] = ajc[j*maxInd+i];
                    ic[i*sizeX + j] = icc[j*maxInd+i];
                    vecX[j] = denseVector[j];
                }
            }
            free(acrc);
            free(accc);
            free(aic);
            free(ajc);
            free(icc);
            free(jcc);
        
            /* if(MNL1 == sizeDenseVector){
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
                    MPI_Send(&val[lastInd*rowSize], (i-lastInd)*rowSize, MPI_DOUBLE, actualProcessus, TAGVAL, MPI_COMM_WORLD);
                    //cout << "Envoie de " << sizeData << " à " << actualProcessus << " (indexProcessus = " << indexProcessus[actualProcessus] << ") length = "<< ptr[i]- ptr[lastInd] << endl;
                    lastInd = i;
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
                    MPI_Recv(&((*result)[indexProcessus[status.MPI_SOURCE]]), length, MPI_DOUBLE, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
                }
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_Bcast(*result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
            } */
        }else{
            my_x = xGrid * (my_rank % ((int)(ceil(1.0*MNL1/xGrid))));
            my_y = 1.0 *yGrid*floor((my_rank*1.0/ceil(MNL1*1.0/xGrid)));
            //sizeX = max(0, min((int)(MNL1), (my_x+1) * xGrid) - my_x * xGrid);
            sizeX = max(0, min((int)(MNL1), my_x + xGrid) - my_x);
            sizeY = max(0, min(maxInd, (my_y + 1) * yGrid) - my_y * yGrid);
            //cout << my_rank << ": my_x : " << my_x << " my_y : "<< my_y << " sizeX : " << sizeX << " sizeY : " << sizeY << endl;
            acr = (double *)malloc(sizeof(double) * sizeX * sizeY);
            ai = (long int *)malloc(sizeof(long int) * sizeX * sizeY);
            jc = (long int *)malloc(sizeof(long int) * sizeY * sizeX);
            acc = (double *)malloc(sizeof(double) * sizeX * sizeY);
            aj = (long int *)malloc(sizeof(long int) * sizeX * sizeY);
            ic = (long int *)malloc(sizeof(long int) * sizeY * sizeX);
            int k;
            //cout << my_rank << " : sizeX : " << sizeX << " sizeY : " << sizeY << endl; 
            for(int i = 0; i < sizeX * sizeY; i++){
                MPI_Recv(&k, 1, MPI_LONG, 0, TAGSIZE, MPI_COMM_WORLD, &status);
                //int ind = floor((k-(my_y*yGrid*MNL1 + my_x))/MNL1)*sizeX + ((k-(my_y*yGrid*MNL1 + my_x)))%sizeX;
                MPI_Recv(&ai[i], 1, MPI_LONG, 0, TAGAI, MPI_COMM_WORLD, &status);
                MPI_Recv(&jc[i], 1, MPI_LONG, 0, TAGJC, MPI_COMM_WORLD, &status);
                MPI_Recv(&acr[i], 1, MPI_DOUBLE, 0, TAGACR, MPI_COMM_WORLD, &status);
                MPI_Recv(&aj[i], 1, MPI_LONG, 0, TAGAJ, MPI_COMM_WORLD, &status);
                MPI_Recv(&ic[i], 1, MPI_LONG, 0, TAGIC, MPI_COMM_WORLD, &status);
                MPI_Recv(&acc[i], 1, MPI_DOUBLE, 0, TAGACC, MPI_COMM_WORLD, &status);
            }
            vecX = (double *)malloc(sizeof(double) * sizeX);
            if(sizeX * sizeY != 0){
                for(int i = 0; i < sizeX; i++){
                    MPI_Recv(&vecX[i], 1, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &status);
                }
            }
            /* cout << my_rank << "ACC : ";
            for(int i =0; i < sizeY * sizeX; i++){
                cout << acc[i] << ' ';
            }
            cout << endl;
            cout << my_rank << "ACR : ";
            for(int i =0; i < sizeY * sizeX; i++){
                cout << acr[i] << ' ';
            }
            cout << endl; */
                
        }
        MPI_Barrier(MPI_COMM_WORLD); 
        //cout << "Tous le monde est là" << endl;
        double * t = (double *)malloc(sizeof(double) * sizeX * sizeY);
        double * tres = (double *)malloc(sizeof(double) * sizeX * sizeY);
        for(int i = 0; i < sizeX * sizeY; i++){
            tres[i] = 0;
            t[i] = 0;
        }
        for(int i = 0; i < sizeY; i++){
            for(int j = 0; j < sizeX; j++){
                if(ai[i*sizeX+j] != -1){
                    //t[i*sizeX + j] += vecX[my_x+j] * acr[i*sizeX + j];
                    t[i*sizeX + j] += vecX[j] * acr[i*sizeX + j];
                    int nombre = (my_x + j) * maxInd + (my_y + i);
                    
                    //cout << "T = " << t[i*sizeX+j] << endl;
                    //cout << my_rank << " : (i " << i << " j "<< j << ") rank " << floor(jc[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(ai[i*sizeX + j]/xGrid) << " Jc : " << jc[i*sizeX + j] << " ai : " << ai[i*sizeX+j] <<endl; 
                    MPI_Isend(&t[i*sizeX + j], 1, MPI_DOUBLE, floor(jc[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(ai[i*sizeX + j]/xGrid), nombre, MPI_COMM_WORLD, &request);
                    //cout << my_rank << " ENVOI A " << floor(jc[i*sizeX + j]/yGrid) * (1.0*MNL1/xGrid) + floor(ai[i*sizeX + j]/xGrid) << " (" <<  nombre << ")" << endl;
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); 
        cout << "Tous le monde est là2" << endl;
        
        for(int i = 0; i < sizeY; i++){
            for(int j = 0; j < sizeX; j++){
                if(aj[i*sizeX+j] != -1){
                    double recep;
                    //cout << my_rank << " ATTEND " << floor(ic[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(aj[i*sizeX + j]/xGrid)<< " (" << aj[i*sizeX+j] * 7+ 3*ic[i*sizeX+j] << ")    " << aj[i*sizeX + j] << ' ' << i << j  << endl;
                    //cout << my_rank << " : rank : " << floor(ic[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(aj[i*sizeX + j]/xGrid) << " ic " << ic[i*sizeX+j] << " aj " << aj[i*sizeX + j] << endl;
                    MPI_Recv(&tres[i*sizeX+j], 1, MPI_DOUBLE, floor(ic[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(aj[i*sizeX + j]/xGrid), aj[i*sizeX+j] * maxInd+ ic[i*sizeX+j] , MPI_COMM_WORLD, &status);
                    //cout << "reception de "<< floor(ic[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(aj[i*sizeX + j]/xGrid) << " par " << my_rank << endl;
                    //cout << my_rank << " ATT OK " << floor(ic[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(aj[i*sizeX + j]/xGrid)<< " (" << aj[i*sizeX+j] * 7+ 3*ic[i*sizeX+j] << ")"<< endl;
                }
            }
        }
        /* for(int i = 0; i < sizeY; i++){
            for(int j = 0; j < sizeX; j++){
                if(ai[i*sizeX+j] != -1){
                    //t[i*sizeX + j] += vecX[my_x+j] * acr[i*sizeX + j];
                    t[i*sizeX + j] += vecX[j] * acr[i*sizeX + j];
                    
                    cout << "T = " << t[i*sizeX+j] << endl;
                    //cout << my_rank << " : (i " << i << " j "<< j << ") rank " << floor(jc[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(ai[i*sizeX + j]/xGrid) << " Jc : " << jc[i*sizeX + j] << " ai : " << ai[i*sizeX+j] <<endl; 
                    MPI_Isend(&t[i*sizeX + j], 1, MPI_DOUBLE, floor(jc[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(ai[i*sizeX + j]/xGrid), TAGEXCHANGE, MPI_COMM_WORLD, &request);
                    cout << my_rank << " ENVOI A " << floor(jc[i*sizeX + j]/yGrid) * (1.0*MNL1/xGrid) + floor(ai[i*sizeX + j]/xGrid) << endl;
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD); 
        cout << "Tous le monde est là2" << endl;
        
        for(int i = 0; i < sizeY; i++){
            for(int j = 0; j < sizeX; j++){
                if(aj[i*sizeX+j] != -1){
                    double recep;
                    cout << my_rank << " ATTEND " << floor(ic[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(aj[i*sizeX + j]/xGrid)<< endl;
                    cout << my_rank << " : rank : " << floor(ic[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(aj[i*sizeX + j]/xGrid) << " ic " << ic[i*sizeX+j] << " aj " << aj[i*sizeX + j] << endl;
                    MPI_Recv(&tres[i*sizeX+j], 1, MPI_DOUBLE, floor(ic[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(aj[i*sizeX + j]/xGrid), TAGEXCHANGE, MPI_COMM_WORLD, &status);
                    cout << "reception de "<< floor(ic[i*sizeX + j]/yGrid) * (MNL1/xGrid) + floor(aj[i*sizeX + j]/xGrid) << " par " << my_rank << endl;
                }
            }
        } */
        MPI_Barrier(MPI_COMM_WORLD);
        cout << "Tous le monde est là3" << endl;
        free(t);
        if(my_y == 0){

            //cout << my_rank << " AVANT : ";
            /* for(int i = 0; i < sizeX * sizeY; i ++){
                cout << tres[i] << " ";
            } */
            //cout << endl;
            //cout << my_rank << " PASSE PAR LA" << endl;
            double * r = (double *)malloc(sizeof(double) * sizeX);
            for(int i = 1; i < ceil(maxInd*1.0/yGrid); i++){
                MPI_Recv(r, sizeX, MPI_DOUBLE, MPI_ANY_SOURCE, TAGREDUCE, MPI_COMM_WORLD, &status);
                for(int j = 0; j < sizeX; j++){
                    tres[j] = tres[j] + r[j]; 
                }
            }
            for(int i = 1; i < sizeY; i++){
                for(int j = 0; j < sizeX; j++){
                    tres[j] = tres[j] + tres[i*sizeX + j]; 
                }
            }
            free(r);
            for(int i = 1; i < ceil(maxInd*1.0/yGrid); i++){
                MPI_Send(tres, sizeX, MPI_DOUBLE, my_x + i*ceil(MNL1/sizeX), TAGSPREAD, MPI_COMM_WORLD);
            }
            for(int i = 1; i < sizeY; i++){
                for(int j = 0; j < sizeX; j++){
                    tres[j+ i*sizeX] = tres[j]; 
                }
            }
            /* cout << my_rank << " : ";
            for(int i = 0; i < sizeX * sizeY; i ++){
                cout << tres[i] << " ";
            }
            cout << endl; 
             */

        }else{
            if(sizeX != 0 && sizeY != 0){
                MPI_Send(tres, sizeX, MPI_DOUBLE, my_x/xGrid, TAGREDUCE, MPI_COMM_WORLD);
                MPI_Recv(tres, sizeX, MPI_DOUBLE, my_x/xGrid, TAGSPREAD, MPI_COMM_WORLD, &status);
                for(int i = 1; i < sizeY; i++){
                    for(int j = 0; j < sizeX; i++){
                        tres[i*sizeX + j] = tres[j];
                    }
                }
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        //Spread all the result vector to all processus
        *result = (double *)malloc(sizeof(double)*MNL1);
        if(my_rank == 0){
            for(int i = 0; i < sizeX; i++){
                (*result)[i] = tres[i];
            }
            double * res = (double *)malloc(sizeof(double) * sizeX);
            for(int i = 1; i < ceil(1.0*MNL1/sizeX); i ++){
                MPI_Recv(res, sizeX, MPI_DOUBLE, MPI_ANY_SOURCE, TAGGET, MPI_COMM_WORLD, &status);
                for(int j = 0; j < sizeX; j++){
                    (*result)[status.MPI_SOURCE*sizeX + j] = res[j];
                    //cout << "RES FROM " << status.MPI_SOURCE << " = " << res[j] << endl;
                }
            }
            free(res);
        }else{
            if(my_y == 0){
                MPI_Send(tres, sizeX, MPI_DOUBLE, 0, TAGGET, MPI_COMM_WORLD);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Bcast(*result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size (size = "<< sizeDenseVector << " MNL = "<< MNL1<< ")" << endl;
        }
    }

}
    
