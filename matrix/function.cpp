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

using namespace std;

vector<double> spmv_mpi(COO *matrix, vector<double> denseVector){
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
    double * result;
    result = (double *)malloc(sizeof(double)*MNL1);
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
        if(MNL1 == denseVector.size()){
            for(long int i = 0; i < MNL1; i++){
                result[i] = 0;
            }
            //Send data to processus
            for(int i = 1; i < p; i++){
                int sizeData = max(0, (int)(min(rowSize, (int)((i + 1) * ceil(rowSize*1.0/p))) - i * ceil(rowSize*1.0/p)));
                MPI_Send(&row[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_LONG, i, TAGROW, MPI_COMM_WORLD);
                MPI_Send(&val[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_DOUBLE, i, TAGVAL, MPI_COMM_WORLD);
                MPI_Send(&col[(int)(i * ceil(rowSize*1.0/p))], sizeData, MPI_LONG, i, TAGCOL, MPI_COMM_WORLD);
            }
            for(long int i = 0; i < min(rowSize, (int)(ceil(rowSize*1.0/p))); i++){
                result[row[i]] += val[i]*denseVector[col[i]];
            }
            free(row);
            free(col);
            free(val);
        }else{
        //If not in the same size, stop the function
            cout << "the vector is not of the right size" << endl;
            return denseVector;
        }
    }else{
        if(MNL1 == denseVector.size()){
            for(long int i = 0; i < MNL1; i++){
                result[i] = 0;
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
                result[row[i]] += val[i]*denseVector[col[i]];
            }

            free(row);
            free(col);
            free(val);
        }else{
            return denseVector;
        }
    }
    
    if(MNL1 == denseVector.size()){
        double * resultReduce;
        resultReduce = (double *)malloc(sizeof(double)*MNL1);
        //MPI_Reduce(result, resultReduce, MNL1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Allreduce(result, resultReduce, MNL1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        RowDouble resultVector(MNL1);
        for(int i = 0; i < MNL1; i++){
            resultVector[i] = resultReduce[i];
        }
        free(result);
        free(resultReduce);
        return resultVector;
    }
    
    

}

vector<double> spmv_mpi(CSR * matrix, vector<double> denseVector){
    int my_rank;
    int p;
    int tag=0;
    long int MNL1 = matrix->getMNL(1);
    int valSize = matrix->getValSize();
    double * result;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    typedef vector<double> RowDouble;
    //RowDouble result(MNL1);
    
    if(my_rank == 0){
        result = (double *)malloc(sizeof(double)*MNL1);
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
        if(MNL1 == denseVector.size()){
            for(long int i = 0; i < MNL1; i++){
                result[i] = 0;
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
                MPI_Send(&lastInd, 1, MPI_LONG, actualProcessus, TAGINDEX, MPI_COMM_WORLD);
                indexProcessus[actualProcessus] = lastInd;
                MPI_Send(&ind[ptr[lastInd]], ptr[i] - ptr[lastInd], MPI_LONG, actualProcessus, TAGIND, MPI_COMM_WORLD);
                MPI_Send(&val[ptr[lastInd]], ptr[i] - ptr[lastInd], MPI_DOUBLE, actualProcessus, TAGVAL, MPI_COMM_WORLD);
                cout << "Envoie de " << sizeData << " à " << actualProcessus << " (indexProcessus = " << indexProcessus[actualProcessus] << ") length = "<< ptr[i]- ptr[lastInd] << endl;
                lastInd = i;
            }
            if(actualProcessus >= p){
                for(int i = lastInd-1; i < MNL1; i++){
                    for(int j = 0; j < ptr[i+1]-ptr[i]; j++){
                        result[i] += denseVector[ind[j+ptr[i]]] * val[j+ptr[i]]; 
                        //result[i] +=  denseVector[ind[j+ptr[i]-1]] * val[j+ptr[i]-1];
                    }
                }
                actualProcessus = actualProcessus - 1;
            }else{
                long int sizeData = -1;
                for(int i = actualProcessus+1; i<p; i++){
                    cout << "Envoie -1 à " << i << endl;;
                    MPI_Send(&sizeData, 1, MPI_LONG, i, TAGSIZE, MPI_COMM_WORLD);
                }
            }
            int length = 0;
            for(long int i = 0; i < actualProcessus; i++){
                MPI_Recv(&length, 1, MPI_LONG, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
                MPI_Recv(&result[indexProcessus[status.MPI_SOURCE]], length, MPI_DOUBLE, status.MPI_SOURCE, tag, MPI_COMM_WORLD, &status);
            }
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }else{
        //If not in the same size, stop the function
            cout << "the vector is not of the right size" << endl;
            return denseVector;
        }
    }else{
        if(MNL1 == denseVector.size()){
            
            
            long int sizePtr;
            MPI_Recv(&sizePtr, 1, MPI_LONG, 0, TAGSIZE, MPI_COMM_WORLD, &status);
            if(sizePtr != -1){
                long int * ptr = (long int *)malloc(sizeof(long int)*(sizePtr+1));
                MPI_Recv(ptr, sizePtr+1, MPI_LONG, 0, TAGPTR, MPI_COMM_WORLD, &status);
                long int sizeData = ptr[sizePtr]-ptr[0];
                long int * ind = (long int *)malloc(sizeof(long int)*sizeData);
                double * val = (double *)malloc(sizeof(double) * sizeData);
                long int startInd;
                MPI_Recv(&startInd,  1, MPI_LONG, 0, TAGINDEX, MPI_COMM_WORLD, &status);
                MPI_Recv(ind, sizeData, MPI_LONG, 0, TAGIND, MPI_COMM_WORLD, &status);
                MPI_Recv(val, sizeData, MPI_DOUBLE, 0, TAGVAL, MPI_COMM_WORLD, &status);
                result = (double *)malloc(sizeof(double)*sizePtr);
                for(long int i = 0; i < sizePtr; i++){
                    result[i] = 0;
                }
                for (long int  i = 0; i < sizePtr; i++){
                    for(long int j = 0; j < ptr[i+1]-ptr[i]; j++){
                        result[i] +=  denseVector[ind[j+ptr[i]-ptr[0]]] * val[j+ptr[i]-ptr[0]];
                    }
                }
                MPI_Send(&sizePtr, 1, MPI_LONG, 0, tag, MPI_COMM_WORLD);
                MPI_Send(result, sizePtr, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD);
                free(result);
                
            }else{
                cout << my_rank << " : Ok, I don't work" << endl;
            }
            result = (double *)malloc(MNL1 * sizeof(double));
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(result, MNL1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        }else{
            return denseVector;
        }
    }    
    
    if(MNL1 == denseVector.size()){
        RowDouble resultVector(MNL1);
        for(int i = 0; i < MNL1; i++){
            resultVector[i] = result[i];
        }
        return resultVector;
    }
}
    
