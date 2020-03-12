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

vector<double> spmv_mpi(int argc, char** argv, COO matrix, vector<double> denseVector){
    int my_rank;
    int p;
    int tag=0;
    long int MNL1 = matrix.getMNL(1);
    int rowSize = matrix.getRowSize();

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    typedef vector<double> RowDouble;
    RowDouble result(MNL1);
    if(my_rank == 0){
        vector<long int> row = matrix.getRow();
        vector<long int> col = matrix.getCol();
        vector<double> val = matrix.getVal();
        
        //Check if denseVector is in the same size
        if(MNL1 == denseVector.size()){
            for(long int i = 0; i < result.size(); i++){
                result[i] = 0;
            }
            //Send data to processus
            for(int i = 1; i < p; i++){
                int sizeData = min(rowSize, (int)((i + 1) * ceil(rowSize/p))) - i * ceil(rowSize/p);
                MPI_Send(&row[i * ceil(rowSize/p)], sizeof(long int) * sizeData, MPI_LONG, i, tag, MPI_COMM_WORLD);
                MPI_Send(&val[i * ceil(rowSize/p)], sizeof(double) * sizeData, MPI_DOUBLE, i, tag, MPI_COMM_WORLD);
                MPI_Send(&col[i * ceil(rowSize/p)], sizeof(long int) * sizeData, MPI_LONG, i, tag, MPI_COMM_WORLD);
            }
            for(long int i = 0; i < min(rowSize, (int)(ceil(rowSize/p))); i++){
                result[row[i]] += val[i]*denseVector[col[i]];
            }
        }else{
        //If not in the same size, stop the function
            for(long int i = 0; i < result.size(); i++){
                result[i] = NULL;
            }
            cout << "the vector is not of the right size" << endl;
        }
    }else{
        if(MNL1 == denseVector.size()){
            for(long int i = 0; i < result.size(); i++){
                result[i] = 0;
            }
            int sizeData = min(rowSize, (int)((my_rank + 1) * ceil(rowSize/p))) - my_rank * ceil(rowSize/p);
            vector<long int> row(sizeData);
            vector<double> val(sizeData);
            vector<long int> col(sizeData);
            MPI_Recv(&row, sizeof(long int) * sizeData, MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&val, sizeof(double) * sizeData, MPI_DOUBLE, 0, tag, MPI_COMM_WORLD, &status);
            MPI_Recv(&col, sizeof(long int) * sizeData, MPI_LONG, 0, tag, MPI_COMM_WORLD, &status);
            for (long int  i = 0; i < sizeData; i++){
                result[row[i]] += val[i]*denseVector[col[i]];
            }
        }
    }
    if(MNL1 == denseVector.size()){
        MPI_Reduce(&result, &result, sizeof(result), MPI_DOUBLE, MPI_SUM, tag, MPI_COMM_WORLD);
    }
    if(my_rank ==0){
        return result;
    }

}