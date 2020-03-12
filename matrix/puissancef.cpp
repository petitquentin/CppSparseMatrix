#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <matrix/matrix.hpp>
#include <mpi.h>
#include <unistd.h>

#define EPSILON 1e-5

//Double vector

void * puissance(ELLf *matrix, double * vecX, int sizeVecX, double ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    /* if(result != NULL){
        free(result);
        result = NULL;
    } */
    int MNL1 = matrix->getMNL(1);
    if(sizeVecX == MNL1){
        double * X = (double *)malloc(sizeof(double) * MNL1);
        double * Y = NULL;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        double error = INFINITY;
        int nbIt = 0;
        double norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            if(p != 1){
                spmv_mpi(matrix, X, sizeVecX, &Y);
            }else{
                matrix->spmv(X, sizeVecX, &Y);
            }
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += Y[i] * Y[i];
            }
            norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/((float)norm);
            }
            error = 0;
            for(int i =0; i < sizeVecX; i++){
                error += (X[i] - Y[i])*(X[i] - Y[i]);
            }
            error = sqrt(error);
            for(int i =0; i < sizeVecX; i++){
                X[i] = Y[i];
            }
            free(Y);
            Y = NULL;
            nbIt++;
        }
        if(my_rank == 0){
            cout << "Nombre d'itérations : " << nbIt << endl;
            cout << "Norme en sortie : " << norm << endl;
            cout << "Erreur en sortie : " << error << endl;
        }
        *result = X;
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size" << endl;
        }
    }
}

void * puissance(COOf *matrix, double * vecX, int sizeVecX, double ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    /* if(result != NULL){
        free(result);
        result = NULL;
    } */
    int MNL1 = matrix->getMNL(1);
    if(sizeVecX == MNL1){
        double * X = (double *)malloc(sizeof(double) * MNL1);
        
        double * Y = NULL;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        double error = INFINITY;
        int nbIt = 0;
        double norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            /* cout << "X = ";
            for(int i = 0; i < sizeVecX; i ++){
                cout << X[i] << ' ';
            }
            cout << endl; */
            if(p != 1){
                spmv_mpi(matrix, X, sizeVecX, &Y);
            }else{
                matrix->spmv(X, sizeVecX, &Y);
            }
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += Y[i] * Y[i];
            }
            norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/norm;
            }
            error = 0;
            for(int i =0; i < sizeVecX; i++){
                error += (X[i] - Y[i])*(X[i] - Y[i]);
            }
            error = sqrt(error);
            for(int i =0; i < sizeVecX; i++){
                X[i] = Y[i];
            }
            free(Y);
            Y = NULL;
            nbIt++;
        }
        if(my_rank == 0){
            cout << "Nombre d'itérations : " << nbIt << endl;
            cout << "Norme en sortie : " << norm << endl;
            cout << "Erreur en sortie : " << error << endl;
        }
        *result = X;
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size" << endl;
        }
    }
}

void * puissance(CSRf *matrix, double * vecX, int sizeVecX, double ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int MNL1 = matrix->getMNL(1);
    /* if(result != NULL){
        free(result);
        result = NULL;
    } */
    if(sizeVecX == MNL1){
        double * X = (double *)malloc(sizeof(double) * MNL1);
        
        double * Y = NULL;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        double error = INFINITY;
        int nbIt = 0;
        double norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            /* cout << "X = ";
            for(int i = 0; i < 8; i ++){
                cout << X[i] << ' ';
            }
            cout << endl; */
            //spmv_mpi(matrix, X, sizeVecX, &Y);
            if(p != 1){
                spmv_mpi(matrix, X, sizeVecX, &Y);
            }else{
                matrix->spmv(X, sizeVecX, &Y);
            }
            //cout << matrix->norm() << " : NORME "<< endl;
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += Y[i] * Y[i];
            }
            norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = 1.0*Y[i]/norm;
            }
            error = 0;
            for(int i =0; i < sizeVecX; i++){
                error += (X[i] - Y[i])*(X[i] - Y[i]);
            }
            error = sqrt(error);
            for(int i =0; i < sizeVecX; i++){
                X[i] = Y[i];
            }
            free(Y);
            Y = NULL;
            nbIt++;
        }
        if(my_rank == 0){
            cout << "Nombre d'itérations : " << nbIt << endl;
            cout << "Norme en sortie : " << norm << endl;
            cout << "Erreur en sortie : " << error << endl;
        }
        *result = X;
        MPI_Barrier(MPI_COMM_WORLD);
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size" << endl;
        }
    }
}

//Float vector

void * puissance(ELLf *matrix, float * vecX, int sizeVecX, float ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    /* if(result != NULL){
        free(result);
        result = NULL;
    } */
    int MNL1 = matrix->getMNL(1);
    if(sizeVecX == MNL1){
        float * X = (float *)malloc(sizeof(float) * MNL1);
        float * Y = NULL;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        float error = INFINITY;
        int nbIt = 0;
        float norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            if(p != 1){
                spmv_mpi(matrix, X, sizeVecX, &Y);
            }else{
                matrix->spmv(X, sizeVecX, &Y);
            }
            
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += Y[i] * Y[i];
            }
            norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/norm;
            }
            error = 0;
            for(int i =0; i < sizeVecX; i++){
                error += (X[i] - Y[i])*(X[i] - Y[i]);
            }
            error = sqrt(error);
            for(int i =0; i < sizeVecX; i++){
                X[i] = Y[i];
            }
            free(Y);
            Y = NULL;
            nbIt++;
        }
        if(my_rank == 0){
            cout << "Nombre d'itérations : " << nbIt << endl;
            cout << "Norme en sortie : " << norm << endl;
            cout << "Erreur en sortie : " << error << endl;
        }
        *result = X;
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size" << endl;
        }
    }
}

void * puissance(COOf *matrix, float * vecX, int sizeVecX, float ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    /* if(result != NULL){
        free(result);
        result = NULL;
    } */
    int MNL1 = matrix->getMNL(1);
    if(sizeVecX == MNL1){
        float * X = (float *)malloc(sizeof(float) * MNL1);
        
        float * Y = NULL;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        float error = INFINITY;
        int nbIt = 0;
        float norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            /* cout << "X = ";
            for(int i = 0; i < sizeVecX; i ++){
                cout << X[i] << ' ';
            }
            cout << endl; */
            if(p != 1){
                spmv_mpi(matrix, X, sizeVecX, &Y);
            }else{
                matrix->spmv(X, sizeVecX, &Y);
            }
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += Y[i] * Y[i];
            }
            norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/norm;
            }
            error = 0;
            for(int i =0; i < sizeVecX; i++){
                error += (X[i] - Y[i])*(X[i] - Y[i]);
            }
            error = sqrt(error);
            for(int i =0; i < sizeVecX; i++){
                X[i] = Y[i];
            }
            free(Y);
            Y = NULL;
            nbIt++;
        }
        if(my_rank == 0){
            cout << "Nombre d'itérations : " << nbIt << endl;
            cout << "Norme en sortie : " << norm << endl;
            cout << "Erreur en sortie : " << error << endl;
        }
        *result = X;
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size" << endl;
        }
    }
}

void * puissance(CSRf *matrix, float * vecX, int sizeVecX, float ** result){
    int my_rank;
    int p;
    int tag=0;
    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    int MNL1 = matrix->getMNL(1);
    /* if(result != NULL){
        free(result);
        result = NULL;
    } */
    if(sizeVecX == MNL1){
        float * X = (float *)malloc(sizeof(float) * MNL1);
        
        float * Y = NULL;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        float error = INFINITY;
        int nbIt = 0;
        float norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            /* cout << "X = ";
            for(int i = 0; i < 8; i ++){
                cout << X[i] << ' ';
            }
            cout << endl; */
            //spmv_mpi(matrix, X, sizeVecX, &Y);
            if(p != 1){
                spmv_mpi(matrix, X, sizeVecX, &Y);
            }else{
                matrix->spmv(X, sizeVecX, &Y);
            }
            //cout << matrix->norm() << " : NORME "<< endl;
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += Y[i] * Y[i];
            }
            norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = 1.0*Y[i]/norm;
            }
            error = 0;
            for(int i =0; i < sizeVecX; i++){
                error += (X[i] - Y[i])*(X[i] - Y[i]);
            }
            error = sqrt(error);
            for(int i =0; i < sizeVecX; i++){
                X[i] = Y[i];
            }
            free(Y);
            Y = NULL;
            nbIt++;
        }
        if(my_rank == 0){
            cout << "Nombre d'itérations : " << nbIt << endl;
            cout << "Norme en sortie : " << norm << endl;
            cout << "Erreur en sortie : " << error << endl;
        }
        *result = X;
        MPI_Barrier(MPI_COMM_WORLD);
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size" << endl;
        }
    }
}