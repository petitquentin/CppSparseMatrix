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

#define EPSILON 1e-6

//Double vector

void * pageRank(ELLf *matrix, double * vecX, int sizeVecX, double ** result){
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
        double somme = 0;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        double error = INFINITY;
        int nbIt = 0;
        double beta = 0.8;
        double norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            spmv_mpi(matrix, X, sizeVecX, &Y);
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += X[i];
            }
            //norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]*beta + (1.0-beta)*norm/sizeVecX;
            }
            somme = 0;
            for(int i = 0; i < sizeVecX; i++){
                somme += Y[i];
            }
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/somme;
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

void * pageRank(CSRf *matrix, double * vecX, int sizeVecX, double ** result){
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
        double somme = 0;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        double error = INFINITY;
        int nbIt = 0;
        double beta = 0.8;
        double norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            spmv_mpi(matrix, X, sizeVecX, &Y);
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += X[i];
            }
            //norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]*beta + (1.0-beta)*norm/sizeVecX;
            }
            somme = 0;
            for(int i = 0; i < sizeVecX; i++){
                somme += Y[i];
            }
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/somme;
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

void * pageRank(COOf *matrix, double * vecX, int sizeVecX, double ** result){
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
        double somme = 0;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        double error = INFINITY;
        int nbIt = 0;
        double beta = 0.8;
        double norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            spmv_mpi(matrix, X, sizeVecX, &Y);
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += X[i];
            }
            //norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]*beta + (1.0-beta)*norm/sizeVecX;
            }
            somme = 0;
            for(int i = 0; i < sizeVecX; i++){
                somme += Y[i];
            }
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/somme;
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

//float vector

void * pageRank(ELLf *matrix, float * vecX, int sizeVecX, float ** result){
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
        double somme = 0;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        double error = INFINITY;
        int nbIt = 0;
        double beta = 0.8;
        double norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            spmv_mpi(matrix, X, sizeVecX, &Y);
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += X[i];
            }
            //norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]*beta + (1.0-beta)*norm/sizeVecX;
            }
            somme = 0;
            for(int i = 0; i < sizeVecX; i++){
                somme += Y[i];
            }
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/somme;
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
        (*result) = X;
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size" << endl;
        }
    }
}

void * pageRank(CSRf *matrix, float * vecX, int sizeVecX, float ** result){
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
        double somme = 0;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        double error = INFINITY;
        int nbIt = 0;
        double beta = 0.8;
        double norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            spmv_mpi(matrix, X, sizeVecX, &Y);
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += X[i];
            }
            //norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]*beta + (1.0-beta)*norm/sizeVecX;
            }
            somme = 0;
            for(int i = 0; i < sizeVecX; i++){
                somme += Y[i];
            }
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/somme;
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
        (*result) = X;
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size" << endl;
        }
    }
}

void * pageRank(COOf *matrix, float * vecX, int sizeVecX, float ** result){
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
        double somme = 0;
        for(int  i = 0; i < sizeVecX; i++){
            X[i] = vecX[i];
        }
        double error = INFINITY;
        int nbIt = 0;
        double beta = 0.8;
        double norm = INFINITY;
        while(error > EPSILON){
            //cout << "Itération " << nbIt << ", erreur actuel : " << error << " (" << norm << ")" << endl;
            spmv_mpi(matrix, X, sizeVecX, &Y);
            norm = 0;
            for(int i = 0; i < sizeVecX; i++){
                norm += X[i];
            }
            //norm = sqrt(norm);
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]*beta + (1.0-beta)*norm/sizeVecX;
            }
            somme = 0;
            for(int i = 0; i < sizeVecX; i++){
                somme += Y[i];
            }
            for(int i = 0; i < sizeVecX; i++){
                Y[i] = Y[i]/somme;
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
        (*result) = X;
    }else{
        if(my_rank == 0){
            cout << "the vector is not of the right size" << endl;
        }
    }
}

