#include <matrix/matrix.hpp>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <unistd.h>
#include <mpi.h>
#include <ctime>
#include <sys/time.h>

#include <chrono>

using namespace std;
    
double my_gettimeofday(){
    struct timeval tmp_time;
    gettimeofday(&tmp_time, NULL);
    return tmp_time.tv_sec + (tmp_time.tv_usec * 1.0e-6L);
}

int main(int argc, char** argv){
    //int size = 17281;
    string file = "data/matrix/rma10.mtx";
    int size = 46835;
    double start, differenceBest;
    int my_rank, p;
    double timer;
    double * result = NULL;
    MPI_Init(&argc, &argv);
    double * X = (double*)malloc(sizeof(double) * size);
    for(int i = 0; i < size; i ++){
        X[i] = pow(-1.0, i) * 1.0/size;
    }
    double * bestResult;
    double bestPrecision = 1e-8;
    int biggestBestResultInd = 0;
    float * resultFloat = NULL;
    float * Xf = (float*)malloc(sizeof(float) * size);
    for(int i = 0; i < size; i ++){
        Xf[i] = pow(-1.0, i) * 1.0/size;
    }

    MPI_Status status;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    if(my_rank == 0){
        cout << "Matrix : " << file << endl;
        cout << "Size : " << size << "x" << size << endl;
        cout << "Nb procs : " << p << endl;
        cout << endl << "----Matrice Double----" << endl;
        cout << "____ELL____" <<endl;
    }
    ELL myMatrix2;
    myMatrix2.initialize(file);
    MPI_Barrier(MPI_COMM_WORLD);

    pageRankBest(&myMatrix2, X, size, &bestResult, bestPrecision);
    double biggestValue = bestResult[0];
    for(int i = 1; i < size; i++){
        if(bestResult[i] > biggestValue){
            biggestValue = bestResult[i];
            biggestBestResultInd = i;
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0){
        start = my_gettimeofday();
    }
    timer = pageRank_time(&myMatrix2, X, size, &result);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS : " <<  my_gettimeofday() - start - timer << " s" << endl;
    }
    //for(int i = 1245; i < 1256; i++){
    //    cout << "result " << i << " : " << result[i] << endl; 
    //}
    differenceBest = 0;
    for(int i = 0; i < size; i++){
        differenceBest += (bestResult[i] - result[i])*(bestResult[i] - result[i]);
    }
    differenceBest = sqrt(differenceBest);
    if(my_rank == 0){
        cout << "Norm with best : " << differenceBest << endl;
        cout << "Difference with biggest value : " << abs(result[biggestBestResultInd] - bestResult[biggestBestResultInd]) << endl;
    }


    free(result);
    result = NULL;

    if(my_rank == 0){
        cout << "____CSR____" <<endl;
    }
    CSR myMatrix3;
    myMatrix3.initialize(file);
    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0){
        start = my_gettimeofday();
    }
    timer = pageRank_time(&myMatrix3, X, size, &result);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS : " <<  my_gettimeofday() - start - timer << " s" << endl;
    }
    differenceBest = 0;
    for(int i = 0; i < size; i++){
        differenceBest += (bestResult[i] - result[i])*(bestResult[i] - result[i]);
    }
    differenceBest = sqrt(differenceBest);
    if(my_rank == 0){
        cout << "Norm with best : " << differenceBest << endl;
        cout << "Difference with biggest value : " << abs(result[biggestBestResultInd] - bestResult[biggestBestResultInd]) << endl;
    }



    free(result);
    result = NULL;

    if(my_rank == 0){
        cout << "____COO____" <<endl;
    }
    COO myMatrix;
    myMatrix.initialize(file);

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    timer = pageRank_time(&myMatrix, X, size, &result);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS : " <<  my_gettimeofday() - start - timer << " s" << endl;
    }
    differenceBest = 0;
    for(int i = 0; i < size; i++){
        differenceBest += (bestResult[i] - result[i])*(bestResult[i] - result[i]);
    }
    differenceBest = sqrt(differenceBest);
    if(my_rank == 0){
        cout << "Norm with best : " << differenceBest << endl;
        cout << "Difference with biggest value : " << abs(result[biggestBestResultInd] - bestResult[biggestBestResultInd]) << endl;
    }

    //Matrix double vector float
    if(my_rank == 0){
        cout << "____ELL____" <<endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    timer = pageRank_time(&myMatrix2, Xf, size, &result);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS : " <<  my_gettimeofday() - start - timer << " s" << endl;
    }
    differenceBest = 0;
    for(int i = 0; i < size; i++){
        differenceBest += (bestResult[i] - result[i])*(bestResult[i] - result[i]);
    }
    differenceBest = sqrt(differenceBest);
    if(my_rank == 0){
        cout << "Norm with best : " << differenceBest << endl;
        cout << "Difference with biggest value : " << abs(result[biggestBestResultInd] - bestResult[biggestBestResultInd]) << endl;
    }


    free(result);
    result = NULL;

    if(my_rank == 0){
        cout << "____CSR____" <<endl;
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    timer = pageRank_time(&myMatrix3, Xf, size, &result);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS : " <<  my_gettimeofday() - start - timer << " s" << endl;
    }
    differenceBest = 0;
    for(int i = 0; i < size; i++){
        differenceBest += (bestResult[i] - result[i])*(bestResult[i] - result[i]);
    }
    differenceBest = sqrt(differenceBest);
    if(my_rank == 0){
        cout << "Norm with best : " << differenceBest << endl;
        cout << "Difference with biggest value : " << abs(result[biggestBestResultInd] - bestResult[biggestBestResultInd]) << endl;
    }


    free(result);
    result = NULL;

    if(my_rank == 0){
        cout << "____COO____" <<endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);

    if(my_rank == 0){
        start = my_gettimeofday();
    }
    timer = pageRank_time(&myMatrix, Xf, size, &result);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS : " <<  my_gettimeofday() - start - timer << " s" << endl;
    }
    differenceBest = 0;
    for(int i = 0; i < size; i++){
        differenceBest += (bestResult[i] - result[i])*(bestResult[i] - result[i]);
    }
    differenceBest = sqrt(differenceBest);
    if(my_rank == 0){
        cout << "Norm with best : " << differenceBest << endl;
        cout << "Difference with biggest value : " << abs(result[biggestBestResultInd] - bestResult[biggestBestResultInd]) << endl;
    }

    myMatrix2.freeData();
    myMatrix.freeData();
    myMatrix3.freeData();

    free(result);
    result = NULL;

    //Float
    if(my_rank == 0){
        cout << endl << "----Matrice Float----" << endl;
        cout << "____ELL____" <<endl;
    }
    ELLf myMatrix2f;
    myMatrix2f.initialize(file);

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    timer = pageRank_time(&myMatrix2f, Xf, size, &resultFloat);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS : " <<  my_gettimeofday() - start - timer << " s" << endl;
    }
    differenceBest = 0;
    for(int i = 0; i < size; i++){
        differenceBest += (bestResult[i] - resultFloat[i])*(bestResult[i] - resultFloat[i]);
    }
    differenceBest = sqrt(differenceBest);
    if(my_rank == 0){
        cout << "Norm with best : " << differenceBest << endl;
        cout << "Difference with biggest value : " << abs(resultFloat[biggestBestResultInd] - bestResult[biggestBestResultInd]) << endl;
    }


    free(resultFloat);
    resultFloat = NULL;

     if(my_rank == 0){
        cout << "____CSR____" <<endl;
    }
    CSRf myMatrix3f;
    myMatrix3f.initialize(file);

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    timer = pageRank_time(&myMatrix3f, Xf, size, &resultFloat);
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS : " <<  my_gettimeofday() - start - timer << " s" << endl;
    }
    differenceBest = 0;
    for(int i = 0; i < size; i++){
        differenceBest += (bestResult[i] - resultFloat[i])*(bestResult[i] - resultFloat[i]);
    }
    differenceBest = sqrt(differenceBest);
    if(my_rank == 0){
        cout << "Norm with best : " << differenceBest << endl;
        cout << "Difference with biggest value : " << abs(resultFloat[biggestBestResultInd] - bestResult[biggestBestResultInd]) << endl;
    }


    free(resultFloat);
    resultFloat = NULL;

    if(my_rank == 0){
        cout << "____COO____" <<endl;
    }
    COOf myMatrixf;
    myMatrixf.initialize(file);

    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        start = my_gettimeofday();
    }
    timer = pageRank_time(&myMatrixf, Xf, size, &resultFloat); 
    MPI_Barrier(MPI_COMM_WORLD);
    if(my_rank == 0){
        cout << "TEMPS : " <<  my_gettimeofday() - start -timer << " s" << endl;
    }
    differenceBest = 0;
    for(int i = 0; i < size; i++){
        differenceBest += (bestResult[i] - resultFloat[i])*(bestResult[i] - resultFloat[i]);
    }
    differenceBest = sqrt(differenceBest);
    if(my_rank == 0){
        cout << "Norm with best : " << differenceBest << endl;
        cout << "Difference with biggest value : " << abs(resultFloat[biggestBestResultInd] - bestResult[biggestBestResultInd]) << endl;
    }

    MPI_Finalize();
    
    return 0;
}
