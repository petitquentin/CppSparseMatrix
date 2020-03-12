#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>

using namespace std;

class Matrix 
{
    public:
        void read_mtx_file(string filePath, vector<double>& val, vector<long int>& row, vector<long int>& col);
        long int getMNL(int num);
        virtual void initialize(string path) = 0;
        virtual void print() = 0;
    protected:
        long int MNL[3];
    private:
};

class COO : public Matrix
{
public:
    void initialize(string path);
    void print();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);
    void spmv(float * denseVector, int sizeDenseVector, double ** result);
    double norm();
    //vector<double> spmv_mpi(vector<double> denseVector);
    int getRowSize();
    long int * getRow();
    long int * getCol();
    double * getVal();
    void freeData();
private:
    double * val = NULL;
    long int * row = NULL;
    long int * col = NULL;
};

class COOf : public Matrix
{
public:
    void initialize(string path);
    void print();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);
    void spmv(float * denseVector, int sizeDenseVector, double ** result);
    void spmv(float * denseVector, int sizeDenseVector, float ** result);
    //vector<double> spmv_mpi(vector<double> denseVector);
    int getRowSize();
    long int * getRow();
    long int * getCol();
    float * getVal();
    float norm();
    void freeData();
private:
    float * val = NULL;
    long int * row = NULL;
    long int * col = NULL;
};

class ELL : public Matrix // Ellpack format
{
public: 
    void initialize(string path);
    void print();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);

    void spmv(float * denseVector, int sizeDenseVector, double ** result);
    double norm();
    //vector<double> spmv_mpi(vector<double> denseVector);
    long int sizeColInd(int num);
    double * getVal();
    long int * getColInd();
    void freeData();
private:
    long int * colInd = NULL;
    double * val = NULL;
    long int maxColInd;

};

class ELLf : public Matrix // Ellpack format
{
public: 
    void initialize(string path);
    void print();
    float norm();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);
    void spmv(float * denseVector, int sizeDenseVector, double ** result);
    void spmv(float * denseVector, int sizeDenseVector, float ** result);
    //vector<double> spmv_mpi(vector<double> denseVector);
    long int sizeColInd(int num);
    float * getVal();
    long int * getColInd();
    void freeData();
private:
    long int * colInd = NULL;
    float * val = NULL;
    long int maxColInd;

};

class CSR : public Matrix //CRS Format
{
public:
    void initialize(string path);
    void print();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);
    void spmv(float * denseVector, int sizeDenseVector, double ** result);
    //vector<double> spmv_mpi(vector<double> denseVector);
    int getValSize();
    long int * getPtr();
    long int * getInd();
    double * getVal();
    double norm();
    void freeData();
private:
    double * val = NULL;
    long int * ind = NULL;
    long int * ptr = NULL;
};

class CSRf : public Matrix //CRS Format
{
public:
    void initialize(string path);
    void print();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);
    void spmv(float * denseVector, int sizeDenseVector, double ** result);
    void spmv(float * denseVector, int sizeDenseVector, float ** result);
    //vector<double> spmv_mpi(vector<double> denseVector);
    int getValSize();
    long int * getPtr();
    long int * getInd();
    float * getVal();
    float norm();
    void freeData();
private:
    float * val = NULL;
    long int * ind = NULL;
    long int * ptr = NULL;
};

class SGP : public Matrix
{
public:
    void initialize(string path);
    void print();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);
    void spmv(float * denseVector, int sizeDenseVector, double ** result);
    void data();
    long int getMaxInd();
    long int * getAi();
    long int * getAj();
    long int * getJc();
    long int * getIc();
    double * getAcc();
    double * getAcr();
    double norm();
protected:
    double * acc = NULL;
    double * acr = NULL;
    long int * ai = NULL;
    long int * aj = NULL;
    long int * jc = NULL;
    long int * ic = NULL;
    int maxInd = -1;
};

//Matrix value : double

void spmv_mpi(COO *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(CSR *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(ELL *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(SGP *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(COO *matrix, float * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(CSR *matrix, float * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(ELL *matrix, float * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(SGP *matrix, float * denseVector, int sizeDenseVector, double ** result);

void spmvs_mpi(COO *matrix, double * denseVectorX, double * denseVectorB,  int sizeVector, double ** result);
void spmvs_mpi(CSR *matrix, double * denseVectorX, double * denseVectorB,  int sizeVector, double ** result);
void spmvs_mpi(ELL *matrix, double * denseVectorX, double * denseVectorB,  int sizeVector, double ** result);
void spmvs_mpi(SGP *matrix, double * denseVectorX, double * denseVectorB,  int sizeVector, double ** result);

//Matrix value : float

void spmv_mpi(COOf *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(CSRf *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(ELLf *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(COOf *matrix, float * denseVector, int sizeDenseVector, float ** result);
void spmv_mpi(CSRf *matrix, float * denseVector, int sizeDenseVector, float ** result);
void spmv_mpi(ELLf *matrix, float * denseVector, int sizeDenseVector, float ** result);


//Méthode de la puissance double
void * puissance(ELL *matrix, double * vecX, int sizeVecX, double ** result);
void * puissance(COO *matrix, double * vecX, int sizeVecX, double ** result);
void * puissance(CSR *matrix, double * vecX, int sizeVecX, double ** result);
void * puissance(ELL *matrix, float * vecX, int sizeVecX, double ** result);
void * puissance(COO *matrix, float * vecX, int sizeVecX, double ** result);
void * puissance(CSR *matrix, float * vecX, int sizeVecX, double ** result);

//Méthode de la puissance float
void * puissance(ELLf *matrix, double * vecX, int sizeVecX, double ** result);
void * puissance(COOf *matrix, double * vecX, int sizeVecX, double ** result);
void * puissance(CSRf *matrix, double * vecX, int sizeVecX, double ** result);
void * puissance(ELLf *matrix, float * vecX, int sizeVecX, float ** result);
void * puissance(COOf *matrix, float * vecX, int sizeVecX, float ** result);
void * puissance(CSRf *matrix, float * vecX, int sizeVecX, float ** result);

//Méthode de la puissance double time
double puissance_time(ELL *matrix, double * vecX, int sizeVecX, double ** result);
double puissance_time(COO *matrix, double * vecX, int sizeVecX, double ** result);
double puissance_time(CSR *matrix, double * vecX, int sizeVecX, double ** result);
double puissance_time(ELL *matrix, float * vecX, int sizeVecX, double ** result);
double puissance_time(COO *matrix, float * vecX, int sizeVecX, double ** result);
double puissance_time(CSR *matrix, float * vecX, int sizeVecX, double ** result);

//Méthode de la puissance float time
double puissance_time(ELLf *matrix, double * vecX, int sizeVecX, double ** result);
double puissance_time(COOf *matrix, double * vecX, int sizeVecX, double ** result);
double puissance_time(CSRf *matrix, double * vecX, int sizeVecX, double ** result);
double puissance_time(ELLf *matrix, float * vecX, int sizeVecX, float ** result);
double puissance_time(COOf *matrix, float * vecX, int sizeVecX, float ** result);
double puissance_time(CSRf *matrix, float * vecX, int sizeVecX, float ** result);

//PageRank Best Precision
void * pageRankBest(ELL *matrix, double * vecX, int sizeVecX, double ** result, double precision);

//PageRank Double
void * pageRank(ELL *matrix, double * vecX, int sizeVecX, double ** result);
void * pageRank(COO *matrix, double * vecX, int sizeVecX, double ** result);
void * pageRank(CSR *matrix, double * vecX, int sizeVecX, double ** result);
void * pageRank(ELL *matrix, float * vecX, int sizeVecX, double ** result);
void * pageRank(COO *matrix, float * vecX, int sizeVecX, double ** result);
void * pageRank(CSR *matrix, float * vecX, int sizeVecX, double ** result);

//pageRank float
void * pageRank(ELLf *matrix, double * vecX, int sizeVecX, double ** result);
void * pageRank(COOf *matrix, double * vecX, int sizeVecX, double ** result);
void * pageRank(CSRf *matrix, double * vecX, int sizeVecX, double ** result);
void * pageRank(ELLf *matrix, float * vecX, int sizeVecX, float ** result);
void * pageRank(COOf *matrix, float * vecX, int sizeVecX, float ** result);
void * pageRank(CSRf *matrix, float * vecX, int sizeVecX, float ** result);

//PageRank Double time
double pageRank_time(ELL *matrix, double * vecX, int sizeVecX, double ** result);
double pageRank_time(COO *matrix, double * vecX, int sizeVecX, double ** result);
double pageRank_time(CSR *matrix, double * vecX, int sizeVecX, double ** result);
double pageRank_time(ELL *matrix, float * vecX, int sizeVecX, double ** result);
double pageRank_time(COO *matrix, float * vecX, int sizeVecX, double ** result);
double pageRank_time(CSR *matrix, float * vecX, int sizeVecX, double ** result);

//pageRank float time
double pageRank_time(ELLf *matrix, double * vecX, int sizeVecX, double ** result);
double pageRank_time(COOf *matrix, double * vecX, int sizeVecX, double ** result);
double pageRank_time(CSRf *matrix, double * vecX, int sizeVecX, double ** result);
double pageRank_time(ELLf *matrix, float * vecX, int sizeVecX, float ** result);
double pageRank_time(COOf *matrix, float * vecX, int sizeVecX, float ** result);
double pageRank_time(CSRf *matrix, float * vecX, int sizeVecX, float ** result);

//SpMV time
double spmv_mpi_time(COO *matrix, double * denseVector, int sizeDenseVector, double ** result);
double spmv_mpi_time(CSR *matrix, double * denseVector, int sizeDenseVector, double ** result);
double spmv_mpi_time(ELL *matrix, double * denseVector, int sizeDenseVector, double ** result);
double spmv_mpi_time(COO *matrix, float * denseVector, int sizeDenseVector, double ** result);
double spmv_mpi_time(CSR *matrix, float * denseVector, int sizeDenseVector, double ** result);
double spmv_mpi_time(ELL *matrix, float * denseVector, int sizeDenseVector, double ** result);

double spmv_mpi_time(COOf *matrix, double * denseVector, int sizeDenseVector, double ** result);
double spmv_mpi_time(CSRf *matrix, double * denseVector, int sizeDenseVector, double ** result);
double spmv_mpi_time(ELLf *matrix, double * denseVector, int sizeDenseVector, double ** result);
double spmv_mpi_time(COOf *matrix, float * denseVector, int sizeDenseVector, float ** result);
double spmv_mpi_time(CSRf *matrix, float * denseVector, int sizeDenseVector, float ** result);
double spmv_mpi_time(ELLf *matrix, float * denseVector, int sizeDenseVector, float ** result);
#endif