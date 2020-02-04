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
    //vector<double> spmv_mpi(vector<double> denseVector);
    int getRowSize();
    long int * getRow();
    long int * getCol();
    double * getVal();
private:
    double * val = NULL;
    long int * row = NULL;
    long int * col = NULL;
};

class ELL : public Matrix // Ellpack format
{
public: 
    void initialize(string path);
    void print();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);
    //vector<double> spmv_mpi(vector<double> denseVector);
    long int sizeColInd(int num);
    double * getVal();
    long int * getColInd();
private:
    long int * colInd = NULL;
    double * val = NULL;
    long int maxColInd;

};

class CSR : public Matrix //CRS Format
{
public:
    void initialize(string path);
    void print();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);
    //vector<double> spmv_mpi(vector<double> denseVector);
    int getValSize();
    long int * getPtr();
    long int * getInd();
    double * getVal();
private:
    double * val = NULL;
    long int * ind = NULL;
    long int * ptr = NULL;
};

class SGP : public Matrix
{
public:
    void initialize(string path);
    void print();
    void spmv(double * denseVector, int sizeDenseVector, double ** result);
    void data();
    long int getMaxInd();
    long int * getAi();
    long int * getAj();
    long int * getJc();
    long int * getIc();
    double * getAcc();
    double * getAcr();
protected:
    double * acc = NULL;
    double * acr = NULL;
    long int * ai = NULL;
    long int * aj = NULL;
    long int * jc = NULL;
    long int * ic = NULL;
    int maxInd = -1;
};

void spmv_mpi(COO *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(CSR *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(ELL *matrix, double * denseVector, int sizeDenseVector, double ** result);
void spmv_mpi(SGP *matrix, double * denseVector, int sizeDenseVector, double ** result);

#endif