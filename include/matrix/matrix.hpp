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
        virtual void initialize(string path) = 0;
    protected:
        long int MNL[3];
    private:
};

class COO : public Matrix
{
    public:
    void initialize(string path);

    private:
    vector<double> val;
    vector<long int> row;
    vector<long int> col;
};

class ELL : public Matrix // Ellpack format
{
public: 
    void initialize(string path);
    void print();
private:
    vector<vector<long int>> colInd;
    vector<vector<double>> val;

};

class CSR : public Matrix
{
    public:
    void initialize(string path);
    void print();
    
private:
    vector<double> val;
    vector<long int> ind;
    vector<long int> ptr;
};

#endif