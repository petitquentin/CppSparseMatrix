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

#endif