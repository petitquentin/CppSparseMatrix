#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <matrix/matrix.hpp>

using namespace std;

void Matrix::read_mtx_file(string filePath, vector<double>& val, vector<long int>& row, vector<long int>& col){
    ifstream reader(filePath);
    if(reader){
        while(reader.peek() == '%'){
            reader.ignore(2048, '\n');
        }
        reader >> MNL[0] >> MNL[1] >> MNL[2];
        for(int i = 0; i < MNL[2]; i++){
            long int r, c;
            double value;
            reader >> r >> c >> value;
            row.push_back(r-1);
            col.push_back(c-1);
            val.push_back(value);
        }
        reader.close();
    }else{
        cout << "impossible to open the file" << endl;
    }
};

long int Matrix::getMNL(int num){
    if(num < 3 and num >=0){
        return MNL[num];
    }else{
        return -1;
    }
}
