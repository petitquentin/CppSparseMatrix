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
        void read_mtx_file(string filePath, vector<double>& val, vector<long int>& row, vector<long int>& col){
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
                cout << "impossible to open the file";
            }

            /* ifstream reader(filePath);
            if(reader){
                while(reader.peek() == '%'){
                    reader.ignore(2048, '\n');
                }
                reader >> MNL[0] >> MNL[1] >> MNL[2];
                int lastRow = 1;
                ptr.push_back(1);
                for(int i = 0; i < MNL[2]; i++){
                    long int row, col;
                    double value;
                    reader >> row >> col >> value;
                    ind.push_back(col);
                    val.push_back(value);
                    if(row > lastRow){
                        lastRow = row;
                        ptr.push_back(ind.size());
                    }
                }
                ptr.push_back(ind.size() + 1);
                reader.close();
            }else{
                cout << "impossible to open the file";
            } */

        }
        virtual void initialize(string path) = 0;
        long int MNL[3];
    private:
        

};

class COO : public Matrix
{
public:
    void initialize(string path)
    {
        val.clear();
        row.clear();
        col.clear();

        read_mtx_file(path, val, row, col);
    };

private:
    vector<double> val;
    vector<long int> row;
    vector<long int> col;
};

class ELL : public Matrix // Ellpack format
{
public: 
    void initialize(string path){
        vector<double> values;
        vector<long int> row;
        vector<long int> col;

        read_mtx_file(path, values, row, col);
        
        int maxInd = 0;
        int nbValuesInLine;
        for (int i = 1; i < MNL[0]; i++){
            nbValuesInLine = count(row.begin(), row.end(), i);
            if(nbValuesInLine > maxInd){
                maxInd = nbValuesInLine;
            }
        }

        //Initialization of colInd 
        typedef vector<long int> RowLongInt;
        typedef vector<double> RowDouble;
        colInd.clear();
        val.clear();
        for(int i = 0; i < MNL[0]; i++){
            RowLongInt rowLongInt(maxInd);
            RowDouble rowDouble(maxInd);
            for(int j = 0; j < maxInd; j++){
                rowLongInt[j] = -1;
                rowDouble[j] = nan("0");
            }
            colInd.push_back(rowLongInt);
            val.push_back(rowDouble)
            
        }



    }
private:
    vector<vector<long int>> colInd;
    vector<vector<double>> val;

};

class CSR : public Matrix
{
public:
    void initialize(string path){
        vector<double> values;
        vector<long int> row;
        vector<long int> col;

        read_mtx_file(path, values, row, col);
        

        //Initialize ptr
        ptr.clear();
        val.clear();
        ind.clear();

        ptr.push_back(0);
        for(long int i = 0; i < MNL[0]; i++){
            ptr.push_back(0);
        }
        for(long int n = 0; n < MNL[2]; n++){
            ptr[row[n]]++;
            val.push_back(0);
            ind.push_back(0);
        }
        for(long int i = 0, cumsum = 0; i< MNL[0]; i++){
            long int temp = ptr[i];
            ptr[i] = cumsum;
            cumsum += temp;
        }
        ptr[ptr.size()-1] = MNL[2];

        for(long int n = 0; n < MNL[2]; n++){
            long int r = row[n];
            long int dest = ptr[r];
            //cout << dest <<endl;
            ind[dest] = col[n];
            val[dest] = values[n];

            ptr[r]++;
        }

        for(long int i = 0, last = 0; i <= MNL[0]; i++){
            long int temp = ptr[i];
            ptr[i] = last;
            last = temp;
        }

    };

    void print(){
        long int c = 0;
        for(long int i = 1; i < ptr.size(); i++){
            long int start = ptr[i-1];
            long int end = ptr[i];
            vector<long int>::const_iterator first = ind.begin() + start;
            vector<long int>::const_iterator last = ind.begin() + end;
            vector<long int> row(first,last);
            for(int j = 0; j < ptr.size(); j++){
                if(count(row.begin(), row.end(), j) == 0)
                    cout << '0' << ' ';
                else{
                    cout << val[c] << '(' << c << ')' << ' ';
                    c++;
                }
            }
            std::cout << std::endl;
            

        }
        //print ind
        cout << "print ind" << endl;
        for(int i = 0; i < ind.size(); i++){
            cout << ind[i] << ' ';
        }
        cout << endl;
        //print val
        cout << "print val" << endl;
        for(int i = 0; i < val.size(); i++){
            cout << val[i] << ' ';
        }
        cout << endl;
        //print ptr
        cout << "print ptr" << endl;
        for(int i = 0; i < ptr.size(); i++){
            cout << ptr[i] << ' ';
        }
        cout << endl;
        cout << endl;
        cout << val.size() << ' ' << ind.size() << ' ' << ptr.size() << ' ' << ptr[ptr.size()-1] << ' ' << endl;
    }
    
private:
    vector<double> val;
    vector<long int> ind;
    vector<long int> ptr;
};

int main(){
    CSR myMatrix;
    myMatrix.initialize("test.mtx");
    myMatrix.print();

    return 0;
}