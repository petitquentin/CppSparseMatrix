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
    int size = 10000;
    double start;
    int elapsed_seconds;

    cout << "____COO____" << endl; 
    COO myMatrix3;
    myMatrix3.initialize("data/smg2s10000.mtx");
    cout << "Norme = " << myMatrix3.norm()<< endl;

    cout << "____CSR____" << endl;
    CSR myMatrix;
    myMatrix.initialize("data/smg2s10000.mtx");
    cout << "Norme = " << myMatrix.norm()<< endl;

    cout << "____ELL____" <<endl;
    ELL myMatrix2;
    myMatrix2.initialize("data/smg2s10000.mtx");
    cout << "Norme = " << myMatrix2.norm()<< endl;

    cout << "____COOf____" << endl; 
    COO myMatrix3f;
    myMatrix3f.initialize("data/smg2s10000.mtx");
    cout << "Norme = " << myMatrix3f.norm()<< endl;

    cout << "____CSRf____" << endl;
    CSR myMatrixf;
    myMatrixf.initialize("data/smg2s10000.mtx");
    cout << "Norme = " << myMatrixf.norm()<< endl;

    cout << "____ELLf____" <<endl;
    ELL myMatrix2f;
    myMatrix2f.initialize("data/smg2s10000.mtx");
    cout << "Norme = " << myMatrix2f.norm()<< endl;

    cout << "____SGP____" <<endl;
    SGP myMatrix1;
    myMatrix1.initialize("data/smg2s10000.mtx");
    cout << "Norme = " << myMatrix1.norm()<< endl;
    
    return 0;
}
