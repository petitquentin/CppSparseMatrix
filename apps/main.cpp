#include <matrix/matrix.hpp>
#include <iostream>
#include <filesystem>
#include <vector>
#include <stdio.h>

#include <cmath>
#include <fstream>
#include <algorithm>
#include <string>
#include <unistd.h>

using namespace std;
    
int main(){
    cout << "Coucou" << endl;
    CSR myMatrix;
    std::cout << "Current path is " << std::filesystem::current_path() << '\n';
    cout << "Coucou1" << endl;
    myMatrix.initialize("test.mtx");

    cout << "Coucou" << endl;
    myMatrix.print();
    cout << "Coucou1" << endl;

    ELL myMatrix2;
    cout << "Coucou" << endl;
    myMatrix2.initialize("test.mtx");
    cout << "Coucou2" << endl;
    myMatrix2.print();


    return 0;
}
