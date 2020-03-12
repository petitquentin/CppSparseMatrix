#include <opencv2/opencv.hpp>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <cmath>
#include <fstream>

using namespace std;

void read_file(string filePath, vector<double>& val, vector<long int>& row, vector<long int>& col, long int * MNL){
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
            if(r <= MNL[0] && c <= MNL[1]){
                row.push_back(r-1);
                col.push_back(c-1);
                val.push_back(value);
            }
        }
        reader.close();
        if(val.size() != MNL[2]){
            MNL[2] = val.size();
        }
    }else{
        cout << "impossible to open the file" << endl;
    }
};


int main(int argc, char **argv)
{
    vector<double> val;
    vector<long int> row;
    vector<long int> col;
    long int * MNL = (long int *)malloc(sizeof(long int) * 3);

    if (argc < 2) {
        cout << "Pas assez d'argument" << endl;
        cout << "use -p to give the path" << endl;
        return 0;
    }

    char *path;

    std::string spectrum = " ";

    for (int i =0; i < argc; i++){

        if (strcasecmp(argv[i],"-p")==0){
                path = argv[i+1] ;
        }
    }
    string pathS(path);
    read_file(pathS, val, row, col, MNL);
    cout << MNL[1] << endl;
    // Rectangle vert 300 par 200
    cv::Mat img(MNL[1],MNL[1],CV_8UC3,cv::Scalar(255,255,255));
    cv::imshow("Essai", img);
    for(int i = 0; i < MNL[2]; i++){
        cv::Rect rect(row[i], col[i], 1, 1);
        cd::rectangle(img, rect, cv
    }
    // Attente appui sur une touche
    cv::waitKey(0);
    return 0;
}



//g++ exemple.cc -lopencv_core -lopencv_highgui
