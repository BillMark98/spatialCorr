#include <iostream>
#include <fstream>
#include <string>
#include <vector> 
using namespace std;

int main() {
    string fileThetaPath = "../data/matlabData/91Threshold/1RthetaFinal.csv";
    string filePhiPath = "../data/matlabData/91Threshold/1RphiFinal.csv";
    string dataPath = "../data/matlabData/91Threshold/1RRSSFinal.csv";
    ifstream inFile;
    inFile.open(fileThetaPath);
    if(!inFile.is_open()) {
        cout << "File not open\n";
        exit(9);
    }
    double val;
    vector<double> theta;
    while(inFile.good()) {
        inFile >> val;
        cout << val << endl;
        theta.push_back(val);
    }
    inFile.close();

    inFile.open(filePhiPath);
    if(!inFile.is_open()) {
        cout << "File not open\n";
        exit(9);
    }
    vector<double> phi;
    cout << "\n\n Phi\n";
    while(inFile.good()) {
        inFile >> val;
        cout << val << endl;
        phi.push_back(val);
    }

    inFile.close();


    inFile.open(dataPath);
    if(!inFile.is_open()) {
        cout << "File not open\n";
        exit(9);
    }
    vector<double> data;
    double min = 0.0;
    cout << "\n\n data\n";
    while(inFile.good()) {
        inFile >> val;
        if(val < min) {
            min = val;
        }
        cout << val << endl;
        data.push_back(val);
    }
    cout << "the minimum : " << min << endl;
    inFile.close();
    return 0;
}