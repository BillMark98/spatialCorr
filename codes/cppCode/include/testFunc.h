#ifndef TEST_FUNC_H
#define TEST_FUNC_H

#include "healpixCorrelation.h"
#include "poissonPoint2d.h"
#include "fileIO.h"
#include "filePath.h"
#include <sstream>

using std::string;
using std::vector;
typedef vector<vector<string> > vecVecString;
typedef vector<vector<double> > vecVecDouble;
typedef vector<vector<int> > vecVecInt;
ostream & outTwoVecCoord(ostream & os, const vec_coord & v1, const vec_coord & v2);
ostream & operator<<(ostream & os, const vec_coord & vec);
ostream & operator<<(ostream & os, const vec_3DCoord & vec);
ostream & operator<<(ostream & os, const std::vector<int> &vec);
template <typename T>
ostream & operator<<(ostream & os, const std::vector<T> &vec);

template <typename T>
ostream & operator<<(ostream & os, const std::vector<T> &vec) {
    for (typename std::vector<T>::const_iterator iter = vec.begin(); iter != vec.end(); iter++) {
        os << *iter << endl;
    }
    return os;
}


template <typename T>
ostream & operator<<(ostream & os, const std::vector<vector<T> > &vec) {
    for (typename std::vector<vector<T> >::const_iterator iter = vec.begin(); iter != vec.end(); iter++) {
        os << *iter << endl;
    }
    return os;
}

void readData(dataVecType & _thetaVec, dataVecType & _phiVec,
    dataVecType & _dataVec, const std::string & thetaFile, const std::string & phiFile,
    const std::string & dataFile);
/**
 * used to read in config files
 * @param filePath, specifies the config file path
 * @param dataFolderParentPath, vector of strings which points to the parent path of the data folder
 * @param v_folderNames vector of vector of strings, saving the data folderNames
 * @param writeDataFolderParentPath: vector of strings, the path to the output
 * @param thresholdVec, vectors saving threshold values vectors
 * @param modeVec, vecVecint saving the modes
 * */
void readConfig(const string & filePath, 
    vector<string> & dataFolderParentPath, vecVecString & v_folderNames, 
    vector<string> & writeDataFolderParentPath, 
    vecVecString & writeDataFolderNames,
    vecVecDouble & thresholdVec,
    vecVecInt & modeVec, const filePathType & loggerDir = "./logs/");


/**
 * used to read in config files
 * @param inFile, ifstream reference
 * @param filePath, specifies the config file path
 * @param dataFolderParentPath, the string which points to the parent path of the data folder
 * @param v_folderNames vector of strings, saving the data folderNames
 * @param writeDataFolderParentPath: string, the path to the output
 * @param thresholdVec, vectors saving threshold values
 * @param modeVec, vector saving the modes
 * */
ifstream& readConfig(ifstream & inFile,const string & filePath, 
    string & dataFolderParentPath, vector<string>& v_folderNames, 
    string & writeDataFolderParentPath, vector<string>& vecWriteFolderNames,
    vector<double>& thresholdVec,
    vector<int> & modeVec, const string & loggerName = "IOLogger");

// lambda daughter is not being used.
void maternTest(const string & dataFilePath, const vector<myDouble> & dTVec, const vector<myDouble> & ThetaBegin,
    const vector<myDouble> & ThetaEnd, 
    const vector<myDouble> & daughterThetaBegin,
    const vector<myDouble> & daughterThetaEnd,
    int runTime, const myDouble & lambdaParent, const myDouble & lambdaDaughter,const myDouble& smallAngle, ofstream & outFile,
        const string & outFilePre, const myDouble & endTheta,
        sizeType _Npoint, sizeType _M);


/**
 * 
 *  @param dataFilePath string that gives the path to the output file
 *  e.g "../data/" + datum + "poissonTest/"
 * */
void poissonTest(const string & dataFolderPath,
    const string & dataFileRelativePath,
    const vector<myDouble> & dTVec, 
    const vector<myDouble> & ThetaBegin,
    const vector<myDouble> & ThetaEnd, 
    int runTime, const myDouble & lambdaRaw,
    ofstream & outFile,
    const string & outFilePre, const myDouble & endTheta);

// /**
//  * @param configFile: string, the path to the configFile
//  * */
void dataTest(const string & configFile, const myDouble & endTheta,
    ofstream & outFile, const string & logFolder = "");

/**
 * 
 *  @param outDataFilePath string that gives the parent folder to the output data folder
 *  e.g "../data/" + datum + "poissonTest/"
 *  @param outDataFileSubPath: str, the folder name
 *  @param fileThetaPath: str,  the path to the csv which saves the theta info
 *  @param filePhiPath: str,  the path to the csv which saves the phi info
 *  @param dataPath: str,  the path to the csv which saves the RSS info
 *  @param dTVec: vector<myDouble>  vector that stores the deltaTheta of the pixels
 *  @param thresholdVec: vector<myDouble> vector that stores the thresholdValue
 *  @param modeVec: vector<int>  mode that specify the data conversion
 *  @param outFile: ostream & used to output result
 *  @param outFilePre: prefix for output data folder
 *  @param logFolder: string, folder for the logger
 * */    

void dataTest(const string & outDataFilePath,
    const string & outDataFileSubPath,  
    const string & fileThetaPath, 
    const string & filePhiPath,
    const string & dataPath, const vector<myDouble> & dTVec,
    const vector<myDouble> & thresholdVec, const vector<int> & modeVec, const myDouble & endTheta,
    ofstream & outFile, const string & outFilePre,
    const string & logFolder = "");
/**
 * main test, wrapping all test functions
 * @param configFile: string, string to the config file
 * @param testName: string, test name
 * */
void mainTest(const string & configFile, const string & testName, const myDouble & endTheta);

#endif