#include "testFunc.h"

// #define SIMPLE_IO
#define MAIN_TEST
int main() {
#ifdef SIMPLE_IO
    string configFile = "/home/panwei/Desktop/bachelorThese/myCode/cplusCode/testConfig/config.txt";
    string dataFolderParentPath;
    vector<string> v_folderNames;
    string writeDataFolderParentPath;
    vector<double> thresholdVec;
    vector<int> modeVec;

    vector<string>  dataFolderParentPaths;
    vecVecString  vecVecFolderNames;
    vector<string>  writeDataFolderParentPaths;
    vecVecString writeDataFolderNames;
    vecVecDouble  vecVecThresholds;
    vecVecInt  vecVecModes;

    readConfig(configFile, dataFolderParentPaths, vecVecFolderNames,
        writeDataFolderParentPaths, 
        vecVecFolderNames, vecVecThresholds, vecVecModes);

    cout << "dataFolderP: " << dataFolderParentPaths << endl;
    cout << "v_folder: " << endl <<  vecVecFolderNames << endl;
    cout << "writeDataFolderParent: " << writeDataFolderParentPaths << endl;
    cout << "thresholdVec" << endl << vecVecThresholds << endl;
    cout << "modeVec" << endl << vecVecModes << endl;
#endif

    // readConfig(inFile, configFile, dataFolderParentPath,
    //     v_folderNames, writeDataFolderParentPath,
    //     thresholdVec, modeVec);
    // cout << "print out data" << endl;
    // cout << "dataFolderP: " << dataFolderParentPath << endl;
    // cout << "v_folder: " << endl <<  v_folderNames << endl;
    // cout << "writeDataFolderParent: " << writeDataFolderParentPath << endl;
    // cout << "thresholdVec" << endl << thresholdVec << endl;
    // cout << "modeVec" << endl << modeVec << endl;
#ifdef MAIN_TEST
    // std::string dataConfig = "/home/panwei/Desktop/bachelorThese/myCode/cplusCode/testConfig/simpleCity.txt";
    // myDouble endTheta = 0.1;
    // mainTest(dataConfig, "dataTest", endTheta);
    
    std::string dataConfig = "/home/panwei/Desktop/bachelorThese/myCode/cplusCode/testConfig/cityData.txt";
    myDouble endTheta = M_PI * 0.6;
    mainTest(dataConfig, "dataTest", endTheta);
#endif
    return 0;
}