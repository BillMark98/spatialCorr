#include "testFunc.h"

// help Variables
int dataTestCounts = 0;

// help function

ifstream & eatBlanks(ifstream & inFile);
/**
 *  help function to get a correct interpreted quoted string
 * */
void getQuotedString(ifstream & inFile, string & dataString, const string & loggerName = "IOLogger");
/**
 * help function get correct quoted string vecs
 * */
void getQuotedVecString(ifstream & inFile, vector<string>& vecString, const string & loggerName = "IOLogger");

/**
 * help function read in correct numerical vecs
 * */
void getDoubleVec(ifstream & inFile, vector<myDouble> & vecDoubles, const string & loggerName = "IOLogger");

/**
 * help function read in correct int
 * */
void getIntVec(ifstream & inFile, vector<int> & vecInts, const string & loggerName = "IOLogger");

ostream & outTwoVecCoord(ostream & os, const vec_coord & v1, const vec_coord & v2) {
    if (v1.size() != v2.size()) {
        cout << "vector not of same length\n";
        cout << "do not print out\n";
        return os;
    }
    for (typeIndex index = 0; index < v1.size(); index++) {
        os << v1[index].first << "\t" << v1[index].second << "\t" << v2[index].second;
        if(index != v1.size() - 1) {
            os << endl;
        }
    }
    return os;
}

ostream & operator<<(ostream & os, const vec_coord & vec) {
    for (unsigned long index = 0;index < vec.size();index++) {
        os << vec[index].first << "\t" << vec[index].second;
        if(index != vec.size() - 1) {
            os << endl;
        }
    }
    return os;
}
ostream & operator<<(ostream & os, const vec_3DCoord & vec) {
    for (unsigned long index = 0;index < vec.size();index++) {
        os << vec[index].getXCoord() << "\t" << setprecision(5) << vec[index].getYCoord() << "\t" << vec[index].getZCoord() << endl;
        // if(index != vec.size() - 1) {
        //     os << endl;
        // }
    }
    return os;
}


ostream & operator<<(ostream & os, const std::vector<int> &vec) {
    for (std::vector<int>::const_iterator iter = vec.begin(); iter != vec.end(); iter++) {
        os << *iter << endl;
    }
    return os;
}


void readData(dataVecType & _thetaVec, dataVecType & _phiVec,
    dataVecType & _dataVec, const std::string & thetaFile, const std::string & phiFile,
    const std::string & dataFile)
{
    ifstream inFile;
    inFile.open(thetaFile);
    if(!inFile.is_open()) {
        cout << "theta File not open\n";
        exit(9);
    }
    double val;
    while(inFile.good()) {
        inFile >> val;
        // cout << val << endl;
        val = 90.0 - val;   // corresponds to the theta in spherical coordinate
        _thetaVec.push_back(val);
    }
    inFile.close();

    inFile.open(phiFile);
    if(!inFile.is_open()) {
        cout << "phi File not open\n";
        exit(9);
    }
    while(inFile.good()) {
        inFile >> val;
        // cout << val << endl;
        _phiVec.push_back(val);
    }

    inFile.close();


    inFile.open(dataFile);
    if(!inFile.is_open()) {
        cout << "data File not open\n";
        exit(9);
    }

    string tempStr;
    while(inFile.good()) {
        inFile >> tempStr;
        if (tempStr == "nan") {
            val = MIN_RSS;
        }
        else {
            val = std::stod(tempStr);
        }
        _dataVec.push_back(val);
    }
    // if (inFile.good()) {

    // }
    inFile.close();
}

void readConfig(const string & filePath, 
    vector<string> & dataFolderParentPaths, 
    vecVecString & vecVecFolderNames, 
    vector<string> & writeDataFolderParentPaths, 
    vecVecString & writeDataFolderNames,
    vecVecDouble & vecVecThresholds,
    vecVecInt & vecVecModes, const filePathType & loggerDir) {

    ifstream inFile;
    inFile.open(filePath);
    if(!inFile.is_open()) {
        cout << "config File not open\n";
        exit(9);
    }
    auto ioLogger = spdlog::basic_logger_mt("IOLogger", loggerDir +  "/IOLogger.txt");
    ioLogger -> set_pattern("%v");
    ioLogger -> info("begin Reading config files");
    // set all vectors to 0 size
    dataFolderParentPaths.clear();
    vecVecFolderNames.clear();
    writeDataFolderParentPaths.clear();
    vecVecThresholds.clear();
    vecVecModes.clear();
    int count = 0;
    while (inFile.good()) {
        string dataFolderParentPath;
        vector<string> v_folderNames;
        string writeDataFolderParentPath;
        vector<string> v_writeDataFolderNames;
        vector<double> thresholdVec;
        vector<int> modeVec;
        readConfig(inFile, filePath, dataFolderParentPath,
        v_folderNames, writeDataFolderParentPath, v_writeDataFolderNames,
        thresholdVec, modeVec, "IOLogger");
        if (dataFolderParentPath.length() > 0) {
            dataFolderParentPaths.push_back(dataFolderParentPath);
            vecVecFolderNames.push_back(v_folderNames);
            writeDataFolderParentPaths.push_back(writeDataFolderParentPath);
            writeDataFolderNames.push_back(v_writeDataFolderNames);
            vecVecThresholds.push_back(thresholdVec);
            vecVecModes.push_back(modeVec);
        }
        else {
            ioLogger -> info("read in end");
            break;
        }
    }
    inFile.close();
}


ifstream& readConfig(ifstream & inFile, const string & filePath,
    string & dataFolderParentPath, vector<string>& v_folderNames, 
    string & writeDataFolderParentPath, 
    vector<string>& vecWriteFolderNames,
    vector<double> &thresholdVec,
    vector<int>& modeVec, const string & loggerName) {
    
    eatBlanks(inFile);
    auto ioLogger = spdlog::get(loggerName);

    if(!inFile.good() && inFile.eof()) {
        ioLogger -> info("EOF reached, file read successfully");
    }
    else {
        ioLogger -> warn("something wrong with config file reading");
    }

    // std::getline(inFile, dataFolderParentPath);
    getQuotedString(inFile, dataFolderParentPath);
    // cout << "dataFolderParentPath: " << dataFolderParentPath << endl;
    ioLogger -> info("dataFolderParentPath: {}", dataFolderParentPath);
    if(dataFolderParentPath.size() == 0) {
        return inFile;
    }
    // first clear the folderNames
    v_folderNames.clear();

    // read a line
    getQuotedVecString(inFile, v_folderNames);
    // string compactString;
    // std::getline(inFile, compactString);
    // // cout << "the compactString : " << compactString << endl;
    // // split the string
    // std::istringstream iss(compactString);
    // do {
    //     string folderName;
    //     iss >> std::quoted(folderName);
    //     // iss >> folderName;
    //     if (folderName.length() >= 1)
    //         v_folderNames.push_back(folderName);
    //     // cout << "folderName : " << folderName << endl;
    // }while(iss);
    ioLogger -> info("folderNames: {}", v_folderNames);
    if (inFile.good()) {
        // std::getline(inFile, writeDataFolderParentPath);
        getQuotedString(inFile, writeDataFolderParentPath);
        // cout << "writeDataFolderParentPath: " << writeDataFolderParentPath << endl;
        ioLogger -> info("writeDataFolderParentPath: {}", writeDataFolderParentPath);
    }
    else {
        // cout << "something is wrong with read\n";
        ioLogger -> error("Something wrong with reading");
        exit(9);
    }

    // clear the vecWriteFolderNames,
    vecWriteFolderNames.clear();

    getQuotedVecString(inFile, vecWriteFolderNames);
    // std::getline(inFile, compactString);
    // // cout << "the compactString : " << compactString << endl;
    // // split the string
    // iss = std::istringstream(compactString);
    // do {
    //     string folderName;
    //     // iss >> std::quoted(folderName);
    //     iss >> folderName;
    //     if (folderName.length() >= 1)
    //         vecWriteFolderNames.push_back(folderName);
    //     // cout << "folderName : " << folderName << endl;
    // }while(iss);

    // clear thresholdVal
    thresholdVec.clear();
    getDoubleVec(inFile, thresholdVec);
    // read in a line
    // std::getline(inFile, compactString);
    // iss = std::istringstream(compactString);
    // do {
    //     string tempStr;
    //     iss >> tempStr;
    //     double val;
    //     if (tempStr == "nan") {
    //         val = MIN_RSS;
    //     }
    //     else {
    //         try {
    //             val = std::stod(tempStr);
    //         }
    //         catch(std::invalid_argument & ) {
    //             if (tempStr != "")
    //                 // cout << "invalid argument, tempStr: " << tempStr << endl;
    //                 ioLogger -> warn("Invalid Argument: {}", tempStr);
    //             continue;
    //         }
    //     }
    //     thresholdVec.push_back(val);
    // }while(iss);
    // cout << "the thresholdVec:\n";
    // cout << thresholdVec << endl;
    ioLogger -> info("thresholdVec: {}", thresholdVec);
    //
    modeVec.clear();
    getIntVec(inFile, modeVec);
    // cout << "the modeVec\n";
    // cout << modeVec << endl;

    ioLogger -> info("modeVec: {}", modeVec);
    ioLogger -> flush();
    return inFile;
}

void maternTest(const string & dataFilePath, const vector<myDouble> & dTVec, const vector<myDouble> & ThetaBegin,
    const vector<myDouble> & ThetaEnd, 
    const vector<myDouble> & daughterThetaBegin,
    const vector<myDouble> & daughterThetaEnd,
    int runTime, const myDouble & lambdaParent, const myDouble & lambdaDaughter,const myDouble& smallAngle, ofstream & outFile,
        const string & outFilePre, const myDouble & endTheta,
        sizeType _Npoint, sizeType _M)
{

    PoissonPoint2d poissonP1;
    ofstream localOutFile;
    for (typeIndex indexBeg = 0; indexBeg < ThetaBegin.size(); indexBeg++) {
        for (typeIndex indexEnd = 0; indexEnd < ThetaEnd.size(); indexEnd++) {
            for (typeIndex d_indexBeg = 0; d_indexBeg < daughterThetaBegin.size(); d_indexBeg++) {
                for (typeIndex d_indexEnd = 0; d_indexEnd < daughterThetaEnd.size(); d_indexEnd++) {
                   for (typeIndex indexdT = 0; indexdT < dTVec.size(); indexdT++) {
                        for (int count = 0; count < runTime; count++) {

                        myDouble thetaBegin = ThetaBegin[indexBeg];
                        myDouble thetaEnd = ThetaEnd[indexEnd];
                        myDouble d_thetaBegin = daughterThetaBegin[d_indexBeg];
                        myDouble d_thetaEnd = daughterThetaEnd[d_indexEnd];
                        vec_3DCoord vecPosition;
                        poissonP1.rPoissonSphereDegree(vecPosition,lambdaParent,thetaEnd - daughterThetaEnd[d_indexEnd], thetaBegin + daughterThetaEnd[d_indexEnd]);
                        countType rawPoint = vecPosition.size();

                        vec_3DCoord vecData;
                    #ifdef MATERN_PARALLEL
                        cout << "the raw Point: " << rawPoint << endl;
                        poissonP1.MaternProcessChooseNPointParallel(vecPosition,vecData, _Npoint,_M, d_thetaEnd, d_thetaBegin);
                    #else
                        poissonP1.MaternProcessNPoint(vec3d2,_Npoint, smallAngle,thetaBegin);
                    #endif

                        unsigned long totalPoint = vecData.size();


    std::stringstream stream;
    stream << "points";
    stream << vecData.size();
    stream << "lambda0_";
    stream << lambdaParent;
    stream << "avgTime_" << AVERAGE_TIME;
    stream << "dT";
    stream << std::fixed << std::setprecision(5) << dTVec[indexdT];
    // std::string str_dT = stream.str();
    stream << "_endT_" << std::fixed << std::setprecision(5) << endTheta;
    stream << "pos_begin_" << std::fixed << std::setprecision(5) << thetaBegin;
    stream << "pos_end_" << std::fixed << std::setprecision(5)<< thetaEnd;
    stream << "d_begin_" << std::fixed << std::setprecision(5) << d_thetaBegin;
    stream << "d_end_" << std::fixed << std::setprecision(5)<< d_thetaEnd;
                        
                        std::string fileName = stream.str();
                        stream << "/runNum_" << count;
                        std::string str_endT = stream.str();

                        string strFileName = dataFilePath + "stringName.txt";
                        outFile.open(strFileName,std::ofstream::out | std::ofstream::app);
                        outFile << str_endT << endl;
                        outFile.close();


                    string path = dataFilePath + fileName;
                    if(!IsPathExist(path)) {
                        if (mkdir(path.c_str(), 0777) != 0 ) {
                            cout << "error creating direc: " << path << endl;
                            continue;
                        }
                    }
                    string fileNamePrefix = dataFilePath + str_endT;

                        string spherePointName = fileNamePrefix + "spherePoint.txt";
                        string healpixMapFileName = fileNamePrefix + "poissonHealpixRandomMap.txt";
                        string pixelHealpixMap = fileNamePrefix + "poissonHealpixRandomPixelMap.txt";
                        // string randomHealpixMapFileName = fileNamePrefix + "Pi_randomHealpixRandom.txt";
                        string peebles1 = fileNamePrefix + "poissonHealpixRandomMode1.txt";
                        string peebles2 = fileNamePrefix + "poissonHealpixRandomMode2.txt";
                        string hamilton = fileNamePrefix + "poissonHealpixRandomMode3.txt";
                        string landy = fileNamePrefix + "poissonHealpixRandomMode4.txt";
                        string randomMap = fileNamePrefix + "randomMap.txt";
                        string calculationInfo = fileNamePrefix + "calculationInfo.txt";
                        outFile.open(spherePointName);
                        outFile << vecData;
                        outFile.close();

                        std::stringstream fileNameTemp;
                        fileNameTemp << "totalPoint_" << std::fixed << std::setprecision(5) << totalPoint;
                        fileNameTemp << "count_" << count << ".txt";
                        string fileOutName = outFilePre + fileNameTemp.str();
                        localOutFile.open(fileOutName);


                        // sphereTwoCorrelation spT(vecData, poissonP1);

                        std::vector<myDouble> healpixMapInfo;
                        std::vector<myDouble> randomMapInfo;

                        vec_coord peebles1VC,peebles2VC,hamiltonVC,landyVC;
                        countType scale = 1;    
                        Healpix_Correlation mapTest(dTVec[indexdT],thetaEnd,thetaBegin);
                        // mapTest.SetPixValue(vecData);
                #ifndef PRINT_POISSON
                        mapTest.pointProcessSweep2PCF(dTVec[indexdT],endTheta,vecData,
                            peebles1VC, peebles2VC, hamiltonVC,landyVC,
                            AVERAGE_TIME, scale,healpixMapInfo,outFile);
                        myDouble dT = dTVec[indexdT];
                        myDouble spannedAngle = mapTest.GetSpannedAngle();
                        int nTime = AVERAGE_TIME;
                        countType map1Points = vecData.size();
                        countType pointCount = mapTest.GetPointsCount();
                        cout << "pointCount: " << map1Points << "\t map1Points: " << map1Points << endl;
                        countType randomCount = map1Points * scale;
                        // Healpix_Correlation Hmap(d,Angle);
                        // Hmap.SetPixValue(vec3d2);
                        // spT.averageHealpixTwoPointCorrNoSet(deltaTheta,endTheta, vC, Hmap,healpixMapInfo,pointCount, spannedAngle, dT, nTime,HEALPIX_RANDOM_POINT_COUNT,1);
                        
                        // spT.sweepAverageHealpix2PCF(dTVec[indexdT],endTheta,vecData,
                        //     thetaBegin, thetaEnd,
                        //     peebles1VC,peebles2VC,
                        //     hamiltonVC,landyVC,healpixMapInfo,pointCount,
                        //     spannedAngle,dT,localOutFile,nTime,randomCount);
                        // mapTest.sweepAverageTwoPCorrelation(dTVec[indexdT],endTheta,
                        //     peebles1VC,peebles2VC,hamiltonVC,landyVC,AVERAGE_TIME,randomCount,
                        //     healpixMapInfo,localOutFile);  
                                                  
                        outFile.open(calculationInfo);
                        outFile << rawPoint << endl;
                        outFile << randomCount << endl;

                        outFile << dT << endl;
                        outFile << nTime << endl;
                        outFile << spannedAngle << endl;
                        outFile << endTheta << endl;
                        outFile << _M << endl;
                        outFile << smallAngle << endl;
                        outFile << _Npoint << endl;
                        outFile << totalPoint << endl;
                        outFile << thetaBegin << endl;
                        outFile << thetaEnd << endl;
                        outFile << d_thetaBegin << endl;
                        outFile << d_thetaEnd << endl;
                    #ifdef signalMode
                        outFile << 1 << endl;
                    #else   
                        outFile << 0 << endl;
                    #endif

                    #ifdef OFFSET_POISSON
                        outFile << 1 << endl;
                    #else
                        outFile << 0 << endl;
                    #endif


                #ifdef MATERN_INFO_OUT
                    cout << "rawPoint: " << rawPoint << endl;
                    cout << "default randomPoint count: " << HEALPIX_RANDOM_POINT_COUNT << endl;
                    cout << "deltaTheta: " << dT << endl;
                    cout << "default avg time: " << AVERAGE_TIME << endl;
                    cout << "spannedAngle: " << spannedAngle << endl;
                    cout << "endTheta: " << endTheta << endl;
                    cout << "_M: " << _M << endl;
                    cout << "smallAngle" << smallAngle << endl;
                    cout << "_Npoint: " << _Npoint << endl;
                    cout << "totalPoint: " << totalPoint << endl;
                #endif

                        outFile.close();

                        outFile.open(peebles1);
                        outFile << peebles1VC;
                        outFile.close();

                        outFile.open(healpixMapFileName);
                        outFile << healpixMapInfo;
                        outFile.close();

                        // outFile.open("../data/randomHealpix.txt");
                        // outFile << randomMapInfo;
                        // outFile.close();

                        // spT.averageHealpixTwoPointCorrNoSet(deltaTheta,endTheta, vC,Hmap,nTime,HEALPIX_RANDOM_POINT_COUNT,2);
                        outFile.open(peebles2);
                        outFile << peebles2VC;
                        outFile.close();

                        // spT.averageHealpixTwoPointCorrNoSet(deltaTheta,endTheta, vC,Hmap,nTime,HEALPIX_RANDOM_POINT_COUNT,3);
                        outFile.open(hamilton);
                        outFile << hamiltonVC;
                        outFile.close();

                        // spT.averageHealpixTwoPointCorrNoSet(deltaTheta,endTheta, vC,Hmap,nTime,HEALPIX_RANDOM_POINT_COUNT,4);
                        outFile.open(landy);
                        outFile << landyVC;
                        outFile.close();

                        // for the angular power spectrum analysis
                        mapTest.outMapAndRandom(healpixMapInfo,randomMapInfo);
                        outFile.open(randomMap);
                        outFile << randomMapInfo;
                        outFile.close();

                        localOutFile.close();

                #else
                        vec_coord poissonPeebles1VC, poissonPeebles2VC,poissonHamiltonVC,poissonLandyVC;

                       mapTest.pointProcessPrintPoissonSweep2PCF(dTVec[indexdT],endTheta,vecData,
                            peebles1VC, peebles2VC, hamiltonVC,landyVC,
                            poissonPeebles1VC,poissonPeebles2VC,poissonHamiltonVC, poissonLandyVC,
                            AVERAGE_TIME, scale,healpixMapInfo,outFile);
                        myDouble dT = dTVec[indexdT];
                        myDouble spannedAngle = mapTest.GetSpannedAngle();
                        int nTime = AVERAGE_TIME;
                        countType map1Points = vecData.size();
                        countType pointCount = mapTest.GetPointsCount();
                        cout << "pointCount: " << map1Points << "\t map1Points: " << map1Points << endl;
                        countType randomCount = map1Points * scale;
                        // Healpix_Correlation Hmap(d,Angle);
                        // Hmap.SetPixValue(vec3d2);
                        // spT.averageHealpixTwoPointCorrNoSet(deltaTheta,endTheta, vC, Hmap,healpixMapInfo,pointCount, spannedAngle, dT, nTime,HEALPIX_RANDOM_POINT_COUNT,1);
                        
                        // spT.sweepAverageHealpix2PCF(dTVec[indexdT],endTheta,vecData,
                        //     thetaBegin, thetaEnd,
                        //     peebles1VC,peebles2VC,
                        //     hamiltonVC,landyVC,healpixMapInfo,pointCount,
                        //     spannedAngle,dT,localOutFile,nTime,randomCount);
                        // mapTest.sweepAverageTwoPCorrelation(dTVec[indexdT],endTheta,
                        //     peebles1VC,peebles2VC,hamiltonVC,landyVC,AVERAGE_TIME,randomCount,
                        //     healpixMapInfo,localOutFile);  
                                                  
                        outFile.open(calculationInfo);
                        outFile << rawPoint << endl;
                        outFile << randomCount << endl;

                        outFile << dT << endl;
                        outFile << nTime << endl;
                        outFile << spannedAngle << endl;
                        outFile << endTheta << endl;
                        outFile << _M << endl;
                        outFile << smallAngle << endl;
                        outFile << _Npoint << endl;
                        outFile << totalPoint << endl;
                        outFile << thetaBegin << endl;
                        outFile << thetaEnd << endl;
                        outFile << d_thetaBegin << endl;
                        outFile << d_thetaEnd << endl;
                    #ifdef signalMode
                        outFile << 1 << endl;
                    #else   
                        outFile << 0 << endl;
                    #endif

                    #ifdef OFFSET_POISSON
                        outFile << 1 << endl;
                    #else
                        outFile << 0 << endl;
                    #endif


                #ifdef MATERN_INFO_OUT
                    cout << "rawPoint: " << rawPoint << endl;
                    cout << "default randomPoint count: " << HEALPIX_RANDOM_POINT_COUNT << endl;
                    cout << "deltaTheta: " << dT << endl;
                    cout << "default avg time: " << AVERAGE_TIME << endl;
                    cout << "spannedAngle: " << spannedAngle << endl;
                    cout << "endTheta: " << endTheta << endl;
                    cout << "_M: " << _M << endl;
                    cout << "smallAngle" << smallAngle << endl;
                    cout << "_Npoint: " << _Npoint << endl;
                    cout << "totalPoint: " << totalPoint << endl;
                #endif

                        outFile.close();

                        outFile.open(peebles1);
                        outTwoVecCoord(outFile,peebles1VC,poissonPeebles1VC);
                        outFile.close();

                        outFile.open(healpixMapFileName);
                        outFile << healpixMapInfo;
                        outFile.close();

                        // outFile.open("../data/randomHealpix.txt");
                        // outFile << randomMapInfo;
                        // outFile.close();

                        // spT.averageHealpixTwoPointCorrNoSet(deltaTheta,endTheta, vC,Hmap,nTime,HEALPIX_RANDOM_POINT_COUNT,2);
                        outFile.open(peebles2);
                        outTwoVecCoord(outFile, peebles2VC,poissonPeebles2VC);
                        outFile.close();

                        // spT.averageHealpixTwoPointCorrNoSet(deltaTheta,endTheta, vC,Hmap,nTime,HEALPIX_RANDOM_POINT_COUNT,3);
                        outFile.open(hamilton);
                        outTwoVecCoord(outFile, hamiltonVC,poissonHamiltonVC);
                        outFile.close();

                        // spT.averageHealpixTwoPointCorrNoSet(deltaTheta,endTheta, vC,Hmap,nTime,HEALPIX_RANDOM_POINT_COUNT,4);
                        outFile.open(landy);
                        outTwoVecCoord(outFile,landyVC,poissonLandyVC);
                        outFile.close();


                #endif
                        }
                   }
                }
            }
        }
    }    
}    


void poissonTest(const string & dataFolderPath,
    const string & dataFileRelativePath,
    const vector<myDouble> & dTVec, 
    const vector<myDouble> & ThetaBegin,
    const vector<myDouble> & ThetaEnd, 
    int runTime, const myDouble & lambdaRaw,
    ofstream & outFile,
    const string & outFilePre, const myDouble & endTheta)
{
    PoissonPoint2d poissonP1;
    // std::string poissonPrefix = "poissonTest/";
    string dataFilePath = dataFolderPath + "/" + dataFileRelativePath;
    string fileNamePrefix = dataFilePath;
    string strFileName = fileNamePrefix + "stringName.txt";
    ofstream localOutFile;
    for (typeIndex indexBeg = 0; indexBeg < ThetaBegin.size(); indexBeg++) {
        for (typeIndex indexEnd = 0; indexEnd < ThetaEnd.size(); indexEnd++) {
            for (typeIndex indexdT = 0; indexdT < dTVec.size(); indexdT++) {
                for (int count = 0; count < runTime; count++) {
                    myDouble thetaBegin = ThetaBegin[indexBeg];
                    myDouble thetaEnd = ThetaEnd[indexEnd];
                    vec_3DCoord vecPoint;
                    poissonP1.rPoissonSphereDegree(vecPoint,lambdaRaw,thetaEnd, thetaBegin);

                    countType map1Points = vecPoint.size();
                    countType scale = 1;
                    countType map2Points = map1Points * scale;
                    countType randomCount = map2Points;
                    int nTime = AVERAGE_TIME;
                    vec_coord peebles1VC;
                    vec_coord peebles2VC;
                    vec_coord hamiltonVC;
                    vec_coord landyVC;

                            // int vecLen = ceil((double)endTheta/dT) + 1;
                            // std::vector<myDouble> bins;
                            // bins.resize(vecLen);
                            // for(int index = 0; index < vecLen; index+= 1) {
                            //     // the end may not be accurate
                            //     bins[index] = index * dT;
                            // }
                            
                            // map3.twoPCorrelation(deltaTheta,endTheta,vC,map4,healpixMapInfo,randomMapInfo,1);
                            // cout << "healpixMapInfo length: " << healpixMapInfo.size() << endl;

                            // auto t2 = std::chrono::high_resolution_clock::now();
                            // std::chrono::duration<double, std::milli> sequentialTime = t2 - t1;
                            // cout << "the time for the first part is: " << sequentialTime.count() << endl;
                            // cout << "using the pi\n";

                            // t1 = std::chrono::high_resolution_clock::now();

                                // countType map2Points = HEALPIX_RANDOM_POINT_COUNT;
        std::stringstream stream;
        stream << "pt";
        stream << vecPoint.size();
        stream << "lda_";
        stream << lambdaRaw;
        stream << "agT_" << AVERAGE_TIME;
        stream << "nr_" << randomCount;
        stream << "dT";
        stream << std::fixed << std::setprecision(5) << dTVec[indexdT];
        // std::string str_dT = stream.str();
        stream << "eT_" << std::fixed << std::setprecision(4) << endTheta;
        stream << "tB_" << std::fixed << std::setprecision(5) << thetaBegin;
        stream << "tE_" << std::fixed << std::setprecision(5) << thetaEnd;

                    std::string fileName = stream.str();
                    stream << "/rN_" << count;
                    std::string str_endT = stream.str();
                    outFile.open(strFileName,std::ofstream::out | std::ofstream::app);
                    outFile << str_endT << endl;
                    outFile.close();

                    // nameHealStr.push_back(str_endT);
                    string path = dataFilePath + fileName;
                    // if(!IsPathExist(path)) {
                    //     if (mkdir(path.c_str(), 0777) != 0 ) {
                    //         cout << "error creating direc: " << path << endl;
                    //         continue;
                    //     }
                    // }
                    if (createPath(0777,dataFolderPath,dataFileRelativePath + fileName) == -1) {
                        cout << "error creating dir_" << path << endl;
                        continue;
                    }    
                    string fileNamePrefix = dataFilePath + str_endT;

            // fileNamePrefix += str_endT;
            string spherePointName = fileNamePrefix +  "spherePoint.txt";
            
            string healpixMapFileName = fileNamePrefix + "poissonHealpixRandomMap.txt";
            string pixelHealpixMap = fileNamePrefix + "poissonHealpixRandomPixelMap.txt";
            // string randomHealpixMapFileName = fileNamePrefix + "Pi_randomHealpixRandom.txt";
            string peebles1 = fileNamePrefix + "poissonHealpixRandomMode1.txt";
            string peebles2 = fileNamePrefix + "poissonHealpixRandomMode2.txt";
            string hamilton = fileNamePrefix + "poissonHealpixRandomMode3.txt";
            string landy = fileNamePrefix + "poissonHealpixRandomMode4.txt";
            string randomMap = fileNamePrefix + "randomMap.txt";
            string calculationInfo = fileNamePrefix + "calculationInfo.txt";

            // cout << spherePointName << endl;
            // if (outFile.is_open()) {
            //     cout << "file opend\n";
            // }
            // else {
            //     break;
            // }
            outFile.open(spherePointName);
            outFile << vecPoint;
            outFile.close();

                    std::stringstream fileNameTemp;
                    fileNameTemp << "map1P_" << std::fixed << std::setprecision(5) << map1Points;
                    fileNameTemp << "count_" << count << ".txt";
                    string fileOutName = outFilePre + fileNameTemp.str();
                    localOutFile.open(fileOutName);

    
                    outFile.open(calculationInfo);
                    outFile << vecPoint.size() << endl;
                    outFile << randomCount << endl;
                    outFile << dTVec[indexdT] << endl;
                    outFile << AVERAGE_TIME << endl;
                    outFile << thetaEnd - thetaBegin << endl;
                    outFile << endTheta << endl;   
                    outFile << thetaBegin << endl;
                    outFile << thetaEnd << endl;             
                    outFile.close();

                    Healpix_Correlation mapTest(dTVec[indexdT],thetaEnd,thetaBegin);
                    mapTest.SetPixValue(vecPoint);
                    // sphereTwoCorrelation spT(vecPoint, poissonP1);

                    // std::vector<unsigned long> healpixMapInfo;
                    // std::vector<unsigned long> randomMapInfo;
                    std::vector<myDouble> healpixMapInfo;
                    std::vector<myDouble> randomMapInfo;

                    
                    // myDouble dT;
                    // myDouble spannedAngle;
                    // countType pointCount;

                    // Healpix_Correlation Hmap(deltaTheta,Angle);
                    // Hmap.SetPixValue(vec3d2);
                    // spT.averageHealpixTwoPointCorrNoSet(deltaTheta,endTheta, vC, Hmap,healpixMapInfo,pointCount, spannedAngle, dT, nTime,HEALPIX_RANDOM_POINT_COUNT,1);
                    // spT.sweepAverageHealpix2PCF(dTVec[indexdT],endTheta,vecPoint,thetaBegin,
                    //     thetaEnd, peebles1VC,peebles2VC,
                    //     hamiltonVC,landyVC,healpixMapInfo,pointCount,
                    //     spannedAngle,dT,localOutFile,AVERAGE_TIME,randomCount);
                    mapTest.sweepAverageTwoPCorrelation(dTVec[indexdT],endTheta,
                        peebles1VC,peebles2VC,hamiltonVC,landyVC,AVERAGE_TIME,map2Points,
                        healpixMapInfo,localOutFile);
                    outFile.open(peebles1);
                    outFile << peebles1VC;
                    outFile.close();

                    outFile.open(healpixMapFileName);
                    outFile << healpixMapInfo;
                    outFile.close();

                    outFile.open(peebles2);
                    outFile << peebles2VC;
                    outFile.close();


                    outFile.open(hamilton);
                    outFile << hamiltonVC;
                    outFile.close();

                    outFile.open(landy);
                    outFile << landyVC;
                    outFile.close();
                    
                    // for the angular power spectrum analysis
                    mapTest.outMapAndRandom(healpixMapInfo,randomMapInfo);
                    outFile.open(randomMap);
                    outFile << randomMapInfo;
                    outFile.close();

                    localOutFile.close();
                }
            }
        }
    }

}

void dataTest(const string & configFile, const myDouble & endTheta,
    ofstream & outFile, const string & logFolder) {
    vector<string>  dataFolderParentPaths;
    vecVecString  vecVecFolderNames;
    vector<string>  writeDataFolderParentPaths;
    vecVecString writeFolderNames;
    vecVecDouble  vecVecThresholds;
    vecVecInt  vecVecModes;
    readConfig(configFile, dataFolderParentPaths, vecVecFolderNames,
        writeDataFolderParentPaths, writeFolderNames,
        vecVecThresholds, vecVecModes);
    
    string thetaSuffix = "thetaFinal.csv";
    string phiSuffix = "phiFinal.csv";
    string rssSuffix = "RSSFinal.csv";   
    vector<myDouble> dTVec = {0.01};
    // myDouble endTheta = M_PI * 0.6;

    for (typeIndex count = 0; count < dataFolderParentPaths.size(); count++) {
        for (typeIndex nameIndex = 0; nameIndex < vecVecFolderNames[count].size(); nameIndex++) {
            // build the paths to files
            string folderName = vecVecFolderNames[count][nameIndex];
            filePathType dataFilePre = dataFolderParentPaths[count];
            string fileThetaPath = dataFilePre + folderName + "/" + thetaSuffix;
            string filePhiPath = dataFilePre + folderName + "/" + phiSuffix;
            string dataPath = dataFilePre + folderName + "/" + rssSuffix;
            vector<double> thresholdVec = vecVecThresholds[count];
            vector<int> modeVec = vecVecModes[count];
            dataTest(writeDataFolderParentPaths[count],
            writeFolderNames[count][nameIndex],
            fileThetaPath, filePhiPath,dataPath,dTVec,thresholdVec,
            modeVec, 
            endTheta, outFile, folderName,
            logFolder);

        }
    }
}

void dataTest(const string & outDataFileRootPath, 
    const string & outDataFileSubPath, 
    const string & fileThetaPath, const string & filePhiPath,
    const string & dataPath, const vector<myDouble> & dTVec,
    const vector<myDouble> & thresholdVec, 
    const vector<int> & modeVec,
    const myDouble & endTheta,
    ofstream & outFile, const string & outFilePre,
    const string & logFolder)
{
    dataTestCounts++;
    std::string logFolderName;
    // build the logger Name
    if (logFolder == "") {
        // use the default
        logFolderName = "test" + std::to_string(dataTestCounts);
    }
    else {
        logFolderName = logFolder;
    }
    filePathType logFolderPath = "./logs/" + logFolderName;
    auto infoLogger = spdlog::basic_logger_mt(logFolderName, logFolderPath+ "/infoLogger.txt");
    dataVecType thetaVec;
    dataVecType phiVec;
    dataVecType dataVec;
    readData(thetaVec,phiVec,dataVec,fileThetaPath,
        filePhiPath, dataPath);
#ifdef DATA_VALUE_SET_DEBUG
    cout << "the data:\n";
    cout << dataVec<< endl;
#endif
    cout << "dataFile len: " << dataVec.size() << endl;
    dataVecType::iterator maxElem = std::max_element(thetaVec.begin(),thetaVec.end());
    dataVecType::iterator minElem = std::min_element(thetaVec.begin(),thetaVec.end());

    cout << "maxElem: " << *maxElem << endl;
    cout << "minElem: " << *minElem << endl;
    myDouble thetaBegin = (*minElem) / 180.0 * M_PI;
    myDouble thetaEnd = (*maxElem) / 180.0 * M_PI;
    string fileNamePrefix = outDataFileRootPath;
    string strFileName = fileNamePrefix + "/stringName.txt";
    ofstream localOutFile;    

    const myDouble phiE = twopi;
    const myDouble phiB = 0.0;
    string logDir = logFolderPath + "/";

    for (typeIndex indexMode = 0; indexMode < modeVec.size(); indexMode++) {
        for (typeIndex indexdT = 0; indexdT < dTVec.size(); indexdT++) {
            for (typeIndex indexTHR = 0; indexTHR < thresholdVec.size(); indexTHR++) {
                int mode = modeVec[indexMode];
                myDouble deltaTheta = dTVec[indexdT];
                myDouble threshold = thresholdVec[indexTHR];
                spContainer vecPoint;
                infoLogger -> debug("need to modify here, set the range of map manually!!");
                Healpix_Correlation hp(deltaTheta,
                thetaEnd, thetaBegin,
                phiE, phiB, logDir);
                hp.SetPixValue(thetaVec,phiVec,dataVec,vecPoint,threshold,mode);
        // typeIndex endIndex = hp.GetBottomPixel();
        // for(typeIndex index = 0; index <= endIndex; index++) {
        //     cout << "index: " << index << " value: " << hp[index] << endl;
        // }
                string fileNamePrefix = outDataFileRootPath;
                countType pointsCount = hp.GetPointsCount();
                myDouble scale = 2;
                countType randomCount = pointsCount * scale;
                myDouble spannedAngle = hp.GetSpannedAngle();

            std::stringstream stream;
            // change that threshold to a string
            int thresholdInt = abs(int(threshold));

            stream << thresholdInt << "/";
            stream << "mode_" << mode << "/";
            stream << "points";
            stream << pointsCount;
            stream << "avgTime_" << AVERAGE_TIME;
            stream << "nr_" << randomCount;
            stream << "thd_" << threshold;
            stream << "dT";
            stream << "tB_"<< thetaBegin;
            stream << "tE_" << thetaEnd;
            stream << std::fixed << std::setprecision(5) << dTVec[indexdT];
            // std::string str_dT = stream.str();
            stream << "_endT_" << std::fixed << std::setprecision(5) << endTheta;
            // stream << "ft_" << FACTOR_SIGVAL;            

                        std::string str_endT = outDataFileSubPath + "/" + stream.str();
                        
                        localOutFile.open(outFilePre + str_endT + ".txt");
                        infoLogger -> info("localOutFile outFilePre: {},  str_endT: {}", outFilePre, str_endT);
                        // nameHealStr.push_back(str_endT);
                        string path = outDataFileRootPath + "/" + str_endT;
                        // if(!IsPathExist(path)) {
                        //     if (mkdir(path.c_str(), 0777) != 0 ) {
                        //         cout << "error creating direc: " << path << endl;
                        //         continue;
                        //     }
                        // }
                        if (createPath(0777, outDataFileRootPath,str_endT) == -1) {
                            cout << "error creating dir_" << path << endl;
                            continue;
                        }    
                        infoLogger -> info("The strFileName: {}", strFileName);
                        infoLogger -> flush();

                        outFile.open(strFileName,std::ofstream::out | std::ofstream::app);
                        outFile << str_endT << endl;
                        outFile.close();
                    string filePathName = path + "/" ;
                string healpixMapFileName = filePathName + "SignalAverageHealpixMap.txt";
                // string randomHealpixMapFileName = fileNamePrefix + "Pi_randomHealpixRandom.txt";
                string peebles1 = filePathName + "SignalAverageHealpixMode1.txt";
                string peebles2 = filePathName + "SignalAverageHealpixMode2.txt";
                string hamilton = filePathName + "SignalAverageHealpixMode3.txt";
                string landy = filePathName + "SignalAverageHealpixMode4.txt";
                string randomMap = filePathName + "randomMap.txt";
                string calculationInfo = filePathName + "calculationInfo.txt";
                string spherePointName = filePathName + "spherePoint.txt";


                outFile.open(spherePointName);
                if (!outFile.is_open()) {
                    cout << "spherePointName: " << spherePointName << " not opened\n";
                    continue;
                }
                outFile << vecPoint;
                outFile.close();

                cout << "The point count is " << pointsCount << endl;

                outFile.open(calculationInfo);
                outFile << pointsCount << endl;
                outFile << randomCount << endl;
                outFile << deltaTheta << endl;
                outFile << AVERAGE_TIME << endl;
                outFile << spannedAngle << endl;
                outFile << endTheta << endl;
                outFile << thetaBegin << endl;
                outFile << thetaEnd << endl;
                outFile << threshold << endl;
                outFile << FACTOR_SIGVAL << endl;
                outFile << outFilePre << endl;
                outFile.close();

                vec_coord peebles1VC, peebles2VC, hamiltonVC,landyVC;
                dataVecType healpixMapInfo;
                dataVecType randomMapInfo;
                hp.sweepAverageTwoPCorrelation(deltaTheta,endTheta,peebles1VC,
                    peebles2VC,hamiltonVC,landyVC,AVERAGE_TIME,randomCount,healpixMapInfo, localOutFile);

                outFile.open(peebles1);
                outFile << peebles1VC;
                outFile.close();

                outFile.open(healpixMapFileName);
                outFile << healpixMapInfo;
                outFile.close();

                outFile.open(peebles2);
                outFile << peebles2VC;
                outFile.close();


                outFile.open(hamilton);
                outFile << hamiltonVC;
                outFile.close();

                outFile.open(landy);
                outFile << landyVC;
                outFile.close();    

                hp.outMapAndRandom(healpixMapInfo,randomMapInfo);
                outFile.open(randomMap);
                outFile << randomMapInfo;
                outFile.close();

                localOutFile.close();
            }
        }
    }
    infoLogger -> flush();
}

void healpixMapTest(const string & dataFilePath, 
    const vector<myDouble> & dTVec, 
    const vector<myDouble> & ThetaBegin,
    const vector<myDouble> & ThetaEnd, 
    const vector<myDouble> & PhiBegin,
    const vector<myDouble> & PhiEnd,
    const vector<countType> & pCount,
    int runTime,
        ofstream & outFile,
    const string & outFilePreCompact, const myDouble & endTheta,
    const myDouble & r2dFactor)
{


    vector<string> nameHealStr;
    ofstream localOutFile;
    for (typeIndex indexdT = 0; indexdT < dTVec.size(); indexdT++) {
        for (typeIndex indexBeg = 0; indexBeg < ThetaBegin.size(); indexBeg++) {
            for (typeIndex indexEnd = 0; indexEnd < ThetaEnd.size(); indexEnd++) {
                for (typeIndex indexPhiBeg = 0; indexPhiBeg < PhiBegin.size(); indexPhiBeg++) { 
                    for (typeIndex indexPhiEnd = 0; indexPhiEnd < PhiEnd.size(); indexPhiEnd++) {
                        for (typeIndex indexPoints = 0; indexPoints < pCount.size(); indexPoints++) {
                            for (int count = 0; count < runTime; count++) {

                                Healpix_Correlation map3(dTVec[indexdT], ThetaEnd[indexEnd],ThetaBegin[indexBeg], PhiEnd[indexPhiEnd], PhiBegin[indexPhiBeg]);
                                map3.SetPixValue(pCount[indexPoints]);
                                // cout << "the order: " << map3.Order() << endl;
                                // cout << "bottompixel: " << map3.GetBottomPixel() << endl;
                                // cout << "spanned: " << map3.GetSpannedAngle() << endl;
                                // cout << "pointCount: " << map3.GetPointsCount() << endl;
                                // cout << "deltaTheta: " << map3.GetDeltaTheta() << endl;
                                int nTime = AVERAGE_TIME;
                                vec_coord peebles1VC;
                                vec_coord peebles2VC;
                                vec_coord hamiltonVC;
                                vec_coord landyVC;
                                std::vector<healpixMapType> healpixMapInfo;
                                std::vector<double> randomMapInfo;
                                // int vecLen = ceil((double)endTheta/dT) + 1;
                                // std::vector<myDouble> bins;
                                // bins.resize(vecLen);
                                // for(int index = 0; index < vecLen; index+= 1) {
                                //     // the end may not be accurate
                                //     bins[index] = index * dT;
                                // }
                                
                                // map3.twoPCorrelation(deltaTheta,endTheta,vC,map4,healpixMapInfo,randomMapInfo,1);
                                // cout << "healpixMapInfo length: " << healpixMapInfo.size() << endl;

                                // auto t2 = std::chrono::high_resolution_clock::now();
                                // std::chrono::duration<double, std::milli> sequentialTime = t2 - t1;
                                // cout << "the time for the first part is: " << sequentialTime.count() << endl;
                                // cout << "using the pi\n";

                                // t1 = std::chrono::high_resolution_clock::now();

                                countType map1Points = map3.GetPointsCount();
                                    // countType map2Points = HEALPIX_RANDOM_POINT_COUNT;
                                countType map2Points = map1Points * r2dFactor;
                    std::stringstream fileNameTemp;
                    fileNameTemp << "map1P_" << std::fixed << std::setprecision(5) << map1Points;
                    fileNameTemp << "count_" << count << ".txt";
                                string fileOutName = outFilePreCompact + fileNameTemp.str();
                                localOutFile.open(fileOutName);

                    std::stringstream stream;
                    stream << "pt";
                    stream << pCount[indexPoints];
                    stream << "agP_";
                    stream << map2Points;
                    stream << "agT_" << AVERAGE_TIME;
                    stream << "dT";
                    stream << std::fixed << std::setprecision(5) << dTVec[indexdT];
                    // std::string str_dT = stream.str();
                    stream << "tB";
                    stream << std::fixed << std::setprecision(5) << ThetaBegin[indexBeg];
                    stream << "tE";
                    stream << std::fixed << std::setprecision(5) << ThetaEnd[indexEnd];
                    stream << "pB_";
                    stream << std::fixed << std::setprecision(5) << PhiBegin[indexPhiBeg];
                    stream << "pE_";
                    stream << std::fixed << std::setprecision(5) << PhiEnd[indexPhiEnd];

                    stream << "eT_" << std::fixed << std::setprecision(5) << endTheta;
                    std::string fileName = stream.str();
                    stream << "/rN_" << count;
                    std::string str_endT = stream.str();
                                
                                    string strFileName = dataFilePath + "stringName.txt";
                                    outFile.open(strFileName,std::ofstream::out | std::ofstream::app);
                                    outFile << str_endT << endl;
                                    outFile.close();



                                string path = dataFilePath + fileName;
                                if(!IsPathExist(path)) {
                                    if (mkdir(path.c_str(), 0777) != 0 ) {
                                        cout << "error creating direc: " << path << endl;
                                        continue;
                                    }
                                }
                                if (createPath(0777, dataFilePath,str_endT) == -1) {
                                    cout << "error creating dir_" << path << endl;
                                    continue;
                                }                                    
                                string fileNamePrefix = dataFilePath + str_endT;
                                string healpixMapFileName = fileNamePrefix + "Pi_HealpixAverageRandomMap.txt";
                                // string randomHealpixMapFileName = fileNamePrefix + "Pi_randomHealpixRandom.txt";
                                string peebles1 = fileNamePrefix + "Pi_HealpixAverageRandomMode1.txt";
                                string peebles2 = fileNamePrefix + "Pi_HealpixAverageRandomMode2.txt";
                                string hamilton = fileNamePrefix + "Pi_HealpixAverageRandomMode3.txt";
                                string landy = fileNamePrefix + "Pi_HealpixAverageRandomMode4.txt";
                                string calculationInfo = fileNamePrefix + "calculationInfo.txt";
                                string randomMap = fileNamePrefix + "randomMap.txt";
                                myDouble dT = map3.GetDeltaTheta();
                                myDouble spannedAngle = map3.GetSpannedAngle();
                                
                                cout << "avergeRR nk : " << map2Points << endl;
                                cout << "calculationInfo: " << calculationInfo << endl;

                
                                outFile.open(calculationInfo);
            if (outFile.is_open()) {
                            cout << "file opend\n";
                        }
                        else {
                            break;
                        }


                                outFile << map1Points << endl;
                                outFile << map2Points << endl;
                                outFile << dT << endl;
                                outFile << AVERAGE_TIME << endl;
                                outFile << spannedAngle << endl;
                                outFile << endTheta << endl;   
                                outFile << ThetaBegin[indexBeg]<< endl;
                                outFile << ThetaEnd[indexEnd] << endl;
                                outFile << PhiBegin[indexPhiBeg]<<endl;
                                outFile << PhiEnd[indexPhiEnd]<< endl;
                                outFile.close();

                                map3.sweepAverageTwoPCorrelation(dTVec[indexdT],endTheta,peebles1VC,
                                    peebles2VC,hamiltonVC,landyVC,nTime,map2Points,healpixMapInfo,localOutFile);
                                
                                outFile.open(peebles1);
                                outFile << peebles1VC;
                                outFile.close();

                                outFile.open(healpixMapFileName);
                                outFile << healpixMapInfo;
                                outFile.close();

                                outFile.open(peebles2);
                                outFile << peebles2VC;
                                outFile.close();


                                outFile.open(hamilton);
                                outFile << hamiltonVC;
                                outFile.close();

                                outFile.open(landy);
                                outFile << landyVC;
                                outFile.close();
                                
                                                    // for the angular power spectrum analysis
                                map3.outMapAndRandom(healpixMapInfo,randomMapInfo);
                                outFile.open(randomMap);
                                outFile << randomMapInfo;
                                outFile.close();
                                localOutFile.close();
                            }
                        }
                    }
                }
            }
        }

    }
}

void mainTest(const string & configFile, const string & testName,
const myDouble & endTheta) {
    // vector<string>  dataFolderParentPaths;
    // vecVecString  vecVecFolderNames;
    // vector<string>  writeDataFolderParentPaths;
    // vecVecString writeFolderNames;
    // vecVecDouble  vecVecThresholds;
    // vecVecInt  vecVecModes;
    // readConfig(configFile, dataFolderParentPaths, vecVecFolderNames,
    //     writeDataFolderParentPaths, writeFolderNames,
    //     vecVecThresholds, vecVecModes);

    ofstream outFile;
    spdlog::set_level(spdlog::level::debug);
    if (testName == "dataTest") {
        dataTest(configFile, endTheta,outFile);
    }
}



ifstream & eatBlanks(ifstream & inFile) {
    char ch;
    inFile >> ch;
    while (ch == '\n') {
        inFile >> ch;
    }
    // put back the char
    inFile.putback(ch);
    return inFile;
}


void getQuotedString(ifstream & inFile, 
string & dataString,
const string & loggerName) {
    string compactString;
    std::getline(inFile, compactString);
    // cout << "the compactString : " << compactString << endl;
    // split the string
    std::istringstream iss(compactString);
    iss >> std::quoted(dataString);
}

void getQuotedVecString(ifstream & inFile, 
    vector<string>& vecString,
    const string & loggerName) {
    // clear if necessary
    vecString.clear();
    string compactString;
    std::getline(inFile, compactString);
    // cout << "the compactString : " << compactString << endl;
    // split the string
    std::istringstream iss(compactString);
    do {
        string folderName;
        iss >> std::quoted(folderName);
        // iss >> folderName;
        if (folderName.length() >= 1)
            vecString.push_back(folderName);
        // cout << "folderName : " << folderName << endl;
    }while(iss);
}


void getDoubleVec(ifstream & inFile, 
    vector<myDouble> & vecDoubles,
    const string & loggerName) {
    vecDoubles.clear();
    string compactString;
    std::getline(inFile, compactString);
    std::istringstream iss(compactString);
    auto ioLogger = spdlog::get(loggerName);

    do {
        string tempStr;
        iss >> tempStr;
        double val;
        if (tempStr == "nan") {
            val = MIN_RSS;
        }
        else {
            try {
                val = std::stod(tempStr);
            }
            catch(std::invalid_argument & ) {
                if (tempStr != "")
                    // cout << "invalid argument, tempStr: " << tempStr << endl;
                    ioLogger -> warn("Invalid Argument: {}", tempStr);
                continue;
            }
        }
        vecDoubles.push_back(val);
    }while(iss);    
}


/**
 * help function read in correct int
 * */
void getIntVec(ifstream & inFile, 
    vector<int> & vecInts,
    const string & loggerName) {

    vecInts.clear();
    string compactString;
    std::getline(inFile, compactString);
    std::istringstream iss(compactString);
    auto ioLogger = spdlog::get(loggerName);
    
    do {
        string tempStr;
        iss >> tempStr;
        int mode;
        try {
            mode = std::stoi(tempStr);
        }
        catch (std::invalid_argument &) {
            if (tempStr != "")
                // cout << "invalid argument: tempStr : " << tempStr << endl;
                ioLogger -> warn("Invalid Argument: {}", tempStr);
            continue;
        }
        vecInts.push_back(mode);
    }while(iss);

}
