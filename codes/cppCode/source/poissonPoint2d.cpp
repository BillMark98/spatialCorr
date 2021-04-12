#include "poissonPoint2d.h"
std::mutex poissonPointMutex;

void PoissonPoint2d::rPoissonSphereDegreeLambdaDefault(vec_3DCoord & vec3d)
{
    rPoissonSphereDegree(vec3d, lambdaDefault);
}

void PoissonPoint2d::rPoissonSphereDegree(vec_3DCoord & vec3d, 
    myDouble lambda0, myDouble thetaEnd, myDouble thetaBegin)
{   
    if (thetaEnd < thetaBegin) {
        cout << "\n!!!!thetaEnd < thetaBegin!!!!!" << endl;
        cout << "cant generate poisson points\n";
        exit(ANGLEBEGIN_LARGER_THAN_ANGLEEND);
    }
    std::default_random_engine generator;
    myDouble area = twopi * (cos(thetaBegin) - cos(thetaEnd)) * square(radius);
    std::poisson_distribution<int> distribution((double)lambda0 * area);
    sizeType number = distribution(generator);
#ifdef POISSON_DEBUG
    cout << "in the rPoissonSphereDegree\n";
    cout << "number of poisson points on a sphere: " << number << endl;
#endif
    rBinomialSphereDegreeN(vec3d,number,thetaEnd,thetaBegin);
}

void PoissonPoint2d::rPoissonSphereDegree(vec_3DCoord & vec3d, 
    myDouble lambda0, myDouble thetaEnd, 
    myDouble thetaBegin, myDouble phiEnd, myDouble phiBegin)
{   
    if (thetaEnd < thetaBegin) {
        cout << "\n!!!!thetaEnd < thetaBegin!!!!!" << endl;
        cout << "cant generate poisson points\n";
        exit(ANGLEBEGIN_LARGER_THAN_ANGLEEND);
    }
    std::default_random_engine generator;
    myDouble area = (phiEnd - phiBegin)* (cos(thetaBegin) - cos(thetaEnd)) * square(radius);
    std::poisson_distribution<int> distribution((double)lambda0 * area);
    sizeType number = distribution(generator);
#ifdef POISSON_DEBUG
    cout << "in the rPoissonSphereDegree\n";
    cout << "number of poisson points on a sphere: " << number << endl;
#endif
    rBinomialSphereDegreeN(vec3d,number,thetaEnd,thetaBegin, phiEnd, phiBegin);
}


void PoissonPoint2d::rPoissonNSphereDegree(vec_3DCoord & vec3d, sizeType _N,
    myDouble _thetaEnd,
    myDouble _thetaBegin,
    myDouble _phiEnd,
    myDouble _phiBegin)
{
    myDouble _lambda = calculateLambda(_N,radius,_thetaEnd, _thetaBegin,_phiEnd, _phiBegin);
    rPoissonSphereDegree(vec3d, _lambda,_thetaEnd, _thetaBegin,_phiEnd, _phiBegin);
}


void PoissonPoint2d::MaternProcessLambda( const vec_3DCoord & vecPosition, vec_3DCoord & vecData, myDouble _lambda, 
    myDouble _dthetaEnd, myDouble _dthetaBegin)
{
    vec3DConstIter endIter = vecPosition.end();
    for(vec3DConstIter iter = vecPosition.begin(); iter != vecPosition.end(); iter++) {
        lambdaPointPoisson(*iter,_lambda,vecData, _dthetaEnd, _dthetaBegin);
    }
}

void PoissonPoint2d::MaternProcessChooseLambda(const vec_3DCoord & vecPosition, vec_3DCoord & vecData,
    myDouble _lambda, sizeType _M, myDouble _dthetaEnd, myDouble _dthetaBegin)
{

    sizeType _len = vecPosition.size();
    if(_len <= _M) {
        MaternProcessLambda(vecPosition, vecData,_lambda,_dthetaEnd, _dthetaBegin);
    }
    else {
        vec_Index vecIndex;
        cout << "To do generateNIndices in MaternProcessChooseLambda, use _dthetaEnd instead\n";
        GenerateNIndices(_len, _M,_dthetaEnd, vecIndex,vecPosition);
        for(typeIndex i = 0; i < _M; i++) {
            typeIndex index = vecIndex[i];
            lambdaPointPoisson(vecPosition[index], _lambda,vecData, _dthetaEnd, _dthetaBegin);
        }
    }
}

void PoissonPoint2d::MaternProcessNPoint(const vec_3DCoord & vecPosition,
    vec_3DCoord & vecData,
    sizeType _N, myDouble _thetaEnd, myDouble _thetaBegin)
{
//     vec3DIter endIter = vec3d.end();


// // #else
//     for(vec3DIter iter = vec3d.begin(); iter != endIter; iter++) {
//         nPointPoisson(*iter,_N,_theta,vec3d);
//     }

    sizeType _len = vecPosition.size();
    for(typeIndex i = 0; i < _len; i++) {
        nPointPoisson(vecPosition[i],vecData,_N, _thetaEnd, _thetaBegin);
    }
}

void PoissonPoint2d::MaternProcessChooseNPoint(const vec_3DCoord & vecPosition, 
    vec_3DCoord & vecData,
    sizeType _N, sizeType _M, myDouble _thetaEnd, myDouble _thetaBegin)
{
    sizeType _len = vecPosition.size();
    if(_len <= _M) {
        MaternProcessNPoint(vecPosition,vecData,_N,_thetaEnd, _thetaBegin);
    }
    else {
        vec_Index vecIndex;
        cout << "GenerateNIndices in MaternProcessChooseNpoint\nstay tuned how far should the cluster points be apart¿\n";
        GenerateNIndices(_len, _M, _thetaEnd,vecIndex, vecPosition);
        for(typeIndex i = 0; i < _M; i++) {
            typeIndex index = vecIndex[i];
            nPointPoisson(vecPosition[index], vecData,  _N, _thetaEnd,_thetaBegin);
        }
    }
}

void PoissonPoint2d::MaternProcessNPointParallel(const vec_3DCoord & vec3d,
    vec_3DCoord & vecData,
    sizeType _N, myDouble _thetaEnd, myDouble _thetaBegin)
{
    int threadCount = MAX_THREAD;
    sizeType _len = vec3d.size();
    sizeType stepSize = ceil(_len / threadCount);
    std::vector<std::thread> vt;
    for(int i = 0; i < threadCount; i++) {
        int _begin, _end;
        _begin = i * stepSize;
        if(i != threadCount - 1) {
            _end = (i + 1) * stepSize;
        }
        else {
            _end = _len;
        }
        vt.push_back(std::thread(std::bind(&PoissonPoint2d::indexParallelPoissonNPoint, 
            *this, _begin, _end, _N,std::cref(vec3d), std::ref(vecData), _thetaEnd, _thetaBegin)));
    }
    for(int i = 0; i < threadCount; i++) {
        vt[i].join();
    }
}            


void PoissonPoint2d::MaternProcessChooseNPointParallel(const vec_3DCoord & vecPosition, 
    vec_3DCoord & vecData, sizeType _N, sizeType _M, myDouble _thetaEnd, myDouble _thetaBegin)
{
    sizeType _len = vecPosition.size();

#ifdef MATERN_DEBUG
    cout << "_len: " << _len << " vecPosition.size(): " << vecPosition.size() << endl;
#endif
#ifndef MATERN_DEBUG
    if(_M >= _len) {
        MaternProcessNPointParallel(vecPosition,vecData, _N, _thetaEnd, _thetaBegin);
    }
    else {
#ifdef MATERN_DEBUG
    cout << "_len: " << _len << " vec3d.size(): " << vec3d.size() << endl;
#endif
        vec_Index vecIndex;
        cout <<" GenerateNIndices in MaternProcessChooseNpointParallel\nstay tuned, how far should the cluster center be apart¿\n";
        GenerateNIndices(_len, _M, _thetaEnd, vecIndex,vecPosition);
        int threadCount = MAX_THREAD;
        sizeType _len = vecIndex.size();
        sizeType stepSize = ceil(_len / threadCount);
        std::vector<std::thread> vt;
        for(int i = 0; i < threadCount; i++) {
            int _begin, _end;
            _begin = i * stepSize;
            if(i != threadCount - 1) {
                _end = (i + 1) * stepSize;
            }
            else {
                _end = _len;
            }
            vt.push_back(std::thread(std::bind(&PoissonPoint2d::vecIndexParallelPoissonNPoint, 
                *this, _begin, _end, std::cref(vecIndex), _N, std::cref(vecPosition), std::ref(vecData), _thetaEnd, _thetaBegin)));
        }
        for(int i = 0; i < threadCount; i++) {
            vt[i].join();
        }
    }
#endif
}            

void PoissonPoint2d::MaternProcess(vec_3DCoord & vecData,myDouble lambdaParent, myDouble lambdaDaughter, myDouble _spanThetaEnd, myDouble _spanThetaBegin,
    myDouble _dthetaEnd, myDouble _dthetaBegin)
{
    if (_spanThetaEnd < _spanThetaBegin || _dthetaEnd < _dthetaBegin || _spanThetaEnd - _dthetaEnd < _spanThetaEnd + _dthetaEnd) {
        cout << "angles invalid, _spanThetaEnd >= _spanThetaBegin, _dthetaEnd >= _dthetaBegin, _spanThetaEnd - _dthetaEnd >= _spanThetaEnd + _dthetaEnd\n";
        cout << "_spanThetaEnd: " << _spanThetaEnd << ", _spanThetaBegin: " << _spanThetaBegin << ", _dthetaEnd: " << _dthetaEnd << ", _dthetaBegin: " << _dthetaBegin << endl;
        exit(ANGLEBEGIN_LARGER_THAN_ANGLEEND);
    }
    vec_3DCoord vecPosition;
    rPoissonSphereDegree(vecPosition, lambdaParent, _spanThetaEnd - _dthetaEnd, _spanThetaBegin + _dthetaEnd);
    MaternProcessLambda(vecPosition,vecData,lambdaDaughter,_dthetaEnd, _dthetaBegin);
}   

void PoissonPoint2d::deleteCluster(vec_3DCoord & vec3d, 
    sizeType _M, myDouble _theta)
{
    sizeType _len = vec3d.size();
    std::srand(time(0));
    typeIndex index = rand() % _len;
    sizeType count = _M;
    while(count > 0 && !(vec3d.empty()) ) {
        vec3d[index].popThetaRing(vec3d, index, _theta);
        index = rand() % vec3d.size();
    }
}    

void PoissonPoint2d::nPointPoisson(const SphericalPoint & p0, 
    vec_3DCoord & vec3d, 
    sizeType _n, 
    myDouble _thetaEnd,
    myDouble _thetaBegin)
{
    vec_3DCoord vecTemp;
    rPoissonNSphereDegree(vecTemp,_n,_thetaEnd,_thetaBegin);
    for(vec3DIter iter = vecTemp.begin(); iter != vecTemp.end(); iter++) {
        iter -> setCoord(p0);
    }
    // insert at the end will be best choice considering the time complexity.
    vec3d.insert(vec3d.end(),std::make_move_iterator(vecTemp.begin()),std::make_move_iterator(vecTemp.end())); 

}

void PoissonPoint2d::indexParallelPoissonNPoint(typeIndex _start, 
    typeIndex _end,  sizeType _N,
    const vec_3DCoord & vecPosition, 
    vec_3DCoord & vecData,
    myDouble _thetaEnd,
    myDouble _thetaBegin)
{
    vec_3DCoord vecTemp;
    for(typeIndex i = _start; i < _end; i++) {
        nPointPoisson(vecPosition[i],vecTemp, _N, _thetaEnd, _thetaBegin);
    }
    std::lock_guard<std::mutex> guard(poissonPointMutex);
    vecData.insert(vecData.end(), std::make_move_iterator(vecTemp.begin()),std::make_move_iterator(vecTemp.end()));
}

void PoissonPoint2d::vecIndexParallelPoissonNPoint(typeIndex _start,
    typeIndex _end, const vec_Index & vecIndex, sizeType _N,
    const vec_3DCoord & vecPosition,
    vec_3DCoord & vecData,
    myDouble _thetaEnd,
    myDouble _thetaBegin)
{
    vec_3DCoord vecTemp;

    for(typeIndex i = _start; i != _end; i++) {
        nPointPoisson(vecPosition[vecIndex[i]], vecTemp, _N, _thetaEnd,_thetaBegin);
    }
    std::lock_guard<std::mutex> guard(poissonPointMutex);
    vecData.insert(vecData.end(), std::make_move_iterator(vecTemp.begin()),std::make_move_iterator(vecTemp.end()));
    #ifdef MATERN_PARALLEL_POISSON_DEBUG
    cout << "******************************\n";
    cout << "check the point not zero in PoissonPoint2d::vecIndexParallelPoissonNPoint\n";
    for(typeIndex i =_start; i != _end; i++) {
        cout << "index: i" << i << endl;
        cout << vec3d[i] << endl;
        cout << "vecIndex[i]: " << vecIndex[i] << endl;
        cout << vec3d[vecIndex[i]] << endl;
    }
#endif
}            

void PoissonPoint2d::lambdaPointPoisson(const SphericalPoint & p0, 
    myDouble _lambda, vec_3DCoord & vec3d, myDouble _thetaEnd, myDouble _thetaBegin)
{
    vec_3DCoord vecTemp;
    rPoissonSphereDegree(vecTemp,_lambda,_thetaEnd, _thetaBegin);
    for(vec3DIter iter = vecTemp.begin(); iter != vecTemp.end(); iter++) {
        iter -> setCoord(p0);   
    }
    vec3d.insert(vec3d.end(),std::make_move_iterator(vecTemp.begin()),
        std::make_move_iterator(vecTemp.end()));
}

void PoissonPoint2d::rBinomialSphereDegreeN(vec_3DCoord & vec3d, 
    sizeType n, myDouble thetaEnd, myDouble thetaBegin,
    myDouble phiEnd, myDouble phiBegin)
{
     // generate random points
    sizeType vect_len = 2 * n;
    vector<myDouble> points(vect_len);
    //Will be used to obtain a seed for the random number engine
    std::random_device rd;  
    //Standard mersenne_twister_engine seeded with rd()
    std::mt19937 gen(rd()); 
    std::uniform_real_distribution<> dis(0.0, 1.0);
    for (typeIndex index = 0; index < vect_len; ++index) {
        // Use dis to transform the random unsigned int generated by gen into a 
        // myDouble in [0, 1). Each call to dis(gen) generates a new random myDouble
        points[index] = dis(gen);
    }
#ifdef RANDOM_SHUFFLE
        std::srand(time(0));
        std::random_shuffle(points.begin(), points.end());
#endif
    vec3d.resize(n);
    myDouble cosBegin = cos(thetaBegin);
    myDouble cosEnd = cos(thetaEnd);
    myDouble phiDiff = phiEnd - phiBegin;
    bool outPoint = false;
    for (typeIndex index = 0; index < n; index++) {
        myDouble phiTemp = points[index * 2] * phiDiff;
        myDouble thetaTemp = acos(cosBegin - points[index * 2 + 1] * (cosBegin- cosEnd));
    if (thetaTemp < thetaBegin - MY_EPSILON || thetaTemp > thetaEnd + MY_EPSILON || phiTemp < phiBegin - MY_EPSILON || phiTemp > phiEnd + MY_EPSILON) {
        cout << "\n!!!!!theta or phi is out of the begin, end range!!!!\n";
        cout << "Index: " << index << endl;
        cout << "Point: " << endl;
        outPoint = true;
    }
        vec3d[index].setRTP(radius,thetaTemp, phiTemp);
        if (outPoint) {
            cout << vec3d[index] << endl;
            outPoint = false;
        }
    }
    
#ifdef POISSON_DEBUG
    cout << "in the rBinomialSphereDegreeN\n";
    cout << "the function successfully called" << endl;
#endif
}


void PoissonPoint2d::setN() {
    myDouble sphere = 2 * pi * square(radius) * (1 - cos(spannedAngle));
    // the mean is lambdaDefault * sphere = N
    N = (int)lambdaDefault * sphere;
}

void PoissonPoint2d::setLambdaDefault() {
    myDouble sphere = 2 * pi * square(radius) * (1 - cos(spannedAngle));
    // the mean is lambdaDefault * sphere = N
    lambdaDefault = N/(sphere);
}

myDouble PoissonPoint2d::calculateLambda(sizeType n, 
    myDouble R, myDouble thetaEnd,myDouble thetaBegin,
    myDouble phiEnd, myDouble phiBegin)
{
    myDouble sphere = (phiEnd - phiBegin)* square(R) * (cos(thetaBegin) - cos(thetaEnd));
    // the mean is lambdaDefault * sphere = N
    return n/(sphere);
}

void PoissonPoint2d::GenerateNIndices(sizeType _N, sizeType _M, 
     myDouble smallAngle, 
    vec_Index &_vecIndex, const vec_3DCoord & vecPosition)
{
    if(_N < VEC_SET_THRESHOLD) {
        cout << "warning!\nin GenerateNIndices, do not consider the chosen point must be at least some distance apart\n";
        vec_Index vecTemp;
        vecTemp.resize(_N);
        for(typeIndex i = 0; i < _N; i++) {
            vecTemp[i] = i;
        }
        std::srand((unsigned) time(0));
        std::random_shuffle(vecTemp.begin(),vecTemp.end());
        vec_Index::iterator iter = vecTemp.begin();
        std::advance(iter, _M - 1);
        _vecIndex.insert(_vecIndex.end(),std::make_move_iterator(vecTemp.begin()), std::make_move_iterator(iter));
    }
    else {
        std::vector<typeIndex> numbers;
        sizeType count = 0;
        while (numbers.size() < _M && count < MAX_LOOP_TIME)
        {
            typeIndex tempNum = rand() % _N;
            bool flagAdmit = true;
            // make sure that cluster center are far enough apart,
            // i.e. d(p1, p2) >= smallAngle
            for(typeIndex _i = 0; _i < numbers.size(); _i++) {
                if(vecPosition[numbers[_i]].withinThetaRing(vecPosition[tempNum],0,smallAngle)) {
                    flagAdmit = false;
                    break;
                }
            }
            if(flagAdmit) {
                numbers.push_back(tempNum);
            }
            count++;
        }
        _vecIndex.insert(_vecIndex.end(),std::make_move_iterator(numbers.begin()), std::make_move_iterator(numbers.end()));
        if(count >= MAX_LOOP_TIME) {
            cout << "PoissonPoint2d::GenerateNIndices has run too many loops, the returned vector index" << endl;
            cout << "may not be in the wanted length" << endl;
        }
    }

}