#include "sphereTwoCorrelation.h"

// sphereTwoCorrelation::sphereTwoCorrelation(spherePointsContainer & v1, myDouble Angle, countType num) : sPContainer(v1),
// fpSpherePoint(v1), fpPoisson(PoissonPoint2d()){
//     // update the pointCount
//     // have to implement before setSphereAngle, because otherwise
//     // the isEmpty function in maxTheta will be evaluated to be true, thus quit the function
//     pointCount = v1.size();
//     if (fabs(Angle - M_PI_2) < MY_EPSILON) {
//         setSphereAngle();
//     }
//     setRadius();
    
// }

sphereTwoCorrelation::sphereTwoCorrelation(spherePointsContainer & v1,PoissonPoint2d & pd, myDouble Angle, bool setAngle) : sPContainer(v1),
fpSpherePoint(v1), fpPoisson(pd),spCont(v1), hpMap(), randomHpMap(){
    // update the pointCount
    // have to implement before setSphereAngle, because otherwise
    // the isEmpty function in maxTheta will be evaluated to be true, thus quit the function
    pointCount = v1.size();
    if(setAngle) {
        setSphereAngle();
    }
    else {
        sphereAngle = Angle;
    }
    setRadius();
    healpixMapSet = false;
    deltaTheta = 0.;
}


void sphereTwoCorrelation::pushPoint(SphericalPoint & p) {
    sPContainer.push_back(p);
}

long sphereTwoCorrelation::DDTheta(myDouble theta, myDouble _dT) const {
    long count = 0L;
    for (spCont_ConstIter i = sPContainer.begin(); i != sPContainer.end(); i++) {
        for (spCont_ConstIter j = i; j != sPContainer.end(); j++) {
            if (*j == *i) {
                continue;
            }
            if (i -> withinThetaRing(*j, theta, _dT)) {
                count++;
            }
        }
    }
    return count;
}

void sphereTwoCorrelation::setSphereAngle() {
    myDouble maxThe = maxTheta();
    myDouble minThe = minTheta();
    myDouble maxPh = maxPhi();
    myDouble minPh = minPhi();
    sphereAngle = (cos(minThe) - cos(maxThe)) * (maxPh - minPh);
}
// assume all the points lie on the sphere
void sphereTwoCorrelation::setRadius() {
    if (isEmpty()) {
        cerr << "sphere data points empty in setRadius" << endl;
        exit(NO_POINTS_ON_SPHERE);
    }
    radius = sPContainer[0].getRadius();
}
myDouble sphereTwoCorrelation::twoPointCorrelation(myDouble theta, 
    myDouble dT, int mode){
    if (mode == 0) {
        FindPair fp(sPContainer);
        long np = fpSpherePoint.findPairs((double)theta, (double)dT);
        // long ddCount = DDTheta(theta, _dT);
        // cout << "theta: " << theta << "\tddCount: " << ddCount << "\t fp: " << np << "is_equal? " << (np == ddCount ? "T" : "F") << endl;
        // here approximately deltaOmega = radius * sin(theta) * _dT * _dT (delta phi = _dT)
        // replace the varying theta with a fixed theta
        myDouble fixTheta = M_PI_4;
        myDouble val = (sphereAngle * np)/(pointCount * (pointCount - 1) * M_PI * radius * sin(fixTheta) * square(dT)) - 1;
        return val;
    }
    // use the estimation 1 + w(theta) = n_p/n_p(t) * (N_t/N)^2 where
    // n_p are the pairs of the given distribution, n_p(t) are that of a poisson point process
    // N_t are the total number of points of a poisson point process
    // N are that of the given data set
    if (mode == 1) {
        unsigned int num = fpPoisson.getNumber();
        // PoissonPoint2d poiPoint(radius, num);
        // FindPair fp(poiPoint);
        long npt = fpPoisson.findPairs((double)theta, (double)dT);
        long np = fpSpherePoint.findPairs((double)theta, (double)dT);
        // cout << "theta: " << theta << "\tnpt: " << npt << "\tnp: " << np << endl;
        // myDouble meta = square<myDouble>((myDouble)num/pointCount);
        // cout << "num: " << num << "\t pointCount: " << pointCount << "\t meta: " << meta << endl;
        // cout << endl;
        myDouble val = (myDouble)np/npt * square<myDouble>(myDouble(num)/pointCount) - 1;
        return val;
    }

    if (mode == 2) {
        unsigned int num = fpPoisson.getNumber();
        // PoissonPoint2d poiPoint(radius, num);
        // FindPair fp(poiPoint);
        long npt = fpPoisson.findPairs((double)theta, (double)dT);
        long np = fpSpherePoint.findPairs((double)theta,(double)dT);
        // cout << "theta: " << theta << "\tnpt: " << npt << "\tnp: " << np << endl;
        // myDouble meta = square<myDouble>((myDouble)num/pointCount);
        // cout << "num: " << num << "\t pointCount: " << pointCount << "\t meta: " << meta << endl;
        // cout << endl;
        myDouble val = (myDouble)np/npt * square<myDouble>(myDouble(num)/pointCount) - 1;
        return val;
    }

    if (mode == 5) {
        if(!isHealpixSet()) {
            hpMap.SetDeltaTheta(dT);
            randomHpMap.SetDeltaTheta(dT);
            healpixMapSet = true;
            hpMap.SetPixValue(spCont);
            randomHpMap.SetPixValue(HEALPIX_RANDOM_POINT_COUNT);
        }

    }
}   

myDouble sphereTwoCorrelation::maxTheta() const {
    if (isEmpty()) {
        cerr << " there are no points, cannot calculate maxTheta" << endl;
        exit(NO_POINTS_ON_SPHERE);
    }
    myDouble maxTheta = sPContainer.begin() -> getTheta();
    for(spCont_ConstIter iter = sPContainer.begin(); iter != sPContainer.end(); iter++) {
        myDouble degree = iter -> getTheta();
        if (degree > maxTheta) {
            maxTheta = degree;
        }
    }
    return maxTheta;
}
myDouble sphereTwoCorrelation::minTheta() const {
    if (isEmpty()) {
        cerr << " there are no points, cannot calculate minTheta" << endl;
        exit(NO_POINTS_ON_SPHERE);
    }
    myDouble minTheta = sPContainer.begin() -> getTheta();
    for(spCont_ConstIter iter = sPContainer.begin(); iter != sPContainer.end(); iter++) {
        myDouble degree = iter -> getTheta();
        if (degree < minTheta) {
            minTheta = degree;
        }
    }
    return minTheta;
}
myDouble sphereTwoCorrelation::maxPhi() const {
    if (isEmpty()) {
        cerr << " there are no points, cannot calculate maxPhi" << endl;
        exit(NO_POINTS_ON_SPHERE);
    }
    myDouble maxPhi = sPContainer.begin() -> getPhi();
    for(spCont_ConstIter iter = sPContainer.begin(); iter != sPContainer.end(); iter++) {
        myDouble degree = iter -> getPhi();
        if (degree > maxPhi) {
            maxPhi = degree;
        }
    }
    return maxPhi;
}
myDouble sphereTwoCorrelation::minPhi() const {
    if (isEmpty()) {
        cerr << " there are no points, cannot calculate minPhi" << endl;
        exit(NO_POINTS_ON_SPHERE);
    }
    myDouble minPhi = sPContainer.begin() -> getPhi();
    for(spCont_ConstIter iter = sPContainer.begin(); iter != sPContainer.end(); iter++) {
        myDouble degree = iter -> getPhi();
        if (degree < minPhi) {
            minPhi = degree;
        }
    }
    return minPhi;
}

countType sphereTwoCorrelation::Neighbors(myDouble theta, spherePointsContainer & v1, int order) {
    Healpix_Correlation map(order, RING);
    // initialize the map structure
    pixelCount pxCount;
    for(spherePointsContainer::const_iterator iter = v1.begin(); iter != v1.end(); iter++) {
        // assume distinct keys
        int pixel = map.vec2pix(vec3((double)iter ->getXCoord(),(double) iter -> getYCoord(), (double) iter -> getZCoord()));
        pixelCount::iterator pxcIter = pxCount.find(pixel);
        if(pxcIter != pxCount.end()) {
            pxcIter -> second += 1;
        }
        else {
            pxCount[pixel] = 1;
        }    
    }


    pixelCount::iterator pxIter;
    countType neighbourCount = 0;
    for(pxIter = pxCount.begin(); pxIter != pxCount.end(); pxIter++) {
        rangeset<int> tempRange;
        map.query_disc_pixel_Internal(pxIter -> first,theta,tempRange);
        vector<int> pixelIndex;
        tempRange.toVector(pixelIndex);
        countType sum = 0;
        for(typeIndex index = 0; index < pixelIndex.size(); index++) {
            sum += pxCount[(int)pixelIndex[index]];
        }
        neighbourCount += sum;
    }
    return neighbourCount;

}

void sphereTwoCorrelation::healpixTwoPointCorr(myDouble dT, myDouble endTheta, spherePointsContainer & v1, vec_coord & vC,int mode) 
{
    Healpix_Correlation map(dT, sphereAngle);
    map.SetPixValue(v1);
    Healpix_Correlation randomMap(dT, sphereAngle);
    randomMap.SetPixValue(HEALPIX_RANDOM_POINT_COUNT);
#ifndef DISCONTINUOUS_SETTING   
    dT = map.GetDeltaTheta();   
#endif    
    // e.g 0,0.5,1,1.5,2,2.2    ceil(2.2/0.5) = 5, starting from 0, so there are 6 end points.
    map.twoPCorrelation(dT,endTheta,vC,randomMap,mode);
}


void sphereTwoCorrelation::healpixTwoPointCorr(myDouble dT, myDouble endTheta, spherePointsContainer & v1, vec_coord & vC, 
    std::vector<myDouble> & healpixMapInfo, std::vector<myDouble> & randomHealpixInfo, int mode) 
{
    
    Healpix_Correlation map(dT,sphereAngle);
    map.SetPixValue(v1);
    
    Healpix_Correlation randomMap(dT, sphereAngle);  
    randomMap.SetPixValue(HEALPIX_RANDOM_POINT_COUNT);

#ifndef DISCONTINUOUS_SETTING   
    dT = map.GetDeltaTheta();   
#endif    
    map.twoPCorrelation(dT,endTheta,vC,randomMap,healpixMapInfo,randomHealpixInfo,mode);
}

void sphereTwoCorrelation::averageHealpixTwoPointCorrNoSet(myDouble dT, myDouble endTheta, 
    vec_coord & vC, Healpix_Correlation & Hmap, 
    int nTime, countType nk, 
    int mode)
{
    Hmap.averageTwoPCorrelation(dT,endTheta,vC,nTime,nk,mode);
}

void sphereTwoCorrelation::averageHealpixTwoPointCorrNoSet(myDouble _dT, 
    myDouble endTheta, vec_coord & vC,
    Healpix_Correlation & Hmap, 
    std::vector<myDouble> & healpixMapInfo,  
    countType & pointCount,
    myDouble & spannedAngle, myDouble & dT, 
    int nTime, 
    countType nk,
    int mode )
{
    Hmap.averageTwoPCorrelation(_dT,endTheta,vC,nTime,nk,healpixMapInfo,mode);
    dT = Hmap.GetDeltaTheta();
    spannedAngle = Hmap.GetSpannedAngle();
    pointCount = Hmap.GetPointsCount();
}

void sphereTwoCorrelation::averageHealpixTwoPointCorr(myDouble _dT, myDouble endTheta, spherePointsContainer & v1, vec_coord & vC, int nTime, countType nk,int mode) {
    Healpix_Correlation map(_dT,sphereAngle);
    map.SetPixValue(v1);
    map.averageTwoPCorrelation(_dT,endTheta,vC,nTime,nk,mode);
}

// for the plotting, additional healpixMap info
void sphereTwoCorrelation::averageHealpixTwoPointCorr(myDouble _dT, myDouble endTheta, spherePointsContainer & v1, vec_coord & vC, std::vector<myDouble> & healpixMapInfo,  int nTime, countType nk,int mode) {
    Healpix_Correlation map(_dT,sphereAngle);
    map.SetPixValue(v1);
    map.averageTwoPCorrelation(_dT,endTheta,vC,nTime,nk,healpixMapInfo,mode);
}

void sphereTwoCorrelation::averageHealpixTwoPointCorr(myDouble _dT, myDouble endTheta, spherePointsContainer & v1, vec_coord & vC, std::vector<myDouble> & healpixMapInfo, countType & pointCount, 
        myDouble & spannedAngle, myDouble & dT, int nTime, countType nk,int mode) {
    Healpix_Correlation map(_dT,sphereAngle);
    map.SetPixValue(v1);
    map.averageTwoPCorrelation(_dT,endTheta,vC,nTime,nk,healpixMapInfo, mode);
    dT = map.GetDeltaTheta();
    spannedAngle = map.GetSpannedAngle();
    pointCount = map.GetPointsCount();
}

void sphereTwoCorrelation::sweepAverageHealpix2PCF(myDouble _dT, myDouble endTheta, 
    spherePointsContainer &v1, vec_coord & peebles1, vec_coord & peebles2,
    vec_coord & hamilton, vec_coord & landy,int nTime, countType nk)
{
    cout << "SOMETHING WRONG with this function, dont use it right now\n";
    exit(DONT_CALL_ME);

    Healpix_Correlation map(_dT, sphereAngle);
    std::vector<myDouble> healpixDummy;
    map.SetPixValue(v1);
    map.sweepAverageTwoPCorrelation(_dT,endTheta,peebles1,peebles2,
        hamilton, landy,nTime,nk,healpixDummy);

}
void sphereTwoCorrelation::sweepAverageHealpix2PCF(myDouble _dT, myDouble endTheta, 
    spherePointsContainer &v1, vec_coord & peebles1, vec_coord & peebles2,
    vec_coord & hamilton, vec_coord & landy,
    std::vector<myDouble> & healpixMapInfo,  countType & pointCount,
    myDouble & spannedAngle, myDouble & dT, ofstream & os, int nTime, 
    countType nk)
{
    cout << "SOMETHING WRONG with this function, dont use it right now\n";
    exit(DONT_CALL_ME);


    Healpix_Correlation map(_dT, sphereAngle);
    map.SetPixValue(v1);
    map.sweepAverageTwoPCorrelation(_dT,endTheta,peebles1,peebles2,
        hamilton, landy,nTime,nk,healpixMapInfo,os);
    dT = map.GetDeltaTheta();
    spannedAngle = map.GetSpannedAngle();
    pointCount = map.GetPointsCount();
#ifdef OFFSET_POISSON
    Healpix_Correlation map2(_dT,sphereAngle);
    map2.SetPixValue(v1.size());
    vec_coord refPeebles1,refPeebles2,refHamilton,refLandy;
    map2.sweepAverageTwoPCorrelation(_dT,endTheta,refPeebles1,refPeebles2,refHamilton,
        refLandy,nTime,nk);
    typeIndex _len = refPeebles1.size();
    for(typeIndex _i = 0; _i < _len; _i++) {
        peebles1[_i].second -= refPeebles1[_i].second;
        peebles2[_i].second -= refPeebles2[_i].second;
        hamilton[_i].second -= refHamilton[_i].second;
        landy[_i].second -= refLandy[_i].second;
    }
#endif
}

void sphereTwoCorrelation::sweepAverageHealpix2PCF(myDouble _dT, myDouble endTheta, 
    spherePointsContainer &v1, myDouble thetaBegin, myDouble thetaEnd, vec_coord & peebles1, vec_coord & peebles2,
    vec_coord & hamilton, vec_coord & landy,
    std::vector<myDouble> & healpixMapInfo,  countType & pointCount,
    myDouble & spannedAngle, myDouble & dT, ofstream & os, int nTime, 
    countType nk)
{

    cout << "SOMETHING WRONG with this function, dont use it right now\n";
    exit(DONT_CALL_ME);

    Healpix_Correlation map(_dT, thetaEnd,thetaBegin);
    map.SetPixValue(v1);
    map.sweepAverageTwoPCorrelation(_dT,endTheta,peebles1,peebles2,
        hamilton, landy,nTime,nk,healpixMapInfo,os);
    dT = map.GetDeltaTheta();
    spannedAngle = map.GetSpannedAngle();
    pointCount = map.GetPointsCount();
#ifdef OFFSET_POISSON
    Healpix_Correlation map2(_dT,sphereAngle);
    map2.SetPixValue(v1.size());
    vec_coord refPeebles1,refPeebles2,refHamilton,refLandy;
    map2.sweepAverageTwoPCorrelation(_dT,endTheta,refPeebles1,refPeebles2,refHamilton,
        refLandy,nTime,nk);
    typeIndex _len = refPeebles1.size();
    for(typeIndex _i = 0; _i < _len; _i++) {
        peebles1[_i].second -= refPeebles1[_i].second;
        peebles2[_i].second -= refPeebles2[_i].second;
        hamilton[_i].second -= refHamilton[_i].second;
        landy[_i].second -= refLandy[_i].second;
    }
#endif

}    

void sphereTwoCorrelation::SetDeltaTheta(myDouble dT) {
    if(deltaTheta == dT) {
        return;
    }
    else {
        deltaTheta = dT;
        healpixMapSet = false;
    }
}

bool sphereTwoCorrelation::isHealpixSet() {
    return healpixMapSet;
}