// #include "sphericalPoint.h" // already in poissonPoint2d
// #include "myDefineConst.h"
// #include "poissonPoint2d.h"
#include "findPair.h"
// #include <vector>
#include "healpixCorrelation.h"
#include <map>

using std::vector;

typedef std::vector<SphericalPoint> spherePointsContainer;
typedef std::vector<SphericalPoint>::const_iterator spCont_ConstIter;
// map whose key is the pixel number and the value is the number of points,
// that lie in that pixel
typedef std::map<int, countType> pixelCount;
#ifndef __CLASS_sphereTwoCorrelation_H
#define __CLASS_shpereTwoCorrelation_H

class sphereTwoCorrelation {
    private:
        spherePointsContainer sPContainer;
        // the sphere solid angle, that all the points span
        myDouble sphereAngle;
        countType pointCount;
        myDouble radius;
        FindPair fpSpherePoint;
        FindPair fpPoisson;
        spherePointsContainer & spCont;

        // for the healpix Method
        Healpix_Correlation hpMap;
        Healpix_Correlation randomHpMap;
        // flag indicating if the Healpix_Correlation hpMap initialized
        bool healpixMapSet;
        myDouble deltaTheta;
    public:
        // sphereTwoCorrelation() { }
        // sphereTwoCorrelation(spherePointsContainer & v1, myDouble Angle = M_PI/2, countType num = POISSON_POINT_COUNT);
        sphereTwoCorrelation(spherePointsContainer & v1, PoissonPoint2d & pd, myDouble Angle = M_PI/2, bool setAngle = false);
        // calculate the sphere solid angle, that all points span
        void setSphereAngle();
        // set the radius
        void setRadius();
        void pushPoint(SphericalPoint & p);
        // count the number of pairs in the range (theta, theta + deltaTheta)
        long DDTheta(myDouble theta, myDouble deltaTheta = DELTA_THETA) const;
        // the twoPoint correlation function at angle theta
        myDouble twoPointCorrelation(myDouble theta, myDouble deltaTheta = DELTA_THETA, int mode = 0);

        // return the number of points
        long getSize() const { return (int) pointCount; }
        // isEmpty
        bool isEmpty() const { return pointCount == 0;}
        // for calculating the boundary
        myDouble maxTheta() const;
        myDouble minTheta() const;
        myDouble maxPhi() const;
        myDouble minPhi() const;

        // the healpix method
        // return the number of neighbors who lie within the disc of theta
        countType Neighbors(myDouble theta,spherePointsContainer & v1, int order);
        // // return the number of neighbors at disc theta, given healpix map
        // countType Neighbors(myDouble theta, )
        // given the deltaTheta, and the pointsContainer v1 and the vector vC, calculate the 
        // the 2pCF using the healpix method, results are in vC
        void healpixTwoPointCorr(myDouble deltaTheta, myDouble endTheta, spherePointsContainer & v1, vec_coord & vC, int mode = 1);
        // for the plotting, return also information about the healpixMap
        // the healpixMapInfo has the content 1. order 2.pixel value
        //  void healpixTwoPointCorr(myDouble deltaTheta, myDouble endTheta, spherePointsContainer & v1, vec_coord & vC, 
        //  std::vector<countType> healpixMapInfo, std::vector<countType> randomHealpixInfo, int mode = 1);

         void healpixTwoPointCorr(myDouble deltaTheta, myDouble endTheta, spherePointsContainer & v1, vec_coord & vC, 
         std::vector<myDouble> & healpixMapInfo, std::vector<myDouble> & randomHealpixInfo, int mode = 1);

        // calculate the averaged 2pCF, based on healpix Map without setting the healpix map value.
        void averageHealpixTwoPointCorrNoSet(myDouble deltaTheta, myDouble endTheta, 
            vec_coord & vC, Healpix_Correlation & Hmap, 
            int nTime = AVERAGE_TIME, countType nk = HEALPIX_RANDOM_POINT_COUNT, 
            int mode = 1);
        // calculate averaged 2pCF given heapix Map and set corresponding variable for the info of the map.
        void averageHealpixTwoPointCorrNoSet(myDouble deltaTheta, myDouble endTheta, vec_coord & vC,
            Healpix_Correlation & Hmap, std::vector<myDouble> & healpixMapInfo,  
            countType & pointCount,
            myDouble & spannedAngle, myDouble & dT, 
            int nTime = AVERAGE_TIME, 
            countType nk = HEALPIX_RANDOM_POINT_COUNT,
            int mode = 1);
        // calculate the healpix version of 2pCF based on averaged DR, and RR, with random points nk and nTime
        void averageHealpixTwoPointCorr(myDouble deltaTheta, myDouble endTheta,
         spherePointsContainer & v1, vec_coord & vC, 
         int nTime = AVERAGE_TIME, countType nk = HEALPIX_RANDOM_POINT_COUNT,
         int mode = 1);
        // for the plotting, additional healpixMap info
        void averageHealpixTwoPointCorr(myDouble deltaTheta, myDouble endTheta, 
            spherePointsContainer & v1, vec_coord & vC, 
            std::vector<myDouble> & healpixMapInfo, 
            int nTime = AVERAGE_TIME, countType nk = HEALPIX_RANDOM_POINT_COUNT,
            int mode = 1);
        // for the plotting, additional healpixMap info, also set the calculatioin info,
        // e.g pointCount, spannedAngle, deltaTheta
        void averageHealpixTwoPointCorr(myDouble deltaTheta, myDouble endTheta, 
            spherePointsContainer & v1, vec_coord & vC,
            std::vector<myDouble> & healpixMapInfo,  countType & pointCount,
            myDouble & spannedAngle, myDouble & dT, int nTime = AVERAGE_TIME, 
            countType nk = HEALPIX_RANDOM_POINT_COUNT,int mode = 1);

        // sweep all estimators based on the healpix method.
        void sweepAverageHealpix2PCF(myDouble deltaTheta, myDouble endTheta, 
            spherePointsContainer &v1, vec_coord & peebles1, vec_coord & peebles2,
            vec_coord & hamilton, vec_coord & landy,int nTime = AVERAGE_TIME, 
            countType nk = HEALPIX_RANDOM_POINT_COUNT);
        // sweep all estimators based on the healpix method. 
        // store also calculation info
        void sweepAverageHealpix2PCF(myDouble deltaTheta, myDouble endTheta, 
            spherePointsContainer &v1, vec_coord & peebles1, vec_coord & peebles2,
            vec_coord & hamilton, vec_coord & landy,
            std::vector<myDouble> & healpixMapInfo,  countType & pointCount,
            myDouble & spannedAngle, myDouble & dT, ofstream & os, int nTime = AVERAGE_TIME, 
            countType nk = HEALPIX_RANDOM_POINT_COUNT);

        // will give the range thetaBegin and thetaEnd
        void sweepAverageHealpix2PCF(myDouble deltaTheta, myDouble endTheta, 
            spherePointsContainer &v1, myDouble thetaBegin, myDouble thetaEnd, vec_coord & peebles1, vec_coord & peebles2,
            vec_coord & hamilton, vec_coord & landy,
            std::vector<myDouble> & healpixMapInfo,  countType & pointCount,
            myDouble & spannedAngle, myDouble & dT, ofstream & os, int nTime = AVERAGE_TIME, 
            countType nk = HEALPIX_RANDOM_POINT_COUNT);        
        // modify the deltaTheta
        void SetDeltaTheta(myDouble dT);
        // judge whether the hpMap is set
        bool isHealpixSet();
        
};
#endif