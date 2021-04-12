

#ifndef HEALPIX_CORRELATION_H
#define HEALPIX_CORRELATION_H


#include <healpix_map.h>
#include <alm.h>
#include <alm_healpix_tools.h>
// #include "myDefineConst.h"
#include "sphericalPoint.h"
#include "binsGenerator.h"
#include <cstdlib>
#include <ctime>
#include <string>
#include <mutex>

#ifdef HEALPIX_DEBUG
#ifndef IOSTREAM_H
#define IOSTREAM_H

#include <iostream>
using std::cout;
using std::endl;

#endif
#endif

#ifdef HEALPIX_PARALLEL_DEBUG
#ifndef IOSTREAM_H
#define IOSTREAM_H

#include <iostream>
using std::cout;
using std::endl;

#endif
#endif
using std::srand;
using std::rand;
using std::time;
using std::vector;
using std::thread;
using std::string;

typedef std::vector<SphericalPoint> spContainer;
typedef std::vector<myDouble> dataVecType;
class Healpix_Correlation : public Healpix_Map<myDouble> {
    protected:
        myDouble spannedAngle;
        int topRing;
        int bottomRing;
        // the last pixel of the map
        int topPixel;
        int bottomPixel;
        bool boundSet;
        myDouble deltaTheta;
        countType pointsCount;
        myDouble hThetaBeg;
        myDouble hThetaEnd;
        myDouble hPhiBeg;
        myDouble hPhiEnd;
        string _logDir;
        // the number of times correlation functins are called
        static int correlationCount;
    public:
        Healpix_Correlation(): Healpix_Map(), spannedAngle(0.0),boundSet(false),pointsCount(0),hThetaBeg(0.0),hThetaEnd(0.0),hPhiBeg(0.0),hPhiEnd(0.0), _logDir("./logs/"){ }

        Healpix_Correlation(int order, Healpix_Ordering_Scheme scheme, 
            const myDouble angle = M_PI_2, const myDouble angleBegin = 0.0,
            const myDouble phiE = twopi, const myDouble phiB = 0.0, const string& logDir = "./logs/") 
            : Healpix_Map(order, scheme), spannedAngle(angle),boundSet(false),pointsCount(0),hThetaBeg(angleBegin),hThetaEnd(angle),
            hPhiBeg(phiB), hPhiEnd(phiE), _logDir("./logs/") { update();}
        /*! Constructs a map with a given \a nside and the ordering
        scheme \a scheme. */
        Healpix_Correlation(int nside, Healpix_Ordering_Scheme scheme, 
            const nside_dummy nd, 
            const myDouble angle = M_PI_2,
            const myDouble angleBegin = 0.0,
            const myDouble phiE = twopi,
            const myDouble phiB = 0.0,
            const string& logDir = "./logs/")
        : Healpix_Map(nside, scheme, nd), spannedAngle(angle),boundSet(false),pointsCount(0),
            hThetaBeg(angleBegin),hThetaEnd(angle),
            hPhiBeg(phiB), hPhiEnd(phiE), _logDir(logDir){ update(); }
        /*! Constructs a map from the contents of \a data and sets the ordering
            scheme to \a Scheme. The size of \a data must be a valid HEALPix
            map size. */
        Healpix_Correlation(const arr<myDouble> &data, Healpix_Ordering_Scheme scheme, 
            const myDouble angle = M_PI_2, const myDouble angleBeg = 0.0,
            const myDouble phiE = twopi, const myDouble phiB = 0.0, const string& logDir = "./logs/")
            : Healpix_Map(data,scheme), spannedAngle(angle), boundSet(false),pointsCount(0),
            hThetaBeg(angleBeg), hThetaEnd(angle), hPhiBeg(phiB), hPhiEnd(phiE), _logDir(logDir){ update(); }
        // given the deltaTheta, automatically modify the binsize
        Healpix_Correlation(const myDouble & dTheta, 
            const myDouble angle = M_PI_2,
            const myDouble angleBegin = 0.0, 
            const myDouble phiE = twopi, 
            const myDouble phiB = 0.0, const string& logDir = "./logs/");
        #ifdef signalMode
        Healpix_Correlation(myDouble theta, myDouble angle = M_PI_2);
        #endif
        // return the rangeset of pixels that lie within the theta-disc (including the boundary) of the given pixel
        void query_disc_pixel_Internal(int pixel, const myDouble & theta,
            rangeset<int> &pixset) const;
        // return the rangeset of pixels that lie within the theta-disc (including the boundary) of the given pixel
        // given the reference hp 
        void query_disc_pixel_Internal(int pixel, const myDouble & theta,
            rangeset<int> &pixset, const Healpix_Correlation & hp) const;
        // return the vector of pixels that lie within the theta-disc
        void query_disc_pixel_Internal(int pixel, const myDouble & theta,
            std::vector<int> & pixset) const;
        // return the vector of pixels that lie within the theta-disc,
        // given the reference to a healpix map hp
        void query_disc_pixel_Internal(int pixel, const myDouble & theta,
            std::vector<int> & pixset, const Healpix_Correlation & hp) const;

        // the vector version of query_disc_pixel_Internal
        // will ignore the pixel ignore
        void query_disc_pixel_InternalVec(int pixel, const myDouble & theta, 
            std::vector<int> & pixset,
            const Healpix_Correlation & hp, int ignore = -1) const;
        // given pixel, theta, hp, calculate the number of points lie in the theta-disc of p directly
        // ignore the pixel ignore
        countType query_disc_count(int pixel, const myDouble & theta,
            const Healpix_Correlation & hp, int ignore = -1) const;
        // calculate the number of points lie in the ring of angle (thetaSmall, thetaLarge) of a pixel
        void query_disc_pixel_NonCumuInternalVec(int pixel, const myDouble & thetaLarge, 
            const myDouble & thetaSmall,
            std::vector<int> & pixset,
            const Healpix_Correlation & hp, std::string fileName = MY_DEFAULT_STRING) const;
        // calculate the number of points lie in the ring of angle (thetaSmall, thetaLarge) of a pixel
        // the pixel is on the (*this) map, while the reference hp, denotes that we consider the neighbor of pixel on that hp map
        // the final neighbor count should be the returned pairs #{(pixel, x), x \in hp and hp[x] != 0} * (*this)[pixel]
        countType query_disc_pixel_NonCumuInternalCount(int pixel, const myDouble & thetaLarge, 
            const myDouble & thetaSmall,
            const Healpix_Correlation & hp) const;            
        // given bins, find the neighbors within each bin for a given pixel, the result stored in to neighborCount
        // assumes neighborCount already has required size
        // besides, if set MULTIPLY_PIXEL_VALUE, will multiply the endresult
        // with (*this)[pixel]
        void query_disc_pixel_Bins(int pixel, const std::vector<myDouble>& bins,
            std::vector<corrType> & neighborCount,
            const Healpix_Correlation & hp);

        void query_disc_pixel_Bins(int pixel, const std::vector<myDouble>& bins,
            std::vector<corrType> & neighborCount,
            const Healpix_Correlation & hp,
            const Healpix_Correlation & origin,
            vector<std::string> fileNames,
            vector<myDouble> thetas);

        void update();
        // update the boundPixels;
        void updateBoundPixel();
        void updateDeltaTheta();

        // set all pixValues to value, sett all other values to zero
        void resetPixValue(int value, int _threadNum);

        // scale all pixValues 
        void scalePixValue(myDouble scale, int _threadNum = MAX_THREAD);
        // set the spanned angle
        void SetAngle(const myDouble & angle);
        // set deltaTheta;
        void SetDeltaTheta(const myDouble & deltaTheta, const myDouble & angle = M_PI_2);
        // given the data point, set the pixel value
        void SetPixValue(const spContainer & vc,int _threadNum = MAX_THREAD);
        // for the random case given total point count, set the value
        // averageTime is nTime
        void SetPixValue(countType pointCount,int nTime = 1, int _threadNum = MAX_THREAD);
        /**
         *  @brief given the experiment data for the signal strength, set the pix value
         *  @param _thetaVec, the vector of theta
         *  @param _phiVec, the vector of phi
         *  @param _dataVec, the signal strength arranged in the order: first fix theta, then sweep phi
         *  @param _threshold, the threshold value for signal strength
         *  @param mode, the mode to convert signal strength to value in healpix map
         * */
        void SetPixValue(const std::vector<myDouble> & _thetaVec,
            const std::vector<myDouble> & _phiVec,
            const std::vector<myDouble> & _dataVec,
            myDouble & _threshold,
            int mode = 1, int _threadNum = MAX_THREAD);

        /**
         * @brief given the experiment data, read in and set pixel value, 
         *  store the data points coord in a vector 
         * 
         * @param _thetaVec, the vector that stores the theta values
         * @param _phiVec, the vector that stores the phi values
         * @param _dataVec, the vector that stores the data values
         * @param vecPoint, the vector that stores the coordinate of points
         * @param _threshold, the threshold value
         * @param _mode: int specify the strategy to convert data to a pixel value
         * if set to 1: simple binary case 
         *           2: still treat each data point as binary, but will accumulate the point
         *           3: use the mW power value if the signal value is above the threshold
         *           4: use the W power value if the signal value is above the threshold
         * */   
        void SetPixValue(const std::vector<myDouble> & _thetaVec,
            const std::vector<myDouble> & _phiVec,
            const std::vector<myDouble> & _dataVec,
            spContainer & vecPoint,
            myDouble & _threshold, 
            int mode = 1, int _threadNum = MAX_THREAD);       

        // in the NeighborVec function and alike, modify the deltaTheta and pixValue
        void ModifySize(myDouble stepTheta);

        // get the pointsCount
        countType GetPointsCount() const { return pointsCount;}
        int GetBottomPixel() const { return bottomPixel;}
        int GetTopPixel() const {return topPixel;}
        // get the DeltaTheta
        myDouble GetDeltaTheta() const { return deltaTheta; }
        // get spannedAngle
        myDouble GetSpannedAngle() const { return spannedAngle; }

        // number of points inside (including the boundary) the theta disc.
        countType Neighbors(const myDouble & theta);
        // number of neighbors of a given pixel in the theta disc
        countType Neighbors(int pixel, const myDouble & theta);

        countType Neighbors(int pixel, const myDouble & theta, std::string filename);

        // number of points inside (including the boundary) the theta disc,
        // using parallel programming, the thread function uses a vector
        countType NeighborsParallelVector(const myDouble & theta,int _threadNum = MAX_THREAD);

        countType NeighborsParallel(const myDouble & theta,int _threadNum = MAX_THREAD);
        
        // //number of neighbors of a given pixel, within the theta disc (including boundary)
        // countType NeighborsPixelParallel(int pixel, const myDouble & theta, int _threadNum = MAX_THREAD);

        // number of points inside (including the boundary) the theta disc
        // with respect to another map ( the DR)
        countType crossNeighbors(const myDouble & theta, Healpix_Correlation & hp);

        countType crossNeighborsParallel(const myDouble & theta, Healpix_Correlation & hp,int _threadNum = MAX_THREAD);
        // given the bins,(uniform), calculate the neighbors within each interval
        // e.g  {0,0.1,0.2,0.3} the the returned is {x,y,z} where x: neighbors in (0,0.1],
        // y : neighbors in (0.1,0.2], z : neighbors in (0.2,0.3]
        std::vector<corrType> NeighborVec(std::vector<myDouble>& bins,int _threadNum = MAX_THREAD);
        // Factor used in averageRR for example
        void NeighborVec(std::vector<myDouble> & bins, std::vector<corrType> & neighborCount, int _threadNum = MAX_THREAD);
        
        // given one pixel, find its neighbors with given bin size, using cumulative method
        void NeighborPixelCumVec(std::vector<myDouble> & bins, int pixel, 
            std::vector<corrType> & neighborCount, 
            std::vector<corrType> & cumulatedCount,
            int threadNum = MAX_THREAD);

        void NeighborPixelCumVec(std::vector<myDouble> & bins, int pixel, 
            std::vector<corrType> & neighborCount, 
            std::vector<corrType> & cumulatedCount,
            vector<std::string> fileNames,
            vector<myDouble> thetas,
            int threadNum = MAX_THREAD);

         // given one pixel, find its neighbors with given bin size, using non cumulative method
        void NeighborPixelNonCumVec(std::vector<myDouble> & bins, int pixel, 
            std::vector<corrType> & neighborCount, int threadNum = MAX_THREAD);

        void NeighborPixelNonCumVec(std::vector<myDouble> & bins, int pixel, 
            std::vector<corrType> & neighborCount, 
            vector<std::string> fileNames,
            vector<myDouble> thetas,            
            int threadNum = MAX_THREAD);
       
        // Non cumulative version to calculate neighbor
        void NeighborNonCumVec(std::vector<myDouble> & bins, 
            std::vector<corrType> & neighborCount,
            const Healpix_Correlation & hp, int Factor, int _threadNum = MAX_THREAD);

        // given the bins,(uniform) calculate the cross neighbors
        std::vector<corrType> crossNeighborVec(std::vector<myDouble>& bins, Healpix_Correlation & hp,int _threadNum = MAX_THREAD);
        void crossNeighborVec(std::vector<myDouble>& bins, std::vector<corrType> & neighborCount, Healpix_Correlation & hp,int _threadNum = MAX_THREAD);

        // averaged version of crossNeighborVec, nk points for randomMap, nTime
        std::vector<corrType> averageCrossNeighborVec(std::vector<myDouble>& bins, countType nk, int nTime,int _threadNum = MAX_THREAD);
        void averageCrossNeighborVec(std::vector<myDouble>& bins, std::vector<corrType> & neighborCount, countType nk, int nTime,int _threadNum = MAX_THREAD);

        // averaged version for calculating RR, nk points, nTime times
        std::vector<corrType> averageRR(std::vector<myDouble>& bins, countType nk, int nTime,int _threadNum = MAX_THREAD);
        void averageRR(std::vector<myDouble>& bins, std::vector<corrType>& RR, countType nk, int nTime,int _threadNum = MAX_THREAD);

        void compactAverageRDR(std::vector<myDouble>& bins, 
            std::vector<corrType>&DD, std::vector<corrType> &DR,
            std::vector<corrType>&RR, countType nk,
            int nTime,int _threadNum = MAX_THREAD);
        // given the deltaTheta, calculate the order
        int deltaTheta2Order(const myDouble & deltaTheta);

        // judge whether two maps are of the same structure?
        bool isStructureIsomorphic(const Healpix_Correlation & hp);

        // help function for parallel programming
        // calculate the number of pairs lie in the the theta disc, center starting from begin to end [begin, end), result
        // stored in vec
        void NeighbourIntervalCountVector(typeIndex begin, typeIndex end, const myDouble & theta, std::vector<corrType> & vec);
        // help function for parallel programming
        // calculate the number of pairs lie in the the theta disc, center starting from begin to end [begin, end), result
        // stored in variable pairCount
        void NeighbourIntervalCount(typeIndex begin, typeIndex end, const myDouble & theta, countType & pairCount);
        void crossNeighbourIntervalCount(typeIndex begin, typeIndex end, const myDouble & theta, const Healpix_Correlation & hp, 
            countType & pairCount);

        void crossNeighborNonCumIntervalCount(typeIndex begin, typeIndex end, 
            const std::vector<myDouble> & bins,
            std::vector<corrType> & neighborCount,
            const Healpix_Correlation & hp);
        // calculate the neighbor of hp, use Factor to update neighborCount
        // neighborCount[index] = (neighborCount[index] * (Factor) + newCouunt)/(Factor + 1),
        // besides, the neighborCount[index] after the calculation is 2 times the actual value
        // so if Factor stands for the final factor, make sure to >> 1, to retrieve the desired result
        void NeighbourNonCumIntervalCount(typeIndex begin, typeIndex end, 
            const std::vector<myDouble> & bins, 
            std::vector<corrType> & neighborCount, 
            const Healpix_Correlation & hp, int Factor = 0);

        // set the value of pixel in [begin, end) as value
        void parallelSetPixValue(typeIndex begin, typeIndex end, countType & nonZeroCount, int value = 0);
        // scale the value of pixel in [begin,end)
        void parallelScalePixValue(typeIndex begin, typeIndex end, myDouble & scale);
        // calculate the 2pCF
        void twoPCorrelation(myDouble dT, myDouble endTheta, vec_coord & vec, Healpix_Correlation & hp, countType nk = HEALPIX_RANDOM_POINT_COUNT,
            int nTime = AVERAGE_TIME, int mode = 1,int _threadNum = MAX_THREAD);  
        // calculate the 2pCF and output the healpix map info
        void twoPCorrelation(myDouble dT, myDouble endTheta, vec_coord & vC, Healpix_Correlation & hp, 
            std::vector<myDouble> & healpixMapInfo, std::vector<myDouble> & randomHealpixInfo, int mode = 1,int _threadNum = MAX_THREAD);
        
        // // for the rangeset_test function
        // // void twoPCorrelationRangeSet(myDouble dT, myDouble endTheta,
        //     vec_coord rawVC, vec_coord rangeModVC, vec_coord vecVC,
        //     Healpix_Correlation & hp, int mode);
        // using an average version, i.e DR will be calculated several times and find the average, times give by nTime, with nk points for the randomMap
        void averageTwoPCorrelation(myDouble dT, myDouble endTheta,
             vec_coord & vC,int nTime, countType nk, 
             int mode = 1,int _threadNum = MAX_THREAD);
        void averageTwoPCorrelation(myDouble dT, myDouble endTheta, 
            vec_coord & vC,int nTime, countType nk, 
            std::vector<myDouble> & healpixMapInfo,int mode = 1,
            int _threadNum = MAX_THREAD);

        // sweep all estimators
        void sweepAverageTwoPCorrelation(myDouble dT, myDouble endTheta,
            vec_coord & peebles1VC, vec_coord & peebles2VC, 
            vec_coord & hamiltonVC, vec_coord & landyVC,
            int nTime, countType nk,
            std::vector<myDouble> & healpixMapInfo,
            int _threadNum = MAX_THREAD);

        // sweep all estimators and output the DD,RR,DR to the os
        void sweepAverageTwoPCorrelation(myDouble dT, myDouble endTheta,
            vec_coord & peebles1VC, vec_coord & peebles2VC, 
            vec_coord & hamiltonVC, vec_coord & landyVC,
            int nTime, countType nk,
            std::vector<myDouble> & healpixMapInfo, ostream & os,
            int _threadNum = MAX_THREAD);

        // the same as above, no info for healpixMap
        void sweepAverageTwoPCorrelation(myDouble dT, myDouble endTheta,
            vec_coord & peebles1VC, vec_coord & peebles2VC, 
            vec_coord & hamiltonVC, vec_coord & landyVC,
            int nTime, countType nk,
            int _threadNum = MAX_THREAD);
        // sweep the estimators, save the DD,DR,DR, specify the log sub Folder

        /**
         * @param logFolder: const std::string &,  the string specify the log sub folder
         * */
        void sweepAverageTwoPCorrelation(myDouble dT, myDouble endTheta,
            vec_coord & peebles1VC, vec_coord & peebles2VC, 
            vec_coord & hamiltonVC, vec_coord & landyVC,
            int nTime, countType nk, std::vector<corrType> & DD,
            std::vector<corrType> & RR,
            std::vector<corrType> & DR, 
            std::vector<myDouble> & bins,
            int _threadNum = MAX_THREAD,
            const std::string & logFolder = "");    
        // calculate the 2pCF for given point process, which are stored in vecData
        // note that the countType is scale not nk
        void pointProcessSweep2PCF(myDouble dT, myDouble endTheta,
            const vec_3DCoord & vecData,
            vec_coord & peebles1VC, vec_coord & peebles2VC, 
            vec_coord & hamiltonVC, vec_coord & landyVC,
            int nTime, countType scale, std::vector<myDouble> & healpixMapInfo, 
            ostream & os,
            int _threadNum = MAX_THREAD);

        void pointProcessPrintPoissonSweep2PCF(myDouble dT, myDouble endTheta,
            const vec_3DCoord & vecData,
            vec_coord & peebles1VC, vec_coord & peebles2VC, 
            vec_coord & hamiltonVC, vec_coord & landyVC,
            vec_coord & poissonPeebles1VC, vec_coord & poissonPeebles2VC,
            vec_coord & poissonHamiltonVC, vec_coord & poissonLandyVC,
            int nTime, countType scale, std::vector<myDouble> & healpixMapInfo, 
            ostream & os,
            int _threadNum = MAX_THREAD);
        // output the healpix map info, also generate a random map with the same number of points, out put to randomMapinfo
        void outMapAndRandom(std::vector<myDouble> & healpixMapInfo, std::vector<myDouble> & randomMapInfo);
        /**
         *  @brief convert the signal strength to value on the healpixMap
         *  @param _strength the signal strength(in dB)
         *  @param _threshold the threshold value
         *  @param _mode the method to be used
         * */
        myDouble signalStrength2Val(myDouble _strength, myDouble _threshold,
            int mode = 1);


        // some test function

        // test if all pixValues are zero?
        bool isAllZero() const;

        // publicise the Healpix_Map function
        // return the ring num above the z-coord z
        int indexRingAbove(double z) const;
        // return the pointing
        pointing publicPix2Ang(int pixel) const;
        void publicGetRingInfoSmall(int ring, int &startpix, int &ringpix, bool &shifted) const;
        void publicPix2ZPhi(int pix, double &z, double &phi) const;
        // convert theta to ring num
        int publicTheta2Ring(const myDouble & theta) const;


        // append the interval [begin,end) to the rangeset
        // degenerate, should not be used, is not accurate as VecAppendNonZeroPoint
        void appendNonZeroRange(rangeset<int> & pixset,
            const Healpix_Correlation & hp, typeIndex begin, typeIndex end, int ignore = -1) const;
        // the equivalent version of appendNonZeroRange, here use vector, so directly add the point in the interval
        // do not add the pixel ignore to the pixVec, used, for example, in the non cum version of query disc
        void VecAppendNonZeroPoint(vector<int> & pixVec, 
            const Healpix_Correlation & hp, typeIndex begin, typeIndex end, int ignore = -1) const;
        // calculate the number of non zero pix-value from begin to end [begin, end)
        countType IntervallNonZeroPoint(const Healpix_Correlation & hp, 
            typeIndex begin, typeIndex end, int ignore = -1) const;

        /**
         * @brief calculate the alm of the spherical harmonics
         * /
         * Alm< xcomplex< myDouble > > & alm,
         * */
        
        // void sphericalHarmonicALM(
        //     const arr< double > & 	weight,
        //     bool add_alm = false, bool useIter = true, int numIter = NUM_ITER
        //     );
        friend ostream & operator<<(ostream & os, const Healpix_Correlation & hp);
        
        // out put the ring info
        void outRingInfo(ostream & os) const;
        // print out the essential information about the map, pointscount, deltaTheta, forexample
        void essentialInfoOut(ostream & os) const;
        /**
         * print the essenatial Info, specify the logFolder
         * */
        string essentialInfoOut(const string & logFolder = "") const;

};


template<typename T>
const T max(T v1, T v2) {
    return (v1 > v2) ? v1 : v2;
}

template<typename T>
const T min(T v1, T v2) {
    return (v1 < v2) ? v1 : v2;
}


#endif