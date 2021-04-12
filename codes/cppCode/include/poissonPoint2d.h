#ifndef __POISSON2D_H
#define __POISSON2D_H

#include <iostream>

// #include "myDefineConst.h"
#include <random>
#include <iterator>
#include "sphericalPoint.h"
#include <set>
#include <algorithm>


// (x,y,z) or (r, theta, phi)s
// typedef struct threeD{
//     // mode 1 : x,y,z 
//     // mode 0 : r, theta, phi
//     myDouble x;
//     myDouble y;
//     myDouble z;
//     int mode;
//     myDouble dist(const threeD & t) const {
//         if (mode == 1) 
//             return sqrt(square(x - t.x) + square(y - t.y) + square(z - t.z));
//         else {
//             myDouble sum = square(x * sin(y) * cos(z) - t.x * sin(t.y) * cos(t.z))
//                         + square(x * sin(y) * sin(z) - t.x * sin(t.y) * sin(t.z))
//                         + square(x * cos(y) - t.x * cos(t.y));
//             return sqrt(sum);

//         }
//     }
// }threeDCoord;

// typedef std::vector<threeDCoord> vec_3DCoord;



class PoissonPoint2d {
    private:
        // radius of the sphere
        myDouble radius;
        // number of points ( in Average)
        sizeType N;
        myDouble lambdaDefault;
        myDouble spannedAngle;
    public:
        // PoissonPoint2d(){}
        PoissonPoint2d(myDouble angle = M_PI_2, myDouble r = 1, sizeType num = POISSON_POINT_COUNT) : radius(r), N(num),spannedAngle(angle){ setLambdaDefault();}

        PoissonPoint2d(myDouble angle, myDouble lambda, myDouble r ): radius(r), lambdaDefault(lambda), spannedAngle(angle) { setN();}

        // return n poisson points process on a circle, with radius defined
        // by radius. mode 1: cartesian coordinate mode 2: polar coordinate
        // vec_coord rPoissonCircle(sizeType n = 100, myDouble rho = 1, int mode = 1);
        // // return n poisson points process on a sphere, with radius defined
        // // by radius. mode 1: cartesian coordinate mode 2: sphere coordinate
        void rPoissonSphere(vec_3DCoord & vec3d, sizeType n = 100, int mode = 1) = delete;
        // // return the poisson points process on a sphere with number N
        void rPoissonSphereDefault(vec_3DCoord & vec3d, int mode = 1);


        void rPoissonSphereDegreeLambdaDefault(vec_3DCoord & vec3d);
        // return n poisson points process on a sphere with center at the north pole
        // and spans a radian of theta degrees, using the projection method
        // vec_3DCoord rPoissonSphereDegreeProjection(myDouble theta,sizeType n = 100, int mode = 1) = delete;
        
        // return the poisson points process on a sphere with lambdaDefault and spannedAngle
        // vec_3DCoord rPoissonSphereDegreeLambdaDefault() = delete;

        // return a poisson points process on a sphere with center at north pole, 
        // spanned degree from thetaBegin to thetaEnd
        // given the density lambda
        // vec_3DCoord rPoissonSphereDegree(myDouble lambda0 = 1, myDouble thetaEnd = M_PI_2, myDouble thetaBegin = 0) = delete;
        // calculate a poisson points process on a sphere with center at north pole, spanned degree thetaBegin to thetaEnd
        // store the points in the vector vec3d
        void rPoissonSphereDegree(vec_3DCoord & vec3d, myDouble lambda0 = 1, myDouble thetaEnd = M_PI_2, myDouble thetaBegin = 0);

        void rPoissonSphereDegree(vec_3DCoord & vec3d, myDouble lambda0, myDouble thetaEnd, myDouble thetaBegin, myDouble phiEnd, myDouble phiBegin);

        // generate approx. n poisson points on a sphere with center at north pole,
        // spanned _thetaBegin to _thetaEnd, results stored in vec3d
        void rPoissonNSphereDegree(vec_3DCoord & vec3d, sizeType _M,
            myDouble _thetaEnd = M_PI_2, myDouble _thetaBegin = 0,
            myDouble _phiEnd = twopi, myDouble _phiBegin = 0.0);
        // return a vector of approx. n poisson points on a sphere with center at north pole,
        // spanned _thetaBegin to _thetaEnd, results stored in vec3d
        // vec_3DCoord rPoissonNSphereDegree(sizeType N, myDouble _thetaEnd = M_PI_2, myDouble _thetaBegin = 0) = delete;
        
         /**
         *  @brief generate matern process based on the poisson process with given intensity
         *  @param vecPosition, the vector that stores Spherical Point
         *  @param vecData, the vector that stores the data
         *  @param _lambda, the intensity
         *  @param _dthetaEnd, the end spannign angle of the cluster
         *  @param _dthetaBegin, the start spanning angle the cluster
         * */
        void MaternProcessLambda(const vec_3DCoord & vecPosition, vec_3DCoord & vecData,  myDouble _lambda, myDouble _dthetaEnd, myDouble _dthetaBegin = 0.0);

        /**
         *  @brief generate matern process based on the poisson process with given intensity, 
         *  and given number of cluster center
         *  @param vecPosition, the vector that stores the position of cluster center
         *  @param vecData, the vector that stores the data points
         *  @param _lambda, the intensity
         *  @param _M, the number of cluster centers, if larger than the given vector size, use the vector size
         *  @param _dthetaEnd, the end spannign angle of the cluster
         *  @param _dthetaBegin, the start spanning angle the cluster
         * */
        void MaternProcessChooseLambda(const vec_3DCoord & vec3d, vec_3DCoord & vecData, myDouble _lambda, sizeType _M, myDouble _dthetaEnd, myDouble _dthetaBegin = 0.0);


        /**
         *  @brief generate matern process based on the poisson process with given number of points in a cluster
         *  @param vecPosition, the vector that stores parent process
         *  @param vecData, the vector that stores the daughter process
         *  @param _N, the expected number within each cluster
         *  @param _thetaEnd, the degree span End
         *  @param _thetaBegin, the degree span Begin
         * */
        void MaternProcessNPoint(const vec_3DCoord & vecPosition, 
            vec_3DCoord & vecData, sizeType _M, 
            myDouble _thetaEnd,
            myDouble _thetaBegin = 0);

        /**
         *  @brief generate matern process based on the poisson process with given number of points in a cluster, 
         *  and given number of cluster center
         *  @param vec3d, the vector that stores parent process
         *  @param vecData, the vec that stores the daughter process
         *  @param _N, the expected number within each cluster
         *  @param _M, the number of cluster centers, if larger than the given vector size, use the vector size
         *  @param _thetaEnd, the degree span End
         *  @param _thetaBegin, the degree span Begin
         * */
        void MaternProcessChooseNPoint(const vec_3DCoord & vecPosition, 
            vec_3DCoord & vecData, sizeType _N, 
            sizeType _M, myDouble _thetaEnd,
            myDouble _thetaBegin = 0);
        /**
         *  @brief generate matern process based on the poisson process with given number of points in a cluster, 
         *  using parallel method
         *  @param vec3d, the vector that stores parent process
         *  @param vecData, the vector that stores the daughter process
         *  @param _N, the expected number within each cluster
         *  @param _dthetaEnd, the degree span End of the daughter process
         *  @param _dthetaBegin, the degree span Begin of the daughter process
         * */
        void MaternProcessNPointParallel(const vec_3DCoord & vec3d, vec_3DCoord & vecData,
            sizeType _M, myDouble _dthetaEnd, myDouble _dthetaBegin = 0);

        /**
         *  @brief generate matern process based on the poisson process with given number of points in a cluster, 
         *  and given number of cluster center using parallel method
         *  @param vecPosition, the vector that stores parent point process
         *  @param vecData, the generated point process
         *  @param vecIndex, the vector that stores the index of cluster center in vec3d,
         *  @param _N, the expected number within each cluster
         *  @param _M, the number of cluster centers, if larger than the given vector size, use the vector size
         *  @param _thetaEnd, the degree spanned End
         *  @param _thetaBegin, the degree spanned Begin
         * */
        void MaternProcessChooseNPointParallel(const vec_3DCoord & vecPosition, vec_3DCoord & vecData,
            sizeType _N, sizeType _M, myDouble _thetaEnd, myDouble _thetaBegin = 0);

        /**
         *  @brief classical method to generate matern process
         *  @param vecData, the vector that stores the generated point process
         *  @param lambdaParent, the intensity of the parent process
         *  @param lambdaDaughter, the intensity of the daughter process
         *  @param _spanThetaEnd, the spanned theta end of the total spanned matern process
         *  @param _spanThetaBegin, the spanned theta begin of the total matern process
         *  @param _dthetaEnd, the end angle of the spanning within each cluser
         *  @param _dthetaBegin, the begin angle of the spanning within each cluster
         * */
        void MaternProcess(vec_3DCoord & vecData,myDouble lambdaParent, myDouble lambdaDaughter, myDouble _spanThetaEnd, myDouble _spanThetaBegin,
            myDouble _dthetaEnd, myDouble dthetaBegin);
        /**
         *  @brief delete randomly chosen clusters 
         *  @param vec3d, the raw vector
         *  @param _M, the number of clusters to be deleted
         *  @param _theta, the spanned angle of each cluster to be deleted
         * */
        void deleteCluster(vec_3DCoord & vec3d, 
            sizeType _M, myDouble _theta);
        /**
         * @brief generate n (approximately) poisson points within theta range of the point p0
         * @param p0, reference to the spherical point
         * @param _N, expected number within each cluster
         * @param vec3d, the results are stored in vec3d
         * @param _thetaEnd, the degree span End
         * @param _thetaBegin, the degree span Begin
         * */

        void nPointPoisson(const SphericalPoint & p0, vec_3DCoord & vec3d, sizeType _N,  myDouble _thetaEnd, myDouble _thetaBegin);
        
        /**
         * @brief return n (approximately) poisson points within _thetaBegin to _thetaEnd of point p0
         * @param p0, reference to the spherical point
         * @param _N, expected number within each cluster
         * @param _thetaEnd, the degree span End
         * @param _thetaBegin, the degree span Begin
         * */ 
        vec_3DCoord nPointPoisson(const SphericalPoint & p0, 
                sizeType _n, myDouble _thetaEnd, myDouble _thetaBegin = 0) = delete;
        /**
         * @brief generate n (approximately) poisson points within theta range of the point p0
         * @param p0, reference to the spherical point
         * @param _lambda, the intensity
         * @param _theta, the degree spanned
         * @param vec3d, the results are stored in vec3d
         * */
        /**
         * @brief generate poisson points within theta range for the points with index [_start, end)
         * @param _start, start of the index
         * @param _end, end of the index (not included)
         * @param _N, expected number within each cluster
         * @param vecPosition, the position point
         * @param vecData, the vectore that stores the points
         * @param _thetaEnd, the degree span End
         * @param _thetaBegin, the degree spanBegin,
         * 
         * */
        void indexParallelPoissonNPoint(typeIndex _start, 
                typeIndex _end,  sizeType _M,
                const vec_3DCoord & vecPosition, vec_3DCoord & vecData, myDouble _thetaEnd,myDouble _thetaBegin = 0);
        /**
         * @brief generate poisson points within theta range [_thetaBegin, _thetaEnd] for the points 
         * with index [vecIndex[_start], vecInde[end])
         * @param _start, start of the index
         * @param _end, end of the index (not included)
         * @param vecIndex, the vectors that stores the index
         * @param _N, expected number within each cluster
         * @param vecPosition, the vector that stores the position
         * @param vecData, the vector that stores the daughter point process
         * @param _theta, the degree span End
         * @param _thetaBegin, the degree span Begin
         * */
        void vecIndexParallelPoissonNPoint(typeIndex _start,
            typeIndex _end, const vec_Index & vecIndex, sizeType _M,
            const vec_3DCoord & vecPosition, vec_3DCoord & vecData,
            myDouble _thetaEnd, 
            myDouble _thetaBegin = 0);
            
        void lambdaPointPoisson(const SphericalPoint & p0, myDouble _lambda, vec_3DCoord & vec3d, myDouble _thetaEnd, myDouble _thetaBegin = 0.0);
        /**
         * @brief return n (approximately) poisson points within theta range of the point p0
         * @param p0, reference to the spherical point
         * @param _lambda, the intensity
         * @param _theta, the degree spanned
         * */
        vec_3DCoord lambdaPointPoisson(SphericalPoint & p0, myDouble _lambda, myDouble _theta) = delete;

        // return n binomial points process on a sphere with center at the north pole, with radius R0 (parameter given to initialize the class)
        // and spans a radian of from thetaBegin to thetaEnd, the returned are (theta, phi)
        vec_3DCoord rBinomialSphereDegreeN(sizeType n, myDouble theta0 = M_PI_2, myDouble thetaBegin = 0) = delete;
        // calculate n binomial points process on a sphere with center at the north pole, with radius R0 (parameter given to initialize the class)
        // and spans a radian from thetaBegin to thetaEnd, the returned are (theta, phi)
        // the results are stored in vec3d
        void rBinomialSphereDegreeN(vec_3DCoord & vec3d, sizeType n, myDouble theta0 = M_PI_2, myDouble thetaBegin = 0, myDouble phiEnd = M_PI * 2, myDouble phiBegin = 0.0);

        // return n binomial points process on the unitsphere with center at the north pole
        // and spans a radian of theta0 degrees, the returned are (theta, phi)
        vec_3DCoord rBinomialUnitSphereDegreeN(sizeType n,myDouble theta0 = M_PI_2) = delete;

        // return the radius
        myDouble getRadius() const { return radius;}

        // set the Radius
        void setRadius(myDouble R) { radius = R;}

        // set the Number 
        void setNumber(sizeType num) { N = num;}

        sizeType getNumber() const { return N;}

        // given the lambdaDefault and spannedAngle, calculate the N
        void setN();
        // given the spannedAngle and average number of points N, caluclate the density lambdaDefault
        void setLambdaDefault();

        // given N, theta, calculate lambda, so that poisson points
        // with density lambda will generate appr. N points within the thetaBegin to thetaEnd
        // spherical ring
        myDouble calculateLambda(sizeType n, myDouble R, myDouble thetaEnd, myDouble thetaBeign = 0, myDouble phiEnd = twopi, myDouble phiBegin = 0.0);

        /**
         *  @brief given _N, generate _M random distinct indexes from [0,_N)
         *  @param _N, the range 
         *  @param vecIndex, the vector that stores the index
         *  @param vecPosition, the vector that stores the position
         * */
        void GenerateNIndices(sizeType _N, sizeType _M,
              myDouble smallAngle, vec_Index & _vecIndex, const vec_3DCoord & vecPosition);
};
#endif