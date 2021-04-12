#ifndef __CLASS_SPHERICAL_H
#define __CLASS_SPHERICAL_H

#include <iostream>
#include <vector>
// #include <cmath>
#include "myDefineConst.h"
#include "spdlog/spdlog.h"
#include "spdlog/sinks/basic_file_sink-inl.h"
#include "spdlog/sinks/stdout_color_sinks.h"
#include "spdlog/fmt/ostr.h"
#include <fstream>
using std::endl;
using std::cin;
using std::cout;
using std::ifstream;
using std::istream;
using std::ostream;
using std::ofstream;

using std::cerr;
using std::cout;
using std::cin;
using std::endl;
using std::ostream;
using std::vector;

class SphericalPoint;

typedef std::pair<myDouble,myDouble> twoDCoord;
typedef std::vector<twoDCoord> vec_coord;
typedef std::vector<SphericalPoint> vec_3DCoord;
typedef vec_3DCoord::const_iterator vec3DConstIter;
typedef vec_3DCoord::iterator vec3DIter;

typedef std::vector<typeIndex> vec_Index;


class SphericalPoint {
    private:
        // the three parameters for the spherical coordianate
        myDouble radius;
        myDouble theta;
        myDouble phi;
        // mode, mode = 0, for spherical coordinate, mode  = 1 for cartesian
        int mode;
        myDouble xCoord;
        myDouble yCoord;
        myDouble zCoord;

    public:
        SphericalPoint(){}
        SphericalPoint(myDouble r, myDouble t, myDouble p, int mode = 0);
        void sphericalToCartesian();
        void cartesianToSpherical();
        myDouble sphericalAngleDist(const SphericalPoint& p) const;
        // calculate the line distance (i.e the distance in R^3) to another point
        myDouble segmentDist(const SphericalPoint& p) const;
        // test if the point p is with the ring (theta, theta + deltaTheta), centered at "this" point.
        bool withinThetaRing(const SphericalPoint& p, myDouble theta1, myDouble deltaTheta) const;
        // test if two points are identical
        bool operator==(const SphericalPoint & p) const;
        // test if two points are on the same sphere
        bool onTheSameSphere(const SphericalPoint & p) const;

        /**
         *  @brief given vectors of Spherical Point, with point P(given index), delete the points
         *  which lie in the theta ring of the point
         *  @param vec3d, the vectors of Spherical Points
         *  @param _index, the index of the point
         *  @param _theta, the degree
         * */
        void popThetaRing(vec_3DCoord & vec3d, 
            typeIndex _index, myDouble _theta);
        // some get functions
        myDouble getRadius() const { return radius; }
        myDouble getTheta() const { return theta; }
        myDouble getPhi() const { return phi; }
        myDouble getXCoord() const { return xCoord; }
        myDouble getYCoord() const { return yCoord; }
        myDouble getZCoord() const { return zCoord; }
        void setXYZ(myDouble x, myDouble y, myDouble z);
        // given the spherical coordinate, set the point data
        void setRTP(myDouble r, myDouble t, myDouble p);
        /**
         * calculate the coordinate when the z-achis is rotated to point
         * point at the point with polar coordinate (radius, _theta0, _phi0)
         * @param _theta0 the theta value of the new z-axis
         * @param _phi0 the phi value of the new z-axis
         * */
        void thetaRotate(myDouble _theta0,myDouble _phi0);
        // set the standard coordinate given the point p
        // i.e the former coordinate is with respect to p,
        // now set the standard coordinate 
        // i.e, (*this) = (1,1,1), p = (1,0,0),(in cartesian) then after calling
        // the coordinate is (2,1,1)
        void setCoord(const SphericalPoint & p);
        // output function
        friend ostream & operator<<(ostream & os, const SphericalPoint & p);
        friend void swap(SphericalPoint & p1, SphericalPoint & p2);
    protected:
        // test if the three numbers given are in the asending order
        bool ifAscend(myDouble x,myDouble y, myDouble z) const { return x <= y && y <= z; }
    
};

#endif