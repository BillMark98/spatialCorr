// #include "myDefineConst.h"
// #include "sphericalPoint.h"
#include "poissonPoint2d.h"
// #include "sphereTwoCorrelation.h"
#include "myAVL.h"


#ifndef __CLASS_FindPair_H
#define __CLASS_FindPair_H

class FindPair {
    private:
        MyAVL distContainer;
        // PoissonPoint2d & poissonRef;
        unsigned long N;
    public:
        // FindPair():poissonRef(nullptr){}
        FindPair(PoissonPoint2d & pd, bool nullData = false);
        FindPair(const std::vector<SphericalPoint> & sp);
        // find the pairs of points within the range (theta, theta + deltaTheta)
        unsigned long findPairs(double theta, double deltaTheta);
        // void setPoisson();

        unsigned long getNumber() const { return N;}
};
#endif