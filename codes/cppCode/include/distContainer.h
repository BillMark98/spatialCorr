// #include "myDefineConst.h"
#ifndef __DISTCONTAINER_H
#define __DISTCONTAINER_H
class DistContainer {

    public:
    DistContainer(){}
    virtual unsigned long getPair(double theta, double deltaTheta) = 0;
    virtual void putValue(double val) = 0;
};
#endif