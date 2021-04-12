#include "sphericalPoint.h"

void mySwap(myDouble & v1, myDouble & v2) {
    myDouble temp = v1;
    v1 = v2;
    v2 = temp;
}

SphericalPoint::SphericalPoint(myDouble r, myDouble t, myDouble p, int mode) {
    if(mode == 0) {
        radius = r;
        theta = t;
        phi = p;
        mode = 0;
        sphericalToCartesian();
    }
    else if(mode == 1) {
        xCoord = r;
        yCoord = t;
        zCoord = p;
        mode = 1;
        cartesianToSpherical();
    }
    else {
        cerr << "mode have to be either 0 or 1" << endl;
        exit(MODEERROR);
    }
}
void SphericalPoint::sphericalToCartesian() {
    myDouble sT = sin(theta);
    xCoord = radius * sT * cos(phi);
    yCoord = radius * sT * sin(phi);
    zCoord = radius * cos(theta);
}
void SphericalPoint::cartesianToSpherical() {
    radius = sqrt(xCoord * xCoord + yCoord * yCoord + zCoord * zCoord);
    theta = acos(zCoord/radius);
    // myDouble xyRadius = radius * sin(theta);
    // phi = acos(xCoord/xyRadius);
    // if(yCoord > 0) {
    //     phi += M_PI;
    // }
    phi = atan2(yCoord,xCoord);
}
myDouble SphericalPoint::sphericalAngleDist(const SphericalPoint& p) const {
    // if radius not equal, indicating not in the same sphere, error
#ifndef signalMode
    if(fabs(radius - p.radius) > MY_EPSILON) {
        cerr << "the two points not in the same sphere, cant calculate angular distance." << endl;
        exit(NOT_IN_THE_SAME_SPHERE);
    }
#else
    if(radius != p.radius) {
        cerr << "the two points not in the same sphere, cant calculate angular distance." << endl;
        exit(NOT_IN_THE_SAME_SPHERE);
    }
#endif
    myDouble segmentLength = segmentDist(p);
    // using the cosine theorem
    myDouble cosVal = 1 - segmentLength * segmentLength/(2 * radius * radius);
    return acos(cosVal);
}
myDouble SphericalPoint::segmentDist(const SphericalPoint& p) const {
    // if radius not equal, indicating not in the same sphere, error
#ifndef signalMode
    if(fabs(radius - p.radius) > MY_EPSILON) {
        cerr << "the two points not in the same sphere, cant calculate angular distance." << endl;
        exit(NOT_IN_THE_SAME_SPHERE);
    }
#else
    if(radius != p.radius) {
        cerr << "the two points not in the same sphere, cant calculate angular distance." << endl;
        exit(NOT_IN_THE_SAME_SPHERE);
    }
#endif

    myDouble dist = (xCoord - p.xCoord) * (xCoord - p.xCoord) +
                    (yCoord - p.yCoord) * (yCoord - p.yCoord) + 
                    (zCoord - p.zCoord) * (zCoord - p.zCoord);
    return sqrt(dist);
}

bool SphericalPoint::withinThetaRing(const SphericalPoint& p, myDouble theta1, myDouble deltaTheta) const {
    myDouble thetaDist = sphericalAngleDist(p);
    if  (ifAscend(theta1, thetaDist, theta1 + deltaTheta)) {
        return true;
    }
    else {
        return false;
    }
}

bool SphericalPoint::operator==(const SphericalPoint & p) const {
#ifndef signalMode
    if (fabs(radius - p.radius) > MY_EPSILON) {
        return false;
    }
    else if (fabs(theta - p.theta) > MY_EPSILON) {
        return false;
    }
    else if (fabs(phi - p.phi) > MY_EPSILON) {
        return false;
    }
#else
    if (radius != p.radius) {
        return false;
    }
    else if (theta != p.theta) {
        return false;
    }
    else if (phi != p.phi) {
        return false;
    }
#endif
    return true;
}

bool SphericalPoint::onTheSameSphere(const SphericalPoint & p) const
{
#ifndef signalMode
    return fabs(radius - p.radius) < MY_EPSILON;
#else
    return radius == p.radius;
#endif
}

void SphericalPoint::popThetaRing(vec_3DCoord & vec3d,
    typeIndex _index, myDouble _theta)
{
    typeIndex _end = vec3d.size();
    for(typeIndex i, legalIndex = _end - 1; i >= 0; i--) {
        if(i != _index) {
            myDouble dT = sphericalAngleDist(vec3d[i]);
            if(dT <= _theta) {
                if(i == legalIndex) {
                    vec3d.pop_back();
                    legalIndex--;
                }
                else {
                    swap(vec3d[i],vec3d[legalIndex--]);
                    vec3d.pop_back();
                }
            }
        }
        else {
            if(!((*this) == vec3d[i])) {
                cout << "error in SphericalPoint::indexThetaRing\n";
                exit(DELETE_POINT_ERROR);
            }
            if(i == legalIndex) {
                vec3d.pop_back();
                legalIndex--;
            }
            else {
                swap(vec3d[i], vec3d[legalIndex--]);
                vec3d.pop_back();
            }
        }
    }
}    

void SphericalPoint::setXYZ(myDouble x, myDouble y, myDouble z) {
    xCoord = x;
    yCoord = y;
    zCoord = z;
    cartesianToSpherical();
}

void SphericalPoint::setRTP(myDouble r, myDouble t, myDouble p) {
    radius = r;
    theta = t;
    phi = p;
    sphericalToCartesian();
}

void SphericalPoint::thetaRotate(myDouble _theta0,myDouble _phi0)
{
/*
    cos(phi0)^2*cos(theta0)+sin(phi0)^2, sin(phi0)*cos(phi0)*(cos(theta0)-1), cos(phi0)sin(phi0)
    cos(phi0)*sin(phi0)*(cos(theta0)-1), 1+(cos(theta0)-1)*sin(phi0)^2,     sin(theta0)sin(phi0)
    -cos(phi0)*sin(thet0)              , -sin(theta0)*sin(phi0)           , cos(theta0)

*/
    myDouble _x = xCoord;
    myDouble _y = yCoord;
    myDouble _z = zCoord;
    myDouble _cosPhi = cos(_phi0);
    myDouble _sinPhi = sin(_phi0);
    myDouble _cosTheta = cos(_theta0);
    myDouble _sinTheta = sin(_theta0);
    xCoord = (square(_cosPhi) * (_cosTheta - 1) + 1) * _x 
                    + (_sinPhi * _cosPhi * (_cosTheta - 1) * _y)
                    + (_cosPhi * _sinTheta) * _z;
    yCoord = (_cosPhi * _sinPhi * (_cosTheta - 1)) * _x
                    + (square((_sinPhi)) * (_cosTheta - 1) + 1) * _y
                    + _sinTheta * _sinPhi * _z;
    zCoord = -_cosPhi * _sinTheta * _x
                    - _sinTheta * _sinPhi * _y
                    + _cosTheta * _z;
    cartesianToSpherical();

}

void SphericalPoint::setCoord(const SphericalPoint & p) 
{
    if(!onTheSameSphere(p)) {
        cerr << "the point p is not on the same sphere, cant setCoord\n";
        cout << " this : " << *this << " p: " << p << endl;
        // exit(NOT_IN_THE_SAME_SPHERE);
    }
    thetaRotate(p.theta, p.phi);
}

ostream & operator<<(ostream & os, const SphericalPoint & p)
{
    os << "(" << p.xCoord << "," << p.yCoord << "," << p.zCoord  << ")" << ", polar:"
        << "(" << p.radius << "," << p.theta << "," << p.phi << ")" << endl;
    return os;
}

void swap(SphericalPoint & p1, SphericalPoint & p2) {
    mySwap(p1.xCoord,p2.xCoord);
    mySwap(p1.yCoord, p2.yCoord);
    mySwap(p1.zCoord, p2.zCoord);
    mySwap(p1.radius,p2.radius);
    mySwap(p1.theta,p2.theta);
    mySwap(p1.phi, p2.phi);
}