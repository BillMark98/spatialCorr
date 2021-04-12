#include "findPair.h"
FindPair::FindPair(PoissonPoint2d & pd, bool nullData) : distContainer() {
    if(!nullData) {
        // use another version
        // vec_3DCoord vec = pd.rPoissonSphereDefault();
#ifdef FINDPAIR_DEBUG
    cout << "in the constructor FindPair for Poisson\n";
#endif

        vec_3DCoord vec;
        pd.rPoissonSphereDegreeLambdaDefault(vec);
        // distContainer.insert(vec);
#ifdef FINDPAIR_DEBUG
    cout << "in the constructor FindPair for Poisson, after the poissonpoints generated";
#endif
        distContainer.readVector(vec);
        N = pd.getNumber();
    }
}

FindPair::FindPair(const std::vector<SphericalPoint> & sp) : distContainer() {
    // distContainer.insert(vec);
    distContainer.readVector(sp);
    N = sp.size();
}

unsigned long FindPair::findPairs(double theta, double deltaTheta){
    return distContainer.getPair(theta, deltaTheta);
}

// void FindPair::setPoisson() {
//     vec_3DCoord vec = poissonRef.rPoissonSphere();
//     // distContainer.insert(vec);
//     distContainer.readVector(vec);
// }