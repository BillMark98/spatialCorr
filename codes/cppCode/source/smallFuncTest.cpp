#include "sphericalPoint.h"
#include <iomanip>
#include <c_utils.h>
#include <sharp_geomhelpers.h>
#include "healpixCorrelation.h"
// #include <sharp_internal.h>
#include <sharp_cxx.h>
ostream & operator<<(ostream & os, sharp_ringpair);
ostream & operator<<(ostream & os, sharp_geom_info *);
int main() 
{
    // SphericalPoint sp(1,2,3);
    int nrings = 2;
    double *theta=RALLOC(double,nrings);
    cout << "The initial: " << endl;
    for (int i = 0; i < nrings; i++) {
        cout << theta[i] <<" ";
    }
    cout << endl;

    theta[0] = (theta == NULL ? 1 : theta[0]);
    cout << "Now: " << endl;
    for (int i = 0; i < nrings; i++) {
        cout << theta[i] <<" ";
    }
    cout << endl;

    int nside = 2;
    int stride = 1;
    nrings = 4 * nside - 1;
    double * weight = nullptr;
    sharp_geom_info * geom_info = 0;
    sharp_make_subset_healpix_geom_info (nside,stride, nrings,
    NULL, weight, &geom_info);
    cout << "The weight\n";
    // for (int i = 0; i < nrings; i++) {
    //     cout << weight[i] << " ";
    // }
    // cout << endl;
    // cout << sp << endl;
    cout << geom_info << endl;
    sharp_destroy_geom_info(geom_info);

    // int ntrans = 1;
    // int spin = 0;
    // sharp_jobtype type = SHARP_MAP2ALM;
    // int theNum = sharp_nv_oracle(type, spin, ntrans);
    // cout << "the Num: " << theNum << endl;

    int lmax = 4;
    int mmax = 4;
    nside_dummy nd;
    Healpix_Correlation hmap(nside,RING,nd);
    Alm<xcomplex<double> > alm(lmax,mmax);
    bool addAlm = 0;
    void *aptr = const_cast<void *>(reinterpret_cast<const void *> (&alm));
    void *mptr=const_cast<void *>(reinterpret_cast<const void *> (&hmap));
    sharp_alm_info *ainfo;
    sharp_geom_info *ginfo;

    sharp_make_subset_healpix_geom_info (nside,stride, nrings,
    NULL, weight, &ginfo);
    cout << "npairs: " << ginfo -> npairs << endl;
    sharp_make_triangular_alm_info (lmax, mmax, 1, &ainfo); 
    cout << "the flag of ainfo: " << ainfo -> flags << endl;
    int flags=cxxjobhelper__<double>::val | (addAlm ? SHARP_ADD : 0);
    cout << "the flags is : " << flags << endl;
        
    map
    flags|=SHARP_USE_WEIGHTS;
    // sharp_nv_oracle
    flags |= 2;
    sharp_execute (SHARP_MAP2ALM,0,&aptr,&mptr,ginfo,ainfo,1,flags,0,0);

    return 0;
}

ostream & operator<<(ostream & os, sharp_ringpair ringP)
{
    os << "first: " << endl;
    /*
double theta, phi0, weight, cth, sth;
  ptrdiff_t ofs;
  int nph, stride;
    */
    os << setw(10) << "theta" << setw(10) << "phi0" << setw(10) << "weight" 
    << setw(10) << "ofs" << setw(10) << "nph" << setw(10) << "stride" << endl;
    os << setw(10) << ringP.r1.theta << setw(10) << ringP.r1.phi0 << setw(10) 
        << ringP.r1.weight
    << setw(10) << ringP.r1.ofs << setw(10) << ringP.r1.nph << setw(10) << ringP.r1.stride << endl;

    os << "second:" << endl;
    os << setw(10) << "theta" << setw(10) << "phi0" << setw(10) << "weight" 
    << setw(10) << "ofs" << setw(10) << "nph" << setw(10) << "stride" << endl;
    os << setw(10) << ringP.r2.theta << setw(10) << ringP.r2.phi0 << setw(10) 
        << ringP.r2.weight
    << setw(10) << ringP.r2.ofs << setw(10) << ringP.r2.nph << setw(10) << ringP.r2.stride << endl;
    return os;
}
ostream & operator<<(ostream & os, sharp_geom_info *gm)
{
    for (int i = 0; i < gm -> npairs; i++) {
        cout << "*******************\n\nindex" << i << endl;
        cout << gm -> pair[i] << endl;
    }
    return os;
}