#include "healpixCorrelation.h"
std::mutex myMutex;

int Healpix_Correlation::correlationCount = 0;

void convertDataVecIndex(const typeIndex & dataIndx, typeIndex & thetaIndex, typeIndex & phiIndex, 
    const sizeType & thetaLen, const sizeType & phiLen);

myDouble convertDeg2Radian(const myDouble & degree);
// assume angle > angleBegins
Healpix_Correlation::Healpix_Correlation(const myDouble & dTheta, 
    const myDouble angleEnd,const myDouble angleBegin, const myDouble phiEnd, const myDouble phiBegin, const string& logDir) : 
    Healpix_Map(),
    boundSet(false), pointsCount(0),
    hThetaBeg(angleBegin),
    hThetaEnd(angleEnd),
    hPhiBeg(phiBegin),
    hPhiEnd(phiEnd),
    _logDir(logDir)
{
#ifdef HEALPIX_DEBUG    
    cout << "In the Healpix_Correlation (double, double) constructor\n";
#endif
    if (angleBegin > angleEnd) {
        cout << "angleBegin > angleEnd!!!!\n";
        cout << "Healpix_Correlation constructin fails\n";
        exit(ANGLEBEGIN_LARGER_THAN_ANGLEEND);
    }
    if (phiBegin > phiEnd) {
        cout << "phiBegin > phiEnd!!!\n";
        cout << "Healpix_Correlation constructin fails\n";
        exit(ANGLEBEGIN_LARGER_THAN_ANGLEEND);        
    }
    spannedAngle = angleEnd - angleBegin;
    int order = deltaTheta2Order(dTheta);
#ifdef HEALPIX_DEBUG    
    cout << "deltaTheta2Order successfully called\n";
#endif
    Healpix_Map::Set(order,RING);
#ifdef HEALPIX_DEBUG    
    cout << "Set successfully called\n";
#endif
    update();
#ifdef HEALPIX_DEBUG    
    cout << "Healpix_Correlation (double, double) successfully called\n\n";
#endif
}

void Healpix_Correlation::query_disc_pixel_Internal(int pixel, const myDouble & theta,
        rangeset<int> &pixset) const {
            // if zero value, do not count anything
    if((*this)[pixel] < MY_EPSILON) {
        return;
    }
#ifdef QUERY_ZERO
    // cout << "theta: " << theta << "fabs(theta)" << fabs(theta) << "  deltaTheta - MY_EPSILON: " << deltaTheta - MY_EPSILON << endl; 
    // cout << "fabs(theta) < (deltaTheta - MY_EPSILON) : " << ((fabs(theta) < (deltaTheta - MY_EPSILON)) ? "true" : "false") << endl;
    
    #ifndef signalMode
        if(theta < (deltaTheta - MY_EPSILON)) {
            pixset.append(pixel,pixel + 1);
            return;
        }
    #else
        if(theta < deltaTheta) {
            pixset.append(pixel,pixel + 1);
            return;
        }
    #endif
#endif
    pointing ptg = pix2ang(pixel);
    myDouble rlat1 = ptg.theta - theta;
    myDouble zmax = cos((double)rlat1);
    // the ring where the point with zmax as z-coord locates
    int irmin = ring_above((double)zmax) + 1;
    // north pole in the disc
    if ((rlat1 <= 0) && (irmin > 1)) {
        int startpix, ringpix;
        bool dummy;
        get_ring_info_small(irmin - 1,startpix,ringpix,dummy);
        pixset.append(0,startpix + ringpix);
    }   
    myDouble rlat2 = ptg.theta + theta;
    myDouble zmin = cos((double)rlat2);
    int irmax = ring_above((double)zmin) + 1;
    // should be min?
    irmax = min<int>(irmax, bottomRing);
#ifdef HEALPIX_DEBUG
    cout << "inside the query_disc_pixel_Internal for rangeset version, begins loop\n";
#endif
    for (int iz = irmin; iz <= irmax; iz+= 1) {
        myDouble z = ring2z(iz);
        myDouble z0 = cos(ptg.theta);
        myDouble invSinT = 1. / sqrt((1 - z0) * (1 + z0));
        myDouble x = (cos((double)theta) - z * z0) * invSinT;
        myDouble ysq = 1 - x * x - z * z;
        if(ysq > 0) {
            myDouble y = sqrt((double)ysq);
            // note dphi >= 0
            myDouble dphi = atan2((double)y,(double)x);
            // intersect at exactly one pixel cell
            if(dphi == 0) {
                // get the pixel number of that cell,
                // since dphi = 0, the pixel has the same phi value 
                // as the ptg
                int pix =  zphi2pix((double)z, ptg.phi);
                // only add the point if it 's nonzero
                // if((*this)[pix] != 0) {
                //     pixset.append(pix,pix+1);
                // }
                pixset.append(pix,pix+1);
            }
            else {
                int ipix1, nr;
                bool shifted;
                // ipix1 is the start pixel of the current ring
                // nr is the number of pixels in the ring
                get_ring_info_small(iz, ipix1,nr,shifted);
                double shift = 0.0;
                if(shifted) {
                    double zdummy;
                    // shift is the phi value of the ipix1
                    pix2zphi(ipix1,zdummy,shift);
                }
                int ipix2 = ipix1 + nr - 1;
                // e.g imagine the circle divided equally into four parts
                // if the [0,pi/2) is the first one. and ptg.phi - dphi = pi/3, and shift = pi/4
                // then ip_lo = 4/(2pi) * pi/12 = 0, since the left point lies in the first cell
                // (note the first cell has pixel number ipix1 + 0)
                // and for 3/4 * pi    ip_lo = 4/(2pi) * pi/2 = 1 yes
                // int ip_lo = nearestInt(nr * inv_twopi * (ptg.phi - dphi - shift));
                // int ip_hi = nearestInt(nr * inv_twopi * (ptg.phi + dphi - shift));
                int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                if (ip_hi >= nr) {
                    ip_lo -= nr;
                    ip_hi -= nr;
                }
                if (ip_lo < 0 ) {
                    pixset.append(ipix1, ipix1 + ip_hi + 1);
                    pixset.append(ipix1 + ip_lo + nr, ipix2 + 1);
                }
                else {
                    pixset.append(ipix1 + ip_lo, ipix1 + ip_hi + 1);
                }
            }
        }
    }
#ifdef HEALPIX_DEBUG
    cout << "the content of pixset is: " << endl;
    cout << pixset << endl;
#endif
// southern cap
    if ((rlat2 > M_PI) && (irmax < bottomRing)) {
        int startpix, ringpix;
        bool dummy;
        get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
    
#ifdef ALL_ONE_DEBUG
    cout << "pixel: " << pixel << endl;
    cout << "startpix: " << startpix << endl;
    cout << "ringpix: " << ringpix << endl;
#endif
        pixset.append(startpix, bottomPixel);
    } 
}

void Healpix_Correlation::query_disc_pixel_Internal(int pixel, const myDouble & theta,
    rangeset<int> &pixset, const Healpix_Correlation & hp) const {
            // if zero value, do not count anything
    if((*this)[pixel] < MY_EPSILON) {
        return;
    }
#ifdef QUERY_ZERO
    // cout << "theta: " << theta << "fabs(theta)" << fabs(theta) << "  deltaTheta - MY_EPSILON: " << deltaTheta - MY_EPSILON << endl; 
    // cout << "fabs(theta) < (deltaTheta - MY_EPSILON) : " << ((fabs(theta) < (deltaTheta - MY_EPSILON)) ? "true" : "false") << endl;
    
    #ifndef signalMode
        if(theta < (deltaTheta - MY_EPSILON)) {
            pixset.append(pixel,pixel + 1);
            return;
        }
    #else
        if(theta < deltaTheta) {
            pixset.append(pixel,pixel + 1);
            return;
        }
    #endif
#endif
    pointing ptg = pix2ang(pixel);
    myDouble rlat1 = ptg.theta - theta;
    myDouble zmax = cos((double)rlat1);
    // the ring where the point with zmax as z-coord locates
    int irmin = ring_above((double)zmax) + 1;
    // north pole in the disc
    if ((rlat1 <= 0) && (irmin > 1)) {
        int startpix, ringpix;
        bool dummy;
        get_ring_info_small(irmin - 1,startpix,ringpix,dummy);
        pixset.append(0,startpix + ringpix);
    }   
    myDouble rlat2 = ptg.theta + theta;
    myDouble zmin = cos((double)rlat2);
    int irmax = ring_above((double)zmin) + 1;
    irmax = min<int>(irmax, bottomRing);
#ifdef HEALPIX_DEBUG
    cout << "inside the query_disc_pixel_Internal for rangeset version, begins loop\n";
#endif
    for (int iz = irmin; iz <= irmax; iz+= 1) {
        myDouble z = ring2z(iz);
        myDouble z0 = cos(ptg.theta);
        myDouble invSinT = 1. / sqrt((1. - z0) * (1. + z0));
        myDouble x = (cos((double)theta) - z * z0) * invSinT;
        myDouble ysq = 1 - x * x - z * z;
        if(ysq > 0) {
            myDouble y = sqrt((double)ysq);
            // note dphi >= 0
            myDouble dphi = atan2((double)y,(double)x);
            // intersect at exactly one pixel cell
            if(dphi == 0) {
                // get the pixel number of that cell,
                // since dphi = 0, the pixel has the same phi value 
                // as the ptg
                int pix =  zphi2pix((double)z, ptg.phi);
                // only add the point if it 's nonzero
                // reference point is (*this)[pix]
                if((*this)[pix] != 0) {
                    pixset.append(pix,pix+1);
                }

            }
            else {
                int ipix1, nr;
                bool shifted;
                // ipix1 is the start pixel of the current ring
                // nr is the number of pixels in the ring
                get_ring_info_small(iz, ipix1,nr,shifted);
                double shift = 0.0;
                if(shifted) {
                    double zdummy;
                    // shift is the phi value of the ipix1
                    pix2zphi(ipix1,zdummy,shift);
                }
                int ipix2 = ipix1 + nr - 1;
                // e.g imagine the circle divided equally into four parts
                // if the [0,pi/2) is the first one. and ptg.phi - dphi = pi/3, and shift = pi/4
                // then ip_lo = 4/(2pi) * pi/12 = 0, since the left point lies in the first cell
                // (note the first cell has pixel number ipix1 + 0)
                // and for 3/4 * pi    ip_lo = 4/(2pi) * pi/2 = 1 yes
                // int ip_lo = nearestInt(nr * inv_twopi * (ptg.phi - dphi - shift));
                // int ip_hi = nearestInt(nr * inv_twopi * (ptg.phi + dphi - shift));
                int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                if (ip_hi >= nr) {
                    ip_lo -= nr;
                    ip_hi -= nr;
                }
                if (ip_lo < 0 ) {       
                    appendNonZeroRange(pixset,hp,ipix1,ipix1 + ip_hi + 1);
                    appendNonZeroRange(pixset,hp,ipix1 + ip_lo + nr,ipix2 + 1);

                }
                else {
                                  
                    appendNonZeroRange(pixset,hp,ipix1 + ip_lo, ipix1 + ip_hi + 1);

                }
            }
        }
    }
#ifdef HEALPIX_DEBUG
    cout << "the content of pixset is: " << endl;
    cout << pixset << endl;
#endif
    // southern cap
    if ((rlat2 > M_PI) && (irmax < bottomRing)) {
        int startpix, ringpix;
        bool dummy;
        get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
    
#ifdef ALL_ONE_DEBUG
    cout << "pixel: " << pixel << endl;
    cout << "startpix: " << startpix << endl;
    cout << "ringpix: " << ringpix << endl;
#endif
        appendNonZeroRange(pixset,hp,startpix, bottomPixel);
    } 
}

void Healpix_Correlation::query_disc_pixel_Internal(int pixel, const myDouble & theta, std::vector<int> & pixset) const 
{
    rangeset<int> rangeTemp;
    query_disc_pixel_Internal(pixel,theta,rangeTemp);
#ifdef HEALPIX_DEBUG
    cout << "in the query_disc_pixel_Internal for the vector version" << endl;
#endif
    pixset = rangeTemp.toVector();
#ifdef HEALPIX_DEBUG
    cout << "the pixset successfully transformed to vector" << endl;
#endif
}

void Healpix_Correlation::query_disc_pixel_Internal(int pixel, const myDouble & theta,
std::vector<int> & pixset, const Healpix_Correlation & hp) const {
    rangeset<int> rangeTemp;
    query_disc_pixel_Internal(pixel,theta,rangeTemp,hp);
#ifdef HEALPIX_DEBUG
    cout << "in the query_disc_pixel_Internal for the vector version" << endl;
#endif
    pixset = rangeTemp.toVector();
#ifdef HEALPIX_DEBUG
    cout << "the pixset successfully transformed to vector" << endl;
#endif
}


void Healpix_Correlation::query_disc_pixel_InternalVec(int pixel,
 const myDouble & theta, std::vector<int> & pixset,
    const Healpix_Correlation & hp, int ignore) const 
{
                // if zero value, do not count anything
    if((*this)[pixel] < MY_EPSILON) {
        return;
    }
#ifdef QUERY_ZERO
    
    #ifndef signalMode
        if(theta < (deltaTheta - MY_EPSILON)) {
            if (pixel != ignore)
                pixset.push_back(pixel);
            return;
        }
    #else
        if(theta < deltaTheta) {
            if (pixel != ignore)
                pixset.push_back(pixel);
            return;
        }
    #endif
#endif
    pointing ptg = pix2ang(pixel);
    myDouble rlat1 = ptg.theta - theta;
    myDouble zmax = cos((double)rlat1);
    // the ring where the point with zmax as z-coord locates
    int irmin = ring_above((double)zmax) + 1;
    // north pole in the disc
    if ((rlat1 <= 0) && (irmin > 1)) {
        int startpix, ringpix;
        bool dummy;


        get_ring_info_small(irmin - 1,startpix,ringpix,dummy);
    
#ifdef ALL_ONE_DEBUG
    cout << "pixel: " << pixel << endl;
    cout << "startpix: " << startpix << endl;
    cout << "ringpix: " << ringpix << endl;
#endif
        VecAppendNonZeroPoint(pixset,hp,0,startpix + ringpix,ignore);
    }   
    myDouble rlat2 = ptg.theta + theta;
    myDouble zmin = cos((double)rlat2);
    int irmax = ring_above((double)zmin) + 1;
    irmax = min<int>(irmax, bottomRing);
#ifdef HEALPIX_DEBUG
    cout << "inside the query_disc_pixel_Internal for rangeset version, begins loop\n";
#endif        
    myDouble z0 = cos(ptg.theta);
    myDouble invSinT = 1. / sqrt((1 - z0) * (1 + z0));

    for (int iz = irmin; iz <= irmax; iz+= 1) {
        myDouble z = ring2z(iz);
        myDouble x = (cos((double)theta) - z * z0) * invSinT;
        myDouble ysq = 1 - x * x - z * z;
        if(ysq > 0) {
            myDouble y = sqrt((double)ysq);
            // note dphi >= 0
            myDouble dphi = atan2((double)y,(double)x);
            // intersect at exactly one pixel cell
            if(dphi == 0) {
                // get the pixel number of that cell,
                // since dphi = 0, the pixel has the same phi value 
                // as the ptg
                int pix =  zphi2pix((double)z, ptg.phi);
                // only add the point if it 's nonzero
                if((*this)[pix] != 0 && (int)pix != ignore) {
                    // VecAppendNonZeroPoint(pixset,hp, pix,pix+1);
                    pixset.push_back(pix);
                }

            }
            else {
                int ipix1, nr;
                bool shifted;
                // ipix1 is the start pixel of the current ring
                // nr is the number of pixels in the ring
                get_ring_info_small(iz, ipix1,nr,shifted);
                double shift = 0.0;
                if(shifted) {
                    double zdummy;
                    // shift is the phi value of the ipix1
                    pix2zphi(ipix1,zdummy,shift);
                }
                int ipix2 = ipix1 + nr - 1;
                // e.g imagine the circle divided equally into four parts
                // if the [0,pi/2) is the first one. and ptg.phi - dphi = pi/3, and shift = pi/4
                // then ip_lo = 4/(2pi) * pi/12 = 0, since the left point lies in the first cell
                // (note the first cell has pixel number ipix1 + 0)
                // and for 3/4 * pi    ip_lo = 4/(2pi) * pi/2 = 1 yes
                // int ip_lo = nearestInt(nr * inv_twopi * (ptg.phi - dphi - shift));
                // int ip_hi = nearestInt(nr * inv_twopi * (ptg.phi + dphi - shift));
                int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                if (ip_hi >= nr) {
                    ip_lo -= nr;
                    ip_hi -= nr;
                }
                if (ip_lo < 0 ) {
#ifdef ALL_ONE_DEBUG
    cout << "in the ip_lo < 0: condition" << endl;
    cout << "pixel: " << pixel << endl;
    cout << "ipix1: " << ipix1 << endl;
    cout << "ip_lo: " << ip_lo << endl;
    cout << "ipix1 + ip_lo: " << ipix1 + ip_lo << endl;
    cout << "ipix1 + ip_lo + n1: " << ipix1 + ip_lo + nr << endl;
    cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
#endif
                    VecAppendNonZeroPoint(pixset,hp,ipix1,ipix1 + ip_hi + 1,ignore);
                    // need to add an addition test, it is possible, due to the round
                    // ip_lo + nr == ip_hi, which means that the whole ring will be covered
                    // and the direct opposite of pixel will be counted twice
                    // e.g if not modify it
                    /*
                            pixel: 0
                            theta: 0.2
                            pixset: 
                            0	
                            1	
                            2	
                            2	
                            3	
                            4	
                            5	
                            6	
                            11	
                    */
                   // we find the 2 pixel is added twice, which is exactly due to the 
                   // fact that the ring 1 (where 2 lies) is complety covered, but due to the rounding
                   // ip_lo + nr and ip_hi intersect
                    if (ip_lo + nr <= ip_hi) {
                        // ip_lo + nr = ip_hi + 1
                        ip_lo = ip_hi + 1 - nr;
                    }
                    VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_lo + nr,ipix2 + 1,ignore);
                }
                else { // ip_lo >= 0
#ifdef ALL_ONE_DEBUG
    cout << "in the ip_lo >= 0: condition" << endl;
    cout << "pixel: " << pixel << endl;
    cout << "ipix1: " << ipix1 << endl;
    cout << "ip_lo: " << ip_lo << endl;
    cout << "ip_hi: " << ip_hi << endl;

    cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
#endif                
                    VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_lo, ipix1 + ip_hi + 1,ignore);
                }
            }
        }
    }
#ifdef HEALPIX_DEBUG
    cout << "the content of pixset is: " << endl;
    cout << pixset << endl;
#endif
    // southern cap
    // it is possible that irmax == bottomRing and the above 
    // code does not cover the case for bottomRing, ysq < 0
    // for example pixel = 4, dT = 0.1
    // theta = 3.1, will lead to ysq < 0
    // even when irmax < bottomRing, it's still possible, 
    // that irmax is not covered,
    // for example pixel = 12, dT = 0.1
    // theta = 3.1, irmax = 30, will lead to ysq < 0
    // irmin <= irmax is tested, there are some extreme cases where
    // irmin > irmax, e.g pixel = 0, dT = 0.02, at theta = 3.16, irmin = 256, irmax = 254
    // already all the points have been covered, but irmax + 1 = bottomRing = 255, so the southpole region
    // will be added twice, but should not test irmin <= irmax, 
    // e.g pixel = 4, dT = 0.1, theta = 3.1, irmax = 30, irmin = 31, but only to ring 30 is added sofar
    // (should actually cover the south pole, i.e, ring 31)
    // but irmin <= bottomRing also wont work, because pixel = 12, dT = 0.08,
    // theta = 3.2, irmin = 63, irmax = 60, ysq > 0, already covered up to ring 62
    // but the code thinks that ysq > 0, so the irmax case already covered (which is not, due to irmin > irmax, so the 
    // for loop in the middle cap wont cover the case), so will add irmax + 1 up to bottomRing, but 
    // the ring 61,62 has already been included, leading to duplicate counts.
    // instead test ysq > 0 && irmin <= irmax to make sure that the ysq > 0 is already covered.
    if ((rlat2 > M_PI)) {  // 
    //     if (irmax < bottomRing) {
    //         int startpix, ringpix;
    //         bool dummy;
    //         myDouble z = ring2z(irmax);
    //         myDouble x = (cos((double)theta) - z * z0) * invSinT;
    //         myDouble ysq = 1 - x * x - z * z;
    //         if (ysq > 0) {
    //             // irmax case already covered
    //             get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
    //             VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);           
    //         }
    //         else {
    //             // add the whole bottomRing
    //             int startpix, ringpix;
    //             bool dummy;
    //             get_ring_info_small(irmax,startpix,ringpix,dummy);
    //             VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);
    //         }            
    // #ifdef ALL_ONE_DEBUG
    //     cout << "pixel: " << pixel << endl;
    //     cout << "startpix: " << startpix << endl;
    //     cout << "ringpix: " << ringpix << endl;
    // #endif
            
    //     }
    //     else { // irmax == bottomRing
    //         myDouble z = ring2z(bottomRing);
    //         myDouble x = (cos((double)theta) - z * z0) * invSinT;
    //         myDouble ysq = 1 - x * x - z * z;
    //         if (ysq > 0) {
    //             // case already covered
    //             return;              
    //         }
    //         else {
    //             // add the whole bottomRing
    //             int startpix, ringpix;
    //             bool dummy;
    //             get_ring_info_small(bottomRing,startpix,ringpix,dummy);
    //             VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);
    //         }
    //     }
        int startpix, ringpix;
        bool dummy;
        myDouble z = ring2z(irmax);
        myDouble x = (cos((double)theta) - z * z0) * invSinT;
        myDouble ysq = 1 - x * x - z * z;
        // ysq > 0 and irmax >= irmin makes sure that irmax has already be covered
        // only need to add the ring under it
        // if however, irmin > irmax, then up to irmin - 1 is fully covered
        if (ysq > 0 && irmax >= irmin) {
            // irmax case already covered
            get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
            VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);           
        }
        else if (irmin > irmax) {  // it means the whole sphere should be covered (up to bound Ring)
        // since up till now, only those rings up to irmin - 1 is covered in the norther cap case
        // fill the blank to the bottomPixel
            get_ring_info_small(irmin,startpix,ringpix,dummy);
            VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);
        }
        else { // ysq < 0 && irmin <= irmax, so dont worry, that part of the region south to the irmax
            // is already covered
            // add the whole bottomRing
            get_ring_info_small(irmax,startpix,ringpix,dummy);
            VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);
        } 
    } 
}

countType Healpix_Correlation::query_disc_count(int pixel, 
    const myDouble & theta,
    const Healpix_Correlation & hp,
    int ignore) const 
{
    
    if((*this)[pixel] < MY_EPSILON) {
        return 0;
    }
#ifdef QUERY_ZERO
    #ifndef signalMode
    if(theta < (deltaTheta - MY_EPSILON)) {
        if (pixel != ignore) {
            return (hp)[pixel];
        }
        else {
            return 0;
        }
    }
    #else
    if(theta < deltaTheta) {
        if (pixel != ignore) {
            return (hp)[pixel];
        }
        else {
            return 0;
        }
    }
    #endif
#endif
    countType pointCount = 0;
    pointing ptg = pix2ang(pixel);
    myDouble rlat1 = ptg.theta - theta;
    myDouble zmax = cos((double)rlat1);
    // the ring where the point with zmax as z-coord locates
    int irmin = ring_above((double)zmax) + 1;
    // north pole in the disc
    if ((rlat1 <= 0) && (irmin > 1)) {
        int startpix, ringpix;
        bool dummy;
        get_ring_info_small(irmin - 1,startpix,ringpix,dummy);
        pointCount += IntervallNonZeroPoint(hp,0,startpix + ringpix,ignore);
    }   
    myDouble rlat2 = ptg.theta + theta;
    myDouble zmin = cos((double)rlat2);
    int irmax = ring_above((double)zmin) + 1;
    irmax = min<int>(irmax, bottomRing);
#ifdef HEALPIX_DEBUG
    cout << "inside the query_disc_pixel_Internal for rangeset version, begins loop\n";
#endif
    myDouble z0 = cos(ptg.theta);
    myDouble invSinT = 1. / sqrt((1 - z0) * (1 + z0));
    for (int iz = irmin; iz <= irmax; iz+= 1) {
        myDouble z = ring2z(iz);
        myDouble x = (cos((double)theta) - z * z0) * invSinT;
        myDouble ysq = 1 - x * x - z * z;
        if(ysq > 0) {
            myDouble y = sqrt((double)ysq);
            // note dphi >= 0
            myDouble dphi = atan2((double)y,(double)x);
            // intersect at exactly one pixel cell
            if(dphi == 0) {
                // get the pixel number of that cell,
                // since dphi = 0, the pixel has the same phi value 
                // as the ptg
                int pix =  zphi2pix((double)z, ptg.phi);
                // only add the point if it 's nonzero
                // reference point is (*this)[pix] but add hp[pix]
                if((hp)[pix] != 0 && (int)pix != ignore) {
                    pointCount += (hp)[pix];
                }

            }
            else {
                int ipix1, nr;
                bool shifted;
                // ipix1 is the start pixel of the current ring
                // nr is the number of pixels in the ring
                get_ring_info_small(iz, ipix1,nr,shifted);
                double shift = 0.0;
                if(shifted) {
                    double zdummy;
                    // shift is the phi value of the ipix1
                    pix2zphi(ipix1,zdummy,shift);
                }
                int ipix2 = ipix1 + nr - 1;
                // e.g imagine the circle divided equally into four parts
                // if the [0,pi/2) is the first one. and ptg.phi - dphi = pi/3, and shift = pi/4
                // then ip_lo = 4/(2pi) * pi/12 = 0, since the left point lies in the first cell
                // (note the first cell has pixel number ipix1 + 0)
                // and for 3/4 * pi    ip_lo = 4/(2pi) * pi/2 = 1 yes
                // int ip_lo = nearestInt(nr * inv_twopi * (ptg.phi - dphi - shift));
                // int ip_hi = nearestInt(nr * inv_twopi * (ptg.phi + dphi - shift));
                int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                if (ip_hi >= nr) {
                    ip_lo -= nr;
                    ip_hi -= nr;
                }
                if (ip_lo < 0 ) {
                    // pixset.append(ipix1, ipix1 + ip_lo + 1);
                    // pixset.append(ipix1 + ip_lo + nr, ipix2 + 1);
                    pointCount += IntervallNonZeroPoint(hp,ipix1,ipix1 + ip_hi + 1,ignore);
                    // see the notes in the VecAppendNonZeroPoint
                    if (ip_lo + nr <= ip_hi) {
                        // ip_lo + nr = ip_hi + 1
                        ip_lo = ip_hi + 1 - nr;
                    }
                    pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_lo + nr,ipix2 + 1,ignore);
                }
                else {
                    // pixset.append(ipix1 + ip_lo, ipix1 + ip_hi + 1);
                    pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_lo, ipix1 + ip_hi + 1,ignore);
                }
            }
        }
    }

#ifdef HEALPIX_DEBUG
    cout << "the content of pixset is: " << endl;
    cout << pixset << endl;
#endif
    // southern cap
    if ((rlat2 > M_PI)) {

        int startpix, ringpix;
        bool dummy;
        myDouble z = ring2z(irmax);
        myDouble x = (cos((double)theta) - z * z0) * invSinT;
        myDouble ysq = 1 - x * x - z * z;
        // ysq > 0 and irmax >= irmin makes sure that irmax has already be covered
        // only need to add the ring under it
        // if however, irmin > irmax, then up to irmin - 1 is fully covered
        if (ysq > 0 && irmax >= irmin) {
            // irmax case already covered
            get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
            pointCount += IntervallNonZeroPoint(hp,startpix, bottomPixel + 1,ignore);           
        }
        else if (irmin > irmax) {  // it means the whole sphere should be covered (up to bound Ring)
        // since up till now, only those rings up to irmin - 1 is covered in the norther cap case
        // fill the blank to the bottomPixel
            get_ring_info_small(irmin,startpix,ringpix,dummy);
            pointCount += IntervallNonZeroPoint(hp,startpix, bottomPixel + 1,ignore);
        }
        else { // ysq < 0 && irmin <= irmax, so dont worry, that part of the region south to the irmax
            // is already covered
            // add the whole bottomRing
            get_ring_info_small(irmax,startpix,ringpix,dummy);
            pointCount += IntervallNonZeroPoint(hp,startpix, bottomPixel + 1,ignore);
        } 
    } 
    return pointCount;
}

void Healpix_Correlation::query_disc_pixel_NonCumuInternalVec(int pixel, const myDouble & thetaLarge, 
            const myDouble & thetaSmall,
            std::vector<int> & pixset,
            const Healpix_Correlation & hp,
            std::string filename) const
{

#define QUERY_DISC_BASIC_OUTPUT() \
    {                                   \
        cout << "pixel: " << pixel << endl; \
        cout << "thetaLarge: " << thetaLarge << endl; \
        cout << "thetaSmall: " << thetaSmall << endl; \
    }


#define QUERY_DISC_OUTPUT() \
    {                       \
        cout << "pixel: " << pixel << endl; \
        cout << "ring: " << iz << endl; \
        cout << "thetaLarge: " << thetaLarge << endl; \
        cout << "thetaSmall: " << thetaSmall << endl; \
        cout << "x,y,z: " << setw(10) << x << setw(10) << y << setw(10) << z << endl;\
        cout << "xPrev, yPrev,z: " << setw(10) << xPrev << setw(10) << yPrev << setw(10) << z << endl;\
        cout << "ptg.phi: " << ptg.phi << endl;\
        cout << "dphi: " << dphi << endl;\
        cout << "dphiPrev: " << dphiPrev << endl;\
    }
        // cout << "x0,y0,z0: " << setw(10) << x0 << setw(10) << y0 << setw(10) << z0 << endl;

#define QUERY_DISC_IP_OUTPUT() \
    {                                       \
        cout << "nr: " << nr << endl;\
        cout << "ipix1: " << ipix1 << endl; \
        cout << "ipix2: " << ipix2 << endl; \
        cout << "ip_lo: " << ip_lo << endl; \
        cout << "ip_loPrev: " <<ip_loPrev << endl; \
        cout << "ip_hi: " << ip_hi << endl; \
        cout << "ip_hiPrev: " << ip_hiPrev << endl; \
        cout << "shift: " << shift << endl;\
    }

                    // if zero value, do not count anything
    if((*this)[pixel] < MY_EPSILON) {
        return;
    }
#ifdef QUERY_ZERO
    
    #ifndef signalMode
        if (thetaSmall < (MY_EPSILON) && thetaLarge < (deltaTheta - MY_EPSILON)) {
#ifdef QUERY_NON_CUMU_DEBUG
        cout << "\nthetaSmall < epsilon && thetaLarge < deltaTheta - epsilon\n";
        cout << "dT too large, or bin width too small\n";
        QUERY_DISC_BASIC_OUTPUT();
        cout << "deltaTheta: " << deltaTheta << endl;
#endif
            pixset.push_back(pixel);
            return;
        }
    #else
        if(thetaSmall < deltaTheta && thetaLarge < deltaTheta) {
            pixset.push_back(pixel);
            return;
        }
    #endif
#endif
    pointing ptg = pix2ang(pixel);
    myDouble rlat1 = ptg.theta - thetaLarge;
    myDouble rlat1Prev = ptg.theta - thetaSmall;
    myDouble zmax = cos((double)rlat1);
    int irmin = ring_above((double)zmax) + 1;


    // this case appears, only when thetaSmall = 0
    // calculate as in the cum case, except, the pixel itself needs to be excluded
    if (thetaSmall < (deltaTheta - MY_EPSILON)) {
        query_disc_pixel_InternalVec(pixel,thetaLarge,pixset,*this,pixel);
#ifdef QUERY_NON_CUMU_DEBUG
    cout << "\nNon cumu version: \n";
    QUERY_DISC_BASIC_OUTPUT();
    cout << "pixset: \n";
    for (int i = 0; i < (int)pixset.size(); i++) {
        cout << pixset[i] << "\t";
    }
    cout << endl;
#endif
        return;
    }
   
    myDouble zmaxPrev = cos((double)rlat1Prev);
    // the ring where the point with zmax as z-coord locates
    
    int irminPrev = ring_above((double)zmaxPrev) + 1;
    myDouble z0 = cos(ptg.theta);
    myDouble invSinT = 1. / sqrt((1 - z0) * (1 + z0));

    // north pole in the disc
    if ((rlat1 <= 0) && (irmin > 1)) {
        int startpix, ringpix;
        bool dummy;


        get_ring_info_small(irmin - 1,startpix,ringpix,dummy);
        // pixel = 4, dT = 0.08, 
        // thetaSmall = 0.08, thetaLarge = 0.16
        // ringpixPrev == 1, rlat1Prev > 0, but thetaSmall already covers 0,3
        if ((rlat1Prev <= 0) || irminPrev == 1) {
#ifdef QUERY_NON_CUMU_DEBUG
    if (irmin < irminPrev) {
        cout << "\n irmin > irminPrev\n";
        QUERY_DISC_BASIC_OUTPUT();
        cout << "irmin: " << irmin << endl;
        cout << "irminPrev: " << irminPrev << endl;
    }
#endif
            
    #ifdef QUERY_NON_CUMU_DEBUG
            int startpixPrev,ringpixPrev;
            get_ring_info_small(irminPrev - 1, startpixPrev,ringpixPrev, dummy);
            if (startpixPrev + ringpixPrev > startpix + ringpix) {
                cout << "startpixPrev + ringpixPrev > startpix + ringpix\n";
                cout << "pixel: " << pixel << endl;
                cout << "startpix: " << startpix << endl;
                cout << "startpixPrev: " << startpixPrev << endl;
                cout << "ringpix:" << ringpix << endl;
                cout << "ringpixPrev: " << ringpixPrev << endl;

            }
    #endif
            // get the corresponding value for the previous theta (thetaSmall) at the ring
            // irminPrev

            // VecAppendNonZeroPoint(pixset, hp, startpixPrev + ringpixPrev, startpix + ringpix);
    #ifdef QUERY_NON_CUMU_DEBUG
            if (irmin -irminPrev > 1) {
                cout << "irmin -irminPrev > 1\n";
                QUERY_DISC_BASIC_OUTPUT();    
                cout << "startpix: " << startpix << endl;
                cout << "startpixPrev: " << startpixPrev << endl;
                cout << "ringpix:" << ringpix << endl;
                cout << "ringpixPrev: " << ringpixPrev << endl;
                cout << "irmin: " << irmin << endl;
                cout << "irminPrev: " << irminPrev << endl;
            }
    #endif
            // since the result of ring_above may vary if the point is above
            // the pixel of that cell or not, it's better to set irmin to irminPrev
            // incase that irminPrev did not actually cover the entire string
            // e.g 
            /*
                dT = 0.1, all one map
                pixel: 0
                thetaLarge: 0.5
                thetaSmall: 0.4  
                will ignore 19 which is at ring 3 = irminPrev
                in this case irminPrev did not fully cover the ring
            */
           for (int iz = irminPrev; iz < irmin; iz++) {
                myDouble z = ring2z(iz);
                myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
                myDouble ysqPrev = 1 - xPrev * xPrev - z * z;
                if (ysqPrev > 0) {
                    myDouble yPrev = sqrt((double)ysqPrev);
                    myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                    int ipix1, nr;
                    bool shifted;
                    // ipix1 is the start pixel of the current ring
                    // nr is the number of pixels in the ring
                    get_ring_info_small(iz, ipix1,nr,shifted);
                    double shift = 0.0;
                    if(shifted) {
                        double zdummy;
                        // shift is the phi value of the ipix1
                        pix2zphi(ipix1,zdummy,shift);
                    }
                    int ipix2 = ipix1 + nr - 1;

                    int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                    int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                    if (ip_hiPrev >= nr) {
                        ip_loPrev -= nr;
                        ip_hiPrev -= nr;
                    }
                    if (ip_loPrev < 0) {
                        if (ip_loPrev + nr > ip_hiPrev) {
                            // the previous theta did not cover the whole ring
                            // fill the blank
                            VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_loPrev + nr);
                        }
                    }
                    else {
                        VecAppendNonZeroPoint(pixset, hp, ipix1,ipix1 + ip_loPrev);
                        VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev + 1, ipix2 + 1);
                    }                    
                }
                else {
                    // it means this ring has already been covered by the previous theta-disc
#ifdef QUERY_NON_CUMU_DEBUG                    
                    cout << "the northern cap, yspPrev < 0\n";
                    cout << "iz: " << iz << endl;
                    cout << "irmin: " << irmin << endl;
                    cout << "irminPrev: " << irminPrev << endl;
                    QUERY_DISC_BASIC_OUTPUT();
#endif
                }
            }
            // need to consider the case irmin separatedly
            // if irmin at the middle earth loop, and ysq < 0, then
            // it will not be considered there
                myDouble z = ring2z(irmin);
                myDouble x = (cos((double)thetaLarge) - z * z0) * invSinT;
                myDouble ysq = 1 - x * x - z * z;
                if (ysq > 0) {
                    // irmin case will be covered below
                    // get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
                    // VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);           
                }
                else {
                    // add the current ring - pixels covered by the previous 
                    myDouble z = ring2z(irmin);
                    myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
                    myDouble ysqPrev = 1 - xPrev * xPrev - z * z;   
    #ifdef QUERY_NON_CUMU_DEBUG
                    if (ysqPrev < 0) {
                        cout << "northern cap, at irmin: ysqPrev < 0\n";
                        QUERY_DISC_BASIC_OUTPUT();
                    }
    #endif
                    if (ysqPrev >= 0){
                        myDouble yPrev = sqrt((double)ysqPrev);
                        myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                        int ipix1, nr;
                        bool shifted;
                        // ipix1 is the start pixel of the current ring
                        // nr is the number of pixels in the ring
                        get_ring_info_small(irmin, ipix1,nr,shifted);
                        double shift = 0.0;
                        if(shifted) {
                            double zdummy;
                            // shift is the phi value of the ipix1
                            pix2zphi(ipix1,zdummy,shift);
                        }
                        int ipix2 = ipix1 + nr - 1;

                        int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                        int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                        if (ip_hiPrev >= nr) {
                            ip_loPrev -= nr;
                            ip_hiPrev -= nr;
                        }
                        if (ip_loPrev < 0) {
                            if (ip_loPrev + nr > ip_hiPrev) {
                                // the previous theta did not cover the whole ring
                                // fill the blank
                                VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_loPrev + nr);
                            }
                        }
                        else {
                            VecAppendNonZeroPoint(pixset, hp, ipix1,ipix1 + ip_loPrev);
                            VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev + 1, ipix2 + 1);
                        }     
                    }


                }
        }
#ifdef ALL_ONE_DEBUG
    cout << "pixel: " << pixel << endl;
    cout << "startpix: " << startpix << endl;
    cout << "ringpix: " << ringpix << endl;
#endif
        // it makes no use if irmin = 1 , because in the next for loop
        // the ysq < 0, still wont process the data
        // so in the above condition, I deleted the  irminPrev > 1
        // to include explicitly the case when irminPrev = 1, so that
        // we could consider it carefully, not losing some points
        // if use irminPrev > 1, this method fails, e.g
        /*
            pixel: 338, dT = 0.1, all one map
            thetaLarge: 1.6
            thetaSmall: 1.5
            pixset: 
            7	9	16	55	68	78	93	122	137	186	201	250	265	314	329	
            where as the correct should be 
            2, 7, 9, 16, 55, 68, 78, 93, 122, 137, 186, 201, 250, 265, 314, 329
        */

        // else {
        //     //  VecAppendNonZeroPoint(pixset,hp,0,startpix + ringpix);
        //     irmin = 1;
        // }
    }   
    myDouble rlat2 = ptg.theta + thetaLarge;
    myDouble rlat2Prev = ptg.theta + thetaSmall;
    myDouble zmin = cos((double)rlat2);
    myDouble zminPrev = cos((double)rlat2Prev);
    int irmax = ring_above((double)zmin) + 1;
    int irmaxPrev = ring_above((double)zminPrev) + 1;
    irmax = min<int>(irmax, bottomRing);
    irmaxPrev = min<int>(irmaxPrev,bottomRing);
#ifdef QUERY_NON_CUMU_DEBUG    
    if (irmaxPrev > irmax) {

        cout << "\n irmaxPrev > irmax\n";
        QUERY_DISC_BASIC_OUTPUT();
        cout << "irmax: " << irmax << endl;
        cout << "irmaxPrev: " << irmaxPrev << endl;

    }
#endif

#ifdef HEALPIX_DEBUG
    cout << "inside the query_disc_pixel_Internal for rangeset version, begins loop\n";
#endif        


    for (int iz = irmin; iz <= irmax; iz+= 1) {
        myDouble z = ring2z(iz);
        myDouble x = (cos((double)thetaLarge) - z * z0) * invSinT;
        myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
        myDouble ysq = 1 - x * x - z * z;
        myDouble ysqPrev = 1 - xPrev * xPrev - z * z;
        if(ysq > 0) {
            myDouble y = sqrt((double)ysq);
            // note dphi >= 0
            myDouble dphi = atan2((double)y,(double)x);
            if (ysqPrev < 0) {
                // intersect at exactly one pixel cell
                if(fabs(dphi) < MY_EPSILON) {
                    // get the pixel number of that cell,
                    // since dphi = 0, the pixel has the same phi value 
                    // as the ptg
                    int pix =  zphi2pix((double)z, ptg.phi);
                    // only add the point if it 's nonzero
                    if((hp)[pix] != 0) {
                        // VecAppendNonZeroPoint(pixset,hp, pix,pix+1);
                        pixset.push_back(pix);
                    }

                }
                else { // ysq > 0, ysqPrev < 0, fabs(dphi) > MY_EMSILON
                    int ipix1, nr;
                    bool shifted;
                    // ipix1 is the start pixel of the current ring
                    // nr is the number of pixels in the ring
                    get_ring_info_small(iz, ipix1,nr,shifted);
                    double shift = 0.0;
                    if(shifted) {
                        double zdummy;
                        // shift is the phi value of the ipix1
                        pix2zphi(ipix1,zdummy,shift);
                    }
                    int ipix2 = ipix1 + nr - 1;
                    // e.g imagine the circle divided equally into four parts
                    // if the [0,pi/2) is the first one. and ptg.phi - dphi = pi/3, and shift = pi/4
                    // then ip_lo = 4/(2pi) * pi/12 = 0, since the left point lies in the first cell
                    // (note the first cell has pixel number ipix1 + 0)
                    // and for 3/4 * pi    ip_lo = 4/(2pi) * pi/2 = 1 yes
                    // int ip_lo = nearestInt(nr * inv_twopi * (ptg.phi - dphi - shift));
                    // int ip_hi = nearestInt(nr * inv_twopi * (ptg.phi + dphi - shift));
                    int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                    int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                    if (ip_hi >= nr) {
                        ip_lo -= nr;
                        ip_hi -= nr;
                    }
                    if (ip_lo < 0 ) {
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo < 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ipix1 + ip_lo: " << ipix1 + ip_lo << endl;
        cout << "ipix1 + ip_lo + n1: " << ipix1 + ip_lo + nr << endl;
        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif
                        VecAppendNonZeroPoint(pixset,hp,ipix1,ipix1 + ip_hi + 1);
                        if (ip_lo + nr <= ip_hi) {
                            // ip_lo + nr = ip_hi + 1
                            ip_lo = ip_hi + 1 - nr;
                        }
                        VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_lo + nr,ipix2 + 1);
                    }
                    else { // // ysq > 0, ysqPrev < 0, fabs(dphi) > MY_EMSILON, ip_lo >= 0

    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo >= 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ip_hi: " << ip_hi << endl;

        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif                
                        VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_lo, ipix1 + ip_hi + 1);
                    }
                }
            }
            else if(iz < irminPrev || iz > irmaxPrev) {
                // intersect at exactly one pixel cell
                if(fabs(dphi) < MY_EPSILON) {
                    // get the pixel number of that cell,
                    // since dphi = 0, the pixel has the same phi value 
                    // as the ptg
                    int pix =  zphi2pix((double)z, ptg.phi);
                    // only add the point if it 's nonzero
                    if((hp)[pix] != 0) {
                        // VecAppendNonZeroPoint(pixset,hp, pix,pix+1);
                        pixset.push_back(pix);
                    }

                }
                else {
                    int ipix1, nr;
                    bool shifted;
                    // ipix1 is the start pixel of the current ring
                    // nr is the number of pixels in the ring
                    get_ring_info_small(iz, ipix1,nr,shifted);
                    double shift = 0.0;
                    if(shifted) {
                        double zdummy;
                        // shift is the phi value of the ipix1
                        pix2zphi(ipix1,zdummy,shift);
                    }
                    int ipix2 = ipix1 + nr - 1;
                    // e.g imagine the circle divided equally into four parts
                    // if the [0,pi/2) is the first one. and ptg.phi - dphi = pi/3, and shift = pi/4
                    // then ip_lo = 4/(2pi) * pi/12 = 0, since the left point lies in the first cell
                    // (note the first cell has pixel number ipix1 + 0)
                    // and for 3/4 * pi    ip_lo = 4/(2pi) * pi/2 = 1 yes
                    // int ip_lo = nearestInt(nr * inv_twopi * (ptg.phi - dphi - shift));
                    // int ip_hi = nearestInt(nr * inv_twopi * (ptg.phi + dphi - shift));
                    int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                    int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                    if (ip_hi >= nr) {
                        ip_lo -= nr;
                        ip_hi -= nr;
                    }
                    if (ip_lo < 0 ) {
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo < 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ipix1 + ip_lo: " << ipix1 + ip_lo << endl;
        cout << "ipix1 + ip_lo + n1: " << ipix1 + ip_lo + nr << endl;
        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif
                        VecAppendNonZeroPoint(pixset,hp,ipix1,ipix1 + ip_hi + 1);
                        if (ip_lo + nr <= ip_hi) {
                            // ip_lo + nr = ip_hi + 1
                            ip_lo = ip_hi + 1 - nr;
                        }
                        VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_lo + nr,ipix2 + 1);
                    }
                    else {
                        // pixset.append(ipix1 + ip_lo, ipix1 + ip_hi + 1);
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo >= 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ip_hi: " << ip_hi << endl;

        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif                
                        VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_lo, ipix1 + ip_hi + 1);
                    }
                }
            }
            else {  // ysqPrev >= 0 && iz also appears when we are dealing with thetaSmall
                    // i.e iz > izminPrev, iz < izmaxPrev

                myDouble yPrev = sqrt((double)ysqPrev);
                myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                //dphi cant be zero in this case
                if (fabs(dphi) < MY_EPSILON) {
            #ifdef QUERY_NON_CUMU_DEBUG
                    if (fabs(dphiPrev) < MY_EPSILON) {
                        cout << "ysq > 0, ysqPrev > 0, dphi ~ 0, dphiPrev ~ 0\n";
                        QUERY_DISC_BASIC_OUTPUT();
                    }
            #endif
                    if (fabs(dphiPrev) >= MY_EPSILON) {
                        cout << "\n\ndphi is zero when ysqPrev > 0\n";
                        cout << "ysq > 0, ysqPrev > 0, dphi ~ 0, dphiPrev > 0!!\n";
            #ifdef QUERY_NON_CUMU_DEBUG
                        QUERY_DISC_OUTPUT();
            #endif   
                        // exit(QUERY_DISC_ERROR);

                        int pix =  zphi2pix((double)z, ptg.phi);
                        // only add the point if it 's nonzero
                        if((hp)[pix] != 0) {
                            // VecAppendNonZeroPoint(pixset,hp, pix,pix+1);
                            pixset.push_back(pix);
                        }
                    }
                }
                else { 
                    // ysqPrev >= 0 && iz > izminPrev, iz < izmaxPrev
                    // fabs(dphi) > MY_EPSILON
                    int ipix1, nr;
                    bool shifted;
                    // ipix1 is the start pixel of the current ring
                    // nr is the number of pixels in the ring
                    get_ring_info_small(iz, ipix1,nr,shifted);
                    double shift = 0.0;
                    if(shifted) {
                        double zdummy;
                        // shift is the phi value of the ipix1
                        pix2zphi(ipix1,zdummy,shift);
                    }
                    int ipix2 = ipix1 + nr - 1;

                    int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                    int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                    int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                    int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                    if (ip_hi >= nr) {
                        ip_lo -= nr;
                        ip_hi -= nr;
                    }
                    if (ip_hiPrev >= nr) {
                        ip_loPrev -= nr;
                        ip_hiPrev -= nr;
                    }
                    if (ip_lo < 0 ) {
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo < 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ipix1 + ip_lo: " << ipix1 + ip_lo << endl;
        cout << "ipix1 + ip_lo + n1: " << ipix1 + ip_lo + nr << endl;
        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif
                        if (ip_loPrev < 0) {
                            if (ip_hiPrev > ip_hi || ip_loPrev < ip_lo) {
        #ifdef QUERY_NON_CUMU_DEBUG
                            cout << "\nip_lo, ip_loPrev < 0, ip_hiPrev > ip_hi or ip_loPrev < ip_lo\n";
                            QUERY_DISC_OUTPUT();
                            QUERY_DISC_IP_OUTPUT();
        #endif               
                            }  
                            else {
                                VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_hi + 1);
                                if (ip_lo + nr <= ip_hi) {
                                    // ip_lo + nr = ip_hi + 1
                                    ip_lo = ip_hi + 1 - nr;
                                }
                                VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_lo + nr, ipix1 + ip_loPrev + nr);
                            }
                        }
                        else { // ip_loPrev >= 0, ip_lo < 0
                            if (ip_hiPrev > ip_hi) {
        #ifdef QUERY_NON_CUMU_DEBUG
                            cout << "\nip_lo < 0, ip_loPrev > 0, ip_hiPrev > ip_hi\n";
                            QUERY_DISC_OUTPUT();
                            QUERY_DISC_IP_OUTPUT();
        #endif                               
                                // this means that  ip_hiPrev did not cover to the ipix2 , i.e, the last pixel at this ring
                                // but ip_hi passes the ipix2, so we need to add this segment from ipix1 + ip_hiPrev + 1 to that passes the ipix2
                                // the same is on the other side, from [ipix1 + ip_lo + nr (since ip_lo < 0), ipix1 + ip_loPrev)
                                // think of the pixel as the point lying in the fourth quadrant, ip_hi passes ipix2 so lies in the 1st-quadrant
                                // ip_hiPrev still lies in the fourth quadrant, will help you to understand.
                                VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev+1, ipix2 + 1);
                                VecAppendNonZeroPoint(pixset, hp, ipix1, ipix1 + ip_hi + 1);
                                if (ip_lo + nr <= ip_hi) {
                                    // ip_lo + nr = ip_hi + 1
                                    ip_lo = ip_hi + 1 - nr;
                                }
                                VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_lo + nr, ipix1 + ip_loPrev);
                            }
                            else { // ip_lo < 0, ip_loPrev >= 0, ip_hiPrev <= ip_hi
        #ifdef QUERY_NON_CUMU_DEBUG
                            cout << "\nip_lo < 0, ip_loPrev > 0, ip_hiPrev < ip_hi\n";
                            QUERY_DISC_OUTPUT();
                            QUERY_DISC_IP_OUTPUT();
        #endif              
                                // this will happen if ip_lo too negative, i.e, the ptg may be in the first quadrant,
                                // ip_loPrev still not pass the ipix1 to the fourth quadrant, but ip_lo does
                                // that means ip_hi, iphiPrev cant subtract nr, in this case, iphiPrev < iphi
                                // need to add three segments, [ipix1 + ip_hiPrev + 1, ipix1 + ip_hi],
                                // [ipix1 + ip_lo + nr, ipix2], and [ipix1, ipix1 + ip_loPrev), again,
                                // Imagine the pixel is in the first quadrant and (pixel + ip_lo (ip_lo negative) mod nr) + ipix1 is in the
                                // fourth quadrant will help
                                VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_hiPrev + 1,ipix1 + ip_hi + 1);
                                if (ip_lo + nr <= ip_hi) {
                                    // ip_lo + nr = ip_hi + 1
                                    ip_lo = ip_hi + 1 - nr;
                                }
                                VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_lo + nr,ipix2 + 1);
                                VecAppendNonZeroPoint(pixset,hp,ipix1, ipix1 + ip_loPrev);
                            }
                            
                        }   
                    }
                    else {  // ip_lo >= 0
                        // pixset.append(ipix1 + ip_lo, ipix1 + ip_hi + 1);
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo >= 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ip_hi: " << ip_hi << endl;

        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif              
                        if (ip_loPrev < 0 || ip_loPrev < ip_lo || ip_hiPrev > ip_hi) {
    #ifdef QUERY_NON_CUMU_DEBUG
                            cout << endl;
                            cout << "ip_loPrev < 0 || ip_loPrev < ip_lo || ip_hiPrev > ip_hi\n";
                            QUERY_DISC_OUTPUT();
                            QUERY_DISC_IP_OUTPUT();
    #endif
                        }
                        VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_lo, ipix1 + ip_loPrev);
                        VecAppendNonZeroPoint(pixset,hp,ipix1 + ip_hiPrev + 1, ipix1 + ip_hi + 1);
                    }
                }
                
            }
        }
    }
#ifdef HEALPIX_DEBUG
    cout << "the content of pixset is: " << endl;
    cout << pixset << endl;
#endif
    // southern cap
    if ((rlat2 > M_PI)) {
#ifdef QUERY_NON_CUMU_DEBUG
        cout << endl;
        cout << "in the southern cap\n";
        QUERY_DISC_BASIC_OUTPUT();
#endif
        // int startpix, ringpix,startpixPrev,ringpixPrev;
        // bool dummy;
        // get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
    
#ifdef ALL_ONE_DEBUG
    cout << "pixel: " << pixel << endl;
    cout << "startpix: " << startpix << endl;
    cout << "ringpix: " << ringpix << endl;
#endif

        // if rlat2Prev <= M_PI, it does not necessarily mean that the previous
        // did not cover pafts of the pole, so we still need to treat them as in the middle part
        // to find the previous cover at each ring
        // if ((rlat2Prev > M_PI)) {
    #ifdef QUERY_NON_CUMU_DEBU
            if (irmaxPrev < irmax) {
                cout << "in the southern cap, irmaxPrev < irmax!\n";
                QUERY_DISC_BASIC_OUTPUT();
            }
    #endif
            // VecAppendNonZeroPoint(pixset,hp,startpix, startpixPrev);
            // int irmaxPrev = min<int>(bottomRing, irmaxPrev);
            // since the behavior at iz == irmax is complicated,
            // the theta-disc may cover the whole ring, or not fully cover it
            // so we must consider irmax separately 
            // and for all ring from [irmax + 1, irmaxPrev], the new theta disc
            // will cover the whole pixels within this range,
            // so it's necessary to consider the case irmaxPrev > irmax, or irmaxPrev == irmax separately 
            // if (irmaxPrev > irmax) {
                for (int iz = irmaxPrev; iz > irmax; iz--) {
                    myDouble z = ring2z(iz);
                    myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
                    myDouble ysqPrev = 1 - xPrev * xPrev - z * z;
                    if (ysqPrev > 0) {
                        myDouble yPrev = sqrt((double)ysqPrev);
                        myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                        int ipix1, nr;
                        bool shifted;
                        // ipix1 is the start pixel of the current ring
                        // nr is the number of pixels in the ring
                        get_ring_info_small(iz, ipix1,nr,shifted);
                        double shift = 0.0;
                        if(shifted) {
                            double zdummy;
                            // shift is the phi value of the ipix1
                            pix2zphi(ipix1,zdummy,shift);
                        }
                        int ipix2 = ipix1 + nr - 1;

                        int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                        int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                        if (ip_hiPrev >= nr) {
                            ip_loPrev -= nr;
                            ip_hiPrev -= nr;
                        }
                        if (ip_loPrev < 0) {
                            if (ip_loPrev + nr > ip_hiPrev) {
                                // the previous theta did not cover the whole ring
                                // fill the blank
                                VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_loPrev + nr);
                            }
                        }
                        else {
                            VecAppendNonZeroPoint(pixset, hp, ipix1,ipix1 + ip_loPrev);
                            VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev + 1, ipix2 + 1);
                        }                    
                    }
                    else {
                        // it means that the ring at iz is now fully covered by the previous theta
                        // so no need to fill the blank
    #ifdef QUERY_NON_CUMU_DEBUG                    
                        cout << "the southern cap, yspPrev < 0\n";
                        cout << "iz: " << iz << endl;
                        cout << "irmax: " << irmax << endl;
                        cout << "irmaxPrev: " << irmaxPrev << endl;
                        QUERY_DISC_BASIC_OUTPUT();
    #endif
                    }                
                }
                // the case irmax 
                myDouble z = ring2z(irmax);
                myDouble x = (cos((double)thetaLarge) - z * z0) * invSinT;
                myDouble ysq = 1 - x * x - z * z;
                if (ysq > 0) {
                    // irmax case already covered
                    // get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
                    // VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);           
                }
                else {
                    // add the current ring - pixels covered by the previous 
                    myDouble z = ring2z(irmax);
                    myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
                    myDouble ysqPrev = 1 - xPrev * xPrev - z * z;   
    #ifdef QUERY_NON_CUMU_DEBUG
                    if (ysqPrev < 0) {
                        cout << "southern cap, irmaxPrev > irmax, at irmax: ysqPrev < 0\n";
                        QUERY_DISC_BASIC_OUTPUT();
                    }
    #endif
                    if (ysqPrev >= 0){
                        myDouble yPrev = sqrt((double)ysqPrev);
                        myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                        int ipix1, nr;
                        bool shifted;
                        // ipix1 is the start pixel of the current ring
                        // nr is the number of pixels in the ring
                        get_ring_info_small(irmax, ipix1,nr,shifted);
                        double shift = 0.0;
                        if(shifted) {
                            double zdummy;
                            // shift is the phi value of the ipix1
                            pix2zphi(ipix1,zdummy,shift);
                        }
                        int ipix2 = ipix1 + nr - 1;

                        int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                        int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                        if (ip_hiPrev >= nr) {
                            ip_loPrev -= nr;
                            ip_hiPrev -= nr;
                        }
                        if (ip_loPrev < 0) {
                            if (ip_loPrev + nr > ip_hiPrev) {
                                // the previous theta did not cover the whole ring
                                // fill the blank
                                VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_loPrev + nr);
                            }
                        }
                        else {
                            VecAppendNonZeroPoint(pixset, hp, ipix1,ipix1 + ip_loPrev);
                            VecAppendNonZeroPoint(pixset, hp, ipix1 + ip_hiPrev + 1, ipix2 + 1);
                        }     
                    }
                } 
    } 
#ifdef QUERY_NON_CUMU_DEBUG
    cout << "\nNon cumu version: \n";
    QUERY_DISC_BASIC_OUTPUT();
    cout << "pixset: \n";
    for (int i = 0; i < (int)pixset.size(); i++) {
        cout << pixset[i] << "\t";
    }
    cout << endl;
#endif

    if (filename != "no") {
        std::ofstream outFile;
        outFile.open(filename);
        for (vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter++) {
            outFile << *iter << endl;
        }
        outFile.close();
    }
}            

countType Healpix_Correlation::query_disc_pixel_NonCumuInternalCount(int pixel, const myDouble & thetaLarge, 
    const myDouble & thetaSmall,
    const Healpix_Correlation & hp) const
{

    if((*this)[pixel] < MY_EPSILON) {
        return 0;
    }
#ifdef QUERY_ZERO
    
    #ifndef signalMode
        if (thetaSmall < (MY_EPSILON) && thetaLarge < (deltaTheta - MY_EPSILON)) {
#ifdef QUERY_NON_CUMU_DEBUG
        cout << "\nthetaSmall < epsilon && thetaLarge < deltaTheta - epsilon\n";
        cout << "dT too large, or bin width too small\n";
        QUERY_DISC_BASIC_OUTPUT();
        cout << "deltaTheta: " << deltaTheta << endl;
#endif
    #ifdef MULTIPLY_PIXEL_VALUE
            // return (hp)[pixel] * (hp)[pixel];
    #else
            return (hp)[pixel];
    #endif
        // means that the disc is small so actually we consider a single cell
        // the neighbors are actually in one "room"
            return (hp)[pixel];
        }
    #else
        if(thetaSmall < deltaTheta && thetaLarge < deltaTheta) {
            pixset.push_back(pixel);
            return;
        }
    #endif
#endif

    // this case appears, only when thetaSmall = 0
    // calculate as in the cum case, except, the pixel itself needs to be excluded
    if (thetaSmall < (deltaTheta - MY_EPSILON)) {
        return query_disc_count(pixel,thetaLarge,hp,pixel);
    }
    countType pointCount = 0;
    pointing ptg = pix2ang(pixel);
    myDouble rlat1 = ptg.theta - thetaLarge;
    myDouble rlat1Prev = ptg.theta - thetaSmall;
    myDouble zmax = cos((double)rlat1);
    int irmin = ring_above((double)zmax) + 1;


   
    myDouble zmaxPrev = cos((double)rlat1Prev);
    // the ring where the point with zmax as z-coord locates
    
    int irminPrev = ring_above((double)zmaxPrev) + 1;
    myDouble z0 = cos(ptg.theta);
    myDouble invSinT = 1. / sqrt((1 - z0) * (1 + z0));

    // north pole in the disc
    if ((rlat1 <= 0) && (irmin > 1)) {
        int startpix, ringpix,startpixPrev,ringpixPrev;
        bool dummy;


        get_ring_info_small(irmin - 1,startpix,ringpix,dummy);
        if ((rlat1Prev <= 0) || irminPrev == 1) {
            get_ring_info_small(irminPrev - 1, startpixPrev,ringpixPrev, dummy);

            // since the result of ring_above may vary if the point is above
            // the pixel of that cell or not, it's better to set irmin to irminPrev
            // incase that irminPrev did not actually cover the entire string
            // e.g 
            /*
                dT = 0.1, all one map
                pixel: 0
                thetaLarge: 0.5
                thetaSmall: 0.4  
                will ignore 19 which is at ring 3 = irminPrev
                in this case irminPrev did not fully cover the ring
            */
           for (int iz = irminPrev; iz < irmin; iz++) {
                myDouble z = ring2z(iz);
                myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
                myDouble ysqPrev = 1 - xPrev * xPrev - z * z;
                if (ysqPrev > 0) {
                    myDouble yPrev = sqrt((double)ysqPrev);
                    myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                    int ipix1, nr;
                    bool shifted;
                    // ipix1 is the start pixel of the current ring
                    // nr is the number of pixels in the ring
                    get_ring_info_small(iz, ipix1,nr,shifted);
                    double shift = 0.0;
                    if(shifted) {
                        double zdummy;
                        // shift is the phi value of the ipix1
                        pix2zphi(ipix1,zdummy,shift);
                    }
                    int ipix2 = ipix1 + nr - 1;

                    int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                    int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                    if (ip_hiPrev >= nr) {
                        ip_loPrev -= nr;
                        ip_hiPrev -= nr;
                    }
                    if (ip_loPrev < 0) {
                        if (ip_loPrev + nr > ip_hiPrev) {
                            // the previous theta did not cover the whole ring
                            // fill the blank
                            pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_loPrev + nr);
                        }
                    }
                    else {
                        pointCount += IntervallNonZeroPoint(hp, ipix1,ipix1 + ip_loPrev);
                        pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev + 1, ipix2 + 1);
                    }                    
                }
                else {
                    // it means this ring has already been covered by the previous theta-disc
#ifdef QUERY_NON_CUMU_DEBUG                    
                    cout << "the northern cap, yspPrev < 0\n";
                    cout << "iz: " << iz << endl;
                    cout << "irmin: " << irmin << endl;
                    cout << "irminPrev: " << irminPrev << endl;
                    QUERY_DISC_BASIC_OUTPUT();
#endif
                }
            }
            // need to consider the case irmin separatedly
            // if irmin at the middle earth loop, and ysq < 0, then
            // it will not be considered there
                myDouble z = ring2z(irmin);
                myDouble x = (cos((double)thetaLarge) - z * z0) * invSinT;
                myDouble ysq = 1 - x * x - z * z;
                if (ysq > 0) {
                    // irmin case will be covered below
                    // get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
                    // VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);           
                }
                else {
                    // add the current ring - pixels covered by the previous 
                    myDouble z = ring2z(irmin);
                    myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
                    myDouble ysqPrev = 1 - xPrev * xPrev - z * z;   
    #ifdef QUERY_NON_CUMU_DEBUG
                    if (ysqPrev < 0) {
                        cout << "northern cap, at irmin: ysqPrev < 0\n";
                        QUERY_DISC_BASIC_OUTPUT();
                    }
    #endif
                    if (ysqPrev >= 0){
                        myDouble yPrev = sqrt((double)ysqPrev);
                        myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                        int ipix1, nr;
                        bool shifted;
                        // ipix1 is the start pixel of the current ring
                        // nr is the number of pixels in the ring
                        get_ring_info_small(irmin, ipix1,nr,shifted);
                        double shift = 0.0;
                        if(shifted) {
                            double zdummy;
                            // shift is the phi value of the ipix1
                            pix2zphi(ipix1,zdummy,shift);
                        }
                        int ipix2 = ipix1 + nr - 1;

                        int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                        int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                        if (ip_hiPrev >= nr) {
                            ip_loPrev -= nr;
                            ip_hiPrev -= nr;
                        }
                        if (ip_loPrev < 0) {
                            if (ip_loPrev + nr > ip_hiPrev) {
                                // the previous theta did not cover the whole ring
                                // fill the blank
                                pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_loPrev + nr);
                            }
                        }
                        else {
                            pointCount += IntervallNonZeroPoint(hp, ipix1,ipix1 + ip_loPrev);
                            pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev + 1, ipix2 + 1);
                        }     
                    }


                }
        }
#ifdef ALL_ONE_DEBUG
    cout << "pixel: " << pixel << endl;
    cout << "startpix: " << startpix << endl;
    cout << "ringpix: " << ringpix << endl;
#endif
    }   
    myDouble rlat2 = ptg.theta + thetaLarge;
    myDouble rlat2Prev = ptg.theta + thetaSmall;
    myDouble zmin = cos((double)rlat2);
    myDouble zminPrev = cos((double)rlat2Prev);
    int irmax = ring_above((double)zmin) + 1;
    int irmaxPrev = ring_above((double)zminPrev) + 1;
    irmax = min<int>(irmax, bottomRing);
    irmaxPrev = min<int>(irmaxPrev,bottomRing);
#ifdef QUERY_NON_CUMU_DEBUG    
    if (irmaxPrev > irmax) {

        cout << "\n irmaxPrev > irmax\n";
        QUERY_DISC_BASIC_OUTPUT();
        cout << "irmax: " << irmax << endl;
        cout << "irmaxPrev: " << irmaxPrev << endl;

    }
#endif

#ifdef HEALPIX_DEBUG
    cout << "inside the query_disc_pixel_Internal for rangeset version, begins loop\n";
#endif        


    for (int iz = irmin; iz <= irmax; iz+= 1) {
        myDouble z = ring2z(iz);
        myDouble x = (cos((double)thetaLarge) - z * z0) * invSinT;
        myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
        myDouble ysq = 1 - x * x - z * z;
        myDouble ysqPrev = 1 - xPrev * xPrev - z * z;
        if(ysq > 0) {
            myDouble y = sqrt((double)ysq);
            // note dphi >= 0
            myDouble dphi = atan2((double)y,(double)x);
            if (ysqPrev < 0) {
                // intersect at exactly one pixel cell
                if(fabs(dphi) < MY_EPSILON) {
                    // get the pixel number of that cell,
                    // since dphi = 0, the pixel has the same phi value 
                    // as the ptg
                    int pix =  zphi2pix((double)z, ptg.phi);
                    // only add the point if it 's nonzero
                    if((hp)[pix] != 0) {
                        // VecAppendNonZeroPoint(pixset,hp, pix,pix+1);
                        pointCount += (hp)[pix];
                    }

                }
                else { // ysq > 0, ysqPrev < 0, fabs(dphi) > MY_EMSILON
                    int ipix1, nr;
                    bool shifted;
                    // ipix1 is the start pixel of the current ring
                    // nr is the number of pixels in the ring
                    get_ring_info_small(iz, ipix1,nr,shifted);
                    double shift = 0.0;
                    if(shifted) {
                        double zdummy;
                        // shift is the phi value of the ipix1
                        pix2zphi(ipix1,zdummy,shift);
                    }
                    int ipix2 = ipix1 + nr - 1;
                    // e.g imagine the circle divided equally into four parts
                    // if the [0,pi/2) is the first one. and ptg.phi - dphi = pi/3, and shift = pi/4
                    // then ip_lo = 4/(2pi) * pi/12 = 0, since the left point lies in the first cell
                    // (note the first cell has pixel number ipix1 + 0)
                    // and for 3/4 * pi    ip_lo = 4/(2pi) * pi/2 = 1 yes
                    // int ip_lo = nearestInt(nr * inv_twopi * (ptg.phi - dphi - shift));
                    // int ip_hi = nearestInt(nr * inv_twopi * (ptg.phi + dphi - shift));
                    int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                    int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                    if (ip_hi >= nr) {
                        ip_lo -= nr;
                        ip_hi -= nr;
                    }
                    if (ip_lo < 0 ) {
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo < 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ipix1 + ip_lo: " << ipix1 + ip_lo << endl;
        cout << "ipix1 + ip_lo + n1: " << ipix1 + ip_lo + nr << endl;
        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif
                        pointCount += IntervallNonZeroPoint(hp,ipix1,ipix1 + ip_hi + 1);
                        if (ip_lo + nr <= ip_hi) {
                            // ip_lo + nr = ip_hi + 1
                            ip_lo = ip_hi + 1 - nr;
                        }
                        pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_lo + nr,ipix2 + 1);
                    }
                    else { // // ysq > 0, ysqPrev < 0, fabs(dphi) > MY_EMSILON, ip_lo >= 0

    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo >= 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ip_hi: " << ip_hi << endl;

        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif                
                        pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_lo, ipix1 + ip_hi + 1);
                    }
                }
            }
            else if(iz < irminPrev || iz > irmaxPrev) {
                // intersect at exactly one pixel cell
                if(fabs(dphi) < MY_EPSILON) {
                    // get the pixel number of that cell,
                    // since dphi = 0, the pixel has the same phi value 
                    // as the ptg
                    int pix =  zphi2pix((double)z, ptg.phi);
                    // only add the point if it 's nonzero
                    if((hp)[pix] != 0) {
                        // VecAppendNonZeroPoint(pixset,hp, pix,pix+1);
                        pointCount += (hp)[pix];
                    }

                }
                else {
                    int ipix1, nr;
                    bool shifted;
                    // ipix1 is the start pixel of the current ring
                    // nr is the number of pixels in the ring
                    get_ring_info_small(iz, ipix1,nr,shifted);
                    double shift = 0.0;
                    if(shifted) {
                        double zdummy;
                        // shift is the phi value of the ipix1
                        pix2zphi(ipix1,zdummy,shift);
                    }
                    int ipix2 = ipix1 + nr - 1;
                    // e.g imagine the circle divided equally into four parts
                    // if the [0,pi/2) is the first one. and ptg.phi - dphi = pi/3, and shift = pi/4
                    // then ip_lo = 4/(2pi) * pi/12 = 0, since the left point lies in the first cell
                    // (note the first cell has pixel number ipix1 + 0)
                    // and for 3/4 * pi    ip_lo = 4/(2pi) * pi/2 = 1 yes
                    // int ip_lo = nearestInt(nr * inv_twopi * (ptg.phi - dphi - shift));
                    // int ip_hi = nearestInt(nr * inv_twopi * (ptg.phi + dphi - shift));
                    int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                    int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                    if (ip_hi >= nr) {
                        ip_lo -= nr;
                        ip_hi -= nr;
                    }
                    if (ip_lo < 0 ) {
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo < 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ipix1 + ip_lo: " << ipix1 + ip_lo << endl;
        cout << "ipix1 + ip_lo + n1: " << ipix1 + ip_lo + nr << endl;
        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif
                        pointCount += IntervallNonZeroPoint(hp,ipix1,ipix1 + ip_hi + 1);
                        if (ip_lo + nr <= ip_hi) {
                            // ip_lo + nr = ip_hi + 1
                            ip_lo = ip_hi + 1 - nr;
                        }
                        pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_lo + nr,ipix2 + 1);
                    }
                    else {
                        // pixset.append(ipix1 + ip_lo, ipix1 + ip_hi + 1);
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo >= 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ip_hi: " << ip_hi << endl;

        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif                
                        pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_lo, ipix1 + ip_hi + 1);
                    }
                }
            }
            else {  // ysqPrev >= 0 && iz also appears when we are dealing with thetaSmall
                    // i.e iz > izminPrev, iz < izmaxPrev

                myDouble yPrev = sqrt((double)ysqPrev);
                myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                //dphi cant be zero in this case
                if (fabs(dphi) < MY_EPSILON) {
                    if (fabs(dphiPrev) >= MY_EPSILON) {
                        cout << "\n\ndphi is zero when ysqPrev > 0\n";
                        cout << "ysq > 0, ysqPrev > 0, dphi ~ 0, dphiPrev > 0!!\n";
                        int pix =  zphi2pix((double)z, ptg.phi);
                        // only add the point if it 's nonzero
                        if((hp)[pix] != 0) {
                            // VecAppendNonZeroPoint(pixset,hp, pix,pix+1);
                            pointCount += (hp)[pix];
                        }
                    }
                }
                else { 
                    // ysqPrev >= 0 && iz > izminPrev, iz < izmaxPrev
                    // fabs(dphi) > MY_EPSILON
                    int ipix1, nr;
                    bool shifted;
                    // ipix1 is the start pixel of the current ring
                    // nr is the number of pixels in the ring
                    get_ring_info_small(iz, ipix1,nr,shifted);
                    double shift = 0.0;
                    if(shifted) {
                        double zdummy;
                        // shift is the phi value of the ipix1
                        pix2zphi(ipix1,zdummy,shift);
                    }
                    int ipix2 = ipix1 + nr - 1;

                    int ip_lo = round(nr * inv_twopi * (ptg.phi - dphi - shift));
                    int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                    int ip_hi = round(nr * inv_twopi * (ptg.phi + dphi - shift));
                    int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                    if (ip_hi >= nr) {
                        ip_lo -= nr;
                        ip_hi -= nr;
                    }
                    if (ip_hiPrev >= nr) {
                        ip_loPrev -= nr;
                        ip_hiPrev -= nr;
                    }
                    if (ip_lo < 0 ) {
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo < 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ipix1 + ip_lo: " << ipix1 + ip_lo << endl;
        cout << "ipix1 + ip_lo + n1: " << ipix1 + ip_lo + nr << endl;
        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif
                        if (ip_loPrev < 0) {
                            if (ip_hiPrev > ip_hi || ip_loPrev < ip_lo) {
        #ifdef QUERY_NON_CUMU_DEBUG
                            cout << "\nip_lo, ip_loPrev < 0, ip_hiPrev > ip_hi or ip_loPrev < ip_lo\n";
                            QUERY_DISC_OUTPUT();
                            QUERY_DISC_IP_OUTPUT();
        #endif               
                            }  
                            else {
                                pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_hi + 1);
                                if (ip_lo + nr <= ip_hi) {
                                    // ip_lo + nr = ip_hi + 1
                                    ip_lo = ip_hi + 1 - nr;
                                }
                                pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_lo + nr, ipix1 + ip_loPrev + nr);
                            }
                        }
                        else { // ip_loPrev >= 0, ip_lo < 0
                            if (ip_hiPrev > ip_hi) {
        #ifdef QUERY_NON_CUMU_DEBUG
                            cout << "\nip_lo < 0, ip_loPrev > 0, ip_hiPrev > ip_hi\n";
                            QUERY_DISC_OUTPUT();
                            QUERY_DISC_IP_OUTPUT();
        #endif                               
                                // this means that  ip_hiPrev did not cover to the ipix2 , i.e, the last pixel at this ring
                                // but ip_hi passes the ipix2, so we need to add this segment from ipix1 + ip_hiPrev + 1 to that passes the ipix2
                                // the same is on the other side, from [ipix1 + ip_lo + nr (since ip_lo < 0), ipix1 + ip_loPrev)
                                // think of the pixel as the point lying in the fourth quadrant, ip_hi passes ipix2 so lies in the 1st-quadrant
                                // ip_hiPrev still lies in the fourth quadrant, will help you to understand.
                                pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev+1, ipix2 + 1);
                                pointCount += IntervallNonZeroPoint(hp, ipix1, ipix1 + ip_hi + 1);
                                if (ip_lo + nr <= ip_hi) {
                                    // ip_lo + nr = ip_hi + 1
                                    ip_lo = ip_hi + 1 - nr;
                                }
                                pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_lo + nr, ipix1 + ip_loPrev);
                            }
                            else { // ip_lo < 0, ip_loPrev >= 0, ip_hiPrev <= ip_hi
        #ifdef QUERY_NON_CUMU_DEBUG
                            cout << "\nip_lo < 0, ip_loPrev > 0, ip_hiPrev < ip_hi\n";
                            QUERY_DISC_OUTPUT();
                            QUERY_DISC_IP_OUTPUT();
        #endif              
                                // this will happen if ip_lo too negative, i.e, the ptg may be in the first quadrant,
                                // ip_loPrev still not pass the ipix1 to the fourth quadrant, but ip_lo does
                                // that means ip_hi, iphiPrev cant subtract nr, in this case, iphiPrev < iphi
                                // need to add three segments, [ipix1 + ip_hiPrev + 1, ipix1 + ip_hi],
                                // [ipix1 + ip_lo + nr, ipix2], and [ipix1, ipix1 + ip_loPrev), again,
                                // Imagine the pixel is in the first quadrant and (pixel + ip_lo (ip_lo negative) mod nr) + ipix1 is in the
                                // fourth quadrant will help
                                pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_hiPrev + 1,ipix1 + ip_hi + 1);
                                if (ip_lo + nr <= ip_hi) {
                                    // ip_lo + nr = ip_hi + 1
                                    ip_lo = ip_hi + 1 - nr;
                                }
                                pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_lo + nr,ipix2 + 1);
                                pointCount += IntervallNonZeroPoint(hp,ipix1, ipix1 + ip_loPrev);
                            }
                            
                        }   
                    }
                    else {  // ip_lo >= 0
                        // pixset.append(ipix1 + ip_lo, ipix1 + ip_hi + 1);
    #ifdef ALL_ONE_DEBUG
        cout << "in the ip_lo >= 0: condition" << endl;
        cout << "pixel: " << pixel << endl;
        cout << "ipix1: " << ipix1 << endl;
        cout << "ip_lo: " << ip_lo << endl;
        cout << "ip_hi: " << ip_hi << endl;

        cout << "ipix1 + ip_hi + 1" << ipix1 + ip_hi + 1 << endl;
    #endif              
                        pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_lo, ipix1 + ip_loPrev);
                        pointCount += IntervallNonZeroPoint(hp,ipix1 + ip_hiPrev + 1, ipix1 + ip_hi + 1);
                    }
                }
                
            }
        }
    }
#ifdef HEALPIX_DEBUG
    cout << "the content of pixset is: " << endl;
    cout << pixset << endl;
#endif
    // southern cap
    if ((rlat2 > M_PI)) {
#ifdef QUERY_NON_CUMU_DEBUG
        cout << endl;
        cout << "in the southern cap\n";
        QUERY_DISC_BASIC_OUTPUT();
#endif
        // int startpix, ringpix,startpixPrev,ringpixPrev;
        // bool dummy;
        // get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
    
#ifdef ALL_ONE_DEBUG
    cout << "pixel: " << pixel << endl;
    cout << "startpix: " << startpix << endl;
    cout << "ringpix: " << ringpix << endl;
#endif

    #ifdef QUERY_NON_CUMU_DEBUG
            // get_ring_info_small(irmaxPrev + 1, startpixPrev, ringpixPrev, dummy);
            // if (startpix > startpixPrev) {

            //     cout << endl;
            //     cout << "(rlat2Prev > M_PI) && (irmaxPrev < bottomRing)\n";
            //     QUERY_DISC_BASIC_OUTPUT();
                
            // }
            if (irmaxPrev < irmax) {
                cout << "in the southern cap, irmaxPrev < irmax!\n";
                QUERY_DISC_BASIC_OUTPUT();
            }
    #endif
            // VecAppendNonZeroPoint(pixset,hp,startpix, startpixPrev);
            // int irmaxPrev = min<int>(bottomRing, irmaxPrev);
            // since the behavior at iz == irmax is complicated,
            // the theta-disc may cover the whole ring, or not fully cover it
            // so we must consider irmax separately 
            // and for all ring from [irmax + 1, irmaxPrev], the new theta disc
            // will cover the whole pixels within this range,
            // so it's necessary to consider the case irmaxPrev > irmax, or irmaxPrev == irmax separately 
        for (int iz = irmaxPrev; iz > irmax; iz--) {
            myDouble z = ring2z(iz);
            myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
            myDouble ysqPrev = 1 - xPrev * xPrev - z * z;
            if (ysqPrev > 0) {
                myDouble yPrev = sqrt((double)ysqPrev);
                myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                int ipix1, nr;
                bool shifted;
                // ipix1 is the start pixel of the current ring
                // nr is the number of pixels in the ring
                get_ring_info_small(iz, ipix1,nr,shifted);
                double shift = 0.0;
                if(shifted) {
                    double zdummy;
                    // shift is the phi value of the ipix1
                    pix2zphi(ipix1,zdummy,shift);
                }
                int ipix2 = ipix1 + nr - 1;

                int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                if (ip_hiPrev >= nr) {
                    ip_loPrev -= nr;
                    ip_hiPrev -= nr;
                }
                if (ip_loPrev < 0) {
                    if (ip_loPrev + nr > ip_hiPrev) {
                        // the previous theta did not cover the whole ring
                        // fill the blank
                        pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_loPrev + nr);
                    }
                }
                else {
                    pointCount += IntervallNonZeroPoint(hp, ipix1,ipix1 + ip_loPrev);
                    pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev + 1, ipix2 + 1);
                }                    
            }
            else {
                // it means that the ring at iz is now fully covered by the previous theta
                // so no need to fill the blank
#ifdef QUERY_NON_CUMU_DEBUG                    
                cout << "the southern cap, yspPrev < 0\n";
                cout << "iz: " << iz << endl;
                cout << "irmax: " << irmax << endl;
                cout << "irmaxPrev: " << irmaxPrev << endl;
                QUERY_DISC_BASIC_OUTPUT();
#endif
            }                
        }
        // the case irmax 
        myDouble z = ring2z(irmax);
        myDouble x = (cos((double)thetaLarge) - z * z0) * invSinT;
        myDouble ysq = 1 - x * x - z * z;
        if (ysq > 0) {
            // irmax case already covered
            // get_ring_info_small(irmax + 1,startpix,ringpix,dummy);
            // VecAppendNonZeroPoint(pixset,hp,startpix, bottomPixel + 1,ignore);           
        }
        else {
            // add the current ring - pixels covered by the previous 
            myDouble z = ring2z(irmax);
            myDouble xPrev = (cos((double)thetaSmall) - z * z0)* invSinT;
            myDouble ysqPrev = 1 - xPrev * xPrev - z * z;   
#ifdef QUERY_NON_CUMU_DEBUG
            if (ysqPrev < 0) {
                cout << "southern cap, irmaxPrev > irmax, at irmax: ysqPrev < 0\n";
                QUERY_DISC_BASIC_OUTPUT();
            }
#endif
            if (ysqPrev >= 0){
                myDouble yPrev = sqrt((double)ysqPrev);
                myDouble dphiPrev = atan2((double)yPrev,(double)xPrev);
                int ipix1, nr;
                bool shifted;
                // ipix1 is the start pixel of the current ring
                // nr is the number of pixels in the ring
                get_ring_info_small(irmax, ipix1,nr,shifted);
                double shift = 0.0;
                if(shifted) {
                    double zdummy;
                    // shift is the phi value of the ipix1
                    pix2zphi(ipix1,zdummy,shift);
                }
                int ipix2 = ipix1 + nr - 1;

                int ip_loPrev = round(nr * inv_twopi * (ptg.phi - dphiPrev - shift));
                int ip_hiPrev = round(nr * inv_twopi * (ptg.phi + dphiPrev - shift));
                if (ip_hiPrev >= nr) {
                    ip_loPrev -= nr;
                    ip_hiPrev -= nr;
                }
                if (ip_loPrev < 0) {
                    if (ip_loPrev + nr > ip_hiPrev) {
                        // the previous theta did not cover the whole ring
                        // fill the blank
                        pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev + 1, ipix1 + ip_loPrev + nr);
                    }
                }
                else {
                    pointCount += IntervallNonZeroPoint(hp, ipix1,ipix1 + ip_loPrev);
                    pointCount += IntervallNonZeroPoint(hp, ipix1 + ip_hiPrev + 1, ipix2 + 1);
                }     
            }
        }
    }

#ifdef MULTIPLY_PIXEL_VALUE
    // pointCount *= (hp)[pixel];
#endif
    return pointCount;
}    

void Healpix_Correlation::query_disc_pixel_Bins(int pixel, 
    const std::vector<myDouble>& bins,
    std::vector<corrType> & neighborCount,
    const Healpix_Correlation & hp)
{
    int len = bins.size();
    vector<int> pixset;
    if ((*this)[pixel] < MY_EPSILON) {
        return;
    }
    for (int index = 1; index < len; index++) {
        countType pair = 0;
#ifdef QUERY_COUNT_DIRECTLY
        pair = query_disc_pixel_NonCumuInternalCount(pixel,bins[index],bins[index-1],hp);
#else
        query_disc_pixel_NonCumuInternalVec(pixel, bins[index],
            bins[index-1],pixset,hp);

    #ifdef OUTPUT_PIX_NEIGHBOR_COUNT
        myDouble theta = bins[index];
        if (fabs(theta - TO_PRINT_THETA) < 0.01)
            cout << "*******\npixel: " << pixel << "\ttheta: " << theta << endl;
    #endif
        for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
            pair += (hp)[*iter];
    #ifdef OUTPUT_PIX_NEIGHBOR_COUNT
        if (fabs(theta - TO_PRINT_THETA) < 0.01)
            cout << *iter << "\t";
    #endif
        }
    #ifdef OUTPUT_PIX_NEIGHBOR_COUNT
        if (fabs(theta - TO_PRINT_THETA) < 0.01)
            cout << endl << "neighbor: " << pair << endl;
    #endif             
        pixset.clear();
#endif

    #ifdef MULTIPLY_PIXEL_VALUE
        pair *= (*this)[pixel];
    #endif
        neighborCount[index - 1] += pair;
    }
}

void Healpix_Correlation::query_disc_pixel_Bins(int pixel, const std::vector<myDouble>& bins,
    std::vector<corrType> & neighborCount,
    const Healpix_Correlation & hp,
    const Healpix_Correlation & origin,
    vector<std::string> fileNames,
    vector<myDouble> thetas)
{

    cout << "in query_disc_pixel_Bins, need to add origin !\n";
    exit(TODO_FUNCTION);
#ifndef QUERY_COUNT_DIRECTLY
    int len = bins.size();
    vector<int> pixset;
    int len1 = fileNames.size();
    int len2 = thetas.size();
    if (len1 != len2) {
        cout << "File size not matching in query_disc_pixel_Bins\n";
        exit(FILE_NUM_NOT_MATCHING);
    }
    bool noMorePlot = false;
    int _thetasIndex = 0;
    myDouble thetaPlot;
    std::string filename;
    if (len1 == 0) {
        noMorePlot = true;
    }
    else {
        thetaPlot = thetas[_thetasIndex];
    }
    for (int index = 1; index < len; index++) {
        countType pair = 0;
        if (!noMorePlot && fabs(thetaPlot - bins[index]) < 0.01) {
            filename = fileNames[_thetasIndex++];
            query_disc_pixel_NonCumuInternalVec(pixel,bins[index],
                bins[index-1], pixset, *this, filename);
            if (_thetasIndex > len1) {
                noMorePlot = true;
            }
            else {
                thetaPlot = thetas[_thetasIndex];
            }
        }
        else {
            query_disc_pixel_NonCumuInternalVec(pixel, bins[index],
                bins[index-1],pixset,*this);
        }

        for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
            pair += (*this)[*iter];
        }
        neighborCount[index - 1] += pair;
        pixset.clear();
    }
#else
    if (!fileNames.empty()) {
        cout << "outputting result to file is not supported for query_disc_pixel_Bins\n";
        cout << "reset the QUERY_COUNT_DIRECTLY to enable the function\n";
    }
    int len = bins.size();
    for (int index = 1; index < len; index++) {
        neighborCount[index - 1] = query_disc_pixel_NonCumuInternalCount(pixel,bins[index],bins[index-1],*this);
    }
#endif
}   

void Healpix_Correlation::update() 
{
    updateBoundPixel();
    updateDeltaTheta();
}

void Healpix_Correlation::updateBoundPixel() 
{
    int ring = pix2ring(zphi2pix(cos((double)hThetaEnd), 0));
    bottomRing = ring;
    ring = pix2ring(zphi2pix(cos((double)hThetaBeg),0));
    topRing = ring;
    int startpix;
    int ringpix;
    bool dummy;
    get_ring_info_small(bottomRing,startpix,ringpix,dummy);
    bottomPixel = startpix + ringpix - 1;
    get_ring_info_small(topRing, startpix, ringpix, dummy);
    topPixel = startpix;
    boundSet = true;
#ifdef HEALPIX_DEBUG
    cout << "In the updateBottomPixel, bottomRing: " << bottomRing << ",  bottomPixel: " << bottomPixel << endl;
#endif
}

void Healpix_Correlation::updateDeltaTheta() {
    // very naive calculation
    // from the north pole ring to the equator there are 2 * N_side rings
    int N_side = Nside();
    if (N_side == 0) {
        cout << "healpix map has no size\n" << endl;
        exit(EMPTY_HEALPIX_MAP);
    }
    deltaTheta = pi/(4 * N_side);
#ifdef HEALPIX_DEBUG
    cout << "In the updateDeltaTheta, the Nside: " << N_side << ",  deltaTheta: " << deltaTheta << endl;
#endif 
}

void Healpix_Correlation::resetPixValue(int value, int _threadNum) 
{
    int nsideTemp = Nside();
    int endPixel = 12 * nsideTemp * nsideTemp - 1;    
    if(hPhiEnd - hPhiBeg > twopi - MY_EPSILON) {
        #ifdef PARALLEL
            pointsCount = 0;
            int threadCount = _threadNum;
            countType stepSize = ceil(bottomPixel/threadCount);
            std::vector<thread> vt;
            std::vector<int> index;
            index.resize(threadCount);

            for (int pix = 0; pix < topPixel; pix += 1) {
                (*this)[pix] = 0;
            }
            for (int i = 0; i < threadCount; i+= 1) {
                typeIndex begin, end;
                begin = i * stepSize;
                if(i != threadCount - 1) {
                    end = (i + 1) * stepSize;
                    index[i] = end;
                }
                else {
                    end = bottomPixel + 1;
                    index[i] = end;
                }
        #ifdef RESET_VALUE_DEBUG
                // cout << "in the caller: \n";
                // cout << "begin: " << begin << endl;
                // cout << "end: " << end << endl;
        #endif
                vt.push_back(thread(std::bind(&Healpix_Correlation::parallelSetPixValue,this,begin,end,std::ref(pointsCount),value)));

        #ifdef JOIN_THREAD_MOVE
                vt[i].join();
        #endif
            }
            // for (int i = 0; i < threadCount; i++) {
            //     typeIndex begin,end;
            //     if (i == 0) {
            //         begin = 0;
            //         end = index[0];
            //     }
            //     else {
            //         begin = index[i-1];
            //         end = index[i];
            //     }
            //     vt.push_back(thread(std::bind(&Healpix_Correlation::parallelSetPixValue,this,begin,end,std::ref(pointsCount),value)));
            // }
            
        #ifndef JOIN_THREAD_MOVE
            for (int i = 0; i < threadCount; i+= 1) {
                vt[i].join();
            }
        #endif
            if (value) {
        #ifdef RESET_VALUE_DEBUG
                    cout << "In the caller: " << endl;
                    cout << "pointsCount: " << pointsCount << endl;
        #endif        
                    // pointsCount = nonZeroCount;
                    
            }
            for (int pix = bottomPixel + 1; pix <= endPixel; pix += 1) {
                (*this)[pix] = 0;
            }
        #else
            for (int pix = 0; pix < topPixel; pix += 1) {
                (*this)[pix] = 0;
            }
            for(int pix = topPixel; pix <= bottomPixel; pix+= 1) {
                (*this)[pix] = value;
                if (fabs(value) > MY_EPSILON) {
                    pointsCount++;
                }
            }
            for (int pix = bottomPixel + 1; pix <= endPixel; pix += 1) {
                (*this)[pix] = 0;
            }
        #endif        
    }
    else {  
        cout << "has phi constraint,, does not support parallel set value\n";
        for (int pix = 0; pix < topPixel; pix += 1) {
            (*this)[pix] = 0;
        }
        for(int pix = topPixel; pix <= bottomPixel; pix+= 1) {
            (*this)[pix] = value;
            if (fabs(value) > MY_EPSILON) {
                pointsCount++;
            }
        }
        for (int pix = bottomPixel + 1; pix <= endPixel; pix += 1) {
            (*this)[pix] = 0;
        }

    }

#ifdef RESET_VALUE_DEBUG
    int checkCount = 0;
    for (int pix = 0; pix <=bottomPixel; pix++) {
        if ((*this)[pix]) {
            checkCount++;
        }
    }
    cout << "\n\ncheckCount: " << checkCount << endl;
#endif
}

void Healpix_Correlation::scalePixValue(myDouble scale, int _threadNum)
{
#ifdef PARALLEL
    int threadCount = _threadNum;
    countType stepSize = ceil(bottomPixel/threadCount);
    std::vector<thread> vt;
    std::vector<int> index;
    index.resize(threadCount);
    for (int i = 0; i < threadCount; i+= 1) {
        typeIndex begin, end;
        begin = i * stepSize;
        if(i != threadCount - 1) {
            end = (i + 1) * stepSize;
            index[i] = end;
        }
        else {
            end = bottomPixel + 1;
            index[i] = end;
        }
#ifdef RESET_VALUE_DEBUG
        // cout << "in the caller: \n";
        // cout << "begin: " << begin << endl;
        // cout << "end: " << end << endl;
#endif
        vt.push_back(thread(std::bind(&Healpix_Correlation::parallelScalePixValue,this,begin,end,std::ref(scale))));

#ifdef JOIN_THREAD_MOVE
        vt[i].join();
#endif
    }
#ifndef JOIN_THREAD_MOVE
    for (int i = 0; i < threadCount; i+= 1) {
        vt[i].join();
    }
#endif    

#else
    for (typeIndex pix = 0; pix <= bottomPixel; pix++) {
        (*this)[pix] *= scale;
    }
#endif
    pointsCount *= scale;
}

void Healpix_Correlation::SetAngle(const myDouble & angle) {
    spannedAngle = angle;
    updateBoundPixel();
}

void Healpix_Correlation::SetDeltaTheta(const myDouble & deltaTheta, const myDouble & angle) {
    int order = deltaTheta2Order(deltaTheta);
    Healpix_Map::Set(order,RING);
    spannedAngle = angle;
    update();
}

void Healpix_Correlation::SetPixValue(const spContainer & vc,int _threadNum) {
    if(!boundSet) {
        cout << "bound not set!\n";
        exit(HEALPIX_BOUND_NOT_SET);
    }
    resetPixValue(0,_threadNum);
    if(!isAllZero()) {
        cout << "initial setting up is wrong, some non zero values\n";
        exit(INITIAL_MAP_NONZERO);
    }
    pointsCount = 0L;
    for (spContainer::const_iterator iter = vc.begin(); iter != vc.end(); iter+= 1) {
        int pixel = vec2pix(vec3((double)iter -> getXCoord(), 
            (double)iter -> getYCoord(), (double)iter -> getZCoord()));
        (*this)[pixel]+= 1;
        pointsCount += 1;
    }
}

void Healpix_Correlation::SetPixValue(countType pointCount,
    int nTime,
    int _threadNum) 
{
    if(!boundSet) {
        cout << "bound not set!\n";
        exit(HEALPIX_BOUND_NOT_SET);
    }
    resetPixValue(0,_threadNum);
    if(!isAllZero()) {
        cout << "the reset is imcomplete, some unzero.\n";
        exit(INITIAL_MAP_NONZERO);
    }
#ifdef TIME_FIXED
    srand(TIME_FIXED);
#else
    srand((unsigned) time(0));
#endif
    pointsCount = 0;
    if (hPhiEnd - hPhiBeg > twopi - MY_EPSILON) {
        int pixRange = bottomPixel - topPixel + 1;
        for (int nt = 0; nt < nTime; nt++) {
            for (countType count = 0; count < pointCount; count+= 1) {
                int pix = rand() % (pixRange) + topPixel;
                (*this)[pix]+= 1;
                pointsCount+= 1;
            }
        }
    }
    else {
        cout << "phi constraint, computation heavier method than no constraint for set pix value\n";
        int pixRange = bottomPixel - topPixel + 1;
        int count = 0;
        int loopTime = 0;
        for (int nt = 0; nt < nTime; nt++) { 
            while(count < pointCount && loopTime < HEALPIX_MAX_LOOP_TIME) {
                // theta constraint
                int pix = rand() % (pixRange) + topPixel;
                // check the phi constraint
                pointing ptg = pix2ang(pix);
                myDouble pixPhi = ptg.phi;
                if (pixPhi > hPhiBeg - MY_EPSILON && pixPhi < hPhiEnd + MY_EPSILON) {
                    (*this)[pix]+= 1;
                    pointsCount+= 1;
                    count++;
                }
                loopTime++;
            }
        }        
    }
    scalePixValue(1.0/nTime,_threadNum);
}

void Healpix_Correlation::SetPixValue(const std::vector<myDouble> & _thetaVec,
    const std::vector<myDouble> & _phiVec,
    const std::vector<myDouble> & _dataVec,
    myDouble & _threshold,
    int mode, int _threadNum)
{
    spContainer dummy;
    SetPixValue(_thetaVec,_phiVec,_dataVec,dummy, _threshold,
        mode, _threadNum);
}

void Healpix_Correlation::SetPixValue(const std::vector<myDouble> & _thetaVec,
    const std::vector<myDouble> & _phiVec,
    const std::vector<myDouble> & _dataVec,
    spContainer & vecPoint,
    myDouble & _threshold, 
    int mode, int _threadNum)
{
    if(!boundSet) {
        cout << "bound not set!\n";
        exit(HEALPIX_BOUND_NOT_SET);
    }
    resetPixValue(0,_threadNum);
    if(!isAllZero()) {
        cout << "the reset is imcomplete, some unzero.\n";
        exit(INITIAL_MAP_NONZERO);
    }
#ifdef SIGNAL_SET_DEBUG
    cout << "\n ****+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1************\n";
    cout << "the SetPixValue function!!!\n";
#endif
    vecPoint.clear();
    // myDouble theta,phi,val;
    sizeType thetaLen = _thetaVec.size();
    sizeType phiLen = _phiVec.size();
    for (typeIndex dataIndex = 0; dataIndex < _dataVec.size(); dataIndex += 1) {
        // #ifdef SIMPLE_THRESHOLD
        myDouble val = signalStrength2Val(_dataVec[dataIndex], _threshold, mode) * FACTOR_SIGVAL;
    #ifdef QUADRATIC
            val *= val;
    #endif
        if (val > MY_EPSILON) {
            typeIndex thetaIndex,phiIndex;
            convertDataVecIndex(dataIndex,thetaIndex,phiIndex,thetaLen,phiLen);
            myDouble thetaTemp = convertDeg2Radian(_thetaVec[thetaIndex]);
    #ifdef DATA_VALUE_SET_DEBUG
            cout << "thetaTemp: " << thetaTemp << "\t degree: " << _thetaVec[thetaIndex] << endl;
    #endif
            myDouble phiTemp = convertDeg2Radian((double)_phiVec[phiIndex]);
            int pixel = zphi2pix(cos(thetaTemp),phiTemp);
            if (mode == SIMPLE_THRESHOLD) {
                if((*this)[pixel] < MY_EPSILON) {
                    (*this)[pixel] = val ;
                    pointsCount += val;

                }                
                else {
                    continue;
                }
            }
            else if (mode == SIMPLE_INCREMENT) {
                (*this)[pixel] += val;
                pointsCount += val;
            }
            else if (mode == SIMPLE_POWER_mW_ADD || mode == SIMPLE_POWER_W_ADD) {
                (*this)[pixel] += val;
                // treat each point as 1
                pointsCount += 1;
            }
           vecPoint.push_back(SphericalPoint(1,thetaTemp, phiTemp));
    #ifdef DATA_VALUE_SET_DEBUG
            SphericalPoint tempP(1,thetaTemp, phiTemp);
            cout << tempP << endl;
    #endif
        }
// #ifdef DATA_VALUE_SET_DEBUG

//     cout << "pointsCount: " << pointsCount << " \t dataIndx: " << dataIndex << endl;
// #endif              
// #endif       
#ifdef SIGNAL_SET_DEBUG
cout << "\n ****************\n";
cout << "pix: " << pixel << endl;
cout << "val : " << val << endl;
cout << "pointsCount: " << pointsCount << endl;        
#endif
    }
}


void Healpix_Correlation::ModifySize(myDouble stepTheta) {
    // to do
    cout << "need to reduce the deltaTheta of healpix" << endl;
    exit(TODO_FUNCTION);
    SetDeltaTheta(stepTheta);
}

countType Healpix_Correlation::Neighbors(const myDouble & theta) {
    countType pair = 0;
#ifdef HEALPIX_DEBUG
    cout << "\n**************** " << endl;
    cout << "function Neighbors with theta: " << theta << endl;
    cout << "the bottomPixel is " << bottomPixel << endl;
#endif
    for(int pix = 0; pix <= bottomPixel; pix+= 1) {
        if((*this)[pix] == 0) {
            continue;
        }
        std::vector<int> pixset;
#ifdef HEALPIX_DEBUG
    cout << "in the first loop of Neighbors, pix: " << pix << endl;
#endif

#ifdef RANGESET_QUERY
    #ifdef APPEND_RANGESET_MODIFY 
        query_disc_pixel_Internal(pix, theta, pixset,*this);
    #else
        query_disc_pixel_Internal(pix, theta, pixset);
    #endif
#else
    #ifdef QUERY_COUNT_DIRECTLY
        pair += query_disc_count(pix,theta,*this);
        continue;
    #else 
        query_disc_pixel_InternalVec(pix,theta,pixset,*this);
    #endif
#endif

#ifdef HEALPIX_DEBUG
    cout << "after the first loop of Neighbors" << endl;
#endif
    #ifdef OUTPUT_PIX_NEIGHBOR_COUNT
        if (fabs(theta - TO_PRINT_THETA) < 0.01)
            cout << "*******\npixel: " << pix << "\ttheta: " << theta << endl;
    #endif
        for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
            pair += (*this)[*iter];
    #ifdef OUTPUT_PIX_NEIGHBOR_COUNT
        if (fabs(theta - TO_PRINT_THETA) < 0.01)
            cout << *iter << "\t";
    #endif
        }
    #ifdef OUTPUT_PIX_NEIGHBOR_COUNT
        if (fabs(theta - TO_PRINT_THETA) < 0.01)
            cout << endl << "neighbor: " << pair << endl;
    #endif        
    }
#ifdef HEALPIX_DEBUG
    cout << "function Neighbors end with pair/2" << pair/2 << endl;
    cout << "\n**************** " << endl;
#endif

#ifndef NEIGHBOR_SCALE_2
    return pair/2;
#else 
    #ifdef PAIR_IS_ODD_TEST
        if (pair % 2) {
            cout << "at theta: " << theta << "\t, cumulated pair:" << pair << "\t odd!\n";
        }
    #endif
    return pair;
#endif
}

countType Healpix_Correlation::Neighbors(int pix, 
    const myDouble & theta)
{
    countType pair = 0;
#ifdef HEALPIX_DEBUG
    cout << "\n**************** " << endl;
    cout << "function Neighbors Pixel with theta: " << theta << endl;
    cout << "the bottomPixel is " << bottomPixel << endl;
#endif

        std::vector<int> pixset;
// #ifdef HEALPIX_DEBUG
//     cout << "in the first loop of Neighbors, pix: " << pix << endl;
// #endif

#ifdef RANGESET_QUERY
    #ifdef APPEND_RANGESET_MODIFY 
        query_disc_pixel_Internal(pix, theta, pixset,*this);
    #else
        query_disc_pixel_Internal(pix, theta, pixset);
    #endif
#else
    #ifdef QUERY_COUNT_DIRECTLY
        pair = query_disc_count(pix,theta,*this);
        return pair;
    #else 
        query_disc_pixel_InternalVec(pix,theta,pixset,*this);
    #endif
#endif

// #ifdef HEALPIX_DEBUG
//     cout << "after the first loop of Neighbors" << endl;
// #endif
        
#ifdef QUERY_NON_CUMU_DEBUG
        cout << "In the Cumu version neighbor\n";
        cout << "pixel: " << pix << endl;
        cout << "theta: " << theta << endl;
        cout << "pixset: " << endl;
#endif
        for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
            pair += (*this)[*iter];
#ifdef QUERY_NON_CUMU_DEBUG
        cout << *iter << "\t" << endl;
#endif
        }

#ifdef QUERY_NON_CUMU_DEBUG
        cout << endl;
#endif
#ifdef HEALPIX_DEBUG
    cout << "function Neighbors end with pair/2" << pair/2 << endl;
    cout << "\n**************** " << endl;
#endif
    return pair;
}

countType Healpix_Correlation::Neighbors(int pix, const myDouble & theta, std::string filename)
{
    countType pair = 0;
#ifdef HEALPIX_DEBUG
    cout << "\n**************** " << endl;
    cout << "function Neighbors Pixel with theta: " << theta << endl;
    cout << "the bottomPixel is " << bottomPixel << endl;
#endif

        std::vector<int> pixset;
// #ifdef HEALPIX_DEBUG
//     cout << "in the first loop of Neighbors, pix: " << pix << endl;
// #endif

#ifdef RANGESET_QUERY
    #ifdef APPEND_RANGESET_MODIFY 
        query_disc_pixel_Internal(pix, theta, pixset,*this);
    #else
        query_disc_pixel_Internal(pix, theta, pixset);
    #endif
#else
    #ifdef QUERY_COUNT_DIRECTLY
        if (filename != MY_DEFAULT_STRING) {
            cout << "outputting result to file is not supported for Neighbors\n";
            cout << "reset the QUERY_COUNT_DIRECTLY to enable the function\n";
        }
        pair = query_disc_count(pix,theta,*this);
        return pair;
    #else 
        query_disc_pixel_InternalVec(pix,theta,pixset,*this);
    #endif
#endif

// #ifdef HEALPIX_DEBUG
//     cout << "after the first loop of Neighbors" << endl;
// #endif
        
#ifdef QUERY_NON_CUMU_DEBUG
        cout << "In the Cumu version neighbor\n";
        cout << "pixel: " << pix << endl;
        cout << "theta: " << theta << endl;
        cout << "pixset: " << endl;
#endif
        for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
            pair += (*this)[*iter];
#ifdef QUERY_NON_CUMU_DEBUG
        cout << *iter << "\t" << endl;
#endif
        }

#ifdef QUERY_NON_CUMU_DEBUG
        cout << endl;
#endif
#ifdef HEALPIX_DEBUG
    cout << "function Neighbors end with pair/2" << pair/2 << endl;
    cout << "\n**************** " << endl;
#endif
    
#ifndef QUERY_COUNT_DIRECTLY
    std::ofstream outFile;
    outFile.open(filename);
    for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
        outFile << *iter << endl;
    }
    outFile.close();
#endif
    return pair;
}

countType Healpix_Correlation::NeighborsParallelVector(const myDouble & theta,int _threadNum) {
    int threadCount = _threadNum;
    countType stepSize = ceil(bottomPixel/threadCount);
    std::vector<corrType> pairCount;
    std::vector<thread> vt;
    for (int i = 0; i < threadCount; i+= 1) {
        typeIndex begin, end;
        begin = i * stepSize;
        if(i != threadCount - 1) {
            end = (i + 1) * stepSize;
        }
        else {
            end = bottomPixel + 1;
        }
        vt.push_back(thread(std::bind(&Healpix_Correlation::NeighbourIntervalCountVector, this, begin, 
            end,std::cref(theta), std::ref(pairCount))));
#ifdef JOIN_THREAD_MOVE
        vt[i].join();
#endif
    }
#ifndef JOIN_THREAD_MOVE
    for (int i = 0; i < threadCount; i+= 1) {
        vt[i].join();
    }
#endif
#ifdef HEALPIX_PARALLEL_DEBUG
    cout << "the length of pairCount: " << pairCount.size() << endl;
    cout << "the content of pairCount\n";
    for(auto l : pairCount) {
        cout << l << ' ' << endl;
    }
#endif 
    countType pair = 0;
    for (auto l : pairCount) {
        pair += l;
    }
    return pair/2;
}

countType Healpix_Correlation::NeighborsParallel(const myDouble & theta, 
    int _threadNum) 
{
    int threadCount = _threadNum;
    countType stepSize = ceil(bottomPixel/threadCount);
    countType pairCount = 0;
    std::vector<thread> vt;
    for (int i = 0; i < threadCount; i+= 1) {
        typeIndex begin, end;
        begin = i * stepSize;
        if(i != threadCount - 1) {
            end = (i + 1) * stepSize;
        }
        else {
            end = bottomPixel + 1;
        }
        vt.push_back(thread(std::bind(&Healpix_Correlation::NeighbourIntervalCount, this, begin, end,
            std::cref(theta), std::ref(pairCount))));
#ifdef JOIN_THREAD_MOVE
        vt[i].join();
#endif
    }
#ifndef JOIN_THREAD_MOVE
    for (int i = 0; i < threadCount; i+= 1) {
        vt[i].join();
    }
#endif
#ifdef HEALPIX_PARALLEL_DEBUG
    cout << "the value of pairCount: " << pairCount << endl;
#endif 

#ifdef NEIGHBOR_SCALE_2
    return pairCount;
#else
    return pairCount/2;
#endif
}

countType Healpix_Correlation::crossNeighbors(const myDouble & theta, Healpix_Correlation & hp) {
    if (!isStructureIsomorphic(hp)) {
        cout << "two healpix map not of similar structure, crossNeighbors impossible" << endl;
        exit(NOT_STRUCTURE_ISOM_HEALPIX);
    }
    countType pair = 0;
    for (int pix = 0; pix <= bottomPixel; pix+= 1) {
        if((*this)[pix] == 0) {
            continue;
        }
        std::vector<int> pixset;
#ifdef RANGESET_QUERY
    #ifdef APPEND_RANGESET_MODIFY 
        query_disc_pixel_Internal(pix, theta, pixset,hp);
    #else
        query_disc_pixel_Internal(pix, theta, pixset);
    #endif
#else
    #ifdef QUERY_COUNT_DIRECTLY
        pair += query_disc_count(pix,theta,hp);
        continue;
    #else 
        query_disc_pixel_InternalVec(pix,theta,pixset,hp);
    #endif
#endif
        for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
            pair += hp[*iter];
        }
    }
    return pair;
}

countType Healpix_Correlation::crossNeighborsParallel(const myDouble & theta, Healpix_Correlation & hp,int _threadNum) {
    int threadCount = _threadNum;
    countType stepSize = ceil(bottomPixel/threadCount);
    countType pairCount = 0;
    std::vector<thread> vt;
    for (int i = 0; i < threadCount; i+= 1) {
        typeIndex begin, end;
        begin = i * stepSize;
        if(i != threadCount - 1) {
            end = (i + 1) * stepSize;
        }
        else {
            end = bottomPixel + 1;
        }
        vt.push_back(thread(std::bind(&Healpix_Correlation::crossNeighbourIntervalCount, this, begin, 
            end,std::cref(theta), std::ref(hp), std::ref(pairCount))));
#ifdef JOIN_THREAD_MOVE
        vt[i].join();
#endif
    }
#ifndef JOIN_THREAD_MOVE
    for (int i = 0; i < threadCount; i+= 1) {
        vt[i].join();
    }
#endif
#ifdef HEALPIX_PARALLEL_DEBUG
    cout << "the value of pairCount: " << pairCount << endl;
#endif 
    return pairCount;
}

std::vector<corrType> Healpix_Correlation::NeighborVec(std::vector<myDouble>& bins,
    int _threadNum) 
{
    // assume uniform bins
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
    std::vector<corrType> vec;
    NeighborVec(bins,vec,_threadNum);
    return vec;
}

void Healpix_Correlation::NeighborVec(std::vector<myDouble>& bins,
    std::vector<corrType> & neighborCount,
    int _threadNum)
{
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
#ifdef HEALPIX_DEBUG
    cout << "\n********************\n";
    cout << "function NeighborVec\n";
    cout << "after the ModifySize in NeighborVec\n";
#endif
    size_t binSize = bins.size();
    if (neighborCount.size() != binSize - 1) 
        neighborCount.resize(binSize - 1);
#ifndef NON_CUM_COUNT
        std::vector<corrType> cumulatedCount;
        cumulatedCount.resize(binSize);
        for(typeIndex index = 0; index < binSize; index+= 1) {
    #ifdef HEALPIX_DEBUG
            cout << "in the first for loop, index = " << index << endl;
    #endif
            // using parallel
    #ifdef PARALLEL
            countType nCount = NeighborsParallel(bins[index],_threadNum);
    #else
            countType nCount = Neighbors(bins[index]);
    #endif
    #ifdef HEALPIX_DEBUG
            cout << "after Neighbors, nCount = " << nCount << endl;
    #endif
            cumulatedCount[index] = nCount;
            if(index != 0) {
                neighborCount[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
    #ifdef NEIGHBOR_SCALE_2
        #ifdef PAIR_IS_ODD_TEST
                if (neighborCount[index - 1] % 2) {
                    cout << "at theta: " << bins[index - 1] << "\tpair:" << neighborCount[index - 1] << "  odd!\n";
                }
        #endif
                neighborCount[index - 1] = neighborCount[index - 1] / 2;
    #endif
            }
        }
    #ifdef HEALPIX_DEBUG
        cout << "Function NeighborVec ends" << endl;
        cout << "\n********************\n";
    #endif
#endif
#ifdef NON_CUM_COUNT
    NeighborNonCumVec(bins,neighborCount,*this,0,_threadNum);
    for (typeIndex i = 0; i < binSize - 1; i++) {
        neighborCount[i] = neighborCount[i] / 2;
    }
#endif
}

void Healpix_Correlation::NeighborPixelCumVec(
    std::vector<myDouble> & bins, 
    int pixel, std::vector<corrType> & neighborCount, 
    std::vector<corrType> & cumulatedCount,
    int threadNum)
{
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
#ifdef HEALPIX_DEBUG
    cout << "\n********************\n";
    cout << "function NeighborVec\n";
    cout << "after the ModifySize in NeighborVec\n";
#endif
    int binSize = bins.size();
    neighborCount.resize(binSize - 1);
    // std::vector<corrType> cumulatedCount;
    cumulatedCount.resize(binSize);
    for(int index = 0; index < binSize; index+= 1) {
#ifdef HEALPIX_DEBUG
        cout << "in the first for loop, index = " << index << endl;
#endif

        countType nCount = Neighbors(pixel,bins[index]);
#ifdef HEALPIX_DEBUG
        cout << "after Neighbors, nCount = " << nCount << endl;
#endif
        cumulatedCount[index] = nCount;
        if(index != 0) {
            neighborCount[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
        }
    }
#ifdef HEALPIX_DEBUG
    cout << "Function NeighborVec ends" << endl;
    cout << "\n********************\n";
#endif
    
}

void Healpix_Correlation::NeighborPixelCumVec(std::vector<myDouble> & bins, int pixel, 
    std::vector<corrType> & neighborCount, 
    std::vector<corrType> & cumulatedCount,
    vector<std::string> fileNames,
    vector<myDouble> thetas,
    int threadNum)
{
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
#ifdef HEALPIX_DEBUG
    cout << "\n********************\n";
    cout << "function NeighborVec\n";
    cout << "after the ModifySize in NeighborVec\n";
#endif
    int binSize = bins.size();
    neighborCount.resize(binSize - 1);
    // std::vector<corrType> cumulatedCount;
    cumulatedCount.resize(binSize);

    int len = fileNames.size();
    int len2 = thetas.size();
    if (len != len2) {
        cout << "File numbers not matching in Neighbors!\n";
        exit(FILE_NUM_NOT_MATCHING);
    }
    bool noMorePlot = false;
    int _thetasIndex = 0;
    myDouble thetaPlot;
    std::string filename;
    if (len == 0) {
        noMorePlot = true;
    }
    else {
        thetaPlot = thetas[_thetasIndex];
    }

    for(int index = 0; index < binSize; index+= 1) {
#ifdef HEALPIX_DEBUG
        cout << "in the first for loop, index = " << index << endl;
#endif
        countType nCount;
        if (!noMorePlot && fabs(thetaPlot - bins[index]) < 0.01) {
            filename = fileNames[_thetasIndex];
            _thetasIndex++;
            nCount = Neighbors(pixel,bins[index], filename);
            if (_thetasIndex >= len) {
                noMorePlot = true;
            }
            else {
                thetaPlot = thetas[_thetasIndex];
            }
        }
        else
            nCount = Neighbors(pixel,bins[index]);
#ifdef HEALPIX_DEBUG
        cout << "after Neighbors, nCount = " << nCount << endl;
#endif
        cumulatedCount[index] = nCount;
        if(index != 0) {
            neighborCount[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
        }
    }
#ifdef HEALPIX_DEBUG
    cout << "Function NeighborVec ends" << endl;
    cout << "\n********************\n";
#endif

}    

void Healpix_Correlation::NeighborPixelNonCumVec(std::vector<myDouble> & bins, int pixel, 
    std::vector<corrType> & neighborCount, int threadNum)
{
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
#ifdef HEALPIX_DEBUG
    cout << "\n********************\n";
    cout << "function NeighborVec\n";
    cout << "after the ModifySize in NeighborVec\n";
#endif
    int binSize = bins.size();
    neighborCount.resize(binSize - 1);

    query_disc_pixel_Bins(pixel,bins,neighborCount,*this);
     
}    

void Healpix_Correlation::NeighborPixelNonCumVec(std::vector<myDouble> & bins, int pixel, 
    std::vector<corrType> & neighborCount, 
    vector<std::string> fileNames,
    vector<myDouble> thetas,            
    int threadNum)
{
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
#ifdef HEALPIX_DEBUG
    cout << "\n********************\n";
    cout << "function NeighborVec\n";
    cout << "after the ModifySize in NeighborVec\n";
#endif
    int binSize = bins.size();
    neighborCount.resize(binSize - 1);

    query_disc_pixel_Bins(pixel,bins,neighborCount,*this,*this,fileNames,thetas);
     
}


void Healpix_Correlation::NeighborNonCumVec(std::vector<myDouble> & bins, 
    std::vector<corrType> & neighborCount,
    const Healpix_Correlation & hp,
    int Factor,
    int _threadNum)
{
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
#ifdef HEALPIX_DEBUG
    cout << "\n********************\n";
    cout << "function NeighborVec\n";
    cout << "after the ModifySize in NeighborVec\n";
#endif
    size_t binSize = bins.size();
    if (neighborCount.size() != binSize - 1)
        neighborCount.resize(binSize - 1);
    std::vector<corrType> TempCount;
    TempCount.resize(binSize - 1);        
    // std::vector<corrType> cumulatedCount;
    // cumulatedCount.resize(binSize);
#ifdef PARALLEL
    int threadCount = _threadNum;
    countType stepSize = ceil(bottomPixel/threadCount);
    std::vector<thread> vt;
    std::vector<int> index;
    index.resize(threadCount);

    for (int i = 0; i < threadCount; i+= 1) {
        typeIndex begin, end;
        begin = i * stepSize;
        if(i != threadCount - 1) {
            end = (i + 1) * stepSize;
            index[i] = end;
        }
        else {
            end = bottomPixel + 1;
            index[i] = end;
        }
        vt.push_back(thread(std::bind(&Healpix_Correlation::NeighbourNonCumIntervalCount, this, begin, end,std::cref(bins),std::ref(TempCount),std::cref(hp),Factor)));
#ifdef JOIN_THREAD_MOVE
        vt[i].join();
#endif
    }

#ifndef JOIN_THREAD_MOVE
    for (int i = 0; i < threadCount; i+= 1) {
        vt[i].join();
    }
#endif

#else

    NeighbourNonCumIntervalCount(0,bottomPixel + 1, bins, TempCount,hp,Factor);
#endif
    // for (int i = 0; i < binSize - 1; i++) {
    //     neighborCount[i] / 2;
    // }
    for (typeIndex index = 0; index < binSize - 1; index++) {
        neighborCount[index] = (neighborCount[index]  * (Factor) + TempCount[index])/(Factor + 1);
    }
}


std::vector<corrType> Healpix_Correlation::crossNeighborVec(std::vector<myDouble>& bins, Healpix_Correlation & hp,int _threadNum) 
{
    
    std::vector<corrType> crossCount;
    crossNeighborVec(bins,crossCount,hp,_threadNum);
    return crossCount;
}

void Healpix_Correlation::crossNeighborVec(std::vector<myDouble>& bins, std::vector<corrType> & neighborCount, Healpix_Correlation & hp,int _threadNum) {
    if (!isStructureIsomorphic(hp)) {
        cout << "two non-structure similar maps in crossNeighborVec" << endl;
        exit(NOT_STRUCTURE_ISOM_HEALPIX);
    }
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
        hp.ModifySize(stepTheta);
    }

    int binSize = bins.size();
    neighborCount.resize(binSize - 1);
    std::vector<corrType> cumulatedCount;
    cumulatedCount.resize(binSize);
    for(int index = 0; index < binSize; index+= 1) {
        // using parallel version
#ifdef PARALLEL
        countType nCount = crossNeighborsParallel(bins[index],hp,_threadNum);
#else
        countType nCount = crossNeighbors(bins[index],hp);
#endif
        cumulatedCount[index] = nCount;
        if(index != 0) {
            neighborCount[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
        }
    }
}

std::vector<corrType> Healpix_Correlation::averageCrossNeighborVec(std::vector<myDouble> &bins, 
    countType nk, int nTime,int _threadNum) {
    std::vector<corrType> neighCount;
    averageCrossNeighborVec(bins,neighCount,nk,nTime,_threadNum);
    return neighCount;
}   

void Healpix_Correlation::averageCrossNeighborVec(std::vector<myDouble>& bins,
     std::vector<corrType> & neighborCount, 
     countType nk, int nTime,int _threadNum) 
{
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
    int binSize = bins.size();
    neighborCount.resize(binSize - 1);
    std::vector<corrType> cumulatedCount;
    cumulatedCount.resize(binSize);
#ifdef SIMILAR_GEOMETRY
    Healpix_Correlation hp(Order(), RING, hThetaEnd,hThetaBeg, hPhiEnd, hPhiBeg);
#else
    Healpix_Correlation hp(Order(), RING, spannedAngle);
#endif

#ifndef AVERAGE_RANDOM_MAP
    for(int term = 0; term < nTime; term+= 1) {
        hp.SetPixValue(nk);
#ifndef NON_CUM_COUNT
        for(int index = 0; index < binSize; index+= 1) {
        // using parallel version
#ifdef PARALLEL
            countType nCount = crossNeighborsParallel(bins[index],hp,_threadNum);
#else
            countType nCount = crossNeighbors(bins[index],hp);
#endif
            cumulatedCount[index] = (cumulatedCount[index]  * (term) + nCount)/(term + 1);
            // the last term
            if (term == nTime - 1) {
                if(index != 0) {
                    neighborCount[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
                }
            }
        }
#else
    NeighborNonCumVec(bins,cumulatedCount,hp,term,_threadNum);
    if (term == nTime - 1) {
        for (int index = 0; index < binSize - 1; index++) {
            // no need to divide by two
            neighborCount[index] = cumulatedCount[index];
        }
    }
#endif
    }

#else

#ifdef AVERAGE_RR_DR
   for(int term = 0; term < nTime; term+= 1) {
        hp.SetPixValue(nk,nTime);
#ifndef NON_CUM_COUNT
        for(int index = 0; index < binSize; index+= 1) {
        // using parallel version
#ifdef PARALLEL
            countType nCount = crossNeighborsParallel(bins[index],hp,_threadNum);
#else
            countType nCount = crossNeighbors(bins[index],hp);
#endif
            cumulatedCount[index] = (cumulatedCount[index]  * (term) + nCount)/(term + 1);
            // the last term
            if (term == nTime - 1) {
                if(index != 0) {
                    neighborCount[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
                }
            }
        }
#else
        NeighborNonCumVec(bins,cumulatedCount,hp,term,_threadNum);
        if (term == nTime - 1) {
            for (int index = 0; index < binSize - 1; index++) {
                // no need to divide by two
                neighborCount[index] = cumulatedCount[index];
            }
        }
   }
#endif
#else
    hp.SetPixValue(nk,nTime);
    #ifndef NON_CUM_COUNT
        #ifdef PARALLEL
            for (int index = 0; index < binSize; index += 1) {
                cumulatedCount[index] = crossNeighborsParallel(bins[index],hp,_threadNum);
        #else
                countType nCount = crossNeighbors(bins[index],hp);
        #endif
                // the last term
                if(index != 0) {
                    neighborCount[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
                }
            }
    #else
        NeighborNonCumVec(bins,cumulatedCount,hp,0,_threadNum);
        for (int index = 0; index < binSize - 1; index++) {
            // no need to divide by two
            neighborCount[index] = cumulatedCount[index];
        }
    #endif
#endif
#endif
}

std::vector<corrType> Healpix_Correlation::averageRR(std::vector<myDouble>& bins, 
    countType nk, int nTime,int _threadNum) {
    std::vector<corrType> RR;
    averageRR(bins, RR, nk ,nTime,_threadNum);
    return RR;
}

void Healpix_Correlation::averageRR(std::vector<myDouble>& bins, 
    std::vector<corrType>& RR, 
    countType nk, int nTime,int _threadNum) 
{
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
    int binSize = bins.size();
    RR.resize(binSize - 1);

    std::vector<corrType> cumulatedCount;
    cumulatedCount.resize(binSize,0);
#ifdef SIMILAR_GEOMETRY
    Healpix_Correlation hp(Order(), RING, hThetaEnd,hThetaBeg, hPhiEnd, hPhiBeg);
#else
    Healpix_Correlation hp(Order(), RING, spannedAngle);
#endif

#ifdef AVERAGE_DEBUG
    cout << "*****************\n";
    cout << "in averageRR, hp created successfully\n";
#endif
#ifndef AVERAGE_RANDOM_MAP
        for(int term = 0; term < nTime; term+= 1) {
            hp.SetPixValue(nk);
            // cout << "wrong setpixValue in hp.SetpixValue in averageRR\n";
            // exit(TODO_FUNCTION);

    #ifdef AVERAGE_DEBUG
            cout << "+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+\n";
            cout << "in the loop, term: " << term << "  SetPixValue successfully";
            cout << "~~~~~~~~~~~~~~~~~~~\n";
    #endif
    #ifndef NON_CUM_COUNT
                for(int index = 0; index < binSize; index+= 1) {
                    // using parallel version
        #ifdef PARALLEL
                    countType nCount = hp.NeighborsParallel(bins[index],_threadNum);
                    #ifdef HEALPIX_ZERO_DEBUG
                        if(index == 0) {
                            cout << "the index = 0 count: " << nCount << endl;
                        }
                        else if(index == 1) {
                            cout << "the index = 1 count: " << nCount << endl;
                        }
                    #endif
        #else
                    countType nCount = hp.Neighbors(bins[index]);
        #endif
                    cumulatedCount[index] = (cumulatedCount[index]  * (term) + nCount)/(term + 1);
                    // the last term
                    if (term == nTime - 1) {
                        if(index != 0) {
                            RR[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
        #ifdef NEIGHBOR_SCALE_2
                            RR[index - 1] = RR[index - 1] / 2;
        #endif
                        }
                    }
                }
    #else
            hp.NeighborNonCumVec(bins,cumulatedCount,hp,term,_threadNum);
            if (term == nTime - 1) {
                for (int index = 0; index < binSize - 1; index++) {
                    RR[index] = cumulatedCount[index] / 2;
                }
            }
    #endif

        }

#else

#ifdef AVERAGE_RR_DR
        for(int term = 0; term < nTime; term+= 1) {
            hp.SetPixValue(nk,nTime);
            // cout << "wrong setpixValue in hp.SetpixValue in averageRR\n";
            // exit(TODO_FUNCTION);

    #ifdef AVERAGE_DEBUG
            cout << "+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+= 1+\n";
            cout << "in the loop, term: " << term << "  SetPixValue successfully";
            cout << "~~~~~~~~~~~~~~~~~~~\n";
    #endif
    #ifndef NON_CUM_COUNT
                for(int index = 0; index < binSize; index+= 1) {
                    // using parallel version
        #ifdef PARALLEL
                    countType nCount = hp.NeighborsParallel(bins[index],_threadNum);
                    #ifdef HEALPIX_ZERO_DEBUG
                        if(index == 0) {
                            cout << "the index = 0 count: " << nCount << endl;
                        }
                        else if(index == 1) {
                            cout << "the index = 1 count: " << nCount << endl;
                        }
                    #endif
        #else
                    countType nCount = hp.Neighbors(bins[index]);
        #endif
                    cumulatedCount[index] = (cumulatedCount[index]  * (term) + nCount)/(term + 1);
                    // the last term
                    if (term == nTime - 1) {
                        if(index != 0) {
                            RR[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
        #ifdef NEIGHBOR_SCALE_2
                            RR[index - 1] = RR[index - 1] / 2;
        #endif
                        }
                    }
                }
    #else
            hp.NeighborNonCumVec(bins,cumulatedCount,hp,term,_threadNum);
            if (term == nTime - 1) {
                for (int index = 0; index < binSize - 1; index++) {
                    RR[index] = cumulatedCount[index] / 2;
                }
            }
    #endif
        }
#else
    hp.SetPixValue(nk,nTime);
    #ifdef AVERAGE_RMAP_DEBUG
        cout << hp << endl;
        cout << hp.GetPointsCount() << endl;
    #endif
    #ifndef NON_CUM_COUNT
        for(int index = 0; index < binSize; index+= 1) {
            // using parallel version
    #ifdef PARALLEL
            countType nCount = hp.NeighborsParallel(bins[index],_threadNum);
            #ifdef HEALPIX_ZERO_DEBUG
                if(index == 0) {
                    cout << "the index = 0 count: " << nCount << endl;
                }
                else if(index == 1) {
                    cout << "the index = 1 count: " << nCount << endl;
                }
            #endif
    #else
            countType nCount = hp.Neighbors(bins[index]);
    #endif
            cumulatedCount[index] = nCount;
            // the last term
            if(index != 0) {
                RR[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
    #ifdef NEIGHBOR_SCALE_2
                RR[index - 1] = RR[index - 1] / 2;
    #endif
            }
        }
    #else
        hp.NeighborNonCumVec(bins,cumulatedCount,hp,0,_threadNum);
        for (int index = 0; index < binSize - 1; index++) {
            RR[index] = cumulatedCount[index] / 2;
        }
    #endif

#endif

#endif
}

void Healpix_Correlation::compactAverageRDR(std::vector<myDouble>& bins, 
    std::vector<corrType>&DD, std::vector<corrType> &DR,
    std::vector<corrType>&RR, countType nk,
    int nTime,int _threadNum)
{
    myDouble stepTheta = bins[1];
    if(stepTheta < deltaTheta) {
        ModifySize(stepTheta);
    }
#ifdef HEALPIX_DEBUG
    cout << "\n********************\n";
    cout << "function NeighborVec\n";
    cout << "after the ModifySize in NeighborVec\n";
#endif
    int binSize = bins.size();
    DD.resize(binSize - 1);
    DR.resize(binSize - 1);

    std::vector<corrType> cumulatedCount;
    cumulatedCount.resize(binSize);
    for(int index = 0; index < binSize; index+= 1) {
#ifdef HEALPIX_DEBUG
        cout << "in the first for loop, index = " << index << endl;
#endif
        // using parallel
#ifdef PARALLEL
        countType nCount = NeighborsParallel(bins[index],_threadNum);
#else
        countType nCount = Neighbors(bins[index]);
#endif
#ifdef HEALPIX_DEBUG
        cout << "after Neighbors, nCount = " << nCount << endl;
#endif
        cumulatedCount[index] = nCount;
        if(index != 0) {
            DD[index - 1] = cumulatedCount[index] - cumulatedCount[index - 1];
        }
    }
#ifdef HEALPIX_DEBUG
    cout << "Function NeighborVec ends" << endl;
    cout << "\n********************\n";
#endif
}

int Healpix_Correlation::deltaTheta2Order(const myDouble & theta) {
    // since deltaTheta = pi/(4 * N_side),
    // 2^(order) = N_side = pi/(4 * deltaTheta)
    int temp;
    #ifdef GEOM_SPACE
    // assume endT = pi which is the max value
        myDouble q = pow(M_PI/theta, 1.0/(NUM_OF_POINTS - 1));
        myDouble dTtemp = theta * (q - 1) /DELTA_THETA_FACTOR;
        // temp = ceil(log2(pi/(4 * theta) *  DELTA_THETA_FACTOR));
        temp = ceil(log2(pi/(4*dTtemp)));
    #else
        temp = ceil(log2(pi/(4 * theta)));
    #endif
#ifdef HEALPIX_DEBUG
    cout << "the given deltaTheta is " << theta << endl;
    cout << "the calculated order is " << temp << endl;
    double dt = pi/(4 * pow(2,temp));
    cout << "the deltaTheta now is " << dt << endl;
#endif
    
    return temp;
}   


bool Healpix_Correlation::isStructureIsomorphic(const Healpix_Correlation & hp) {
    if(spannedAngle != hp.spannedAngle) {
        return false;
    }
    if (Order() != hp.Order()) {
        return false;
    }
    if (bottomPixel != hp.bottomPixel) {
        return false;
    }
    if (deltaTheta != hp.deltaTheta) {
        return false;
    }
    return true;
}


void Healpix_Correlation::NeighbourIntervalCountVector(typeIndex begin, 
    typeIndex end, const myDouble & theta, std::vector<corrType> & vec) 
{
      countType pair = 0;
#ifdef HEALPIX_DEBUG
    cout << "\n**************** " << endl;
    cout << "function Neighbors with theta: " << theta << endl;
    cout << "the bottomPixel is " << bottomPixel << endl;
#endif
    for(typeIndex pix = begin; pix < end; pix+= 1) {
        std::vector<int> pixset;
#ifdef HEALPIX_DEBUG
    cout << "in the first loop of Neighbors, pix: " << pix << endl;
#endif
        if((*this)[pix] == 0) {
            continue;
        }

#ifdef RANGESET_QUERY
    #ifdef APPEND_RANGESET_MODIFY 
        query_disc_pixel_Internal(pix, theta, std::ref(pixset),std::cref(*this));
    #else
        query_disc_pixel_Internal(pix, theta, std::ref(pixset));
    #endif
#else
        query_disc_pixel_InternalVec(pix,theta,std::ref(pixset),std::cref(*this));
#endif
#ifdef HEALPIX_DEBUG
    cout << "after the first loop of Neighbors" << endl;
#endif
        for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
            pair += (*this)[*iter];
        }
    }
#ifdef HEALPIX_DEBUG
    cout << "function Neighbors end with pair/2" << pair/2 << endl;
    cout << "\n**************** " << endl;
#endif
    std::lock_guard<std::mutex> guard(myMutex);
    vec.push_back(pair);
 }

void Healpix_Correlation::NeighbourIntervalCount(typeIndex begin, typeIndex end,
  const myDouble & theta, countType & pairCount) 
{
    countType pair = 0;
#ifdef HEALPIX_DEBUG
    cout << "\n**************** " << endl;
    cout << "function Neighbors with theta: " << theta << endl;
    cout << "the bottomPixel is " << bottomPixel << endl;
#endif
    for(typeIndex pix = begin; pix < end; pix+= 1) {
        countType temp = 0;
#ifdef HEALPIX_DEBUG
    cout << "in the first loop of Neighbors, pix: " << pix << endl;
#endif
        if((*this)[pix] == 0) {
            continue;
        }
#ifdef RANGESET_QUERY
    #ifdef APPEND_RANGESET_MODIFY 
        query_disc_pixel_Internal(pix, theta, std::ref(pixset),std::cref(*this));
    #else
        query_disc_pixel_Internal(pix, theta, std::ref(pixset));
    #endif
#else
     #ifdef QUERY_COUNT_DIRECTLY
        temp = query_disc_count(pix,theta,std::cref(*this));
        #ifdef MULTIPLY_PIXEL_VALUE
            temp *= (*this)[pix];
        #endif
        pair += temp;
        continue;
    #else 
        std::vector<int> pixset;
        query_disc_pixel_InternalVec(pix,theta,std::ref(pixset),std::cref(*this));
        for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
            pair += (*this)[*iter];
        }
    #endif
#endif
#ifdef HEALPIX_DEBUG
    cout << "after the first loop of Neighbors" << endl;
#endif

    }
#ifdef HEALPIX_DEBUG
    cout << "function Neighbors end with pair/2" << pair/2 << endl;
    cout << "\n**************** " << endl;
#endif
    std::lock_guard<std::mutex> guard(myMutex);
    pairCount += pair;
 }


 void Healpix_Correlation::crossNeighbourIntervalCount(typeIndex begin, 
    typeIndex end, const myDouble & theta,
    const Healpix_Correlation & hp, countType & pairCount) 
{
    countType pair = 0;
     // assume structure isomorphic
#ifdef HEALPIX_DEBUG
    cout << "\n**************** " << endl;
    cout << "function Neighbors with theta: " << theta << endl;
    cout << "the bottomPixel is " << bottomPixel << endl;
#endif
    for(typeIndex pix = begin; pix < end; pix+= 1) {
        if((*this)[pix] == 0) {
            continue;
        }
        countType temp = 0;
#ifdef HEALPIX_DEBUG
    cout << "in the first loop of Neighbors, pix: " << pix << endl;
#endif
        
#ifdef RANGESET_QUERY
    #ifdef APPEND_RANGESET_MODIFY 
        query_disc_pixel_Internal(pix, theta, std::ref(pixset),std::cref(hp));
    #else
        query_disc_pixel_Internal(pix, theta,std::ref(pixset));
    #endif
#else
    #ifdef QUERY_COUNT_DIRECTLY
        temp = query_disc_count(pix,theta,std::cref(hp));
        #ifdef MULTIPLY_PIXEL_VALUE
            temp *= (*this)[pix];
        #endif
        pair += temp;
        continue;
    #else 
        std::vector<int> pixset;
        query_disc_pixel_InternalVec(pix,theta,std::ref(pixset),std::cref(hp));
        for(std::vector<int>::iterator iter = pixset.begin(); iter != pixset.end(); iter+= 1) {
            pair += hp[*iter];
        }
    #endif
#endif
        
#ifdef HEALPIX_DEBUG
    cout << "after the first loop of Neighbors" << endl;
#endif

    }
#ifdef HEALPIX_DEBUG
    cout << "function Neighbors end with pair/2" << pair/2 << endl;
    cout << "\n**************** " << endl;
#endif
    std::lock_guard<std::mutex> guard(myMutex);
    pairCount += pair;
 }

// void Healpix_Correlation::crossNeighborNonCumIntervalCount(typeIndex begin, typeIndex end, 
//     const std::vector<myDouble> & bins,
//     std::vector<corrType> & neighborCount,
//     Healpix_Correlation & hp)
// {
//     std::vector<corrType> TempNeighborCount;
//     size_t len = neighborCount.size();
//     TempNeighborCount.resize(len);
//     for (typeIndex pix = begin; pix < end; pix++) {
//         query_disc_pixel
//     }
// }

void Healpix_Correlation::NeighbourNonCumIntervalCount(typeIndex begin, 
    typeIndex end, const std::vector<myDouble> & bins, 
    std::vector<corrType> & neighborCount, 
    const Healpix_Correlation & hp,
    int Factor)
{
    std::vector<corrType> TempNeighborCount;
    int len = neighborCount.size();
    TempNeighborCount.resize(len);
    for (typeIndex pix = begin; pix < end; pix++) {
        query_disc_pixel_Bins(pix, bins, TempNeighborCount,hp);
    }
    std::lock_guard<std::mutex> guard(myMutex);
    for (int i = 0; i < len; i++) {
        neighborCount[i] += TempNeighborCount[i];
    }

}            

void Healpix_Correlation::parallelSetPixValue(typeIndex begin, typeIndex end, countType & nonZeroCount, int value) 
{
    int count = 0;
    for (typeIndex pix = begin; pix < end; pix+= 1) {
        (*this)[pix] = value;
        if (value) {
            // nonZeroCount++;
            count++;
            // pointsCount++;
        }
    }
    std::lock_guard<std::mutex> guard(myMutex);
        pointsCount += count;
#ifdef RESET_VALUE_DEBUG
        cout << "\nbegin: " << begin << endl;
        cout << "\nend: " << end << endl;
        cout << "the stage variable\n";
        // cout << "nonZeroCount: " << nonZeroCount << endl;
        cout << "pointsCount: " << pointsCount << endl;
#endif 
}

void Healpix_Correlation::parallelScalePixValue(typeIndex begin, typeIndex end, myDouble & scale)
{
    for (typeIndex pix = begin; pix < end; pix += 1) {
        (*this)[pix] *= scale;
    }
}

void Healpix_Correlation::twoPCorrelation(myDouble dT, 
    myDouble endTheta, vec_coord & vC, 
    Healpix_Correlation & hp, 
    countType nk,
    int nTime,
    int mode,int _threadNum) 
{
#ifndef DISCONTINUOUS_SETTING 
    dT = deltaTheta;
#endif
    int vecLen = ceil((double)endTheta/dT) + 1;
    std::vector<myDouble> bins;
    bins.resize(vecLen);
    for(int index = 0; index < vecLen; index+= 1) {
        // the end may not be accurate
        bins[index] = index * dT;
    }
    std::vector<corrType> DD;
    std::vector<corrType> RR;
    std::vector<corrType> DR;
    NeighborVec(bins,DD,_threadNum);
    // hp.NeighborVec(bins,RR,_threadNum);
    averageRR(bins,RR,nk,nTime,_threadNum);
    // crossNeighborVec(bins,DR,hp,_threadNum);
    vC.resize(vecLen - 1);
    // Peebles
    myDouble nr = nk;
    myDouble nd = pointsCount;
    if (mode == 1) {
        
        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = nr * (nr - 1)/(nd * (nd - 1)) * DD[index]/RR[index] - 1;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    #ifndef SIGNAL_DATA_DRROUTPUT
        return;
    #endif
    }
    // Peebles 2
    averageCrossNeighborVec(bins,DR,nk,nTime,_threadNum);
    if (mode == 2) {
        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = 2 * nr / (nd - 1) * DD[index]/DR[index] - 1;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    }
    // Hamiliton
    if (mode == 3) {
        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = (myDouble)DD[index] * RR[index]/(DR[index] * DR[index]) * 4 * nd * nr /((nd - 1) * (nr - 1))- 1 + HAMILTON_BIAS;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    }
    // Landy
    if (mode == 4) {
        myDouble const1 = (nr * (nr - 1))/(nd * (nd - 1));
        myDouble const2 = (nr - 1) / nd;

        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = (DD[index] * const1 - const2 *  DR[index])/ RR[index] + 1 + LANDY_BIAS;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    }
#ifdef SIGNAL_DATA_DRROUTPUT
    cout << setw(10) << "bins" << "\t"  << setw(10) << "DD" << setw(15) << "RR" << setw(15) << "DR" <<  setw(20) << "ksi "<< endl;
    cout.precision(4);
    cout.setf(std::ios_base::fixed,std::ios_base::floatfield);
    for(int index = 0; index < vecLen - 1; index+= 1) {
        cout << setw(10) << bins[index] << "\t" <<  setw(10) << DD[index] << setw(15) << 
            RR[index] << setw(15) << DR[index] <<  setw(20) << vC[index].second << endl;
    }
#endif
}

void Healpix_Correlation::twoPCorrelation(myDouble dT, 
    myDouble endTheta, vec_coord & vC, 
    Healpix_Correlation & hp, std::vector<myDouble> & healpixMapInfo, 
    std::vector<myDouble> & randomHealpixInfo, 
    int mode,int _threadNum) 
{
#ifndef DISCONTINUOUS_SETTING 
    dT = deltaTheta;
#endif
    int vecLen = ceil((double)endTheta/dT) + 1;
    std::vector<myDouble> bins;
    bins.resize(vecLen);
    for(int index = 0; index < vecLen; index+= 1) {
        // the end may not be accurate
        bins[index] = index * dT;
    }
    std::vector<corrType> DD;
    std::vector<corrType> RR;
    std::vector<corrType> DR;
    NeighborVec(bins,DD,_threadNum);
    hp.NeighborVec(bins,RR,_threadNum);
    crossNeighborVec(bins,DR,hp,_threadNum);
    vC.resize(vecLen - 1);
    myDouble nr = hp.pointsCount;
    myDouble nd = pointsCount;
    if (mode == 1) {
        
        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = nr * (nr - 1)/(nd * (nd - 1)) * DD[index]/RR[index] - 1;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    }
    // Peebles 2
    if (mode == 2) {
        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = 2 * nr / (nd - 1) * DD[index]/DR[index] - 1;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    }
    // Hamiliton
    if (mode == 3) {
        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = (myDouble)DD[index] * RR[index]/(DR[index] * DR[index]) * 4 * nd * nr /((nd - 1) * (nr - 1))- 1 + HAMILTON_BIAS;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    }
    // Landy
    if (mode == 4) {
        myDouble const1 = (nr * (nr - 1))/(nd * (nd - 1));
        myDouble const2 = (nr - 1) / nd;

        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = (DD[index] * const1 - const2 *  DR[index])/ RR[index] + 1 + LANDY_BIAS;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    }

    // for the healpixMap information
    healpixMapInfo.resize(bottomPixel + 3);
    randomHealpixInfo.resize(bottomPixel + 2);

    healpixMapInfo[0] = (Order());
    healpixMapInfo[1] = pointsCount;
    randomHealpixInfo[0] = (hp.Order());
    for(int pix = 0; pix <= bottomPixel; pix+= 1) {
        healpixMapInfo[pix + 2] = (*this)[pix];
        randomHealpixInfo[pix + 1] = hp[pix];
    }
}


void Healpix_Correlation::averageTwoPCorrelation(myDouble dT, 
    myDouble endTheta, vec_coord & vC,
    int nTime, countType nk, 
    int mode,int _threadNum)
{
#ifndef DISCONTINUOUS_SETTING 
    dT = deltaTheta;
#endif
    int vecLen = ceil((double)endTheta/dT) + 1;
    std::vector<myDouble> bins;
    bins.resize(vecLen);
#ifdef AVERAGE_DEBUG
    cout << "************\n";
    cout << "in averageTwoPCorrelation,vecLen: " << vecLen << endl;
    cout << "=================\n";
#endif
    for(int index = 0; index < vecLen; index+= 1) {
        // the end may not be accurate
        bins[index] = index * dT;
    }
#ifdef AVERAGE_DEBUG
    cout << "************\n";
    cout << "in averageTwoPCorrelation, bins set up ok\n";
    cout << "=================\n";
#endif
    std::vector<corrType> DD;
    std::vector<corrType> RR;
    std::vector<corrType> DR;
    NeighborVec(bins,DD,_threadNum);
#ifdef RANGESET_DEBUG
    cout << "\n********************\nthe DD: " << endl;
    std::for_each(DD.begin(), DD.end(), [](const double & i){ cout << i << " ";});
    cout << endl;
#endif
#ifdef AVERAGE_DEBUG
    cout << "************\n";
    cout << "in averageTwoPCorrelation, DD set up ok\n";
    cout << "=================\n";
#endif
    averageRR(bins,RR,nk,nTime,_threadNum);
#ifdef RANGESET_DEBUG
    cout << "\n********************\nthe RR: " << endl;
    std::for_each(RR.begin(), RR.end(), [](const double & i){ cout << i << " ";});
    cout << endl;
#endif
#ifdef AVERAGE_DEBUG
    cout << "************\n";
    cout << "in averageTwoPCorrelation, RR set up ok\n";
    cout << "=================\n";
#endif
    vC.resize(vecLen - 1);
    myDouble nr = nk;
    myDouble nd = pointsCount;
#ifdef AVERAGE_DEBUG
    cout << "************\n";
    cout << "in averageTwoPCorrelation, initial set up ok\n";
    cout << "=================\n";
#endif
    if (mode == 1) {
        
        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = nr * (nr - 1)/(nd * (nd - 1)) * DD[index]/RR[index] - 1;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    }
    // Peebles 2
    else {
        averageCrossNeighborVec(bins,DR,nk,nTime,_threadNum);
#ifdef RANGESET_DEBUG
    cout << "\n********************\nthe DR: " << endl;
    std::for_each(DR.begin(), DR.end(), [](const double & i){ cout << i << " ";});
    cout << endl;
#endif
        if (mode == 2) {
            for (int index = 0; index < vecLen - 1; index+= 1) {
                myDouble ksi = 2 * nr / (nd - 1) * DD[index]/DR[index] - 1;
                myDouble var = bins[index + 1];
                twoDCoord tempPair(var,ksi);
                vC[index] = tempPair;
            }
        }
        // Hamiliton
        if (mode == 3) {
            for (int index = 0; index < vecLen - 1; index+= 1) {
                myDouble ksi = (myDouble)DD[index] * RR[index]/(DR[index] * DR[index]) * 4 * nd * nr /((nd - 1) * (nr - 1))- 1 + HAMILTON_BIAS;
                myDouble var = bins[index + 1];
                twoDCoord tempPair(var,ksi);
                vC[index] = tempPair;
            }
        }
        // Landy
        if (mode == 4) {
            myDouble const1 = (nr * (nr - 1))/(nd * (nd - 1));
            myDouble const2 = (nr - 1) / nd;

            for (int index = 0; index < vecLen - 1; index+= 1) {
                myDouble ksi = (DD[index] * const1 - const2 *  DR[index])/ RR[index] + 1 + LANDY_BIAS;
                myDouble var = bins[index + 1];
                twoDCoord tempPair(var,ksi);
                vC[index] = tempPair;
            }
        }
    }
#ifdef SIGNAL_DATA_DRROUTPUT

    essentialInfoOut(cout);
    cout << "averageTime: " << nTime << endl;
    cout << "nr : " << nk << endl;
    cout << setw(10) << "bins" << "\t"  << setw(15) << "DD" << setw(15) << "RR" << setw(15) << "DR" << endl;
    cout.precision(4);
    cout.setf(std::ios_base::fixed,std::ios_base::floatfield);
    for(int index = 0; index < vecLen - 1; index+= 1) {
        cout << setw(10) << bins[index] << "\t" <<  setw(15) << DD[index] << setw(15) << 
            RR[index] << setw(15) << DR[index] << endl;
    }
#endif
}
void Healpix_Correlation::averageTwoPCorrelation(myDouble dT, 
    myDouble endTheta, vec_coord & vC,int nTime, 
    countType nk, std::vector<myDouble> & healpixMapInfo, 
    int mode,int _threadNum) 
{
#ifndef DISCONTINUOUS_SETTING 
    dT = deltaTheta;
#endif
    int vecLen = ceil((double)endTheta/dT) + 1;
    std::vector<myDouble> bins;
    bins.resize(vecLen);
    for(int index = 0; index < vecLen; index+= 1) {
        // the end may not be accurate
        bins[index] = index * dT;
    }
    std::vector<corrType> DD;
    std::vector<corrType> RR;
    std::vector<corrType> DR;
    NeighborVec(bins,DD,_threadNum);
    averageRR(bins,RR,nk,nTime,_threadNum);
    vC.resize(vecLen - 1);
    myDouble nr = nk;
    myDouble nd = pointsCount;
    if (mode == 1) {
        
        for (int index = 0; index < vecLen - 1; index+= 1) {
            myDouble ksi = nr * (nr - 1)/(nd * (nd - 1)) * DD[index]/RR[index] - 1;
            myDouble var = bins[index + 1];
            twoDCoord tempPair(var,ksi);
            vC[index] = tempPair;
        }
    }
    // Peebles 2
    else {
        averageCrossNeighborVec(bins,DR,nk,nTime,_threadNum);
        if (mode == 2) {
            for (int index = 0; index < vecLen - 1; index+= 1) {
                myDouble ksi = 2 * nr / (nd - 1) * DD[index]/DR[index] - 1;
                myDouble var = bins[index + 1];
                twoDCoord tempPair(var,ksi);
                vC[index] = tempPair;
            }
        }
        // Hamiliton
        if (mode == 3) {
            for (int index = 0; index < vecLen - 1; index+= 1) {
                myDouble ksi = (myDouble)DD[index] * RR[index]/(DR[index] * DR[index]) * 4 * nd * nr /((nd - 1) * (nr - 1))- 1 + HAMILTON_BIAS;
                myDouble var = bins[index + 1];
                twoDCoord tempPair(var,ksi);
                vC[index] = tempPair;
            }
        }
        // Landy
        if (mode == 4) {
            myDouble const1 = (nr * (nr - 1))/(nd * (nd - 1));
            myDouble const2 = (nr - 1) / nd;

            for (int index = 0; index < vecLen - 1; index+= 1) {
                myDouble ksi = (DD[index] * const1 - const2 *  DR[index])/ RR[index] + 1 + LANDY_BIAS;
                myDouble var = bins[index + 1];
                twoDCoord tempPair(var,ksi);
                vC[index] = tempPair;
            }
        }
    }

    // to do here
    // exit(TODO_FUNCTION);
    // for the healpixMap information
    healpixMapInfo.resize(bottomPixel + 3);
    // randomHealpixInfo.resize(bottomPixel + 2);

    healpixMapInfo[0] = (Order());
    healpixMapInfo[1] = pointsCount;
    // randomHealpixInfo[0] = (hp.Order());
    for(int pix = 0; pix <= bottomPixel; pix+= 1) {
        healpixMapInfo[pix + 2] = (*this)[pix];
        // randomHealpixInfo[pix + 1] = hp[pix];
    }
}

void Healpix_Correlation::sweepAverageTwoPCorrelation(myDouble dT, myDouble endTheta,
            vec_coord & peebles1VC, vec_coord & peebles2VC, 
            vec_coord & hamiltonVC, vec_coord & landyVC,
            int nTime, countType nk,
            std::vector<myDouble> & healpixMapInfo,int _threadNum)
{
    sweepAverageTwoPCorrelation(dT,endTheta,peebles1VC,peebles2VC,
        hamiltonVC,landyVC,nTime,nk,_threadNum);
    // to do here
    // exit(TODO_FUNCTION);
    // for the healpixMap information
    healpixMapInfo.resize(bottomPixel + 3);
    // randomHealpixInfo.resize(bottomPixel + 2);

    healpixMapInfo[0] = (Order());
    healpixMapInfo[1] = pointsCount;
    // randomHealpixInfo[0] = (hp.Order());
    for(int pix = 0; pix <= bottomPixel; pix+= 1) {
        healpixMapInfo[pix + 2] = (*this)[pix];
        // randomHealpixInfo[pix + 2] = hp[pix];
    }
}

void Healpix_Correlation::sweepAverageTwoPCorrelation(myDouble dT, myDouble endTheta,
    vec_coord & peebles1VC, vec_coord & peebles2VC, 
    vec_coord & hamiltonVC, vec_coord & landyVC,
    int nTime, countType nk,
    std::vector<myDouble> & healpixMapInfo, ostream & os,
    int _threadNum)
{
    std::vector<corrType> DD;
    std::vector<corrType> RR;
    std::vector<corrType> DR;    
    std::vector<myDouble> bins;
    sweepAverageTwoPCorrelation(dT, endTheta, peebles1VC,
        peebles2VC, hamiltonVC, landyVC,nTime, nk, DD,
        RR, DR, bins, _threadNum);
    int ddVecLen = DD.size();
    essentialInfoOut(os);
    os << "averageTime: " << nTime << endl;
    os << "nr : " << nk << endl;
    os << setw(10) << "bins" << "\t"  << setw(15) << "DD" << setw(15) << "RR" << setw(15) << "DR" << endl;
    os.precision(4);
    os.setf(std::ios_base::fixed,std::ios_base::floatfield);
    for(int index = 0; index < ddVecLen; index+= 1) {
        os << setw(10) << bins[index] << "\t" <<  setw(15) << DD[index] << setw(15) << 
            RR[index] << setw(15) << DR[index] << endl;
    }    
    healpixMapInfo.resize(bottomPixel + 3);
    // randomHealpixInfo.resize(bottomPixel + 2);

    healpixMapInfo[0] = (Order());
    healpixMapInfo[1] = pointsCount;
    // randomHealpixInfo[0] = (hp.Order());
    for(int pix = 0; pix <= bottomPixel; pix+= 1) {
        healpixMapInfo[pix + 2] = (*this)[pix];
        // randomHealpixInfo[pix + 2] = hp[pix];
    }
}    

void Healpix_Correlation::sweepAverageTwoPCorrelation(myDouble dT, myDouble endTheta,
    vec_coord & peebles1VC, vec_coord & peebles2VC, 
    vec_coord & hamiltonVC, vec_coord & landyVC,
    int nTime, countType nk,
    int _threadNum)
{
    std::vector<corrType> DD;
    std::vector<corrType> RR;
    std::vector<corrType> DR;
    std::vector<myDouble> bins;
    sweepAverageTwoPCorrelation(dT, endTheta, peebles1VC,
        peebles2VC, hamiltonVC, landyVC,nTime, nk, DD,
        RR, DR,bins,_threadNum);
}

void Healpix_Correlation::sweepAverageTwoPCorrelation(myDouble dT, myDouble endTheta,
    vec_coord & peebles1VC, vec_coord & peebles2VC, 
    vec_coord & hamiltonVC, vec_coord & landyVC,
    int nTime, countType nk, std::vector<corrType> & DD,
    std::vector<corrType> & RR,
    std::vector<corrType> & DR, std::vector<myDouble> & bins,
    int _threadNum, const std::string & logFolder)
{
#ifndef DISCONTINUOUS_SETTING 
    dT = deltaTheta;
#endif
    int vecLen;
#ifdef LINSPACE
    vecLen = ceil((double)endTheta/dT) + 1;
    bins.resize(vecLen);
    for(int index = 0; index < vecLen; index+= 1) {
        // the end may not be accurate
        bins[index] = index * dT;
    }
#endif
#ifdef GEOM_SPACE
    BinsGenerator binGen;
    vecLen = NUM_OF_POINTS;
    binGen.geomspace(bins,dT,endTheta,vecLen);
#endif

    NeighborVec(bins,DD,_threadNum);
    averageRR(bins,RR,nk,nTime,_threadNum);
    peebles1VC.resize(vecLen - 1);
    peebles2VC.resize(vecLen - 1);
    hamiltonVC.resize(vecLen - 1);
    landyVC.resize(vecLen - 1);
    myDouble nr = nk;
    myDouble nd = pointsCount;

    for (int index = 0; index < vecLen - 1; index+= 1) {
        myDouble ksi = nr * (nr - 1)/(nd * (nd - 1)) * DD[index]/RR[index] - 1;
        myDouble var = bins[index + 1];
        twoDCoord tempPair(var,ksi);
        peebles1VC[index] = tempPair;
    }
// Peebles 2
    averageCrossNeighborVec(bins,DR,nk,nTime,_threadNum);
#ifdef OUTPUT_DDRR

    cout << setw(9) << "bins:" << setw(10) << "DD" << setw(10) << "DR" << setw(10) << "RR\n";
    for(int i = 0; i < vecLen - 1; i+= 1) {
        cout << setw(9) << setprecision(4) << bins[i + 1] << setw(10) 
            << setprecision(4) << DD[i] << setw(10) << setprecision(4) << DR[i]
            << setw(10)<< setprecision(4) << RR[i] << endl;
    }
#endif
    for (int index = 0; index < vecLen - 1; index+= 1) {
        myDouble ksi = 2 * nr / (nd - 1) * DD[index]/DR[index] - 1;
        myDouble var = bins[index + 1];
        twoDCoord tempPair(var,ksi);
        peebles2VC[index] = tempPair;
    }

// Hamiliton

    for (int index = 0; index < vecLen - 1; index+= 1) {
        myDouble ksi = (myDouble)DD[index] * RR[index]/(DR[index] * DR[index]) * 4 * nd * nr /((nd - 1) * (nr - 1))- 1 + HAMILTON_BIAS;
        myDouble var = bins[index + 1];
        twoDCoord tempPair(var,ksi);
        hamiltonVC[index] = tempPair;
    }

// Landy

    myDouble const1 = (nr * (nr - 1))/(nd * (nd - 1));
    myDouble const2 = (nr - 1) / nd;

    for (int index = 0; index < vecLen - 1; index+= 1) {
        myDouble ksi = (DD[index] * const1 - const2 *  DR[index])/ RR[index] + 1 + LANDY_BIAS;
        myDouble var = bins[index + 1];
        twoDCoord tempPair(var,ksi);
        landyVC[index] = tempPair;
    }

#ifdef SIGNAL_DATA_DRROUTPUT
    // add correlationCount
    correlationCount++;
// try to use spdlog
    string loggerName = essentialInfoOut("");
    auto calcuLogger = spdlog::get(loggerName);
    calcuLogger -> info("averageTime: {}", nTime);
    calcuLogger -> info("nr: {}", nk);

    // essentialInfoOut(cout);
    // cout << "averageTime: " << nTime << endl;
    // cout << "nr : " << nk << endl;
    // cout << setw(10) << "bins" << "\t"  << setw(15) << "DD" << setw(15) << "RR" << setw(15) << "DR" << endl;
    // cout.precision(4);
    // cout.setf(std::ios_base::fixed,std::ios_base::floatfield);
    // for(int index = 0; index < vecLen - 1; index+= 1) {
    //     cout << setw(10) << bins[index] << "\t" <<  setw(15) << DD[index] << setw(15) << 
    //         RR[index] << setw(15) << DR[index] << endl;
    // }
    string drLoggerName = loggerName + "DRLogger";
    auto drLogger = spdlog::basic_logger_mt(drLoggerName, _logDir + loggerName +  "/DR.txt");
    // use ostream
    string fileName = _logDir + loggerName +  "/DR.txt";
    std::ofstream outFile;
    outFile.open(fileName);
    outFile << setw(10) << "bins" << "\t"  << setw(15) << "DD" << setw(15) << "RR" << setw(15) << "DR" << endl;
    for(int index = 0; index < vecLen - 1; index+= 1) {
        outFile << setw(10) << bins[index] << "\t" <<  setw(15) << DD[index] << setw(15) << 
            RR[index] << setw(15) << DR[index] << endl;
    }
    drLogger -> info("Test logger");
#endif


}    

void Healpix_Correlation::pointProcessSweep2PCF(myDouble dT, myDouble endTheta,
    const vec_3DCoord & vecData,
    vec_coord & peebles1VC, vec_coord & peebles2VC, 
    vec_coord & hamiltonVC, vec_coord & landyVC,
    int nTime, countType scale, std::vector<myDouble> & healpixMapInfo, 
    ostream & os,
    int _threadNum)
{
    SetPixValue(vecData);
    countType nk = pointsCount * scale;
    sweepAverageTwoPCorrelation(dT, endTheta,peebles1VC,
    peebles2VC, hamiltonVC, landyVC, nTime, nk,
    healpixMapInfo,os,_threadNum);


#ifdef OFFSET_POISSON
    Healpix_Correlation map2(dT,hThetaEnd,hThetaBeg);
    map2.SetPixValue(pointsCount);
    vec_coord refPeebles1,refPeebles2,refHamilton,refLandy;
    map2.sweepAverageTwoPCorrelation(dT,endTheta,refPeebles1,refPeebles2,refHamilton,
        refLandy,nTime,nk);
    typeIndex _len = refPeebles1.size();
    for(typeIndex _i = 0; _i < _len; _i++) {
        peebles1VC[_i].second -= refPeebles1[_i].second;
        peebles2VC[_i].second -= refPeebles2[_i].second;
        hamiltonVC[_i].second -= refHamilton[_i].second;
        landyVC[_i].second -= refLandy[_i].second;
    }
#endif

}    

void Healpix_Correlation::pointProcessPrintPoissonSweep2PCF(myDouble dT, myDouble endTheta,
    const vec_3DCoord & vecData,
    vec_coord & peebles1VC, vec_coord & peebles2VC, 
    vec_coord & hamiltonVC, vec_coord & landyVC,
    vec_coord & poissonPeebles1VC, vec_coord & poissonPeebles2VC,
    vec_coord & poissonHamiltonVC, vec_coord & poissonLandyVC,
    int nTime, countType scale, std::vector<myDouble> & healpixMapInfo, 
    ostream & os,
    int _threadNum)
{
    SetPixValue(vecData);
    countType nk = pointsCount * scale;
    sweepAverageTwoPCorrelation(dT, endTheta,peebles1VC,
    peebles2VC, hamiltonVC, landyVC, nTime, nk,
    healpixMapInfo,os,_threadNum);


    Healpix_Correlation map2(dT,hThetaEnd,hThetaBeg);
    map2.SetPixValue(pointsCount);
    vec_coord refPeebles1,refPeebles2,refHamilton,refLandy;
    map2.sweepAverageTwoPCorrelation(dT,endTheta,poissonPeebles1VC,
        poissonPeebles2VC,poissonHamiltonVC,
        poissonLandyVC,nTime,nk);

}    

void Healpix_Correlation::outMapAndRandom(std::vector<myDouble> & healpixMapInfo, std::vector<myDouble> & randomMapInfo)
{

    randomMapInfo.resize(bottomPixel + 3);
    Healpix_Correlation map2(deltaTheta, hThetaEnd,hThetaBeg,hPhiEnd,hPhiBeg);
    map2.SetPixValue(pointsCount);
    if (Order() != map2.Order()) {
        cout << "The generated random map in outMapAndRandom is not structure isomorphic to *this!\n";
        exit(NOT_STRUCTURE_ISOM_HEALPIX);
    }

    randomMapInfo[0] = (Order());
    randomMapInfo[1] = pointsCount;
    for (int pix = 0; pix <= bottomPixel; pix += 1) {
        randomMapInfo[pix + 2] = map2[pix];
    }
}

myDouble Healpix_Correlation::signalStrength2Val(myDouble _strength, myDouble _threshold,
    int mode)
{   
    if(mode == SIMPLE_THRESHOLD || mode == SIMPLE_INCREMENT) {
        return (_strength > _threshold) ? 1 : 0;
    }
    else if(mode == SIMPLE_POWER_mW_ADD) {
        // convert to mW
        return (_strength > _threshold) ?  pow(10, _strength/10) : 0;
    }
    else if(mode == SIMPLE_POWER_W_ADD) {
        // convert to W
        return (_strength > _threshold) ?  pow(10, _strength/10) * 1000 : 0;
    }
    else {
        cout << "no other mode supported yet\n";
        exit(NO_SUPPORT);
    }
}

bool Healpix_Correlation::isAllZero() const {
    bool allZeroFlag = true;
    int nsideTemp = Nside();
    int endPixel = 12 * nsideTemp * nsideTemp - 1;
    for(int pix = 0; pix <= endPixel; pix+= 1) {
        if(fabs((*this)[pix]) > MY_EPSILON) {
            cout << "pix: " << pix << " is not zero, with value: " << (*this)[pix] << endl;
            allZeroFlag = false;
        }
    }
    return allZeroFlag;
}

int Healpix_Correlation::indexRingAbove(double z) const {
    return ring_above(z);
}

pointing Healpix_Correlation::publicPix2Ang(int pixel) const
{
    return pix2ang(pixel);
}

void Healpix_Correlation::publicGetRingInfoSmall(int ring, int &startpix, int &ringpix, bool &shifted) const
{
    get_ring_info_small(ring,startpix,ringpix,shifted);
}

void Healpix_Correlation::publicPix2ZPhi(int pix, double &z, double &phi) const
{
    pix2zphi(pix,z,phi);
}

int Healpix_Correlation::publicTheta2Ring(const myDouble & theta) const
{
    return pix2ring(zphi2pix(cos(theta),0));
}

void Healpix_Correlation::appendNonZeroRange(rangeset<int> & pixset, 
   const Healpix_Correlation & hp,  typeIndex begin, typeIndex end, int ignore) const {
    bool continuousFlag = false;
    // the begin and end of the current continuous interval [currLeft, currRIght)
    int currLeft = begin, currRight = begin;
    if (ignore >= 0) {
        cout << "appendNonZeroRange does not support ignore\n";
        exit(DO_NOT_SUPPORT_IGNORE);
    }
    for(typeIndex pix = begin; pix < end; pix+= 1) {
        if((hp)[pix] != 0) {
            if(continuousFlag) {
                currRight+= 1;
            }
            else {
                continuousFlag = true;
                currLeft = pix;
                currRight = pix + 1;
            }
        }
        else {
            // the continuous intervall ends here
            // append this interval and set the flag to false
            if(continuousFlag) {
                pixset.append(currLeft,currRight);
                continuousFlag = false;
                // reset the currLeft, currRight
                currLeft = pix + 1;
                currRight = pix + 1;
            }
        }
    }
    // it is possible, that the last continuous intervall lasts to the end
    // and unable to append that interval because it's out of the loop
    if(continuousFlag) {
        pixset.append(currLeft,currRight);
    }
}

void Healpix_Correlation::VecAppendNonZeroPoint(vector<int> & pixVec,
const Healpix_Correlation & hp, typeIndex begin, typeIndex end, int ignore) const
{
    // need to modify the end, and begin
    begin = max<int>(0,begin);
    end = min<int>(end,bottomPixel + 1);
    for(typeIndex pix = begin; pix < end; pix+= 1) {
        if((hp)[pix] != 0 && (int)pix != ignore){
            // if(pixVec.empty()) {
            //     pixVec.resize(1);
            // }
            pixVec.push_back(pix);
        }
    }
}

countType Healpix_Correlation::IntervallNonZeroPoint(const Healpix_Correlation & hp, 
        typeIndex begin, typeIndex end, int ignore) const 
{
    countType pCount = 0;
    for(typeIndex pix = begin; pix < end; pix+= 1) {
        if((hp)[pix] != 0 && (int)pix != ignore) {
            pCount += (hp)[pix];
        }
    }
    return pCount;
}

ostream & operator<<(ostream & os, 
    const Healpix_Correlation & hp)
{ 
    os << hp.Order() << endl;
    os << hp.GetPointsCount() << endl;
    sizeType Bottom = hp.GetBottomPixel();
    for(typeIndex i = 0; i <= Bottom; i+= 1) {
        os << hp[i] << endl;
    }
    return os;
}

void Healpix_Correlation::outRingInfo(ostream & os) const
{
    int startpix,ringpix;
    bool shifted;
    for (int iz = 1; iz <= bottomRing; iz++) {
        get_ring_info_small(iz,startpix,ringpix,shifted);
        os << "ring: "<<  iz << "\tstart: " << startpix << "\tend: " << startpix + ringpix - 1
            << "  shifted: " << shifted  << endl;
    }
}

void Healpix_Correlation::essentialInfoOut(ostream & outFile) const
{
    outFile << "pointsCount: " << pointsCount << endl;
    outFile << "deltaTheta: " << deltaTheta << endl;
    outFile << "spannedAngle: " << spannedAngle << endl;
    outFile << "bottomRing: " << bottomRing << endl;
    outFile << "bottomPixel: " << bottomPixel << endl;
}


string Healpix_Correlation::essentialInfoOut(const string & logFolder) const
{
    std::string logFolderName;
    if (logFolder == "") {
        // use the default
        logFolderName = "default" + std::to_string(correlationCount);
    }
    else {
        logFolderName = logFolder;
    }
    auto file_logger = spdlog::basic_logger_mt(logFolderName, _logDir + logFolderName + "/calcuInfo.txt");
    // register logger
    file_logger -> set_pattern("%v");
    file_logger ->info("pointsCount: {}",pointsCount);
    file_logger ->info("deltaTheta: {}", deltaTheta);
    file_logger -> info("spannedAngle: {}",spannedAngle);
    file_logger -> info("bottomRing: {}", bottomRing);
    file_logger -> info("bottomPixel: {}",bottomPixel);

    return logFolderName;
}

void convertDataVecIndex(const typeIndex & dataIndex, typeIndex & thetaIndex, typeIndex & phiIndex, 
    const sizeType & thetaLen, const sizeType & phiLen)
{
    // need to know that the theta is sorted in ascending order,
    // so the tempthetaIndex is what makes the formula dataIndex = phiLen * tempthetaIndex + phiIndex
    phiIndex = dataIndex % phiLen;
    typeIndex tempthetaIndex = dataIndex / phiLen;
    thetaIndex = thetaLen - 1 - tempthetaIndex;
    if (thetaIndex >= thetaLen || phiIndex >= phiLen || thetaIndex < 0 || phiIndex < 0) {
        cout << "Invalid index in convertDataVecIndex\n";
        cout << "thetaIndex: " << thetaIndex << "\t thetaLen: " << thetaLen << endl;
        cout << "phiIndex: " << phiIndex << "\t phiLen: " << phiLen << endl;
        exit(INDEX_OUT_OF_BOUND);
    }
}

myDouble convertDeg2Radian(const myDouble & degree)
{
    return degree/180.0 * M_PI;
}
