#ifndef __MY_DEFINE_CONST
#define __MY_DEFINE_CONST

#define MODEERROR 31
#define NOT_IN_THE_SAME_SPHERE 32
#define NO_POINTS_ON_SPHERE 33
#define TREENODE_EMPTY 34
#define FILE_OPEN_ERROR 35
#define HEALPIX_BOUND_NOT_SET 36
#define NOT_STRUCTURE_ISOM_HEALPIX 37
#define EMPTY_HEALPIX_MAP 38
#define DELETE_POINT_ERROR 39
#define QUERY_DISC_ERROR 40
#define DO_NOT_SUPPORT_IGNORE 41   // used for the ignore certain pixel when query theta Disc
#define FILE_NUM_NOT_MATCHING 42   // used for Neighbors, outputting desired vectors
#define DIVIDE_BY_ZERO 43 // for the case dividing by zero
#define ANGLEBEGIN_LARGER_THAN_ANGLEEND 44
#define DONT_CALL_ME 45
#define NO_SUPPORT 46
#define INITIAL_MAP_NONZERO 47
#define INDEX_OUT_OF_BOUND 48
#define TODO_FUNCTION 99

#define DELTA_THETA 0.2
#define MY_EPSILON 1e-5
#define POISSON_POINT_COUNT 200
#define HEALPIX_RANDOM_POINT_COUNT 500
#define MY_DEFAULT_STRING "no"
#define PARALLEL
#define MAX_THREAD 6
// indicating that given deltaTheta (for calculation of 2pCF), use the given one
// instead of using the one of the healpix map
#define DISCONTINUOUS_SETTING 
// if set use the geom spaced bins
// #define GEOM_SPACE
#ifdef GEOM_SPACE
    #define NUM_OF_POINTS 25
    #define DELTA_THETA_FACTOR 2
#endif

#define LINSPACE


// #define AVERAGE_DEBUG
#define LANDY_BIAS -0.0
#define HAMILTON_BIAS -0.0
#define AVERAGE_TIME 10

// if set, will square the val in the SetPixValue for the analysis of data
// #define QUADRATIC
#define FACTOR_SIGVAL 1


#define MAX_LOOP_TIME 10000
// in the Healpix_Correlation::SetPixValue, used if has a constraint on phi
#define HEALPIX_MAX_LOOP_TIME 1000000
#define MIN_RSS -5000
// the number of iteration to calculate alm used in Healpix_Correlation::sphericalHarmonicALM
#define NUM_ITER 2

// if set, first average all the random set map to calculate further
// see Healpix_Corrleation::averageRR
#define AVERAGE_RANDOM_MAP

// if set use the similar geometry for calculating RR, DR
#define SIMILAR_GEOMETRY
// if set, will calculate the RR, DR from several realizations and average them
// #define AVERAGE_RR_DR

// #define AVERAGE_RMAP_DEBUG
// if set , will multiply the (hp)[pix]
#define MULTIPLY_PIXEL_VALUE
#define ALM_TRY

// if set, the neighbor returned is scaled by 2, need to divide it at
// the upper level of the calling function
#define NEIGHBOR_SCALE_2
// if set, test if the calculated neighborpair is odd
// #define PAIR_IS_ODD_TEST

// if defined, use the rangeset based query_disc_internal
// #define RANGESET_QUERY
// if defined, use the modified version of appending element to a rangeset
// #define APPEND_RANGESET_MODIFY 

// if set, the query function will directly gives back the pointCount
// instead of the vector which stores the pixels
#define QUERY_COUNT_DIRECTLY
// #define RANGESET_TEST

// if set, use the non-cum version of count
#define NON_CUM_COUNT

// for the noncumulative version of query disc debug
// #define QUERY_NON_CUMU_DEBUG

// for average debug:
// #define OUTPUT_PIX_NEIGHBOR_COUNT
#define TO_PRINT_THETA 0.2

// for the resetValue of Healpix_Correlation debug
// #define RESET_VALUE_DEBUG

// #define DATA_VALUE_SET_DEBUG

// for fixing the issues curve at zero very large
// #define HEALPIX_ZERO_DEBUG
// set the time seed for random num generator
// #define TIME_FIXED 17

// #define POISSON_DEBUG
// #define FINDPAIR_DEBUG
// #define HEALPIX_DEBUG
// #define HEALPIX_PARALLEL_DEBUG

// #define MATERN_PARALLEL_POISSON_DEBUG

// if set, will output information about the generated matern process
#define MATERN_INFO_OUT

// #define ALL_ONE_DEBUG

#define QUERY_ZERO
// #define JOIN_THREAD_MOVE

// #define SIGNAL_SET_DEBUG

// #define OUTPUT_DDRR
#ifdef OUTPUT_DDRR
    #include <iomanip>
    using std::setw;
    using std::setprecision;
#endif



// for the implementation in PoissonPoint2d::GenerateNIndices
#define VEC_SET_THRESHOLD 3

// random shuffle the generated random number
#define RANDOM_SHUFFLE
// #define RANGESET_DEBUG
#ifdef RANGESET_DEBUG
    #include <algorithm>
#endif

#define MATERN_PARALLEL

// if defined, will subtract a simulated poisson process during
// calculation of 2pcf
// #define OFFSET_POISSON
#ifndef OFFSET_POISSON
    // #define PRINT_POISSON
#endif
// #define MATERN_DEBUG

// which means that the conversion from signal strength to pixel value
// is simply calculated by testing if the value is above that threshold
#define SIMPLE_THRESHOLD 1

// if set, the value of the pixel will increase by one each time a data point
// lies in that cell
#define SIMPLE_INCREMENT 2

// if set, the value of the pixel will be the sum of the corresponding power in mW, 
// if it exceeds the threshold, otherwise set to 0
#define SIMPLE_POWER_mW_ADD 3

// if set, the value of the pixel will be the sum of the corresponding power in W, 
// if it exceeds the threshold, otherwise set to 0
#define SIMPLE_POWER_W_ADD 4


#include <thread>
#include <mutex>
#include <functional>







#include "mathDouble.h"
#include <iostream> 
#include <fstream>

using std::endl;
using std::cout;
const double pi = M_PI;
const double twopi = 2 * M_PI;
const double inv_twopi=1.0/twopi;
typedef size_t typeIndex;
typedef unsigned long sizeType;

// #define signalMode
#define SIGNAL_DATA_DRROUTPUT
#ifdef SIGNAL_DATA_DRROUTPUT
    #include <iomanip>
    using std::setw;
    using std::setprecision;
#endif
#ifdef signalMode
    typedef MathDouble countType;
    typedef MathDouble myDouble;
    typedef MathDouble corrType;
#else
    typedef double countType;
    typedef double myDouble;
    typedef double corrType;
    typedef myDouble healpixMapType;
#endif


#include <string>
typedef std::string filePathType;

// not const T & since will return reference to a temporary
template <typename T>
const T square(T v1) {
    // cout << "in square: " << v1 << endl;
    return v1 * v1;
}



#endif