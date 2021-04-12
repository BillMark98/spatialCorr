#include "myDefineConst.h"
#include <vector>
class BinsGenerator {
    public:
        // create a bin from base^(start), to base^(end), with logarithmically num spaced points
        void logspace(std::vector<myDouble> & bins, myDouble start, myDouble end, sizeType num, myDouble base);
        // create linear spaced bins, from start to end, total num points
        void linspace(std::vector<myDouble> & bins, myDouble start, myDouble end, sizeType num);
        // create geometric series, a_{1} = start, a_{num} = end
        void geomspace(std::vector<myDouble> & bins, myDouble start, myDouble end, sizeType num);
};