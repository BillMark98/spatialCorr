#include "binsGenerator.h"

void BinsGenerator::logspace(std::vector<myDouble> & bins, myDouble start, myDouble end, sizeType num, myDouble base)
{
    bins.resize(num);
    myDouble step = (end - start)/(num - 1);
    bins[0] = pow(base,start);
    myDouble q = pow(base,step);
    for (typeIndex index = 1; index < num; index++) {
        bins[index] = bins[index - 1] * q;
    }
}
// create linear spaced bins, from start to end, total num points
void BinsGenerator::linspace(std::vector<myDouble> & bins, myDouble start, myDouble end, sizeType num)
{
    bins.resize(num);
    myDouble step = (end - start)/(num - 1);
    bins[0] = start;
    for (typeIndex index = 1; index < num; index++) {
        bins[index] = bins[index - 1] + step;
    }
}
// create geometric series, a_{1} = start, a_{num} = end
void BinsGenerator::geomspace(std::vector<myDouble> & bins, myDouble start, myDouble end, sizeType num)
{
    bins.resize(num);
    if (fabs(start) < MY_EPSILON) {
        std::cout << "start to small, cant generate geomspace\n";
        exit(DIVIDE_BY_ZERO);
    }
    myDouble q = pow((end/start), 1.0/(num-1));
    bins[0] = start;
    for (typeIndex index = 1; index < num; index++) {
        bins[index] = bins[index - 1] * q;
    }
}
