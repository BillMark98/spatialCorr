#include "mathDouble.h"

bool MathDouble::operator < (const MathDouble & v1) const
{
    return data < v1.data - RANGE_ERROR;
}

bool MathDouble::operator==(const MathDouble & v1) const
{
    return fabs(data - v1.data) < RANGE_ERROR;
}

bool MathDouble::operator!=(const MathDouble & v1) const
{
    return fabs(data - v1.data) > RANGE_ERROR;
}

bool MathDouble::operator>(const MathDouble & v1) const
{
    return data > v1.data + RANGE_ERROR;
}

bool MathDouble::operator <= (const MathDouble & v1) const
{
    return ((*this) < v1 || (*this) == v1);
}

bool MathDouble::operator >= (const MathDouble & v1) const
{
    return ((*this) > v1 || (*this) == v1);
}

MathDouble MathDouble::operator+(const MathDouble & v1) const
{
    MathDouble temp(data + v1.data);
    return temp;
}

MathDouble MathDouble::operator-(const MathDouble & v1) const
{
    MathDouble temp(data - v1.data);
    return temp;
}

MathDouble MathDouble::operator*(const MathDouble & v1) const
{
    MathDouble temp(data * v1.data);
    return temp;
}

MathDouble MathDouble::operator/(const MathDouble & v1) const
{
    double val;
    if(v1.data != 0) {
        val = data / v1.data;
    }
    else {
        val = MY_INFINITE;
    }
    return MathDouble(val);
}

void MathDouble::operator += (const MathDouble & v1)
{
    data += v1.data;
}
void MathDouble::operator -= (const MathDouble & v1)
{
    data -= v1.data;
}
void MathDouble::operator *= (const MathDouble & v1)
{
    data *= v1.data;
}

void MathDouble::operator /= (const MathDouble & v1)
{
    if(v1 != 0) {
        data /= v1.data;
    }
    else {
        data = MY_INFINITE;
    }
}

MathDouble MathDouble::operator+(double d1) const
{
    MathDouble temp(data + d1);
    return temp;
}

MathDouble MathDouble::operator-(double d1) const
{
    MathDouble temp(data - d1);
    return temp;
}

MathDouble MathDouble::operator*(double d1) const
{
    MathDouble temp(data * d1);
    return temp;
}

MathDouble MathDouble::operator/(double d1) const
{
    double val;
    if(fabs(d1) < RANGE_ERROR) {
        val = MY_INFINITE;
    }
    else {
        val = data / d1;
    }
    return MathDouble(val);
}

void MathDouble::operator += (double d1)
{
    data += d1;
}

void MathDouble::operator -= (double d1)
{
    data -= d1;
}

void MathDouble::operator *= (double d1)
{
    data *= d1;
}

void MathDouble::operator /= (double d1)
{
    if(fabs(d1) < RANGE_ERROR) {
        data = MY_INFINITE;
    }
    else {
        data /= d1;
    }
    
}

void MathDouble::operator++()
{
    data += 1;
}

void MathDouble::operator--()
{
    data -= 1;
}
MathDouble::operator double() const
{
    return data;
}

MathDouble::operator int() const
{
    return int(data);
}


double operator+(double d1, const MathDouble & v1)
{
    return d1 + v1.data;
}

double operator-(double d1, const MathDouble & v1)
{
    return d1 - v1.data;
}

double operator*(double d1, const MathDouble & v1)
{
    return d1 * v1.data;
}

double operator/(double d1, const MathDouble & v1)
{
    if(v1 != 0) 
        return d1 / v1.data;
    else 
        return MY_INFINITE;
}

double operator-(const MathDouble & v1)
{
    return -v1.data;
}

double cos(const MathDouble & v1)
{
    return std::cos(v1.data);
}

double sin(const MathDouble & v1)
{
    return std::sin(v1.data);
}

double sqrt(const MathDouble & v1)
{
    return std::sqrt(v1.data);
}

double atan2(const MathDouble & v1, const MathDouble & v2)
{
    return std::atan2(v1.data,v2.data);
}

double acos(const MathDouble & v1)
{
    return std::acos(v1.data);
}

double fabs(const MathDouble & v1)
{
    return std::fabs(v1.data);
}

std::ostream & operator<<(std::ostream & os, 
    const MathDouble & v1)
{
    os << v1.data;
    return os;
}    

std::istream & operator>>(std::istream & is, 
    const MathDouble & v1)
{
    is >> v1.data;
    return is;
}



