#ifndef _MATHDOUBLE_H
#define _MATHDOUBLE_H

#define RANGE_ERROR 1e-6
#include <cmath>
#include <ostream>
#define MY_INFINITE 1e10
class MathDouble {
    private:
        double data;
    public:
        MathDouble(double _num = 0.) : data(_num) {}
        bool operator < (const MathDouble & v1) const;
        bool operator==(const MathDouble & v1) const;
        bool operator!=(const MathDouble & v1) const;
        bool operator>(const MathDouble & v1) const;
        bool operator <= (const MathDouble & v1) const;
        bool operator >= (const MathDouble & v1) const;
        MathDouble operator+(const MathDouble & v1) const;
        MathDouble operator-(const MathDouble & v1) const;
        MathDouble operator*(const MathDouble & v1) const;
        MathDouble operator/(const MathDouble & v1) const;

        void operator += (const MathDouble & v1);
        void operator -= (const MathDouble & v1);
        void operator *= (const MathDouble & v1);
        void operator /= (const MathDouble & v1);

        MathDouble operator+(double d1) const;
        MathDouble operator-(double d1) const;
        MathDouble operator*(double d1) const;
        MathDouble operator/(double d1) const;

        void operator += (double d1);
        void operator -= (double d1);
        void operator *= (double d1);
        void operator /= (double d1);

        void operator++();
        void operator--();
        explicit operator double() const;
        explicit operator int() const;

        friend double operator+(double, const MathDouble & v1);
        friend double operator-(double, const MathDouble & v1);
        friend double operator*(double, const MathDouble & v1);
        friend double operator/(double, const MathDouble & v1);

        friend double operator-(const MathDouble & v1);


        friend double cos(const MathDouble & v1);
        friend double sin(const MathDouble & v1);
        friend double sqrt(const MathDouble & v1);
        friend double atan2(const MathDouble & v1, const MathDouble & v2);
        friend double acos(const MathDouble & v1);

        friend double fabs(const MathDouble & v1);
        
        friend std::ostream & operator<<(std::ostream & os, const MathDouble & v1);
        friend std::istream & operator>>(std::istream & is, const MathDouble & v1);

};
#endif