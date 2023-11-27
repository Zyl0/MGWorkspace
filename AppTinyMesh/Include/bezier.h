#ifndef BEZIER_H
#define BEZIER_H

#include "mesh.h"
#include "mathematics.h"

#include <stdlib.h>
#include <vector>


static size_t pascalTriangleHeight = 0;
static std::vector<double> pascalTriangle;
void reservePascalTriangle(size_t n);

static std::vector<size_t> factorials;

void reserveFactorials(size_t max);
inline size_t factorial(size_t N) {return factorials[N];}

inline double pascalCni(size_t N, size_t I) 
{
    return pascalTriangle[I + ((N*(N+1)) >> 1)];
}

inline double factorialCni(size_t N, size_t I)
{
    return (((double)factorial(N))/(factorial(I) * factorial(N - I)));
}

inline double Cni(size_t N, size_t I)
{
    //Using pascal triangle
    return pascalCni(N,I);

    //Using precomputed factorials
    //return factorialCni(N, I);
}

inline double B(double U, size_t I, size_t N) 
{
    return (Cni(N, I) * pow(U, I) * pow(1 - U, N - I));
}

class BezierCurve
{
private:
    std::vector<Vector> controlPoints;
public:
    BezierCurve(const std::vector<Vector>& points);
    Vector SamplePoint(double U) const;
};

class BezierSurface
{
private:

    size_t n, m;/* data */
    std::vector<Vector> controlPoints;

    inline const Vector &cp(size_t I, size_t J) const {return controlPoints[I * n + J];}
public:
    Vector SamplePoint(double U, double V) const;
    Vector SampleNormal(double U, double V) const;
    BezierSurface(size_t n, size_t m, const std::vector<Vector> &controlPoints);
    Mesh Polygonize(size_t verticeCountU, size_t verticeCountV) const;
};

class Revolution
{
private:
    Vector A, B;
    BezierCurve curve;
    
    inline Vector Q(double U) const {return ((1 - U) * A + U * B);}
    Vector SamplePoint(double U, double theta) const;
public:
    Revolution(Vector A, Vector B, const std::vector<Vector> &curve);
    Mesh Polygonize(size_t axisSampleCount, size_t circleSampleCount) const;
};

#endif // BEZIER_H
