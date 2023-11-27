#include "bezier.h"

#include <string>
#include <iostream>

#include "mat.h"

//! Cross product.
inline Vector crossProduct(const Vector& u, const Vector& v)
{
    return u / v;
    //return Vector(u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]);
}

inline Vector RotateVectAroundAxis(Vector Vect, Vector axis, float angle) 
{
    return Vect * cos(angle) + crossProduct(axis, Vect) * sin(angle) + axis * dotProd(axis, Vect) * (1 - cos(angle));
}

Vector BezierSurface::SamplePoint(double U, double V) const
{
    Vector p(0.);
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < n; j++)
        {
            Vector t =  cp(i, j);
            double bUV = B(U, i, m - 1) * B(V, j, n - 1);
            p += bUV * t;
        }
    
    return p;
}

Vector BezierSurface::SampleNormal(double U, double V) const
{
    Vector pdU(0.);
    for (size_t i = 0; i < (m - 1); i++)
        for (size_t j = 0; j < n; j++)
            pdU += B(U, i, m - 2) * B(V, j, n - 1) * m * (cp(i + 1, j) - cp(i, j));

    
    Vector pdV(0.);
    for (size_t i = 0; i < m; i++)
        for (size_t j = 0; j < (n - 1); j++)
            pdV += B(U, i, m - 1) * B(V, j, n - 2) * m * (cp(i, j + 1) - cp(i, j));
        
    Vector N = pdU / pdV;
    Normalize(N);
    return N;
}

BezierSurface::BezierSurface(size_t _n, size_t _m, const std::vector<Vector> &_controlPoints) :
    n(_n),
    m(_m),
    controlPoints(_controlPoints)
{
    size_t f = n > m ? n : m;
    reservePascalTriangle(f);
    reserveFactorials(f);
}

Mesh BezierSurface::Polygonize(size_t verticeCountU, size_t verticeCountV) const
{
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<int> indexes;

    for (size_t i = 0; i < verticeCountU; i++)
    {
        double U = ((double)i)/(verticeCountU - 1);
        for (size_t j = 0; j < verticeCountV; j++)
        {
            double V = ((double)j)/(verticeCountV - 1);
            vertices.push_back(SamplePoint(U,V));
            normals.push_back(SampleNormal(U,V));
        }
    }

    for (size_t i = 0; i < verticeCountU - 1; i++)
    {
        for (size_t j = 0; j < verticeCountV - 1; j++)
        {
            indexes.push_back(i       * verticeCountV + j);
            indexes.push_back((i + 1) * verticeCountV + j);
            indexes.push_back((i + 1) * verticeCountV + j + 1);
            indexes.push_back(i       * verticeCountV + j);
            indexes.push_back((i + 1) * verticeCountV + j + 1);
            indexes.push_back(i       * verticeCountV + j + 1);
        }
        
    }

    return Mesh(vertices, normals, indexes, indexes);
}

void reservePascalTriangle(size_t n)
{
    if(n < pascalTriangleHeight) 
        return;
        
    if(pascalTriangleHeight == 0) 
    {
        pascalTriangle.push_back(1.0);
        pascalTriangleHeight = 1;
    }

    for(unsigned int i = pascalTriangleHeight; i <= n; ++i)
    {
        pascalTriangle.push_back(1.0);
        for(unsigned int j = 1; j < i; ++j)
        {
            unsigned int id = j + (i * (i-1)) / 2;
            pascalTriangle.push_back(pascalTriangle[id] + pascalTriangle[id - 1]);
        }
        pascalTriangle.push_back(1);
    }

    pascalTriangleHeight = n;
}

void reserveFactorials(size_t max)
{
    if(factorials.size() > max) 
        return;

    if(factorials.size() == 0)
        factorials.push_back(1.0);

    for (size_t i = factorials.size(); i <= max; i++)
    {
        size_t val = 1;
        for (size_t j = 2; j <= i; j++)
        {
            val *= j;
        }
        factorials.push_back(val);
    }
}

BezierCurve::BezierCurve(const std::vector<Vector> &points) :
    controlPoints(points)
{
}

Vector BezierCurve::SamplePoint(double U) const
{
    Vector p(0.);
    for (size_t i = 0, n = controlPoints.size(); i < n; i++)
    {
        Vector t =  controlPoints[i];
        double bUV = B(U, i, n - 1);
        p += bUV * t;
    }
    
    return p;
}

Vector Revolution::SamplePoint(double U, double theta) const
{
    Vector p = curve.SamplePoint(U);
    
    return Vector();
}

Revolution::Revolution(Vector _A, Vector _B, const std::vector<Vector> &_curve) :
    A(_A),
    B(_B),
    curve(_curve)
{
    reservePascalTriangle(_curve.size());
    reserveFactorials(_curve.size());
}

Mesh Revolution::Polygonize(size_t axisSampleCount, size_t circleSampleCount) const
{
    std::vector<Vector> vertices;
    std::vector<Vector> normals;
    std::vector<int> indexes;

    std::vector<Transform> rotationMatrixes;

    Vector axis = B - A;
    Normalize(axis);
    for (size_t j = 0; j < circleSampleCount; j++)
    {
        double theta = 2.0 * PI * (((double)j) / (circleSampleCount - 1));

        rotationMatrixes.push_back(RotationRad(axis, theta));
    }

    for (size_t i = 0; i < axisSampleCount; i++)
    {
        double U = ((double)i)/(circleSampleCount - 1);
        Vector ro = curve.SamplePoint(U);
        Vector q = Q(U);
        for (size_t j = 0; j < circleSampleCount; j++)
        {
            //Vector p = RotateVectAroundAxis(ro, axis, theta);
            Vector p = rotationMatrixes[j](ro);

            vertices.push_back(p + q);
            Normalize(p);
            normals.push_back(p);
        }
        
    }

    for (size_t i = 0; i < axisSampleCount - 1; i++)
    {
        for (size_t j = 0; j < circleSampleCount; j++)
        {
            Vector p   = vertices[i       * circleSampleCount + j];
            Vector ppj = (j == (circleSampleCount - 1) ?
                         vertices[i       * circleSampleCount + 0] :
                         vertices[i       * circleSampleCount + j + 1]);
            Vector ppi = vertices[(i + 1) * circleSampleCount + j];
            Vector pdi = ppi - p;
            Vector pdj = (j == (circleSampleCount - 1) ? p - ppj : ppj - p);
            Vector n = pdi / pdj;
            Normalize(n);
            normals[i       * circleSampleCount + j] = n * -1. ; //must be inverted
        }

    }

    for (size_t i = 0; i < axisSampleCount - 1; i++)
    {
        for (size_t j = 0; j < circleSampleCount - 1; j++)
        {
            indexes.push_back(i       * circleSampleCount + j);
            indexes.push_back((i + 1) * circleSampleCount + j);
            indexes.push_back((i + 1) * circleSampleCount + j + 1);
            indexes.push_back(i       * circleSampleCount + j);
            indexes.push_back((i + 1) * circleSampleCount + j + 1);
            indexes.push_back(i       * circleSampleCount + j + 1);
        }
        
    }
    return Mesh(vertices, normals, indexes, indexes);
}
