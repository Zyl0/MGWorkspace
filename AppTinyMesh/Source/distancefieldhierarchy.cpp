#include "distancefieldhierarchy.h"

#include "math.h"


double G(double a, double b)
{
    double r =  b - a;
    double h = Math::Max(r - abs(a - b), 0) / r;
    return (1./6.) * r * (h*h*h);  
}

double HDFUnion::Value(const Vector &v) const
{
    return Math::Min(getLeftSon()->Value(v), getRightSon()->Value(v));
}

double HDFIntersection::Value(const Vector &v) const
{
    return Math::Max(getLeftSon()->Value(v), getRightSon()->Value(v));
}

double HDFDiff::Value(const Vector &v) const
{
    return Math::Max(getLeftSon()->Value(v), -(getRightSon()->Value(v)));
}

double HDFBlend::Value(const Vector &v) const
{
    double a =getLeftSon()->Value(v), b = getRightSon()->Value(v);
    
    double g = G(a, b);

    return Math::Min(a, b) - g;
}

double HDFSmouthUnion::Value(const Vector &v) const
{
    double a =getLeftSon()->Value(v), b = -(getRightSon()->Value(v));
    
    double g = G(a, b);

    return Math::Max(a, b) + g;
}
