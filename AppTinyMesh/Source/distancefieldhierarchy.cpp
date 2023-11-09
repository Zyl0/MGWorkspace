#include "distancefieldhierarchy.h"

#include "math.h"


double G(double a, double b, double r = 0.25)
{
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
    
    double g = G(a, b, radius);

    return Math::Min(a, b) - g;
}

double HDFSmouthUnion::Value(const Vector &v) const
{
    double a =getLeftSon()->Value(v), b = -(getRightSon()->Value(v));
    
    double g = G(a, b, radius);

    return Math::Max(a, b) + g;
}

double HDFTransform::Value(const Vector &v) const
{
    vec4 vh = invTransform(vec4(v[0], v[1], v[2], 1));
    return getLeftSon()->Value(Vector(vh.x/vh.w, vh.y/vh.w, vh.z/vh.w));
}

AnalyticScalarField *HierarchalDistanceField::PopLeftSon()
{
    AnalyticScalarField * oldLS = leftSon;
    leftSon = nullptr;
    return oldLS;
}
