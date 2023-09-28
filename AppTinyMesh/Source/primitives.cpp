#include "primitives.h"

double SDFSphere::Value(const Vector &v) const
{
    Vector dV = v - center;
    double dist = SquaredNorm(dV);
    
    return 0.0;
}
