#include "primitives.h"
#include "math.h"

double SDFSphere::Value(const Vector &v) const
{
    return Norm(v - center) - radius;
}

double SDFCapsule::Value(const Vector &v) const
{
    double deltaPA = Norm(v - A);
    double deltaPB = Norm(v - B);
    double deltaAB = Norm(B - A);

    Vector vPA = v - A;
    Vector vPB = v - B;

    Vector u = (B - A) / deltaAB;
    double l = dotProd(vPA, u);

    double dist = 0;

    if (l < 0)
    {
        dist = SquaredNorm(vPA);
    }
    else if (l < deltaAB)
    {
        dist = sqrt(SquaredNorm(vPA) - l*l);
    }
    else
    {
        dist = SquaredNorm(vPB);
    }
    

    return dist - radius;
}

double SDFBox::Value(const Vector &v) const
{
    Vector q = Abs(v) - D;
    double dist = Math::Min(Math::Max(Math::Max(q[0], q[1]), q[2]), 0) + Norm(Vector::Max(q, Vector(0)));
    return dist - radius;
}

double SDFTore::Value(const Vector &v) const
{
    Vector v2 = v - center;
    const double x = Norm(Vector(v2[0], 0, v2[2])) - outerRadius;
    return Norm(Vector(x, v2[1], 0)) - inerRadius;
}