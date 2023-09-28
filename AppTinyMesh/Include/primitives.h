#ifndef PRIMITIVES_H
#define PRIMITIVES_H

#include "implicits.h"

class SDFSphere : public AnalyticScalarField
{
private:
    Vector center;
    double radius;
public:
    SDFSphere(double Radius = 1.0, Vector Center = Vector(0)) :
        radius(Radius),
        center(Center)
    {}

    double Value(const Vector&) const;
};


#endif // PRIMITIVES_H
