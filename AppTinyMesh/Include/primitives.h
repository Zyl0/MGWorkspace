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
        center(Center), radius(Radius)
    {}

    double Value(const Vector&) const override;
};

class SDFCapsule : public AnalyticScalarField
{
private:
    Vector A;
    Vector B;
    double radius;
public:
    SDFCapsule(double Radius = 1.0, Vector A = Vector(0), Vector B = Vector(0)) :
        A(A),
        B(B),
        radius(Radius)
    {}

    double Value(const Vector&) const override;
};

class SDFBox : public AnalyticScalarField
{
private:
    Vector A;
    Vector B;
    Vector C;
    Vector D;
    double radius;
public:
    SDFBox(double Radius = 1.0, Vector A = Vector(0), Vector B = Vector(0)) :
        A(A),
        B(B),
        radius(Radius)
    {
        C = (A + B) / 2;
        D = (B - A) / 2;
    }

    double Value(const Vector&) const override;
};

class SDFTore : public AnalyticScalarField
{
private:
    Vector center;
    double inerRadius;
    double outerRadius;
public:
    SDFTore(double inerRadius = 0.75, double outerRadius = 1.0, Vector center = Vector(0)) :
        center(center),
        inerRadius(inerRadius),
        outerRadius(outerRadius)
    {
    }

    double Value(const Vector&) const override;
};

#endif // PRIMITIVES_H
