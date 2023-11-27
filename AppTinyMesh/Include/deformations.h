#ifndef DEFORMATIONS_H
#define DEFORMATIONS_H

#include "mathematics.h"
#include "mesh.h"

class Twist
{
private:
    Vector axis;
    float frequency;
    
    Vector Warp(const Vector& point);

public:
    Twist(float _frequency, Vector _axis) : frequency(_frequency), axis(_axis)
    {
        Normalize(axis);
    }

    Mesh WarpMesh(Mesh& mesh);
};


#endif
