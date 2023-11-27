#include "deformations.h"

#include "mat.h"

Vector Twist::Warp(const Vector& point) {
    float alpha = (2 * PI * point[2]) / frequency;

    Transform rotator = RotationRad(axis, alpha);

    return rotator(point);
}


Mesh Twist::WarpMesh(Mesh &mesh)
{
    Mesh m(mesh);
    for (int i = 0; i < mesh.Vertexes(); ++i) 
    {
        m.setVertex(i, Warp(mesh.Vertex(i)));
    }

    return m;
}
