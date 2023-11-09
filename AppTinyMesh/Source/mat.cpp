
#include <cassert>
#include <cstdio>
#include <cmath>
#include <algorithm>

#include "mat.h"


double radians( const double deg )
{
    return ((float) M_PI  / 180.f) * deg;
}

double degrees( const double rad )
{
    return (180.f / (float) M_PI) * rad;
}

Transform::Transform (
    const double t00, const double t01, const double t02, const double t03,
    const double t10, const double t11, const double t12, const double t13,
    const double t20, const double t21, const double t22, const double t23,
    const double t30, const double t31, const double t32, const double t33 )
{
    m[0][0]= t00; m[0][1]= t01; m[0][2]= t02; m[0][3]= t03;
    m[1][0]= t10; m[1][1]= t11; m[1][2]= t12; m[1][3]= t13;
    m[2][0]= t20; m[2][1]= t21; m[2][2]= t22; m[2][3]= t23;
    m[3][0]= t30; m[3][1]= t31; m[3][2]= t32; m[3][3]= t33;
}

Transform& Transform::column( const unsigned id, const double t0, const double t1, const double t2, const double t3 )
{
    m[0][id]= t0;
    m[1][id]= t1;
    m[2][id]= t2;
    m[3][id]= t3;
    return *this;
}

vec4 Transform::column( const unsigned id ) const
{
    assert(id < 4);
    return vec4(m[0][id], m[1][id], m[2][id], m[3][id]);
}

vec4 Transform::column( const unsigned id )
{
    assert(id < 4);
    return vec4(m[0][id], m[1][id], m[2][id], m[3][id]);
}


Transform& Transform::row( const unsigned id, const double t0, const double t1, const double t2, const double t3 )
{
    m[id][0]= t0;
    m[id][1]= t1;
    m[id][2]= t2;
    m[id][3]= t3;
    return *this;
}

vec4 Transform::row( const unsigned id ) const
{
    assert(id < 4);
    return vec4(m[id][0], m[id][1], m[id][2], m[id][3]);
}

vec4 Transform::row( const unsigned id )
{
    assert(id < 4);
    return vec4(m[id][0], m[id][1], m[id][2], m[id][3]);
}


Transform& Transform::column_major( const double matrix[16] ) 
{
    for(int i= 0; i < 4; i++)
        column(i, matrix[4*i], matrix[4*i+1], matrix[4*i+2], matrix[4*i+3]);
    return *this;
}

Transform& Transform::row_major( const double matrix[16] )
{
    for(int i= 0; i < 4; i++)
        row(i, matrix[4*i], matrix[4*i+1], matrix[4*i+2], matrix[4*i+3]);
    return *this;
}

Transform::Transform( const Vector& x, const Vector& y, const Vector& z, const Vector& w )
{
    m[0][0] = x[0];	m[0][1] = y[0];	m[0][2] = z[0];	m[0][3] = w[0];
    m[1][0] = x[1];	m[1][1] = y[1];	m[1][2] = z[1];	m[1][3] = w[1];
    m[2][0] = x[2];	m[2][1] = y[2];	m[2][2] = z[2];	m[2][3] = w[2];
    m[3][0] = 0;	m[3][1] = 0;	m[3][2] = 0;	m[3][3] = 1;
}

Transform::Transform( const vec4& x, const vec4& y, const vec4& z, const vec4& w )
{
    m[0][0] = x.x;	m[0][1] = y.x;	m[0][2] = z.x;	m[0][3] = w.x;
    m[1][0] = x.y;	m[1][1] = y.y;	m[1][2] = z.y;	m[1][3] = w.y;
    m[2][0] = x.z;	m[2][1] = y.z;	m[2][2] = z.z;	m[2][3] = w.z;
    m[3][0] = x.w;  m[3][1] = y.w;	m[3][2] = z.w;	m[3][3] = w.w;
}

Vector Transform::operator[] ( const unsigned c ) const
{
    assert(c < 4);
    return Vector(m[0][c], m[1][c], m[2][c]);
}


//! renvoie le point transforme.
Point Transform::operator() ( const Point& p ) const
{
    double x= p.x;
    double y= p.y;
    double z= p.z;

    double xt= m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3];        // dot(vec4(m[0]), vec4(p, 1))
    double yt= m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3];        // dot(vec4(m[1]), vec4(p, 1))
    double zt= m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3];        // dot(vec4(m[2]), vec4(p, 1))
    double wt= m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3];        // dot(vec4(m[3]), vec4(p, 1))

    assert(wt != 0);
    double w= 1.f / wt;
    if(wt == 1.f)
        return Point(xt, yt, zt);
    else
        return Point(xt*w, yt*w, zt*w);
}

//! renvoie le vecteur transforme.
Vector Transform::operator() ( const Vector& v ) const
{
    double x= v[0];
    double y= v[1];
    double z= v[2];

    double xt= m[0][0] * x + m[0][1] * y + m[0][2] * z;                  // dot(vec4(m[0]), vec4(v, 0))
    double yt= m[1][0] * x + m[1][1] * y + m[1][2] * z;                  // dot(vec4(m[1]), vec4(v, 0))
    double zt= m[2][0] * x + m[2][1] * y + m[2][2] * z;                  // dot(vec4(m[2]), vec4(v, 0))
    // dot(vec4(m[3]), vec4(v, 0)) == dot(vec4(0, 0, 0, 1), vec4(v, 0)) == 0 par definition

    return Vector(xt, yt, zt);
}

//! renvoie le point/vecteur homogene transforme.
vec4 Transform::operator() ( const vec4& v ) const
{
    double x= v.x;
    double y= v.y;
    double z= v.z;
    double w= v.w;

    double xt= m[0][0] * x + m[0][1] * y + m[0][2] * z + m[0][3] * w;    // dot(vec4(m[0]), v)
    double yt= m[1][0] * x + m[1][1] * y + m[1][2] * z + m[1][3] * w;    // dot(vec4(m[1]), v)
    double zt= m[2][0] * x + m[2][1] * y + m[2][2] * z + m[2][3] * w;    // dot(vec4(m[2]), v)
    double wt= m[3][0] * x + m[3][1] * y + m[3][2] * z + m[3][3] * w;    // dot(vec4(m[3]), v)

    return vec4(xt, yt, zt, wt);
}

//! renvoie la transposee de la matrice.
Transform Transform::transpose( ) const
{
    return Transform(
        m[0][0], m[1][0], m[2][0], m[3][0],
        m[0][1], m[1][1], m[2][1], m[3][1],
        m[0][2], m[1][2], m[2][2], m[3][2],
        m[0][3], m[1][3], m[2][3], m[3][3]);
}


Transform Transform::operator() ( const Transform& b ) const
{
    return compose_transform(*this, b);
}

//! renvoie la transformation a appliquer aux normales d'un objet transforme par la matrice m.
Transform Transform::normal( ) const
{
    return inverse().transpose();
}


Transform Identity( )
{
    return Transform();
}

Transform Transpose( const Transform& m )
{
    return m.transpose();
}

Transform Inverse( const Transform& m )
{
    return m.inverse();
}

Transform Normal( const Transform& m )
{
    return m.normal();
}

Transform Scale( const double x, const double y, const double z )
{
    return Transform(
        x, 0, 0, 0,
        0, y, 0, 0,
        0, 0, z, 0,
        0, 0, 0, 1);
}

Transform Translation( const Vector& v )
{
    return Transform(
        1, 0, 0, v[0],
        0, 1, 0, v[1],
        0, 0, 1, v[2],
        0, 0, 0, 1);
}

Transform Translation( const double x, const double y, const double z )
{
    return Translation( Vector(x, y, z) );
}

Transform RotationX( const double a )
{
    double sin_t= sinf(radians(a));
    double cos_t= cosf(radians(a));

    return Transform(
        1,     0,      0, 0,
        0, cos_t, -sin_t, 0,
        0, sin_t,  cos_t, 0,
        0,     0,      0, 1 );
}

Transform RotationY( const double a )
{
    double sin_t= sinf(radians(a));
    double cos_t= cosf(radians(a));

    return Transform(
         cos_t, 0, sin_t, 0,
             0, 1,     0, 0,
        -sin_t, 0, cos_t, 0,
             0, 0,     0, 1 );
}

Transform RotationZ( const double a )
{
    double sin_t= sinf(radians(a));
    double cos_t= cosf(radians(a));

    return Transform(
        cos_t, -sin_t, 0, 0,
        sin_t,  cos_t, 0, 0,
            0,      0, 1, 0,
            0,      0, 0, 1 );
}

Transform Rotation( const Vector& axis, const double angle )
{
    Vector a = axis;
    Normalize(a);
    double s= sinf(radians(angle));
    double c= cosf(radians(angle));

    return Transform(
        a[0] * a[0] + (1 - a[0] * a[0] ) * c,
        a[0] * a[1] * (1 - c  ) - a[2] * s,
        a[0] * a[2] * (1 - c  ) + a[1] * s,
        0,
        
        a[0] * a[1] * (1 - c  ) + a[2] * s,
        a[1] * a[1] + (1 - a[1] * a[1] ) * c,
        a[1] * a[2] * (1 - c  ) - a[0] * s,
        0,
        
        a[0] * a[2] * (1 - c  ) - a[1] * s,
        a[1] * a[2] * (1 - c  ) + a[0] * s,
        a[2] * a[2] + (1 - a[2] * a[2] ) * c,
        0,
        
        0, 0, 0, 1);
}


Transform Rotation( const Vector& u, const Vector& v )
{
    Vector a = u;
    Normalize(a);
    Vector b= v;
    Normalize(b);
    
    Vector w= a / b;      // rotation autour de w, un vecteur perpendiculaire a u et v
    double s= Norm(w); // sin theta
    double c= dotProd(a, b); // cos theta
    
    // si u et v sont colineaires, pas d'axe de rotation, renvoyer +1 ou -1
    if(s < float(0.00001))
        //! \todo ajuster epsilon, ou trouver une autre formulation non degeneree...
        return Scale(std::copysign(c, 1));
    
    // normalise l'axe de rotation
    w= w / s;
    
    // meme matrice de rotation qu'au dessus , cf Rotation(axis, angle), l'axe est le vecteur w, s et c sont le sinus et le cosinus de l'angle
    return Transform(
        w[0] * w[0] + (1 - w[0] * w[0] ) * c,
        w[0] * w[1] * (1 - c  ) - w[2] * s,
        w[0] * w[2] * (1 - c  ) + w[1] * s,
        0,
        
        w[0] * w[1] * (1 - c  ) + w[2] * s,
        w[1] * w[1] + (1 - w[1] * w[1] ) * c,
        w[1] * w[2] * (1 - c  ) - w[0] * s,
        0,
        
        w[0] * w[2] * (1 - c  ) - w[1] * s,
        w[1] * w[2] * (1 - c  ) + w[0] * s,
        w[2] * w[2] + (1 - w[2] * w[2] ) * c,
        0,
        
        0, 0, 0, 1);
}


Transform Perspective( const double fov, const double aspect, const double znear, const double zfar )
{
    // perspective, openGL version
    double itan= 1 / tanf(radians(fov) * 0.5f);
    double id= 1 / (znear - zfar);

    return Transform(
        itan/aspect,    0,               0,                 0,
                  0, itan,               0,                 0,
                  0,    0, (zfar+znear)*id, 2.f*zfar*znear*id,
                  0,    0,              -1,                 0);
}


Transform Ortho( const double left, const double right, const double bottom, const double top, const double znear, const double zfar )
{
    double tx= - (right + left) / (right - left);
    double ty= - (top + bottom) / (top - bottom);
    double tz= - (zfar + znear) / (zfar - znear);
   
    return Transform(
        2.f / (right - left),                    0,                     0, tx,
                           0, 2.f / (top - bottom),                     0, ty,
        0,                                       0, -2.f / (zfar - znear), tz,
        0,                                       0,                     0, 1);
}


Transform Viewport( const double width, const double height )
{
    double w= width / 2.f;
    double h= height / 2.f;

    return Transform(
        w, 0,   0,   w,
        0, h,   0,   h,
        0, 0, .5f, .5f,
        0, 0,   0,   1);
}

Transform Lookat( const Point& from, const Point& to, const Vector& up )
{
    Vector dir = Vector(to.x - from.x, to.y - from.y, to.z - from.z);
    Normalize(dir);
    Vector nup = up;
    Normalize(nup);
    Vector right = dir / nup;
    Normalize(right);
    Vector newUp = right / dir;
    Normalize(newUp);

    Transform m(
        right[0], newUp[0], -dir[0], from.x,
        right[1], newUp[1], -dir[1], from.y,
        right[2], newUp[2], -dir[2], from.z,
        0,       0,        0,     1);

    return m.inverse();
}

Transform compose_transform( const Transform& a, const Transform& b )
{
    Transform m;
    for(int i = 0; i < 4; i++)
        for(int j = 0; j < 4; j++)
            m.m[i][j]= a.m[i][0] * b.m[0][j] + a.m[i][1] * b.m[1][j] + a.m[i][2] * b.m[2][j] + a.m[i][3] * b.m[3][j];

    return m;
}

Transform operator* ( const Transform& a, const Transform& b )
{
    return compose_transform(a, b);
}

Transform Transform::inverse( ) const
{
    Transform minv= *this;

    int indxc[4], indxr[4];
    int ipiv[4] = { 0, 0, 0, 0 };

    for (int i = 0; i < 4; i++) {
        int irow = -1, icol = -1;
        double big = 0.f;

        // Choose pivot
        for (int j = 0; j < 4; j++) {
            if (ipiv[j] != 1) {
                for (int k = 0; k < 4; k++) {
                    if (ipiv[k] == 0) {
                        if (fabsf(minv.m[j][k]) >= big) {
                            big = std::abs(minv.m[j][k]);
                            irow = j;
                            icol = k;
                        }
                    }
                    else if (ipiv[k] > 1)
                        printf("singular matrix in Transform::inverse()\n");
                }
            }
        }

        assert(irow >= 0 && irow < 4);
        assert(icol >= 0 && icol < 4);

        ++ipiv[icol];
        // Swap rows _irow_ and _icol_ for pivot
        if (irow != icol) {
            for (int k = 0; k < 4; ++k)
                std::swap(minv.m[irow][k], minv.m[icol][k]);
        }

        indxr[i] = irow;
        indxc[i] = icol;
        if (minv.m[icol][icol] == 0.)
            printf("singular matrix in Transform::inverse()\n");

        // Set $m[icol][icol]$ to one by scaling row _icol_ appropriately
        double pivinv = 1.f / minv.m[icol][icol];
        minv.m[icol][icol] = 1.f;
        for (int j = 0; j < 4; j++)
            minv.m[icol][j] *= pivinv;

        // Subtract this row from others to zero out their columns
        for (int j = 0; j < 4; j++) {
            if (j != icol) {
                double save = minv.m[j][icol];
                minv.m[j][icol] = 0;
                for (int k = 0; k < 4; k++)
                    minv.m[j][k] -= minv.m[icol][k]*save;
            }
        }
    }

    // Swap columns to reflect permutation
    for (int j = 3; j >= 0; j--) {
        if (indxr[j] != indxc[j]) {
            for (int k = 0; k < 4; k++)
                std::swap(minv.m[k][indxr[j]], minv.m[k][indxc[j]]);
        }
    }

    return minv;
}
