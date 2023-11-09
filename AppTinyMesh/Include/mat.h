
#ifndef _MAT_H
#define _MAT_H

#include "mathematics.h"
#include "vec.h"


//! \addtogroup math manipulations de points, vecteur, matrices, transformations
///@{

//! \file
//! transformation de points et vecteurs

//! conversion en radians.
double radians( const double deg );
//! conversion en degres.
double degrees( const double rad );

//! representation d'une transformation, une matrice 4x4, organisee par ligne / row major.
struct Transform
{
    //! constructeur./home/wellow/Desktop/lyon1/M2 SI3D/SI3FWorkspace/src/gKit/vec.cpp /home/wellow/Desktop/lyon1/M2 SI3D/SI3FWorkspace/src/gKit/vec.h
    Transform (
        const double t00= 1, const double t01= 0, const double t02= 0, const double t03= 0,
        const double t10= 0, const double t11= 1, const double t12= 0, const double t13= 0,
        const double t20= 0, const double t21= 0, const double t22= 1, const double t23= 0,
        const double t30= 0, const double t31= 0, const double t32= 0, const double t33= 1 );
    
    //! constructeur a partir de 4 Vector colonnes, met (0, 0, 0, 1) dans la derniere ligne.
    Transform( const Vector& x, const Vector& y, const Vector& z, const Vector& w );
    //! constructeur a partir de 4 colonnes
    Transform( const vec4& x, const vec4& y, const vec4& z, const vec4& w );
    
    //! initialise une colonne de la matrice a partir de 4 doubles.
    Transform& column( const unsigned id, const double t0, const double t1, const double t2, const double t3 );
    //!renvoie une colonne.
    vec4 column( const unsigned id ) const;
    //!renvoie une colonne.
    vec4 column( const unsigned id );
    
    //! initialise une ligne de la matrice.
    Transform& row( const unsigned id, const double t0, const double t1, const double t2, const double t3 );
    //!renvoie une ligne.
    vec4 row( const unsigned id ) const;
    //!renvoie une ligne.
    vec4 row( const unsigned id );
    
    //! initialise la matrice avec 16 doubles organises par colonne.
    Transform& column_major( const double matrix[16] );
    
    //! initialise la matrice avec 16 doubles organises par ligne.
    Transform& row_major( const double matrix[16] );
    
    //! renvoie le Vector colonne c de la matrice
    Vector operator[] ( const unsigned c ) const;

    //! renvoie le point transforme.
    Point operator() ( const Point& p ) const;
    //! renvoie le vecteur transforme.
    Vector operator() ( const Vector& v ) const;
    //! renvoie le point/vecteur homogene transforme.
    vec4 operator() ( const vec4& v ) const;
    
    //! renvoie la composition de la transformation this et b, t = this * b. permet de transformer un point sans "ambiguite" Point q= a(b(c(p)));
    Transform operator() ( const Transform& b ) const;
    
    //! renvoie la transposee de la matrice.
    Transform transpose( ) const;
    //! renvoie l'inverse de la matrice.
    Transform inverse( ) const;
    //! renvoie la transformation a appliquer aux normales d'un objet transforme par la matrice m.
    Transform normal( ) const;  
    
    //! renvoie l'adresse de la premiere valeur de la matrice.
    const double *data( ) const { return &m[0][0]; }
    
    double m[4][4];
};

//! construit la transformation identite.
Transform Identity( );

//! renvoie la transposee de la matrice.
Transform Transpose( const Transform& m );
//! renvoie l'inverse de la matrice.
Transform Inverse( const Transform& m );
//! renvoie la transformation a appliquer aux normales d'un objet transforme par la matrice m.
Transform Normal( const Transform& m );

//! renvoie la matrice representant une mise a l'echelle / etirement.
Transform Scale( const double x, const double y, const double z );
inline Transform Scale( const double s ) { return Scale(s, s, s); }

//! renvoie la matrice representant une translation par un vecteur.
Transform Translation( const Vector& v );
//! renvoie la matrice representant une translation par un vecteur x y z.
Transform Translation( const double x, const double y, const double z );

//! renvoie la matrice representation une rotation de angle degree autour de l'axe X.
Transform RotationX( const double angle );
//! renvoie la matrice representation une rotation de a degree autour de l'axe Y.
Transform RotationY( const double angle );
//! renvoie la matrice representation une rotation de angle degree autour de l'axe Z.
Transform RotationZ( const double angle );
//! renvoie la matrice representation une rotation de angle degree autour de l'axe axis.
Transform Rotation( const Vector& axis, const double angle );

//! renvoie la matrice de rotation entre u et v.
Transform Rotation( const Vector&u, const Vector& v );

//! renvoie la matrice representant une transformation viewport.
Transform Viewport( const double width, const double height );
//! renvoie la matrice representant une transformation projection perspective.
Transform Perspective( const double fov, const double aspect, const double znear, const double zfar );
//! renvoie la matrice representant une transformation orthographique, passage d'un cube []x[]x[] vers [-1 1]x[-1 1]x[-1 1].
Transform Ortho( const double left, const double right, const double bottom, const double top, const double znear, const double zfar );
//! renvoie la matrice representant le placement et l'orientation d'une camera pour observer le point to.
Transform Lookat( const Point& from, const Point& to, const Vector& up );

//! renvoie la composition des transformations a et b, t= a * b.
Transform compose_transform( const Transform& a, const Transform& b );
//! renvoie la composition des transformations a et b, t = a * b.
Transform operator* ( const Transform& a, const Transform& b );

#include <iostream>

inline std::ostream& operator<<(std::ostream& o, const Transform& t)
{
    o << t.m[0][0] << " " << t.m[0][1] << " " << t.m[0][2] << " " << t.m[0][3] << " " << std::endl;
    o << t.m[1][0] << " " << t.m[1][1] << " " << t.m[1][2] << " " << t.m[1][3] << " " << std::endl;
    o << t.m[2][0] << " " << t.m[2][1] << " " << t.m[2][2] << " " << t.m[2][3] << " " << std::endl;
    o << t.m[3][0] << " " << t.m[3][1] << " " << t.m[3][2] << " " << t.m[3][3] << " " << std::endl;
    return o;
}


///@}
#endif
