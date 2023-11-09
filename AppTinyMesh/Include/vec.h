
#ifndef _VEC_H
#define _VEC_H


//! \addtogroup math
///@{

//! \file
//! operations sur points et vecteurs

//! declarations anticipees.
struct vec2;
struct vec3;
struct vec4;
struct Point;

//! representation d'un point 3d.
struct Point
{
    //! constructeur par defaut.
    Point( ) : x(0), y(0), z(0) {}
    explicit Point( const float _x, const float _y, const float _z ) : x(_x), y(_y), z(_z) {}

    //! cree un point a partir des coordonnees du vecteur generique (v.x, v.y, v.z).
    Point( const vec2& v, const float z );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    Point( const vec3& v );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    Point( const vec4& v );   // l'implementation se trouve en fin de fichier, la structure vec3 n'est pas encore connue.
    //! cree un point a partir des coordonnes du vecteur (v.x, v.y, v.z).
    
    //! renvoie la ieme composante du point.
    float operator() ( const unsigned int i ) const; // l'implementation se trouve en fin de fichier
    float& operator() ( const unsigned int i ); // l'implementation se trouve en fin de fichier
    
    float x, y, z;
};

//! renvoie le point origine (0, 0, 0)
Point Origin( );

//! renvoie la distance etre 2 points.
float distance( const Point& a, const Point& b );
//! renvoie le carre de la distance etre 2 points.
float distance2( const Point& a, const Point& b );

//! renvoie le milieu du segment ab.
Point center( const Point& a, const Point& b );

//! renvoie la plus petite composante de chaque point. x, y, z= min(a.x, b.x), min(a.y, b.y), min(a.z, b.z).
Point min( const Point& a, const Point& b );
//! renvoie la plus grande composante de chaque point. x, y, z= max(a.x, b.x), max(a.y, b.y), max(a.z, b.z).
Point max( const Point& a, const Point& b );


//! renvoie le "point" a + b.
Point operator+ ( const Point& a, const Point& b );

//! renvoie le "point" a + b.
Point operator- ( const Point& a, const Point& b );

//! renvoie le "point" k*a;
Point operator* ( const float k, const Point& a );
//! renvoie le "point" a*k;
Point operator* ( const Point& a, const float k );
//! renvoie le "point" v/k;
Point operator/ ( const Point& a, const float k );

//! vecteur generique, utilitaire.
struct vec2
{
    //! constructeur par defaut.
    vec2( ) : x(0), y(0) {}
    explicit vec2( const float _x, const float _y ) : x(_x), y(_y) {}
    
    //! renvoie la ieme composante du vecteur.
    float operator() ( const unsigned int i ) const { return (&x)[i]; }
    float& operator() ( const unsigned int i ) { return (&x)[i]; }

    float x, y;
};


//! vecteur generique, utilitaire.
struct vec3
{
    //! constructeur par defaut.
    vec3( ) : x(0), y(0), z(0) {}
    explicit vec3( const float _x, const float _y, const float _z ) : x(_x), y(_y), z(_z) {}
    //! constructeur par defaut.
    vec3( const vec2& a, const float _z ) : x(a.x), y(a.y), z(_z) {}

    //! cree un vecteur generique a partir des coordonnees du point a.
    vec3( const Point& a );    // l'implementation se trouve en fin de fichier.

    //! renvoie la ieme composante du vecteur.
    float operator() ( const unsigned int i ) const { return (&x)[i]; }
    float& operator() ( const unsigned int i ) { return (&x)[i]; }
    
    operator float*() { return &x; }
    
    float x, y, z;
};


//! vecteur generique 4d, ou 3d homogene, utilitaire.
struct vec4
{
    //! constructeur par defaut.
    vec4( ) : x(0), y(0), z(0), w(0) {}
    explicit vec4( const float _x, const float _y, const float _z, const float _w ) : x(_x), y(_y), z(_z), w(_w) {}
    //! constructeur par defaut.
    vec4( const vec2& v, const float _z= 0, const float _w= 0 ) : x(v.x), y(v.y), z(_z), w(_w) {}
    //! constructeur par defaut.
    vec4( const vec3& v, const float _w= 0 ) : x(v.x), y(v.y), z(v.z), w(_w) {}

    //! cree un vecteur generique a partir des coordonnees du point a, (a.x, a.y, a.z, 1).
    vec4( const Point& a );    // l'implementation se trouve en fin de fichier.
    
    //! renvoie la ieme composante du vecteur.
    float operator() ( const unsigned int i ) const { return (&x)[i]; }
    float& operator() ( const unsigned int i ) { return (&x)[i]; }

    operator float*() { return &x; }

    float x, y, z, w;
};


// implementation des constructeurs explicites.
inline Point::Point( const vec2& v, const float z ) : x(v.x), y(v.y), z(z) {}
inline Point::Point( const vec3& v ) : x(v.x), y(v.y), z(v.z) {}
inline Point::Point( const vec4& v ) : x(v.x), y(v.y), z(v.z) {}

inline vec3::vec3( const Point& a ) : x(a.x), y(a.y), z(a.z) {}

inline vec4::vec4( const Point& a ) : x(a.x), y(a.y), z(a.z), w(1.f) {}

//
inline float Point::operator( ) ( const unsigned int i ) const { return (&x)[i]; }

inline float& Point::operator( ) ( const unsigned int i ) { return (&x)[i]; }

//
#include <iostream>

inline std::ostream& operator<<(std::ostream& o, const Point& p)
{
    o<<"p("<<p.x<<","<<p.y<<","<<p.z<<")";
    return o;
}

///@}
#endif
