
#include <algorithm>
#include <cmath>

#include "vec.h"


Point Origin( )
{
    return Point(0, 0, 0);
}


float distance( const Point& a, const Point& b )
{
    return std::sqrt(distance(a, b));
}

float distance2( const Point& a, const Point& b )
{
    Point v = a - b;
    return v.x * v.x + v.y * v.y + v.z * v.z;
}

Point center( const Point& a, const Point& b )
{
    return Point((a.x + b.x) / 2, (a.y + b.y) / 2, (a.z + b.z) / 2);
}


Point min( const Point& a, const Point& b )
{ 
    return Point( std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z) ); 
}

Point max( const Point& a, const Point& b ) 
{ 
    return Point( std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z) ); 
}

Point operator* ( const float k, const Point& a )
{
    return Point(k * a.x, k * a.y, k * a.z);
}

Point operator* ( const Point& a, const float k )
{
    return k * a;
}

Point operator/ ( const Point& a, const float k )
{ 
    float kk= 1.f / k; 
    return kk * a; 
}

Point operator+ ( const Point& a, const Point& b )
{
    return Point(a.x + b.x, a.y + b.y, a.z + b.z);
}

Point operator- (const Point &a, const Point &b)
{
    return Point(a.x - b.x, a.y - b.y, a.z - b.z);
}
