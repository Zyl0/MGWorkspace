#ifndef DISTANCEFIELDHIERARCHY_H
#define DISTANCEFIELDHIERARCHY_H

#include "implicits.h"

class HierarchalDistanceField : public AnalyticScalarField
{
private:
    AnalyticScalarField *leftSon;
    AnalyticScalarField *rightSon;
protected:
    inline AnalyticScalarField *getLeftSon() const {return leftSon;}
    inline AnalyticScalarField *getRightSon() const {return rightSon;}
public:
    HierarchalDistanceField(AnalyticScalarField *leftSon = nullptr, AnalyticScalarField *rightSon = nullptr) :
        leftSon(leftSon),
        rightSon(rightSon)
    {}
};

class HDFUnion : public HierarchalDistanceField
{
private:
    /* data */
public:
    HDFUnion(AnalyticScalarField *leftSon, AnalyticScalarField *rightSon) :
        HierarchalDistanceField(leftSon, rightSon)
    {}

    double Value(const Vector&) const;
};

class HDFIntersection : public HierarchalDistanceField
{
private:
    /* data */
public:
    HDFIntersection(AnalyticScalarField *leftSon, AnalyticScalarField *rightSon) :
        HierarchalDistanceField(leftSon, rightSon)
    {}

    double Value(const Vector&) const;
};

class HDFDiff: public HierarchalDistanceField
{
private:
    /* data */
public:
    HDFDiff(AnalyticScalarField *leftSon, AnalyticScalarField *rightSon) :
        HierarchalDistanceField(leftSon, rightSon)
    {}

    double Value(const Vector&) const;
};

class HDFBlend: public HierarchalDistanceField
{
private:
    /* data */
public:
    HDFBlend(AnalyticScalarField *leftSon, AnalyticScalarField *rightSon) :
        HierarchalDistanceField(leftSon, rightSon)
    {}

    double Value(const Vector&) const;
};

class HDFSmouthUnion: public HierarchalDistanceField
{
private:
    /* data */
public:
    HDFSmouthUnion(AnalyticScalarField *leftSon, AnalyticScalarField *rightSon) :
        HierarchalDistanceField(leftSon, rightSon)
    {}

    double Value(const Vector&) const;
};


#endif // DISTANCEFIELDHIERARCHY_H
