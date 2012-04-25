/*
 * Region.h
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef REGION_H_
#define REGION_H_

#include <vector>
#include <math.h>
#include <algorithm>
#include <omp.h>
#include "Material.h"
#include "Isotope.h"
#include "Neutron.h"
#include "Surface.h"
#include "Binner.h"

/* Pre-define Surface classes so compiler knows it exists */
class Surface;
class XPlane;
class YPlane;
class ZPlane;
class XCircle;
class YCircle;
class ZCircle;
class Sphere;
class OpenXCylinder;
class OpenYCylinder;
class OpenZCylinder;


typedef enum varianceReduceTypes {
	IMPLICIT_CAPTURE,
	FORCED_COLLISION,
	GEOMETRY_SPLITTING
} varianceReduceType;

	/**
 * The Region class represents a single dimensional region
 * bounded by two planes, or Surface class objects. The region
 * contains a vector of neutrons which live within it and is filled
 * by a set of Isotope class objects, each with a different number
 * density. The Region class contains all of the physics for moving
 * colliding neutrons
 */
class Region
{
protected:
	char* _region_name;

	/* Vector of pointers to the Surfaces which bound this Region */
	std::vector<Surface*> _boundaries;

	/* Region, if any, contained by this Region */
	Region* _interior_region;

	/* Vector of pointers to neutrons contained by this Region */
	std::vector<neutron*> _neutrons;

	/* Binners for tallying */
	std::vector<Binner*> _bin_sets;

	Material* _material;
	float _volume;

	/* Variance reduction */
	bool _use_implicit_capture;
	bool _use_forced_collision;
	float _weight_low;
	float _weight_avg;
public:
	Region();
	virtual ~Region();

    char* getRegionName();
    Material* getMaterial();
    float getVolume();
    int getNumNeutrons();

    void setRegionName(char* region_name);
    void setMaterial(Material* material);
    void setVolume(float volume);
    void addBinner(Binner* bins);
    void setInteriorRegion(Region* region);
    void addNeutron(neutron* neutron);

    void useImplicitCapture(float weight_low, float weight_avg);
    void useForcedCollision(float weight_low, float weight_avg);

    bool playRussianRoulette(neutron* neutron);
    void clearBinners();
    Surface* computeMinSurfDist(neutron* neutron, float* min_dist);
    bool inInteriorRegion(float x, float y, float z);
    void moveNeutrons();

    virtual bool contains(float x, float y, float z) =0;
    virtual bool onBoundary(neutron* neutron) =0;
};

class Parallelepiped : public Region {
private:
	XPlane* _x_left_surface;
	XPlane* _x_right_surface;
	YPlane* _y_left_surface;
	YPlane* _y_right_surface;
	ZPlane* _z_left_surface;
	ZPlane* _z_right_surface;
public:
	Parallelepiped();
	virtual ~Parallelepiped();
	void setXLeftSurface(XPlane* plane);
	void setXRightSurface(XPlane* plane);
	void setYLeftSurface(YPlane* plane);
	void setYRightSurface(YPlane* plane);
	void setZLeftSurface(ZPlane* plane);
	void setZRightSurface(ZPlane* plane);

	bool contains(float x, float y, float z);
    bool onBoundary(neutron* neutron);
};


class XCylinder : public Region {
private:
	XCircle* _left_circle;
	XCircle* _right_circle;
	OpenXCylinder* _cylinder;
public:
	XCylinder();
	virtual ~XCylinder();
	void setXLeftCircle(XCircle* circle);
	void setXRightCircle(XCircle* circle);
	void setOpenXCylinder(OpenXCylinder* cylinder);
	bool contains(float x, float y, float z);
    bool onBoundary(neutron* neutron);
};


class YCylinder : public Region {
private:
	YCircle* _left_circle;
	YCircle* _right_circle;
	OpenYCylinder* _cylinder;
public:
	YCylinder();
	virtual ~YCylinder();
	void setYLeftCircle(YCircle* circle);
	void setYRightCircle(YCircle* circle);
	void setOpenYCylinder(OpenYCylinder* cylinder);
	bool contains(float x, float y, float z);
    bool onBoundary(neutron* neutron);
};


class ZCylinder : public Region {
private:
	ZCircle* _left_circle;
	ZCircle* _right_circle;
	OpenZCylinder* _cylinder;
public:
	ZCylinder();
	virtual ~ZCylinder();
	void setZLeftCircle(ZCircle* circle);
	void setZRightCircle(ZCircle* circle);
	void setOpenZCylinder(OpenZCylinder* cylinder);
	bool contains(float x, float y, float z);
    bool onBoundary(neutron* neutron);
};


class Spheroid : public Region {
private:
	Sphere* _sphere;
public:
	Spheroid();
	virtual ~Spheroid();
	void setSphere(Sphere* sphere);
	bool contains(float x, float y, float z);
	bool onBoundary(neutron* neutron);
};


class SphericalShell : public Region {
private:
	Sphere* _inner;
	Sphere* _outer;
public:
	SphericalShell();
	virtual ~SphericalShell();
	void setInnerSphere(Sphere* sphere);
	void setOuterSphere(Sphere* sphere);
	bool contains(float x, float y, float z);
	bool onBoundary(neutron* neutron);
};


#endif /* REGION_H_ */
