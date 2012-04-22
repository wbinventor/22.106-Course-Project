/*
 * Surface.h
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef SURFACE_H_
#define SURFACE_H_

#include <vector>
#include "Region.h"
#include "Neutron.h"

class Region;

#define PI_OVER_TWO 1.57079633
#define THREE_PI_OVER_TWO 4.71238898
#define TWO_PI 6.28318531
#define TINY_MOVE 1E-3

/* Pre-define the Region1D class so the compiler knows it exists */
class Region1D;

/* The types of boundaries */
typedef enum boundaryTypes {
	VACUUM,
	INTERFACE
} boundaryType;


/**
 * The Surface class represents a plane in two dimensions that is perpendicular
 * to the x-axis. The Surface knows the Region1D class objects which border
 * it and it knows whether it has reflective or vacuum boundary conditions, or
 * whether it is an interface between two regions. The Surface contains a vector
 * of neutrons which are on the Surface, and it either adds them to a bordering
 * Region1D or kills them depending on what type of boundary it is
 */
class Surface
{
protected:
	boundaryType _boundary_type;
	Region* _right_region;
	Region* _left_region;
	std::vector<neutron*> _neutrons;
public:
	Surface();
	virtual ~Surface();

    boundaryType getBoundaryType() const;
    Region* getRightRegion() const;
    Region* getLeftRegion() const;

    void setBoundaryType(boundaryType type);
    void setLeftRegion(Region* left_region);
    void setRightRegion(Region* right_region);

    virtual void addNeutron(neutron* neutron) =0;
    virtual float computeDistance(neutron* neutron) =0;
    virtual bool onSurface(neutron* neutron) =0;
    virtual bool onSurface(float x, float y, float z) =0;
    virtual void moveNeutrons() =0;
};


class XPlane: public Surface {
private:
	float _x;
public:
	XPlane();
	virtual ~XPlane();
	float getX();
	void setX(float x);
    void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};


class YPlane: public Surface {
private:
	float _y;
public:
	YPlane();
	virtual ~YPlane();
	float getY();
	void setY(float y);
    void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};


class ZPlane: public Surface {
private:
	float _z;
public:
	ZPlane();
	virtual ~ZPlane();
	float getZ();
	void setZ(float z);
    void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};


class Circle: public Surface {
protected:
	float _r, _r_squared;
	float _x0, _y0, _z0;
public:
	Circle();
	virtual ~Circle();
	float getX0();
	float getY0();
	float getZ0();
	float getRadius();
	void setX0(float x0);
	void setY0(float y0);
	void setZ0(float z0);
	void setRadius(float r);
    virtual void addNeutron(neutron* neutron) =0;
    virtual float computeDistance(neutron* neutron) =0;
    virtual bool onSurface(neutron* neutron) =0;
    virtual bool onSurface(float x, float y, float z) =0;
    virtual void moveNeutrons() =0;
};


class XCircle: public Circle {
public:
	XCircle();
	virtual ~XCircle();
    virtual void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};


class YCircle: public Circle {
public:
	YCircle();
	virtual ~YCircle();
    void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};


class ZCircle: public Circle {
public:
	ZCircle();
	virtual ~ZCircle();
    void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};


class Sphere: public Surface {
private:
	float _r, _r_squared;
	float _x0, _y0, _z0;
public:
	Sphere();
	virtual ~Sphere();
	float getX0();
	float getY0();
	float getZ0();
	float getRadius();
	void setX0(float x0);
	void setY0(float y0);
	void setZ0(float z0);
	void setRadius(float r);
    void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};

class OpenCylinder: public Surface {
protected:
	float _r, _r_squared;
	float _x0, _y0, _z0;
public:
	OpenCylinder();
	virtual ~OpenCylinder();
	float getX0();
	float getY0();
	float getZ0();
	float getRadius();
	void setX0(float x0);
	void setY0(float y0);
	void setZ0(float z0);
	void setRadius(float r);
    virtual void addNeutron(neutron* neutron) =0;
    virtual float computeDistance(neutron* neutron) =0;
    virtual bool onSurface(neutron* neutron) =0;
    virtual bool onSurface(float x, float y, float z) =0;
    virtual void moveNeutrons() =0;
};

class OpenXCylinder: public OpenCylinder {
private:
	float _x_left, _x_right;
public:
	OpenXCylinder();
	virtual ~OpenXCylinder();
	float getXLeft();
	float getXRight();
	void setXLeft(float x_left);
	void setXRight(float x_right);
    void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};

class OpenYCylinder: public OpenCylinder {
private:
	float _y_left, _y_right;
public:
	OpenYCylinder();
	virtual ~OpenYCylinder();
	float getYLeft();
	float getYRight();
	void setYLeft(float y_left);
	void setYRight(float y_right);
    void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};

class OpenZCylinder: public OpenCylinder {
private:
	float _z_left, _z_right;
public:
	OpenZCylinder();
	virtual ~OpenZCylinder();
	float getZLeft();
	float getZRight();
	void setZLeft(float z_left);
	void setZRight(float z_right);
    void addNeutron(neutron* neutron);
    float computeDistance(neutron* neutron);
    bool onSurface(neutron* neutron);
    bool onSurface(float x, float y, float z);
    void moveNeutrons();
};
#endif /* SURFACE_H_ */
