/*
 * Surface.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Surface.h"


/**
 * Surface constructor sets default values
 */
Surface::Surface() {
	_boundary_type = VACUUM;		/* Default boundary type */
	_left_region = NULL;
	_right_region = NULL;
}


/**
 * Surface destructor
 */
Surface::~Surface() { }



/**
 * Returns this Surface's boundary type (REFLECTIVE, VACUUM, or INTERFACE)
 * @return Surface boundary type
 */
boundaryType Surface::getBoundaryType() const {
    return _boundary_type;
}


/**
 * Returns a pointer to the Region which borders this Surface on the left
 * @return a pointer to the Region
 */
Region* Surface::getLeftRegion() const {
    return _left_region;
}


/**
 * Returns a pointer to the Region which borders this Surface on the right
 * @return a pointer to the Region
 */
Region* Surface::getRightRegion() const {
    return _right_region;
}


/**
 * Sets the boundary type for this Surface (REFLECTIVE, VACUUM, or INTERFACE)
 * @param type the boundary type
 */
void Surface::setBoundaryType(boundaryType type) {
    _boundary_type = type;
}


/**
 * Sets the Region which borders this Surface to the left
 * @param left_region pointer to a Region class object
 */
void Surface::setLeftRegion(Region* left_region) {
    _left_region = left_region;
}


/**
 * Sets the Region which borders this Surface to the right
 * @param right_region pointer to a Region class object
 */
void Surface::setRightRegion(Region* right_region) {
    _right_region = right_region;
}


/******************************************************************************
 *****************************   XPlane   *************************************
 *****************************************************************************/

XPlane::XPlane() {
	_x = 0.0;
};

XPlane::~XPlane() { };

float XPlane::getX() {
	return _x;
}

void XPlane::setX(float x) {
	_x = x;
}


/**
 * Adds a neutron to this XPlane
 * @param neutron a pointer to a neutron struct
 */
void XPlane::addNeutron(neutron* neutron) {

	/* Must compute intersection point with XPlane */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;

	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float t = (_x - x) / vx;

	/* Update neutron's position to intersection point */
	neutron->_x = _x;
	neutron->_y = vy * t + y;
	neutron->_z = vz * t + z;

	_neutrons.push_back(neutron);
}


float XPlane::computeDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
	float x = neutron->_x;
	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));

	if ((x < _x && vx > 0.0) || (x > _x && vx < 0.0)) {

		float y = neutron->_y;
		float z = neutron->_z;

		float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
		float vz = neutron->_mu;

		float t = (_x - x) / vx;
		float y_inters = vy * t + y;
		float z_inters = vz * t + z;
		dist = sqrt((_x - x) * (_x - x) +
					(y_inters - y) * (y_inters - y) +
					(z_inters - z) * (z_inters - z));
	}

	return dist;
}


/**
 * Checks whether a certain x values is on the XPlane
 * (within numerical error)
 * @param x the value to check
 * @return if on the XPlane (true), otherwise false
 */
bool XPlane::onSurface(neutron* neutron) {
	if (fabs(_x - neutron->_x) < 1E-6)
		return true;

	return false;
}


/**
 * Checks whether a certain x values is on the XPlane
 * (within numerical error)
 * @param x the value to check
 * @return if on the XPlane (true), otherwise false
 */
bool XPlane::onSurface(float x, float y, float z) {
	if (fabs(_x - x) < 1E-6)
		return true;

	return false;
}


/**
 * Moves neutrons to either the left or right bordering Region
 * based on the neutron's trajectory. If this Surface has VACUUM
 * boundary conditions, it instead kills the neutron
 */
void XPlane::moveNeutrons() {

	log_printf(DEBUG, "Inside XPlane moveNeutrons method with x = %f", _x);

	neutron* curr;
	float phi;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);
		phi = curr->_phi;

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*TINY_MOVE;
			curr->_y += vy*TINY_MOVE;
			curr->_z += vz*TINY_MOVE;

			/* Figure out which region to put neutron in */
			if (phi > 0.0 && phi <= M_PI)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}



/******************************************************************************
 *****************************   YPlane   *************************************
 *****************************************************************************/

YPlane::YPlane() {
	_y = 0.0;
};

YPlane::~YPlane() { };


float YPlane::getY() {
	return _y;
}

void YPlane::setY(float y) {
	_y = y;
}


/**
 * Adds a neutron to this YPlane
 * @param neutron a pointer to a neutron struct
 */
void YPlane::addNeutron(neutron* neutron) {

	/* Must compute intersection point with YPlane */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;

	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float t = (_y - y) / vy;

	/* Update neutron's position to intersection point */
	neutron->_y = _y;
	neutron->_x = vx * t + x;
	neutron->_z = vz * t + z;

	_neutrons.push_back(neutron);
}


float YPlane::computeDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
	float y = neutron->_y;
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));

	if ((y < _y && vy > 0.0) || (y > _y && vy < 0.0)) {

		float x = neutron->_x;
		float z = neutron->_z;

		float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
		float vz = neutron->_mu;

		float t = (_y - y) / vy;
		float x_inters = vx * t + x;
		float z_inters = vz * t + z;
		dist = sqrt((x_inters - x) * (x_inters - x) +
					(_y - y) * (_y - y) +
					(z_inters - z) * (z_inters - z));
	}

	return dist;
}


/**
 * Checks whether a certain y values is on the YPlane
 * (within numerical error)
 * @param x the value to check
 * @return if on the YPlane (true), otherwise false
 */
bool YPlane::onSurface(neutron* neutron) {
	if (fabs(_y - neutron->_y) < 1E-6)
		return true;

	return false;
}

/**
 * Checks whether a certain y values is on the YPlane
 * (within numerical error)
 * @param x the value to check
 * @return if on the YPlane (true), otherwise false
 */
bool YPlane::onSurface(float x, float y, float z) {
	if (fabs(_y - y) < 1E-6)
		return true;

	return false;
}

/**
 * Moves neutrons to either the left or right bordering Region
 * based on the neutron's trajectory. If this Surface has VACUUM
 * boundary conditions, it instead kills the neutron
 */
void YPlane::moveNeutrons() {

	log_printf(DEBUG, "Inside YPlane moveNeutrons method with y = %f", _y);

	neutron* curr;
	float phi;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);
		phi = curr->_phi;

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*TINY_MOVE;
			curr->_y += vy*TINY_MOVE;
			curr->_z += vz*TINY_MOVE;

			/* Figure out which region to put neutron in */
			if (phi > PI_OVER_TWO && phi <= THREE_PI_OVER_TWO)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}


/******************************************************************************
 *****************************   ZPlane   *************************************
 *****************************************************************************/

ZPlane::ZPlane() {
	_z = 0.0;
};

ZPlane::~ZPlane() { };


float ZPlane::getZ() {
	return _z;
}

void ZPlane::setZ(float z) {
	_z = z;
}


/**
 * Adds a neutron to this ZPlane
 * @param neutron a pointer to a neutron struct
 */
void ZPlane::addNeutron(neutron* neutron) {

	/* Must compute intersection point with YPlane */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;

	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float t = (_z - z) / vz;

	/* Update neutron's position to intersection point */
	neutron->_z = _z;
	neutron->_y = vy * t + y;
	neutron->_x = vx * t + x;

	_neutrons.push_back(neutron);
}


float ZPlane::computeDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
	float z = neutron->_z;
	float vz = neutron->_mu;

	if ((z < _z && vz > 0.0) || (z > _z && vz < 0.0)) {

		float x = neutron->_x;
		float y = neutron->_y;

		float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
		float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));

		float t = (_z - z) / vz;
		float x_inters = vx * t + x;
		float y_inters = vy * t + y;
		dist = sqrt((x_inters - x) * (x_inters - x) +
					(y_inters - y) * (y_inters - y) +
					(_z - z) * (_z - z));
	}

	return dist;
}


/**
 * Checks whether a certain z values is on the ZPlane
 * (within numerical error)
 * @param z the value to check
 * @return if on the ZPlane (true), otherwise false
 */
bool ZPlane::onSurface(neutron* neutron) {
	if (fabs(_z - neutron->_z) < 1E-6)
		return true;

	return false;
}


/**
 * Checks whether a certain z values is on the ZPlane
 * (within numerical error)
 * @param z the value to check
 * @return if on the ZPlane (true), otherwise false
 */
bool ZPlane::onSurface(float x, float y, float z) {
	if (fabs(_z - z) < 1E-6)
		return true;

	return false;
}


/**
 * Moves neutrons to either the left or right bordering Region
 * based on the neutron's trajectory. If this Surface has VACUUM
 * boundary conditions, it instead kills the neutron
 */
void ZPlane::moveNeutrons() {

	log_printf(DEBUG, "Inside ZPlane moveNeutrons method with z = %f", _z);

	neutron* curr;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*TINY_MOVE;
			curr->_y += vy*TINY_MOVE;
			curr->_z += vz*TINY_MOVE;

			/* Figure out which region to put neutron in */
			if (curr->_mu < 0.0)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}


/******************************************************************************
 *****************************   Circle   *************************************
 *****************************************************************************/

Circle::Circle() {
	_x0 = 0;
	_y0 = 0;
	_z0 = 0;
	_r = 0;
	_r_squared = 0;
};


Circle::~Circle() { };


float Circle::getX0() {
	return _x0;
}


float Circle::getY0() {
	return _y0;
}


float Circle::getZ0() {
	return _z0;
}


float Circle::getRadius() {
	return _r;
}


void Circle::setX0(float x0) {
	_x0 = x0;
}


void Circle::setY0(float y0) {
	_y0 = y0;
}


void Circle::setZ0(float z0) {
	_z0 = z0;
}


void Circle::setRadius(float r) {
	_r = r;
	_r_squared = r*r;
}


/******************************************************************************
 *****************************   XCircle   ************************************
 *****************************************************************************/

XCircle::XCircle() { };

XCircle::~XCircle() { };


void XCircle::addNeutron(neutron* neutron) {

	/* Compute intersetion point of neutron with XCircle */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;

	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float t = (_x0 - x) / vx;
	float y_inters = vy * t + y;
	float z_inters = vz * t + z;

	/* Update neutron's position to the intersection point */
	neutron->_x = _x0;
	neutron->_y = y_inters;
	neutron->_z = z_inters;

	_neutrons.push_back(neutron);
}

float XCircle::computeDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
	float x = neutron->_x;
	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));

	if ((x < _x0 && vx > 0.0) || (x > _x0 && vx < 0.0)) {

		float y = neutron->_y;
		float z = neutron->_z;

		float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
		float vz = neutron->_mu;

		float t = (_x0 - x) / vx;
		float y_inters = vy * t + y;
		float z_inters = vz * t + z;
		float r_squared = (y_inters-_y0)*(y_inters-_y0) +
							(z_inters-_z0)*(z_inters-_z0);

		/* Check if the intersection point is inside the circle */
		if (r_squared <= _r_squared) {
			dist = sqrt((_x0 - x) * (_x0 - x) +
					(y_inters - y) * (y_inters - y) +
					(z_inters - z) * (z_inters - z));
		}
	}

	return dist;
}


bool XCircle::onSurface(neutron* neutron) {

	if (fabs(_x0 - neutron->_x) > 1E-8)
		return false;

	float r_squared = (neutron->_y - _y0) * (neutron->_y - _y0) +
						(neutron->_z - _z0) * (neutron->_z - _z0);

	if (r_squared <= _r_squared)
		return true;
	else
		return false;
}


bool XCircle::onSurface(float x, float y, float z) {

	if (fabs(_x0 - x) > 1E-8)
		return false;

	float r_squared = (y - _y0) * (y - _y0) +
						(z - _z0) * (z - _z0);

	if (r_squared <= _r_squared)
		return true;
	else
		return false;
}


void XCircle::moveNeutrons() {
	log_printf(DEBUG, "Inside XCircle moveNeutrons method");

	neutron* curr;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*TINY_MOVE;
			curr->_y += vy*TINY_MOVE;
			curr->_z += vz*TINY_MOVE;

			/* Figure out which region to put neutron in */
			if (vx < 0.0)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}


/******************************************************************************
 *****************************   YCircle   ************************************
 *****************************************************************************/

YCircle::YCircle() { };

YCircle::~YCircle() { };


void YCircle::addNeutron(neutron* neutron) {

	/* Compute intersetion point of neutron with YCircle */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;

	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float t = (_y0 - y) / vy;
	float x_inters = vx * t + x;
	float z_inters = vz * t + z;

	/* Update neutron's position to the intersection point */
	neutron->_y = _y0;
	neutron->_x = x_inters;
	neutron->_z = z_inters;

	_neutrons.push_back(neutron);
}


float YCircle::computeDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
	float y = neutron->_y;
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));

	if ((y < _y0 && vy > 0.0) || (y > _y0 && vy < 0.0)) {

		float x = neutron->_x;
		float z = neutron->_z;

		float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
		float vz = neutron->_mu;

		float t = (_y0 - y) / vy;
		float x_inters = vx * t + x;
		float z_inters = vz * t + z;
		float r_squared = (x_inters-_x0)*(x_inters-_x0) +
							(z_inters-_z0)*(z_inters-_z0);

		/* Check if the intersection point is inside the circle */
		if (r_squared <= _r_squared)
			dist = sqrt((_y0 - y) * (_y0 - y) +
					(x_inters - x) * (x_inters - x) +
					(z_inters - z) * (z_inters - z));
	}

	return dist;
}



bool YCircle::onSurface(neutron* neutron) {

	if (fabs(_y0 - neutron->_y) > 1E-8)
		return false;

	float r_squared = (neutron->_x - _x0) * (neutron->_x - _x0) +
						(neutron->_z - _z0) * (neutron->_z - _z0);

	if (r_squared <= _r_squared)
		return true;
	else
		return false;
}


bool YCircle::onSurface(float x, float y, float z) {

	if (fabs(_y0 - y) > 1E-8)
		return false;

	float r_squared = (x - _x0) * (x - _x0) +
						(z - _z0) * (z - _z0);

	if (r_squared <= _r_squared)
		return true;
	else
		return false;
}


void YCircle::moveNeutrons() {
	log_printf(DEBUG, "Inside YCircle moveNeutrons method for y = %f", _y0);

	neutron* curr;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*TINY_MOVE;
			curr->_y += vy*TINY_MOVE;
			curr->_z += vz*TINY_MOVE;

			if (vy < 0.0)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}


/******************************************************************************
 *****************************   ZCircle   ************************************
 *****************************************************************************/

ZCircle::ZCircle() { };

ZCircle::~ZCircle() { };


void ZCircle::addNeutron(neutron* neutron) {

	/* Compute intersetion point of neutron with ZCircle */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;

	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float t = (_z0 - z) / vz;
	float x_inters = vx * t + x;
	float y_inters = vy * t + y;

	/* Update neutron's position to the intersection point */
	neutron->_z = _z0;
	neutron->_x = x_inters;
	neutron->_y = y_inters;

	_neutrons.push_back(neutron);
}


float ZCircle::computeDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
	float z = neutron->_z;
	float vz = neutron->_mu;

	if ((z < _z0 && vz > 0.0) || (z > _z0 && vz < 0.0)) {

		float x = neutron->_x;
		float y = neutron->_y;

		float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
		float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));

		float t = (_z0 - z) / vz;
		float x_inters = vx * t + x;
		float y_inters = vy * t + y;
		float r_squared = (x_inters-_x0)*(x_inters-_x0) +
							(y_inters-_y0)*(y_inters-_y0);

		/* Check if the intersection point is inside the circle */
		if (r_squared <= _r_squared)
			dist = sqrt((_z0 - z) * (_z0 - z) +
					(x_inters - x) * (x_inters - x) +
					(y_inters - y) * (y_inters - y));
	}

	return dist;
}



bool ZCircle::onSurface(neutron* neutron) {

	if (fabs(_z0 - neutron->_z) > 1E-8)
		return false;

	float r_squared = (neutron->_x - _x0) * (neutron->_x - _x0) +
						(neutron->_y - _y0) * (neutron->_y - _y0);

	if (r_squared <= _r_squared)
		return true;
	else
		return false;
}


bool ZCircle::onSurface(float x, float y, float z) {

	if (fabs(_z0 - z) > 1E-8)
		return false;

	float r_squared = (x - _x0) * (x - _x0) +
						(y - _y0) * (y - _y0);

	if (r_squared <= _r_squared)
		return true;
	else
		return false;
}


void ZCircle::moveNeutrons() {
	log_printf(DEBUG, "Inside ZCircle moveNeutrons method");

	neutron* curr;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*TINY_MOVE;
			curr->_y += vy*TINY_MOVE;
			curr->_z += vz*TINY_MOVE;

			/* Figure out which region to put neutron in */
			if (vz < 0.0)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}



/******************************************************************************
 ******************************   Sphere   ************************************
 *****************************************************************************/

Sphere::Sphere() {
	_x0 = 0;
	_y0 = 0;
	_z0 = 0;
	_r = 0;
	_r_squared = 0;
};


Sphere::~Sphere() { };


float Sphere::getX0() {
	return _x0;
}


float Sphere::getY0() {
	return _y0;
}


float Sphere::getZ0() {
	return _z0;
}


float Sphere::getRadius() {
	return _r;
}


void Sphere::setX0(float x0) {
	_x0 = x0;
}


void Sphere::setY0(float y0) {
	_y0 = y0;
}


void Sphere::setZ0(float z0) {
	_z0 = z0;
}


void Sphere::setRadius(float r) {
	_r = r;
	_r_squared = r*r;
}

void Sphere::addNeutron(neutron* neutron) {

	/* Compute intersetion point of neutron with Sphere */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;

	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float a = vx*vx + vy*vy + vz*vz;
	float b = 2*vx*(x-_x0) + 2*vy*(y-_y0) + 2*vz*(z-_z0);
	float c = _x0*_x0 + _y0*_y0 + _z0*_z0 + x*x + y*y + z*z -
				2*(_x0*x + _y0*y + _z0*z) - _r_squared;
	float discr = b*b - 4*a*c;

	/* There is not an intersection point */
	if (discr < 0) {
		log_printf(NORMAL, "Cannot add neutron to sphere since its "
				"trajectory does not intersect the surface");
		delete neutron;
	}


	/* There is one intersection point */
	else if (discr == 0) {
		float t = -b / (2.0*a);
		float new_x = x + vx*t;
		float new_y = y + vy*t;
		float new_z = z + vz*t;
		float new_r = sqrt((new_x - _x0)*(new_x - _x0) +
						(new_y - _y0)*(new_y-_y0) +
						(new_z - _z0)*(new_z - _z0));
		if (t >= 0 && fabs(new_r - _r) < 1E-6) {
			neutron->_x = new_x;
			neutron->_y = new_y;
			neutron->_z = new_z;
			_neutrons.push_back(neutron);
		}
		else {
			log_printf(NORMAL, "Cannot add neutron to sphere since its "
					"trajectory does not intersect the surface");
			delete neutron;
		}

	}

	/* There are two intersection points */
	else {
		float t1 = (-b + sqrt(discr)) / (2.0*a);
		float t2 = (-b - sqrt(discr)) / (2.0*a);
		float new_x1 = x + vx*t1;
		float new_y1 = y + vy*t1;
		float new_z1 = z + vz*t1;
		float new_x2 = x + vx*t2;
		float new_y2 = y + vy*t2;
		float new_z2 = z + vz*t2;
		float test_dist1 = std::numeric_limits<float>::infinity();
		float test_dist2 = std::numeric_limits<float>::infinity();

		float new_r1 = sqrt((new_x1 - _x0)*(new_x1 - _x0) +
						(new_y1 - _y0)*(new_y1-_y0) +
						(new_z1 - _z0)*(new_z1 - _z0));
		float new_r2 = sqrt((new_x2 - _x0)*(new_x2 - _x0) +
						(new_y2 - _y0)*(new_y2-_y0) +
						(new_z2 - _z0)*(new_z2 - _z0));

		if (t1 >= 0 && fabs(new_r1 - _r) < 1E-4) {
			test_dist1 = sqrt((new_y1 - y)*(new_y1 - y) +
						(new_z1 - z)*(new_z1 - z) +
						(new_x1 - x)*(new_x1 - x));
		}

		if (t2 >= 0 && fabs(new_r2 - _r) < 1E-4) {
			test_dist2 = sqrt((new_y2 - y)*(new_y2 - y) +
						(new_z2- z)*(new_z2 - z) +
						(new_x2 - x)*(new_x2 - x));
		}

		if (test_dist1 < test_dist2) {
			neutron->_x = new_x1;
			neutron->_y = new_y1;
			neutron->_z = new_z1;
			_neutrons.push_back(neutron);
		}
		else {
			neutron->_x = new_x2;
			neutron->_y = new_y2;
			neutron->_z = new_z2;
			_neutrons.push_back(neutron);
		}
	}

	return;
}


float Sphere::computeDistance(neutron* neutron) {

	/* Compute intersetion point of neutron with Sphere */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;

	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();

	/* Formula taken from:
	 *http://www.ccs.neu.edu/home/fell/CSU540/programs/RayTracingFormulas.htm
	 */
	float a = vx*vx + vy*vy + vz*vz;
	float b = 2*vx*(x-_x0) + 2*vy*(y-_y0) + 2*vz*(z-_z0);
	float c = _x0*_x0 + _y0*_y0 + _z0*_z0 + x*x + y*y + z*z -
				2*(_x0*x + _y0*y + _z0*z) - _r_squared;
	float discr = b*b - 4*a*c;

	/* There is not an intersection point */
	if (discr < 0)
		return dist;

	/* There is one intersection point */
	else if (discr == 0) {
		float t = -b / (2.0*a);
		float new_x = x + vx*t;
		float new_y = y + vy*t;
		float new_z = z + vz*t;
		float new_r = sqrt((new_x - _x0)*(new_x - _x0) +
						(new_y - _y0)*(new_y-_y0) +
						(new_z - _z0)*(new_z - _z0));
		if (t >= 0 && fabs(new_r - _r) < 1E-6) {
			float new_z = z + vz*t;
			float new_y = y + vy*t;
			dist = sqrt((new_x - x)*(new_x - x) +
						(new_y - y)*(new_y - y) +
						(new_z - z)*(new_z - z));
		}

		return dist;
	}

	/* There are two intersection points */
	else {
		float t1 = (-b + sqrt(discr)) / (2.0*a);
		float t2 = (-b - sqrt(discr)) / (2.0*a);
		float new_x1 = x + vx*t1;
		float new_y1 = y + vy*t1;
		float new_z1 = z + vz*t1;
		float new_x2 = x + vx*t2;
		float new_y2 = y + vy*t2;
		float new_z2 = z + vz*t2;
		float test_dist1 = std::numeric_limits<float>::infinity();
		float test_dist2 = std::numeric_limits<float>::infinity();

		float new_r1 = sqrt((new_x1 - _x0)*(new_x1 - _x0) +
						(new_y1 - _y0)*(new_y1-_y0) +
						(new_z1 - _z0)*(new_z1 - _z0));
		float new_r2 = sqrt((new_x2 - _x0)*(new_x2 - _x0) +
						(new_y2 - _y0)*(new_y2-_y0) +
						(new_z2 - _z0)*(new_z2 - _z0));

		if (t1 >= 0 && fabs(new_r1 - _r) < 1E-4) {
			test_dist1 = sqrt((new_y1 - y)*(new_y1 - y) +
						(new_z1 - z)*(new_z1 - z) +
						(new_x1 - x)*(new_x1 - x));
		}

		if (t2 >= 0 && fabs(new_r2 - _r) < 1E-4) {
			test_dist2 = sqrt((new_y2 - y)*(new_y2 - y) +
						(new_z2- z)*(new_z2 - z) +
						(new_x2 - x)*(new_x2 - x));
		}

		if (test_dist1 < test_dist2)
			return test_dist1;
		else
			return test_dist2;
	}

	return dist;
}


bool Sphere::onSurface(neutron* neutron) {

	float r_squared = (neutron->_x - _x0) * (neutron->_x - _x0) +
						(neutron->_y - _y0) * (neutron->_y - _y0) +
						(neutron->_z - _z0) * (neutron->_z - _z0);
	float r = sqrt(r_squared);

	if (fabs(r - _r) < 1E-6)
		return true;
	else
		return false;
}


bool Sphere::onSurface(float x, float y, float z) {

	float r_squared = (x - _x0) * (x - _x0) +
						(y - _y0) * (y - _y0) +
						(z - _z0) * (z - _z0);
	float r = sqrt(r_squared);

	if (fabs(r - _r) < 1E-6)
		return true;
	else
		return false;
}

void Sphere::moveNeutrons() {
	log_printf(DEBUG, "Inside Sphere moveNeutrons method");

	neutron* curr;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*1.0*TINY_MOVE;
			curr->_y += vy*1.0*TINY_MOVE;
			curr->_z += vz*1.0*TINY_MOVE;

			float new_r_squared = (curr->_x - _x0)*(curr->_x - _x0) +
									(curr->_y - _y0)*(curr->_y - _y0) +
									(curr->_z - _z0)*(curr->_z - _z0);

			/* Figure out which region to put neutron in */
			if (new_r_squared < _r_squared)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}


/******************************************************************************
 ***************************   OpenCylinder  **********************************
 *****************************************************************************/

OpenCylinder::OpenCylinder() {
	_x0 = 0;
	_y0 = 0;
	_z0 = 0;
	_r = 0;
	_r_squared = 0;
};


OpenCylinder::~OpenCylinder() { };


float OpenCylinder::getX0() {
	return _x0;
}


float OpenCylinder::getY0() {
	return _y0;
}


float OpenCylinder::getZ0() {
	return _z0;
}


float OpenCylinder::getRadius() {
	return _r;
}


void OpenCylinder::setX0(float x0) {
	_x0 = x0;
}


void OpenCylinder::setY0(float y0) {
	_y0 = y0;
}


void OpenCylinder::setZ0(float z0) {
	_z0 = z0;
}


void OpenCylinder::setRadius(float r) {
	_r = r;
	_r_squared = r*r;
}


/******************************************************************************
 **************************   OpenXCylinder   *********************************
 *****************************************************************************/

OpenXCylinder::OpenXCylinder() { };

OpenXCylinder::~OpenXCylinder() { };


float OpenXCylinder::getXLeft() {
	return _x_left;
}


float OpenXCylinder::getXRight() {
	return _x_right;
}


void OpenXCylinder::setXLeft(float x_left) {
	_x_left = x_left;
}


void OpenXCylinder::setXRight(float x_right) {
	_x_right = x_right;
}


void OpenXCylinder::addNeutron(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;
	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	double a = vz*vz + vy*vy;
	double b = 2.0*z*vz - 2.0*_z0*vz + 2.0*y*vy - 2.0*_y0*vy;
	double c = z*z + _z0*_z0 - 2.0*_z0*z + y*y + _y0*_y0 - 2.0*_y0*y - _r_squared;

	double discr = b*b - 4.0*a*c;

	/* There is not an intersection point */
	if (discr < 0) {
		log_printf(WARNING, "Cannot add neutron to sphere since its "
				"trajectory does not intersect the surface");
		delete neutron;
	}

	/* There is one intersection point */
	else if (discr == 0) {
		float t = -b / (2.0*a);
		float new_x = x + vx*t;
		if (t >= 0 && new_x >= _x_left && new_x <= _x_right) {
			neutron->_x = new_x;
			neutron->_y += vy*t;
			neutron->_z += vz*t;
			_neutrons.push_back(neutron);
		}
		else {
			log_printf(WARNING, "Cannot add neutron to sphere since its "
					"trajectory does not intersect the surface");
			delete neutron;
		}
	}

	/* There are two intersection points */
	else {
		double t1 = (-b + sqrt(discr)) / (2.0*a);
		double t2 = (-b - sqrt(discr)) / (2.0*a);
		float new_x1 = x + vx*t1;
		float new_x2 = x + vx*t2;
		float new_y1, new_z1, new_y2, new_z2;
		float test_dist1 = std::numeric_limits<float>::infinity();
		float test_dist2 = std::numeric_limits<float>::infinity();

		if (t1 >= 0 && new_x1 >= _x_left && new_x1 <= _x_right) {
			new_z1 = z + vz*t1;
			new_y1 = y + vy*t1;
			test_dist1 = sqrt((new_y1 - y)*(new_y1 - y) +
						(new_z1 - z)*(new_z1 - z) +
						(new_x1 - x)*(new_x1 - x));
		}

		if (t2 >= 0 && new_x2 >= _x_left && new_x2 <= _x_right) {
			new_z2 = z + vz*t2;
			new_y2 = y + vy*t2;
			test_dist2 = sqrt((new_y2 - y)*(new_y2 - y) +
						(new_z2- z)*(new_z2 - z) +
						(new_x2 - x)*(new_x2 - x));
		}

		if (test_dist1 < test_dist2) {
			neutron->_x += vx*t1;
			neutron->_y += vy*t1;
			neutron->_z += vz*t1;
			_neutrons.push_back(neutron);
		}
		else {
			neutron->_x += vx*t2;
			neutron->_y += vy*t2;
			neutron->_z += vz*t2;
			_neutrons.push_back(neutron);
		}
	}

	return;
}


float OpenXCylinder::computeDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;
	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float a = vz*vz + vy*vy;
	float b = 2.0*z*vz - 2.0*_z0*vz + 2.0*y*vy - 2.0*_y0*vy;
	float c = z*z + _z0*_z0 - 2.0*_z0*z + y*y + _y0*_y0 - 2.0*_y0*y - _r_squared;

	float discr = b*b - 4.0*a*c;

	/* There is not an intersection point */
	if (discr < 0)
		return dist;

	/* There is one intersection point */
	else if (discr == 0) {
		float t = -b / (2.0*a);
		float new_x = x + vx*t;
		if (t >= 0 && new_x >= _x_left && new_x <= _x_right) {
			float new_z = z + vz*t;
			float new_y = y + vy*t;
			dist = sqrt((new_x - x)*(new_x - x) +
						(new_y - y)*(new_y - y) +
						(new_z - z)*(new_z - z));
		}

		return dist;
	}

	/* There are two intersection points */
	else {
		float t1 = (-b + sqrt(discr)) / (2.0*a);
		float t2 = (-b - sqrt(discr)) / (2.0*a);
		float new_x1 = x + vx*t1;
		float new_x2 = x + vx*t2;
		float new_y, new_z;
		float test_dist1 = std::numeric_limits<float>::infinity();
		float test_dist2 = std::numeric_limits<float>::infinity();

		if (t1 >= 0 && new_x1 >= _x_left && new_x1 <= _x_right) {
			new_z = z + vz*t1;
			new_y = y + vy*t1;
			test_dist1 = sqrt((new_y - y)*(new_y - y) +
						(new_z - z)*(new_z - z) +
						(new_x1 - x)*(new_x1 - x));
		}

		if (t2 >= 0 && new_x2 >= _x_left && new_x2 <= _x_right) {
			new_z = z + vz*t2;
			new_y = y + vy*t2;
			test_dist2 = sqrt((new_y - y)*(new_y - y) +
						(new_z- z)*(new_z - z) +
						(new_x2 - x)*(new_x2 - x));
		}

		if (test_dist1 < test_dist2)
			return test_dist1;
		else
			return test_dist2;
	}
}



bool OpenXCylinder::onSurface(neutron* neutron) {

	if (neutron->_x < _x_left || neutron->_x > _x_right)
		return false;

	float r_squared = (neutron->_y - _y0) * (neutron->_y - _y0) +
						(neutron->_z - _z0) * (neutron->_z - _z0);
	float r = sqrt(r_squared);

	if (fabs(r - _r) < 1E-6)
		return true;
	else
		return false;
}


bool OpenXCylinder::onSurface(float x, float y, float z) {

	if (x < _x_left || x > _x_right)
		return false;

	float r_squared = (y - _y0) * (y - _y0) +
						(z - _z0) * (z - _z0);
	float r = sqrt(r_squared);

	if (fabs(r - _r) < 1E-6)
		return true;
	else
		return false;
}


void OpenXCylinder::moveNeutrons() {
	log_printf(DEBUG, "Inside OpenXCylinder moveNeutrons method for y0 = %f", _y0);

	neutron* curr;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*TINY_MOVE;
			curr->_y += vy*TINY_MOVE;
			curr->_z += vz*TINY_MOVE;

			float new_r_squared = (curr->_y - _y0)*(curr->_y - _y0) +
									(curr->_z - _z0)*(curr->_z - _z0);

			log_printf(DEBUG, "new_r_squared = %f, _r_squared = %f", new_r_squared, _r_squared);

			/* Figure out which region to put neutron in */
			if (new_r_squared < _r_squared)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}


/******************************************************************************
 **************************   OpenYCylinder   *********************************
 *****************************************************************************/

OpenYCylinder::OpenYCylinder() { };

OpenYCylinder::~OpenYCylinder() { };


float OpenYCylinder::getYLeft() {
	return _y_left;
}


float OpenYCylinder::getYRight() {
	return _y_right;
}


void OpenYCylinder::setYLeft(float y_left) {
	_y_left = y_left;
}


void OpenYCylinder::setYRight(float y_right) {
	_y_right = y_right;
}

void OpenYCylinder::addNeutron(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;
	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float a = vx*vx + vz*vz;
	float b = 2.0*x*vx - 2.0*_x0*vx + 2.0*z*vz - 2.0*_z0*vz;
	float c = x*x + _x0*_x0 - 2.0*_x0*x + z*z + _z0*_z0 - 2.0*_z0*z - _r_squared;

	float discr = b*b - 4.0*a*c;

	/* There is not an intersection point */
	if (discr < 0) {
		log_printf(WARNING, "Cannot add neutron to sphere since its "
				"trajectory does not intersect the surface");
		delete neutron;
	}

	/* There is one intersection point */
	else if (discr == 0) {
		float t = -b / (2.0*a);
		float new_y = y + vy*t;
		if (t >= 0 && new_y >= _y_left && new_y <= _y_right) {
			neutron->_x += vx*t;
			neutron->_y = new_y;
			neutron->_z += vz*t;
			_neutrons.push_back(neutron);
		}
		else {
			log_printf(WARNING, "Cannot add neutron to sphere since its "
					"trajectory does not intersect the surface");
			delete neutron;
		}
	}

	/* There are two intersection points */
	else {
		float t1 = (-b + sqrt(discr)) / (2.0*a);
		float t2 = (-b - sqrt(discr)) / (2.0*a);
		float new_y1 = y + vy*t1;
		float new_y2 = y + vy*t2;
		float new_x1, new_z1, new_x2, new_z2;
		float test_dist1 = std::numeric_limits<float>::infinity();
		float test_dist2 = std::numeric_limits<float>::infinity();

		if (t1 >= 0 && new_y1 >= _y_left && new_y1 <= _y_right) {
			new_z1 = z + vz*t1;
			new_x1 = x + vx*t1;
			test_dist1 = sqrt((new_x1 - x)*(new_x1 - x) +
						(new_z1 - z)*(new_z1 - z) +
						(new_y1 - y)*(new_y1 - y));
		}

		if (t2 >= 0 && new_y2 >= _y_left && new_y2 <= _y_right) {
			new_z2 = z + vz*t2;
			new_x2 = x + vx*t2;
			test_dist2 = sqrt((new_x2 - x)*(new_x2 - x) +
						(new_z2- z)*(new_z2 - z) +
						(new_y2 - y)*(new_y2 - y));
		}

		if (test_dist1 < test_dist2) {
			neutron->_x += vx*t1;
			neutron->_y += vy*t1;
			neutron->_z += vz*t1;
			_neutrons.push_back(neutron);
		}
		else {
			neutron->_x += vx*t2;
			neutron->_y += vy*t2;
			neutron->_z += vz*t2;
			_neutrons.push_back(neutron);
		}
	}

	return;
}


float OpenYCylinder::computeDistance(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;
	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float a = vx*vx + vz*vz;
	float b = 2.0*x*vx - 2.0*_x0*vx + 2.0*z*vz - 2.0*_z0*vz;
	float c = x*x + _x0*_x0 - 2.0*_x0*x + z*z + _z0*_z0 - 2.0*_z0*z - _r_squared;

	float discr = b*b - 4.0*a*c;

	/* There is not an intersection point */
	if (discr < 0)
		return dist;

	/* There is one intersection point */
	else if (discr == 0) {
		float t = -b / (2.0*a);
		float new_y = y + vy*t;
		if (t >= 0 && new_y >= _y_left && new_y <= _y_right) {
			float new_z = z + vz*t;
			float new_x = x + vx*t;
			dist = sqrt((new_x - x)*(new_x - x) +
						(new_y - y)*(new_y - y) +
						(new_z - z)*(new_z - z));
		}

		return dist;
	}

	/* There are two intersection points */
	else {
		float t1 = (-b + sqrt(discr)) / (2.0*a);
		float t2 = (-b - sqrt(discr)) / (2.0*a);
		float new_y1 = y + vy*t1;
		float new_y2 = y + vy*t2;
		float new_x, new_z;
		float test_dist1 = std::numeric_limits<float>::infinity();
		float test_dist2 = std::numeric_limits<float>::infinity();

		if (t1 >= 0 && new_y1 >= _y_left && new_y1 <= _y_right) {
			new_z = z + vz*t1;
			new_x = x + vx*t1;
			test_dist1 = sqrt((new_x - x)*(new_x - x) +
						(new_z - z)*(new_z - z) +
						(new_y1 - y)*(new_y1 - y));
		}

		if (t2 >= 0 && new_y2 >= _y_left && new_y2 <= _y_right) {
			new_z = z + vz*t2;
			new_x = x + vx*t2;
			test_dist2 = sqrt((new_x - x)*(new_x - x) +
						(new_z- z)*(new_z - z) +
						(new_y2 - y)*(new_y2 - y));
		}

		if (test_dist1 < test_dist2)
			return test_dist1;
		else
			return test_dist2;
	}
}



bool OpenYCylinder::onSurface(neutron* neutron) {

	if (neutron->_y < _y_left || neutron->_y > _y_right)
		return false;

	float r_squared = (neutron->_x - _x0) * (neutron->_x - _x0) +
						(neutron->_z - _z0) * (neutron->_z - _z0);
	float r = sqrt(r_squared);

	if (fabs(r - _r) < 1E-6)
		return true;
	else
		return false;
}


bool OpenYCylinder::onSurface(float x, float y, float z) {

	if (y < _y_left || y > _y_right)
		return false;

	float r_squared = (x - _x0) * (x - _x0) +
						(z - _z0) * (z - _z0);
	float r = sqrt(r_squared);

	if (fabs(r - _r) < 1E-6)
		return true;
	else
		return false;
}


void OpenYCylinder::moveNeutrons() {
	log_printf(DEBUG, "Inside OpenYCylinder moveNeutrons method");

	neutron* curr;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*TINY_MOVE;
			curr->_y += vy*TINY_MOVE;
			curr->_z += vz*TINY_MOVE;

			float new_r_squared = (curr->_x - _x0)*(curr->_x - _x0) +
									(curr->_z - _z0)*(curr->_z - _z0);

			/* Figure out which region to put neutron in */
			if (new_r_squared < _r_squared)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}


/******************************************************************************
 **************************   OpenZCylinder   *********************************
 *****************************************************************************/

OpenZCylinder::OpenZCylinder() { };

OpenZCylinder::~OpenZCylinder() { };


float OpenZCylinder::getZLeft() {
	return _z_left;
}


float OpenZCylinder::getZRight() {
	return _z_right;
}


void OpenZCylinder::setZLeft(float z_left) {
	_z_left = z_left;
}


void OpenZCylinder::setZRight(float z_right) {
	_z_right = z_right;
}

void OpenZCylinder::addNeutron(neutron* neutron) {

	/* Set dist to infinity to begin with */
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;
	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

	float a = vx*vx + vy*vy;
	float b = 2.0*x*vx - 2.0*_x0*vx + 2.0*y*vy - 2.0*_y0*vy;
	float c = x*x + _x0*_x0 - 2.0*_x0*x + y*y + _y0*_y0 - 2.0*_y0*y - _r_squared;

	float discr = b*b - 4.0*a*c;

	/* There is not an intersection point */
	if (discr < 0) {
		log_printf(WARNING, "Cannot add neutron to sphere since its "
				"trajectory does not intersect the surface");
		delete neutron;
	}

	/* There is one intersection point */
	else if (discr == 0) {
		float t = -b / (2.0*a);
		float new_z = z + vz*t;
		if (t >= 0 && new_z >= _z_left && new_z <= _z_right) {
			neutron->_x += vx*t;
			neutron->_y += vy*t;
			neutron->_z = new_z;
			_neutrons.push_back(neutron);
		}
		else {
			log_printf(WARNING, "Cannot add neutron to sphere since its "
					"trajectory does not intersect the surface");
			delete neutron;
		}
	}

	/* There are two intersection points */
	else {
		float t1 = (-b + sqrt(discr)) / (2.0*a);
		float t2 = (-b - sqrt(discr)) / (2.0*a);
		float new_z1 = z + vz*t1;
		float new_z2 = z + vz*t2;
		float new_x1, new_y1, new_x2, new_y2;
		float test_dist1 = std::numeric_limits<float>::infinity();
		float test_dist2 = std::numeric_limits<float>::infinity();

		if (t1 >= 0 && new_z1 >= _z_left && new_z1 <= _z_right) {
			new_y1 = y + vy*t1;
			new_x1 = x + vx*t1;
			test_dist1 = sqrt((new_x1 - x)*(new_x1 - x) +
						(new_y1 - y)*(new_y1 - y) +
						(new_z1 - z)*(new_z1 - z));
		}

		if (t2 >= 0 && new_z2 >= _z_left && new_z2 <= _z_right) {
			new_y2 = y + vy*t2;
			new_x2 = x + vx*t2;
			test_dist2 = sqrt((new_x2 - x)*(new_x2 - x) +
						(new_y2 - y)*(new_y2 - y) +
						(new_z2 - z)*(new_z2 - z));
		}

		if (test_dist1 < test_dist2) {
			neutron->_x += vx*t1;
			neutron->_y += vy*t1;
			neutron->_z += vz*t1;
			_neutrons.push_back(neutron);
		}
		else {
			neutron->_x += vx*t2;
			neutron->_y += vy*t2;
			neutron->_z += vz*t2;
			_neutrons.push_back(neutron);
		}
	}
	return;
}


float OpenZCylinder::computeDistance(neutron* neutron) {

//	log_printf(NORMAL, "Computing distance to z cylinder");

	/* Set dist to infinity to begin with */
	float dist = std::numeric_limits<float>::infinity();
	float x = neutron->_x;
	float y = neutron->_y;
	float z = neutron->_z;
	float vx = cos(neutron->_phi) * sin(acos(neutron->_mu));
	float vy = sin(neutron->_phi) * sin(acos(neutron->_mu));
	float vz = neutron->_mu;

//	log_printf(NORMAL, "x = %f, y = %f, z = %f, vx = %f, vy = %f, vz = %f", x, y, z, vx, vy ,vz);

	float a = vx*vx + vy*vy;
	float b = 2.0*x*vx - 2.0*_x0*vx + 2.0*y*vy - 2.0*_y0*vy;
	float c = x*x + _x0*_x0 - 2.0*_x0*x + y*y + _y0*_y0 - 2.0*_y0*y - _r_squared;

	float discr = b*b - 4.0*a*c;

//	log_printf(NORMAL, "discr = %f", discr);

	/* There is not an intersection point */
	if (discr < 0)
		return dist;

	/* There is one intersection point */
	else if (discr == 0) {
		float t = -b / (2.0*a);
		float new_z = z + vz*t;
		if (t >= 0 && new_z >= _z_left && new_z <= _z_right) {
			float new_y = y + vy*t;
			float new_x = x + vx*t;
			dist = sqrt((new_x - x)*(new_x - x) +
						(new_y - y)*(new_y - y) +
						(new_z - z)*(new_z - z));
		}

		return dist;
	}

	/* There are two intersection points */
	else {
		float t1 = (-b + sqrt(discr)) / (2.0*a);
		float t2 = (-b - sqrt(discr)) / (2.0*a);
		float new_z1 = z + vz*t1;
		float new_z2 = z + vz*t2;
		float new_x, new_y;
		float test_dist1 = std::numeric_limits<float>::infinity();
		float test_dist2 = std::numeric_limits<float>::infinity();

//		log_printf(NORMAL, "new_z1 = %f, new_z2 = %f, t1 = %f, t2 = %f", new_z1, new_z2, t1, t2);

		if (t1 >= 0 && new_z1 >= _z_left && new_z1 <= _z_right) {
			new_y = y + vy*t1;
			new_x = x + vx*t1;
			test_dist1 = sqrt((new_x - x)*(new_x - x) +
						(new_y - y)*(new_y - y) +
						(new_z1 - z)*(new_z1 - z));

//			log_printf(NORMAL, "new_y = %f, new_x = %f, test_dist1 = %f", new_x, new_y, test_dist1);
		}

		if (t2 >= 0 && new_z2 >= _z_left && new_z2 <= _z_right) {
			new_y = y + vy*t2;
			new_x = x + vx*t2;
			test_dist2 = sqrt((new_x - x)*(new_x - x) +
						(new_y- y)*(new_y - y) +
						(new_z2 - z)*(new_z2 - z));
		}

		if (test_dist1 < test_dist2)
			return test_dist1;
		else
			return test_dist2;
	}
}



bool OpenZCylinder::onSurface(neutron* neutron) {

	if (neutron->_z < _z_left || neutron->_z > _z_right)
		return false;

	float r_squared = (neutron->_x - _x0) * (neutron->_x - _x0) +
						(neutron->_y - _y0) * (neutron->_y - _y0);
	float r = sqrt(r_squared);

	if (fabs(r - _r) < 1E-6)
		return true;
	else
		return false;
}


bool OpenZCylinder::onSurface(float x, float y, float z) {

	if (z < _z_left || z > _z_right)
		return false;

	float r_squared = (x - _x0) * (x - _x0) +
						(y - _y0) * (y - _y0);
	float r = sqrt(r_squared);

	if (fabs(r - _r) < 1E-6)
		return true;
	else
		return false;
}


void OpenZCylinder::moveNeutrons() {
	log_printf(DEBUG, "Inside OpenZCylinder moveNeutrons method");

	neutron* curr;

	std::vector<neutron*>::iterator iter;
	for (iter = _neutrons.begin(); iter != _neutrons.end(); ++iter) {

		curr = (*iter);

		log_printf(DEBUG, "Moving neutron %s", neutronToString(curr).c_str());

		/* If the surface is a vacuum, kill neutron */
		if (_boundary_type == VACUUM) {
			iter = _neutrons.erase(iter);
			--iter;
			delete curr;
		}

		/* If the surface is an interface between two regions */
		else {

			float vx = cos(curr->_phi) * sin(acos(curr->_mu));
			float vy = sin(curr->_phi) * sin(acos(curr->_mu));
			float vz = curr->_mu;

			curr->_x += vx*TINY_MOVE;
			curr->_y += vy*TINY_MOVE;
			curr->_z += vz*TINY_MOVE;

			float new_r_squared = (curr->_x - _x0)*(curr->_x - _x0) +
									(curr->_y - _y0)*(curr->_y - _y0);

			/* Figure out which region to put neutron in */
			if (new_r_squared < _r_squared)
				_left_region->addNeutron(curr);
			else
				_right_region->addNeutron(curr);
		}
	}

	/* Clear all neutrons from this surface */
	_neutrons.clear();
}

