/*
 * Region.cpp
 *
 *  Created on: Mar 15, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Region.h"


/**
 * Region constructor sets default values
 */
Region::Region() {

	/* Default region name */
	_region_name = (char*)"";

	/* Default volume */
	_volume = 0.0;

	/* By default region does not use variance reduction */
	_use_implicit_capture = false;
	_use_forced_collision = false;
	_weight_low = 0.0;
	_weight_avg = 0.0;

	_interior_region = NULL;
}


/**
 * Region destructor
 */
Region::~Region() {
	_boundaries.clear();
}


/**
 * Returns the volume of this Region (cm^3)
 * @return the Region's volume
 */
float Region::getVolume() {
	return _volume;
}


/**
 * Returns a pointer to the Material filling this Region
 * @return a pointer to a Material
 */
Material* Region::getMaterial() {
	return _material;
}


/**
 * Return the name, if any, of this region as specified by
 * the user
 * @return a character array representing this region's name
 */
char* Region::getRegionName() {
	return _region_name;
}

/**
 * Returns the number of neutrons inside of this region
 * @return the number of neutrons contained by this region
 */
int Region::getNumNeutrons() {
	return _neutrons.size();
}


/**
 * Sets this Region's name as specified by a user
 * @param region_name the name of this Region
 */
void Region::setRegionName(char* region_name) {
	_region_name = region_name;
}


/**
 * Set the Material filling this Region
 * @param material a pointer to a Material
 */
void Region::setMaterial(Material* material) {
	_material = material;
}


/**
 * Sets this Region's volume (cm^3)
 * @param volume the volume of this Region
 */
void Region::setVolume(float volume) {
	_volume = volume;
}


/**
 * Adds a new Binner to this region for tallying
 * @param bins a pointer to a Binner class object
 */
void Region::addBinner(Binner* bins) {
	_bin_sets.push_back(bins);
}


/**
 * Adds a new region that is contained by this region
 * @param region a Region contained by this Region
 */
void Region::setInteriorRegion(Region* region) {
	_interior_region = region;
}


/**
 * Adds a new neutron to this region
 * @param neutron a pointer to a neutron
 */
void Region::addNeutron(neutron* neutron) {

	if (!contains(neutron->_x, neutron->_y, neutron->_z))
		log_printf(ERROR, "Cannot add a neutron %s to region %s"
				"since it does not contain it",
				neutronToString(neutron).c_str(), _region_name);

	_neutrons.push_back(neutron);
}


/**
 * Sets this Region to use implicit capture variance reduction
 * @param weight_cutoff the weight cutoff below which we use Russian
 * Roulette to decide which neutrons survive
 */
void Region::useImplicitCapture(float weight_low, float weight_avg) {
	_use_implicit_capture = true;
	_weight_low = weight_low;
	_weight_avg = weight_avg;
	return;
}


/**
 * Sets this Region to use forced collision variance reduction
 * @param weight_cutof the weight cutoff below which we use Russian
 * Roulette to decide which neutrons survive
 */
void Region::useForcedCollision(float weight_low, float weight_avg) {
	_use_forced_collision = true;
	_weight_low = weight_low;
	_weight_avg = weight_avg;
	return;
}


/**
 * This method plays Russian Roulette with a neutron by checking if it
 * has a weight below the low weight and if so, it kills the neutron
 * with a probability equal to the neutron's weight divided by weight
 * average
 * @param neutron a pointer to the neutron of interest
 * @return a boolean to kill (true) or not (false)
 */
bool Region::playRussianRoulette(neutron* neutron) {

	bool kill_neutron = false;

	log_printf(DEBUG, "Playing russian roulette in region %s with neutron "
			"with weight = %f against cutoff= %f", _region_name,
			neutron->_weight, _weight_low);

	/* Only play Roulette if the weight of the neutron is below the
	 * cutoff */
	if (neutron->_weight < _weight_low) {
		float test = float(rand()) / RAND_MAX;

		if (test > (neutron->_weight / _weight_avg))
			kill_neutron = true;
	}


	return kill_neutron;
}


/**
 * Clear this region's vector of Binner class object pointers
 */
void Region::clearBinners() {

	/* First iterate over each binner and accumulate tallies */
	std::vector<Binner*>::iterator iter;
	Binner* curr_bin;

	for (iter = _bin_sets.begin(); iter != _bin_sets.end(); ++iter) {
		curr_bin = *iter;
		curr_bin->processTallyAccumulators();
	}

	/* Then remove the binners from this Region */
	_bin_sets.clear();
}


Surface* Region::computeMinSurfDist(neutron* neutron, float min_dist) {

//	log_printf(NORMAL, "Computing min surf dist for region %s. neutron x = %f, "
//			"y = %f, z = %f", _region_name, neutron->_x, neutron->_y, neutron->_z);

	/* Find distance to nearest wall along neutron's trajectory */
	min_dist = std::numeric_limits<float>::infinity();
	float curr_dist = 0.0;;
	Surface* nearest = NULL;

	/* Loop over all bordering surfaces */
	std::vector<Surface*>::iterator iter;

	for (iter = _boundaries.begin(); iter != _boundaries.end(); ++iter) {

		/* Compute distance to this Surface */
		curr_dist = (*iter)->computeDistance(neutron);

//		log_printf(NORMAL, "computed curr dist = %1.8f", curr_dist);

		if (curr_dist < min_dist) {
//			log_printf(NORMAL, "updated min dist");
			min_dist = curr_dist;
			nearest = (*iter);
		}
	}

//	log_printf(NORMAL, "Returning min dist = %f", min_dist);

	return nearest;
}


bool Region::inInteriorRegion(float x, float y, float z) {

	if (_interior_region == NULL)
		return false;
	else if (_interior_region->contains(x, y, z))
		return true;
	else
		return false;
}


void Region::moveNeutrons() {

	log_printf(DEBUG, "Inside moveNeutrons method for region %s with "
			"%d neutrons", _region_name, getNumNeutrons());

	if (!_material->isRescaled())
		log_printf(ERROR, "Region %s is unable to move neutrons since it's "
				"Material %s has not been rescaled",
				_region_name, _material->getMaterialName().c_str());

	int energy_index;
	float sigma_t;
	float sigma_a;
	float path_length;
	float theta;
	float new_x, new_y, new_z;
	neutron* curr;
	Isotope* isotope;
	float min_dist = 0.0;
	Surface* nearest = NULL;
	collisionType collision_type;
	std::vector<neutron*>::iterator iter1;
	std::vector<Binner*>::iterator iter2;

	/**************************************************************************
	 ************************  MONTE CARLO KERNEL  ****************************
	 *************************************************************************/
	/* Iterate over all neutrons that are inside this Region */
	for (iter1 = _neutrons.begin(); iter1 != _neutrons.end(); ++iter1) {

		log_printf(DEBUG, "Region %s contains %d neutrons",
										_region_name, _neutrons.size());

		curr = (*iter1);

		log_printf(DEBUG, "%s", neutronToString(curr).c_str());

		/* Compute total sigma_t from all isotopes */
		energy_index = _material->getEnergyGridIndex(curr->_energy);
		sigma_t = _material->getTotalMacroXS(energy_index);

		/* Compute path length and the new x,y,z coordinates of this neutron */
		path_length = -log(float(rand()) / RAND_MAX) / sigma_t;

		theta = acos(curr->_mu);
		new_x = curr->_x + path_length * cos(curr->_phi) * sin(theta);
		new_y = curr->_y + path_length * sin(curr->_phi) * sin(theta);
		new_z = curr->_z + path_length * curr->_mu;

		log_printf(DEBUG, "sigma_t = %f, path_length = %f, new_x = %f, "
				"new_y = %f, new_z = %f", sigma_t, path_length,
											new_x, new_y, new_z);

		/* The neutron collided within this region */
		if (contains(new_x, new_y, new_z)) {

			/* Figure out which isotope the neutron collided in */
			isotope = _material->sampleIsotope(energy_index);

			/* Figure out the collision type */
			collision_type = isotope->getCollisionType(energy_index);

			log_printf(DEBUG, "Neutron made %d collision type in isotope: %s",
							collision_type, isotope->getIsotopeType().c_str());

			/* Update the neutron time */
			updateNeutronTime(curr, path_length);

			/* Update the neutron position */
			curr->_x = new_x;
			curr->_y = new_y;
			curr->_z = new_z;

			/* Tally neutron for each binner */
			for (iter2 = _bin_sets.begin(); iter2 != _bin_sets.end(); ++iter2)
				(*iter2)->weightedTally(curr, sigma_t, energy_index,
													_material, isotope);


			/******************************************************************
			 **********************  VARIANCE REDUCTION  **********************
			 *****************************************************************/
			/* Implicit capture variance reduction */
			if (_use_implicit_capture) {

				log_printf(DEBUG, "Applying implicit capture to neutron with "
						"weight = %f", curr->_weight);

				/* If we are using implicit capture, force the collision to be
				 * a scattering collision */
				if (_use_implicit_capture) {

					/* Sample a isotope with a scattering cross-section */
					while (isotope->getScatterXS(curr->_energy) == 0)
						isotope = _material->sampleIsotope(energy_index);

					while (collision_type == CAPTURE)
						collision_type =
									isotope->getCollisionType(energy_index);
				}

				/* Get the total absorption cross-section for this Region */
				sigma_a = _material->getAbsorbMacroXS(energy_index);

				/* Reduce weight by sigma_a / sigma_t */
				curr->_weight *= (1.0 - (sigma_a / sigma_t));

				log_printf(DEBUG, "New weight = %f", curr->_weight);
			}

			/* Forced collision variance reduction */
			else if (_use_forced_collision) {

				log_printf(DEBUG, "Applying forced collision to neutron with "
						"weight = %f", curr->_weight);

				/* Find distance to nearest wall along neutron's trajectory */
				nearest = computeMinSurfDist(curr, min_dist);

				/* Compute exponential term in forced collision prob dist */
				float px = exp(-sigma_t * min_dist);

				/* Create a new uncollided neutron that hits next surface
				 * in the region along the collided neutron's trajectory */
				neutron* new_neutron = initializeNewNeutron();
				new_neutron->_energy = curr->_energy;
				new_neutron->_mu = curr->_mu;
				new_neutron->_phi = curr->_phi;
				new_neutron->_thread_num = curr->_thread_num;
				new_neutron->_weight = curr->_weight * px;

				/* Figure out which surface to put uncollided neutron on */
				nearest->addNeutron(new_neutron);

				/* Update collided neutron's properties */
				path_length =  log(1.0 - float(rand()/RAND_MAX) *
											(1.0 - px)) / sigma_t;
				theta = acos(curr->_mu);

				new_neutron->_x = curr->_x + path_length *
												cos(curr->_phi) * sin(theta);
				new_neutron->_y = curr->_y + path_length *
												sin(curr->_phi) * sin(theta);
				new_neutron->_z = curr->_z + path_length * curr->_mu;

				curr->_weight *= (1.0 - px);

				updateNeutronTime(new_neutron, path_length);
			}

			/******************************************************************
			 ***********************  COLLISIONS  *****************************
			 *****************************************************************/
			/* Check if collision type was absorption and if so kill neutron */
			if (collision_type == CAPTURE && !_use_implicit_capture) {

				log_printf(DEBUG, "Capture type collision");

				iter1 = _neutrons.erase(iter1);
				--iter1;
				delete curr;

				continue;
			}

			/* Check if collision type was fission and if so kill neutron */
			else if (collision_type == FISSION && !_use_implicit_capture) {

				log_printf(DEBUG, "Fission type collision");

				iter1 = _neutrons.erase(iter1);
				--iter1;
				delete curr;
				continue;
			}

			/* Check if collision was general scattering */
			else if (collision_type == SCATTER) {

				log_printf(DEBUG, "Normal scatter type collision");

				float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
				curr->_phi = phi;

				/* Isotropic in lab */
				if (isotope->getScatterAngleType() == ISOTROPIC_LAB) {
					float A = float(isotope->getA());
					float energy = curr->_energy * (1.0 -
						(1.0 - isotope->getAlpha())*(float(rand())/RAND_MAX));
					float mu_cm = ((2 * energy / curr->_energy) -
									isotope->getAlpha() - 1.0) / (1.0 -
														isotope->getAlpha());
					curr->_mu = (1.0 + A*mu_cm) / sqrt(A*A + 1.0 + 2*mu_cm*A);

					if (curr->_energy < 4 && isotope->usesThermalScattering())
						curr->_energy =
						isotope->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy = energy;
				}

				/* Isotropic in CM */
				else {
					float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
					int A = isotope->getA();
					float mu_l = (1.0 + A*mu_cm) / (sqrt(A*A + 2*A*mu_cm+1.0));
					curr->_mu = curr->_mu * mu_l +
							sqrt(1.0 - curr->_mu * curr->_mu) *
							sqrt(1.0 - mu_l * mu_l) * sin(phi);

					if (curr->_energy < 4 && isotope->usesThermalScattering())
						curr->_energy =
						isotope->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy *= (A*A + 2*A*mu_cm + 1.0) /
															((A+1.0)*(A+1.0));
				}

				log_printf(DEBUG, "Updated mu= %f, phi = %f, energy = %f",
									curr->_mu, curr->_phi, curr->_energy);
			}

			/* Check if collision type was inelastic scattering */
			else if (collision_type == INELASTIC) {

				log_printf(DEBUG, "Inelastic scatter type collision");

				float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
				curr->_phi = phi;

				/* Isotropic in lab */
				if (isotope->getInelasticAngleType() == ISOTROPIC_LAB)
					curr->_mu = (float(rand()) / RAND_MAX) * 2.0 - 1.0;

				/* Isotropic in CM */
				else {
					float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
					int A = isotope->getA();
					float mu_l = (1.0 + A*mu_cm) / (sqrt(A*A + 2*A*mu_cm+1.0));

					curr->_mu = curr->_mu * mu_l +
							sqrt(1.0 - curr->_mu * curr->_mu) *
							sqrt(1.0 - mu_l * mu_l) * sin(phi);
				}

				/* Update the neutron energy */
				curr->_energy =
						isotope->getInelasticScatterEnergy(curr->_energy);

				log_printf(DEBUG, "Updated mu= %f, phi = %f, energy = %f",
									curr->_mu, curr->_phi, curr->_energy);
			}

			/* Check if collision type was elastic scattering */
			else if (collision_type == ELASTIC) {

				log_printf(DEBUG, "Elastic scatter type collision");

				float phi = (float(rand()) / RAND_MAX) * 2.0 * M_PI;
				curr->_phi = phi;

				/* Isotropic in lab */
				if (isotope->getElasticAngleType() == ISOTROPIC_LAB) {
					float A = float(isotope->getA());
					float energy = curr->_energy * (1.0 -
						(1.0 - isotope->getAlpha())*(float(rand())/RAND_MAX));
					float mu_cm = ((2 * energy / curr->_energy) -
									isotope->getAlpha() - 1.0) / (1.0 -
														isotope->getAlpha());
					curr->_mu = (1.0 + A*mu_cm) / sqrt(A*A + 1.0 + 2*mu_cm*A);

					if (curr->_energy < 4 && isotope->usesThermalScattering())
						curr->_energy =
						isotope->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy = energy;
				}

				/* Isotropic in CM */
				else {
					float mu_cm = (float(rand()) / RAND_MAX) * 2.0 - 1.0;
					int A = isotope->getA();
					float mu_l = (1.0 + A*mu_cm)/(sqrt(A*A + 2*A*mu_cm + 1.0));
					curr->_mu = curr->_mu * mu_l +
							sqrt(1.0 - curr->_mu * curr->_mu) *
							sqrt(1.0 - mu_l * mu_l) * sin(phi);

					if (curr->_energy < 4 && isotope->usesThermalScattering())
						curr->_energy =
						isotope->getThermalScatteringEnergy(curr->_energy);
					else
						curr->_energy *= (A*A + 2*A*mu_cm + 1.0) /
															((A+1.0)*(A+1.0));
				}

				log_printf(DEBUG, "Updated mu= %f, phi = %f, energy = %f",
									curr->_mu, curr->_phi, curr->_energy);
			}

			/* Play Russian Roulette with neutron */
			if (curr->_weight < _weight_low && _use_implicit_capture) {
				if (playRussianRoulette(curr)) {
					log_printf(DEBUG, "Killing neutron");
					iter1 = _neutrons.erase(iter1);
					delete curr;
					--iter1;
				}
				else
					curr->_weight = _weight_avg;
			}
		}

		/******************************************************************
		 ********************  SURFACE INTERSECTIONS  *********************
		 *****************************************************************/
		else {

			/* Check if this neutron is contained by an interior region */
			if (inInteriorRegion(new_x, new_y, new_z))
				nearest = _interior_region->computeMinSurfDist(curr, min_dist);

			/* Find distance to nearest wall along neutron's trajectory */
			else
				nearest = computeMinSurfDist(curr, min_dist);

			updateNeutronTime(curr, min_dist);
			nearest->addNeutron(curr);
			iter1 = _neutrons.erase(iter1);
			--iter1;
		}
	}

	return;
}

/******************************************************************************
 ****************************  Parallelepiped  ********************************
 *****************************************************************************/

Parallelepiped::Parallelepiped() { };


Parallelepiped::~Parallelepiped() { };

void Parallelepiped::setXLeftSurface(XPlane* plane) {
	_x_left_surface = plane;
	_boundaries.push_back(plane);
}

void Parallelepiped::setXRightSurface(XPlane* plane) {
	_x_right_surface = plane;
	_boundaries.push_back(plane);
}

void Parallelepiped::setYLeftSurface(YPlane* plane) {
	_y_left_surface = plane;
	_boundaries.push_back(plane);
}

void Parallelepiped::setYRightSurface(YPlane* plane) {
	_y_right_surface = plane;
	_boundaries.push_back(plane);
}

void Parallelepiped::setZLeftSurface(ZPlane* plane) {
	_z_left_surface = plane;
	_boundaries.push_back(plane);
}

void Parallelepiped::setZRightSurface(ZPlane* plane) {
	_z_right_surface = plane;
	_boundaries.push_back(plane);
}

bool Parallelepiped::contains(float x, float y, float z) {

	if (inInteriorRegion(x, y, z))
		return false;
	else if (x < _x_left_surface->getX() || x > _x_right_surface->getX())
		return false;
	else if (y < _y_left_surface->getY() || y > _y_right_surface->getY())
		return false;
	else if (z < _z_left_surface->getZ() || z > _z_right_surface->getZ())
		return false;
	else if (_x_left_surface->onSurface(x, y, z))
		return false;
	else if (_x_right_surface->onSurface(x, y, z))
		return false;
	else if (_y_left_surface->onSurface(x, y, z))
		return false;
	else if (_y_right_surface->onSurface(x, y, z))
		return false;
	else if (_z_left_surface->onSurface(x, y, z))
		return false;
	else if (_z_right_surface->onSurface(x, y, z))
		return false;
	else
		return true;
}


bool Parallelepiped::onBoundary(neutron* neutron) {

	return (_x_left_surface->onSurface(neutron) ||
			_x_right_surface->onSurface(neutron) ||
			_y_left_surface->onSurface(neutron) ||
			_y_right_surface->onSurface(neutron) ||
			_z_left_surface->onSurface(neutron) ||
			_z_right_surface->onSurface(neutron));
}



/******************************************************************************
 *******************************  XCylinder  **********************************
 *****************************************************************************/

XCylinder::XCylinder() { };


XCylinder::~XCylinder() { };

void XCylinder::setXLeftCircle(XCircle* circle) {
	_left_circle = circle;
	_boundaries.push_back(circle);
}

void XCylinder::setXRightCircle(XCircle* circle) {
	_right_circle = circle;
	_boundaries.push_back(circle);
}

void XCylinder::setOpenXCylinder(OpenXCylinder* cylinder) {
	_cylinder = cylinder;
	_boundaries.push_back(cylinder);
}

bool XCylinder::contains(float x, float y, float z) {

	if (inInteriorRegion(x, y, z))
		return false;
	else if (x < _left_circle->getX0() || x > _right_circle->getX0())
		return false;
	else if (_left_circle->onSurface(x, y, z) ||
			_right_circle->onSurface(x, y, z) ||
					_cylinder->onSurface(x, y, z)) {
		return false;
	}

	float y0 = _left_circle->getY0();
	float z0 = _left_circle->getZ0();
	float r = sqrt((y0 - y)*(y0 - y) + (z0 - z)*(z0 - z));

	if (r >= _left_circle->getRadius())
		return false;
	else
		return true;
}


bool XCylinder::onBoundary(neutron* neutron) {

	return (_left_circle->onSurface(neutron) ||
			_right_circle->onSurface(neutron) ||
			_cylinder->onSurface(neutron));
}


/******************************************************************************
 *******************************  YCylinder  **********************************
 *****************************************************************************/

YCylinder::YCylinder() { };


YCylinder::~YCylinder() { };

void YCylinder::setYLeftCircle(YCircle* circle) {
	_left_circle = circle;
	_boundaries.push_back(circle);
}

void YCylinder::setYRightCircle(YCircle* circle) {
	_right_circle = circle;
	_boundaries.push_back(circle);
}

void YCylinder::setOpenYCylinder(OpenYCylinder* cylinder) {
	_cylinder = cylinder;
	_boundaries.push_back(cylinder);
}

bool YCylinder::contains(float x, float y, float z) {

	if (inInteriorRegion(x, y, z))
		return false;
	else if (y < _left_circle->getY0() || y > _right_circle->getY0())
		return false;
	else if (_left_circle->onSurface(x, y, z) ||
			_right_circle->onSurface(x, y, z) ||
					_cylinder->onSurface(x, y, z)) {
		return false;
	}

	float x0 = _left_circle->getX0();
	float z0 = _left_circle->getZ0();
	float r = sqrt((x0 - x)*(x0 - x) + (z0 - z)*(z0 - z));

	if (r > _left_circle->getRadius())
		return false;

	return true;
}


bool YCylinder::onBoundary(neutron* neutron) {

	return (_left_circle->onSurface(neutron) ||
			_right_circle->onSurface(neutron) ||
			_cylinder->onSurface(neutron));
}


/******************************************************************************
 *******************************  ZCylinder  **********************************
 *****************************************************************************/

ZCylinder::ZCylinder() { };


ZCylinder::~ZCylinder() { };

void ZCylinder::setZLeftCircle(ZCircle* circle) {
	_left_circle = circle;
	_boundaries.push_back(circle);
}

void ZCylinder::setZRightCircle(ZCircle* circle) {
	_right_circle = circle;
	_boundaries.push_back(circle);
}

void ZCylinder::setOpenZCylinder(OpenZCylinder* cylinder) {
	_cylinder = cylinder;
	_boundaries.push_back(cylinder);
}

bool ZCylinder::contains(float x, float y, float z) {

	if (inInteriorRegion(x, y, z))
		return false;
	else if (z < _left_circle->getZ0() || z > _right_circle->getZ0())
		return false;
	else if (_left_circle->onSurface(x, y, z) ||
			_right_circle->onSurface(x, y, z) ||
					_cylinder->onSurface(x, y, z))
		return false;

	float x0 = _left_circle->getX0();
	float y0 = _left_circle->getY0();
	float r = sqrt((x0 - x)*(x0 - x) + (y0 - y)*(y0 - y));

	if (r > _left_circle->getRadius())
		return false;

	return true;
}


bool ZCylinder::onBoundary(neutron* neutron) {

	return (_left_circle->onSurface(neutron) ||
			_right_circle->onSurface(neutron) ||
			_cylinder->onSurface(neutron));
}



/******************************************************************************
 ********************************  Spheroid  **********************************
 *****************************************************************************/

Spheroid::Spheroid() { };

Spheroid::~Spheroid() { };

void Spheroid::setSphere(Sphere* sphere) {
	_sphere = sphere;
	_boundaries.push_back(_sphere);
}


bool Spheroid::contains(float x, float y, float z) {

	if (inInteriorRegion(x, y, z))
		return false;

	float x0 = _sphere->getX0();
	float y0 = _sphere->getY0();
	float z0 = _sphere->getZ0();

	float r = sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0));

	if (inInteriorRegion(x, y, z))
		return false;
	else if (r <= _sphere->getRadius())
		return true;
	else
		return false;
}


bool Spheroid::onBoundary(neutron* neutron) {
	return _sphere->onSurface(neutron);
}


/******************************************************************************
 ****************************  Spherical Shell  *******************************
 *****************************************************************************/

SphericalShell::SphericalShell() { };


SphericalShell::~SphericalShell() { }


void SphericalShell::setInnerSphere(Sphere* sphere) {
	_inner = sphere;
	_boundaries.push_back(sphere);
}

void SphericalShell::setOuterSphere(Sphere* sphere) {
	_outer = sphere;
	_boundaries.push_back(_outer);
}


bool SphericalShell::contains(float x, float y, float z) {

	if (inInteriorRegion(x, y, z))
		return false;

	float x0 = _inner->getX0();
	float y0 = _inner->getY0();
	float z0 = _inner->getZ0();

	float r = sqrt((x - x0)*(x - x0) + (y - y0)*(y - y0) + (z - z0)*(z - z0));

	if (inInteriorRegion(x, y, z))
		return false;
	else if (r >= _inner->getRadius() && r <= _outer->getRadius())
		return true;
	else
		return false;
}


bool SphericalShell::onBoundary(neutron* neutron) {
	return(_inner->onSurface(neutron) || _outer->onSurface(neutron));
}
