/*
 * Neutron.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Neutron.h"

/**
 * Creates a new neutron with a set of default attribute values
 * @return a pointer to the new neutron
 */
neutron* initializeNewNeutron() {

	/* Allocate memory for new neutron struct */
	neutron* new_neutron = new neutron;

	/* Assign default values to neutron attributes */
	new_neutron->_x = 0.0;
	new_neutron->_y = 0.0;
	new_neutron->_z = 0.0;
	new_neutron->_energy = 0.0;
	new_neutron->_phi = 0.0;
	new_neutron->_mu = 0.0;
	new_neutron->_weight = 1.0;
	new_neutron->_energy = 0.0;
	new_neutron->_time = 0.0;
	new_neutron->_thread_num = 0;
	return new_neutron;
}


void updateNeutronTime(neutron* neutron, float path_length) {

	/* Compute velocity in m/s */
//	double velocity = sqrt(2.0 * neutron->_energy / M_N) * C;
	double velocity = sqrt(2.0 * neutron->_energy / M_N) * C * C;

	/* Update neutron's time */
	neutron->_time += path_length / velocity;

	return;
}


void updateNeutronTime(neutron* neutron, float x, float y, float z) {

	float x0 = neutron->_x;
	float y0 = neutron->_y;
	float z0 = neutron->_z;

	/* Compute path length */
	float path_length = sqrt((x0 - x)*(x0 - x) +
							 (y0 - y)*(y0 - y) +
							 (z0 - z)*(z0 - z));

	/* Compute velocity in m/s */
	double velocity = sqrt(2.0 * neutron->_energy / M_N) * C;

	/* Update neutron's time */
	neutron->_time += path_length / velocity;

	return;
}



std::string neutronToString(neutron* neutron) {

	std::stringstream result;

	result << "neutron x = " << neutron->_x << " y = " << neutron->_y;
	result << " z = " << neutron->_z << " energy = " << neutron->_energy;
	result << " phi = " << neutron->_phi << " mu = " << neutron->_mu;
	result << " weight = " << neutron->_weight << " time = " << neutron->_time;
	result << " thread = " << neutron->_thread_num;

	return result.str();
}
