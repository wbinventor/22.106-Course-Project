/*
 * Neutron.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef NEUTRON_H_
#define NEUTRON_H_

#include <string>
#include <sstream>
#include <math.h>

#define C 3E2				/* Speed of light (m/us) */
#define M_N 939.565E6		/* Mass of neutron in eV */


/* Structure to represent a neutron */
struct neutron {
	float _x, _y, _z;
	float _mu, _phi;
	float _energy;
	float _weight;
	double _time;
	int _thread_num;
};


neutron* initializeNewNeutron();
void updateNeutronTime(neutron* neutron, float path_length);
void updateNeutronTime(neutron* neutron, float new_x, float new_y, float new_z);
std::string neutronToString(neutron* neutron);

#endif /* NEUTRON_H_ */
