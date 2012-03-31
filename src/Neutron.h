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
std::string neutronToString(neutron* neutron);

#endif /* NEUTRON_H_ */
