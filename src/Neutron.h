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

/* Structure to represent a neutron */
struct neutron {
	float _x, _y, _z;
	float _mu, _phi;
	float _energy;
	float _weight;
	bool _infuel;
	int _thread_num;
};


neutron* initializeNewNeutron();

#endif /* NEUTRON_H_ */
