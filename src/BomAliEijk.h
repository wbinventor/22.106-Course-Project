/*
 * BomAliEijk.h
 *
 *  Created on: Apr 26, 2012
 *      Author: wboyd
 */

#ifndef BomAliEijk_H_
#define BomAliEijk_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include <omp.h>
#include "log.h"
#include "gnuplot.h"
#include "Timer.h"
#include "BatchBinSet.h"
#include "Options.h"
#include "Isotope.h"
#include "Material.h"
#include "Neutron.h"
#include "Region.h"


void BomAliEijk(Options* options);

#endif /* BOMALIEIJK_H_ */

