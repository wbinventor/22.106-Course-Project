/*
 * pset3.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sstream>
#include <math.h>
#include <omp.h>
#include "log.h"
#include "gnuplot.h"
#include "Timer.h"
#include "Options.h"
#include "BatchBinSet.h"
#include "Isotope.h"
#include "Material.h"
#include "Neutron.h"
#include "Region.h"
#include "BomAliEijk.h"
#include "DatemaBomEijk.h"
#include "BaysoySubasi.h"
#include "testregions.h"


int main(int argc, const char **argv) {

	Options options(argc, argv);
	Timer timer;

	log_setlevel(options.getVerbosity());

	if (options.testRegions())
		testRegions(&options);

	timer.start();
//	BaysoySubasi(&options);

	if (options.bomAliEijk())
		BomAliEijk(&options);

	if (options.datemaBomEijk())
		DatemaBomEijk(&options);
	timer.stop();
	timer.recordSplit("Time");

	timer.printSplits();

	return 0;
}
