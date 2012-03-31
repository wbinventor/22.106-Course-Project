/*
 * binner.h
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#ifndef BINNER_H_
#define BINNER_H_

#include <limits>
#include <math.h>
#include <set>
#include <utility>
#include <omp.h>
#include "log.h"
#include "arraycreator.h"

#define NUM_THREADS 4


/* Bin spacing types */
typedef enum binTypes {
	EQUAL,
	LOGARITHMIC,
	OTHER
} binType;


/* Type of tallies */
typedef enum tallyTypes {
	FLUX_SPATIAL,
	FLUX_ENERGY,
	CAPTURE_RATE_SPATIAL,
	CAPTURE_RATE_ENERGY,
	ABSORPTION_RATE_SPATIAL,
	ABSORPTION_RATE_ENERGY,
	ELASTIC_RATE_SPATIAL,
	ELASTIC_RATE_ENERGY,
	FISSION_RATE_SPATIAL,
	FISSION_RATE_ENERGY,
	TRANSPORT_RATE_SPATIAL,
	TRANSPORT_RATE_ENERGY,
	COLLISION_RATE_SPATIAL,
	COLLISION_RATE_ENERGY,
	DIFFUSION_RATE_SPATIAL,
	DIFFUSION_RATE_ENERGY
} tallyType;


/**
 * This class represents a set of bins which are defined by a set of values
 * defining the edges between bins. This class holds the edges, the centers
 * between bins. It also allows for tallies to be made within each bin.
 */
class Binner{
private:
	char* _name;
	int _num_bins;
	int _num_threads;
	float* _edges;
	double* _centers;
	double* _tallies;
	double* _tallies_squared;
	double* _tally_acc;
	std::set<int> _acc_indices;

	bool _statistics_compute;
	double* _bin_mu;
	double* _bin_variance;
	double* _bin_std_dev;
	double* _bin_rel_err;

	float _bin_delta;
	binType _bin_type;
	tallyType _tally_type;
	char* _isotopes;

	omp_lock_t _lock;
public:
	Binner();
	virtual ~Binner();
	char* getBinnerName();
	int getNumBins();
	float* getBinEdges();
	double* getBinCenters();
	float getBinDelta();
	float getBinDelta(float sample);
	binType getBinType();
	tallyType getTallyType();
	double* getTallies();
	double getTally(int bin_index);
	double* getTalliesSquared();
	double getTallySquared(int bin_index);
	double getMaxTally();
	double getMinTally();
	double* getBinMu();
	double* getBinVariance();
	double* getBinStdDev();
	double* getBinRelativeError();
	int getBinIndex(float sample);
	char* getIsotopes();

	void setBinnerName(char* name);
	void setNumThreads(int num_threads);
	void setTallyType(tallyType type);
	void setBinEdges(float* edges, int num_edges);
	void setIsotopes(char* isotopes);

	void generateBinEdges(float start, float end, int num_bins, binType type);
	void generateBinCenters();

	void tally(int thread_num, float sample);
	void weightedTally(int thread_num, float sample, float weight);
	void processTallyAccumulators();
	void normalizeTallies();
	void normalizeTallies(float scale_factor);
	void computeScaledHistoryStatistics(float scale_factor);
	void outputHistoryStatistics(const char* filename);
};

#endif /* BINNER_H_ */
