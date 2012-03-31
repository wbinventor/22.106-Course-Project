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
#include "Neutron.h"
#include "Material.h"
#include "Isotope.h"

class Isotope;
class Material;

/* Bin spacing types */
typedef enum binTypes {
	LINEAR,
	LOGARITHMIC,
	OTHER
} binType;


/* Tally type */
typedef enum tallyDomainTypes {
	X,
	Y,
	Z,
	ENERGY,
	TIME
} tallyDomainType;


/**
 * This class represents a set of bins which are defined by a set of values
 * defining the edges between bins. This class holds the edges, the centers
 * between bins. It also allows for tallies to be made within each bin.
 */
class Binner{
protected:
	char* _name;
	int _num_bins;
	int _num_threads;
	float* _edges;
	double* _centers;
	double* _tallies;
	double* _tallies_squared;
	double* _tally_acc;

	bool _statistics_compute;
	double* _bin_mu;
	double* _bin_variance;
	double* _bin_std_dev;
	double* _bin_rel_err;

	float _bin_delta;
	binType _bin_type;
	tallyDomainType _tally_domain_type;
	char* _isotopes;

	omp_lock_t _lock;
public:
	Binner();
	virtual ~Binner();
	char* getBinnerName();
	binType getBinType();
	tallyDomainType getTallyDomainType();
	int getNumBins();
	float* getBinEdges();
	double* getBinCenters();
	float getBinDelta();
	float getBinDelta(neutron* neutron);
	double* getTallies();
	double getTally(int bin_index);
	double* getTalliesSquared();
	double getTallySquared(int bin_index);
	double* getBinMu();
	double* getBinVariance();
	double* getBinStdDev();
	double* getBinRelativeError();
	int getBinIndex(neutron* neutron);
	int getBinIndex(float sample);
	char* getIsotopes();

	void setBinnerName(char* name);
	void setNumThreads(int num_threads);
	void setTallyDomainType(tallyDomainType type);
	void setBinEdges(float* edges, int num_edges);
	void setIsotopes(char* isotopes);

	void generateBinEdges(float start, float end, int num_bins, binType type);
	void generateBinCenters();

	void tally(int thread_num, neutron* neutron);
	void tally(int thread_num, float sample);
	virtual void weightedTally(neutron* neutron, float sigma_t,
			int energy_index, Material* material, Isotope* isotope) =0;
	void processTallyAccumulators();
	void computeScaledHistoryStatistics(float scale_factor);
	void outputHistoryStatistics(const char* filename);
};


class FluxBinner: public Binner {
public:
	FluxBinner();
	virtual ~FluxBinner();
	void weightedTally(neutron* neutron, float sigma_t,
						int energy_index, Material* material, Isotope* isotope);
};

class CaptureRateBinner: public Binner {
public:
	CaptureRateBinner();
	virtual ~CaptureRateBinner();
	void weightedTally(neutron* neutron, float sigma_t,
						int energy_index, Material* material, Isotope* isotope);
};


class AbsorptionRateBinner: public Binner {
public:
	AbsorptionRateBinner();
	virtual ~AbsorptionRateBinner();
	void weightedTally(neutron* neutron, float sigma_t,
						int energy_index, Material* material, Isotope* isotope);
};


class ElasticRateBinner: public Binner {
public:
	ElasticRateBinner();
	virtual ~ElasticRateBinner();
	void weightedTally(neutron* neutron, float sigma_t,
						int energy_index, Material* material, Isotope* isotope);
};


class InelasticRateBinner: public Binner {
public:
	InelasticRateBinner();
	virtual ~InelasticRateBinner();
	void weightedTally(neutron* neutron, float sigma_t,
						int energy_index, Material* material, Isotope* isotope);
};


class FissionRateBinner: public Binner {
public:
	FissionRateBinner();
	virtual ~FissionRateBinner();
	void weightedTally(neutron* neutron, float sigma_t,
						int energy_index, Material* material, Isotope* isotope);
};


class TransportRateBinner: public Binner {
public:
	TransportRateBinner();
	virtual ~TransportRateBinner();
	void weightedTally(neutron* neutron, float sigma_t,
						int energy_index, Material* material, Isotope* isotope);
};


class CollisionRateBinner: public Binner {
public:
	CollisionRateBinner();
	virtual ~CollisionRateBinner();
	void weightedTally(neutron* neutron, float sigma_t,
						int energy_index, Material* material, Isotope* isotope);
};


class DiffusionRateBinner: public Binner {
public:
	DiffusionRateBinner();
	virtual ~DiffusionRateBinner();
	void weightedTally(neutron* neutron, float sigma_t,
						int energy_index, Material* material, Isotope* isotope);
};


#endif /* BINNER_H_ */
