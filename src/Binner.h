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


/* Binner types (subclasses of Binner) */
typedef enum binnerTypes {
	FLUX,
	CAPTURE_RATE,
	ABSORPTION_RATE,
	ELASTIC_RATE,
	INELASTIC_RATE,
	FISSION_RATE,
	TRANSPORT_RATE,
	COLLISION_RATE,
	DIFFUSION_RATE

} binnerType;


/* Bin spacing types */
typedef enum spacingTypes {
	LINEAR,
	LOGARITHMIC,
	OTHER
} spacingType;


/* Tally type */
typedef enum tallyDomainTypes {
	X,
	Y,
	Z,
	R_X,
	R_Y,
	R_Z,
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
	float* _edges;
	float* _centers;
	double* _tallies;

	float _bin_delta;
	binnerType _binner_type;
	spacingType _spacing_type;
	tallyDomainType _tally_domain_type;
	char* _isotopes;

	omp_lock_t _lock;

public:
	Binner();
	virtual ~Binner();
	char* getBinnerName();
	binnerType getBinnerType();
	spacingType getBinType();
	tallyDomainType getTallyDomainType();
	int getNumBins();
	float* getBinEdges();
	float* getBinCenters();
	float getBinDelta();
	float getBinDelta(neutron* neutron);
	double* getTallies();
	double getTally(int bin_index);
	int getBinIndex(neutron* neutron);
	int getBinIndex(float sample);
	char* getIsotopes();

	void setBinnerName(char* name);
	void setNumThreads(int num_threads);
	void setBinnerType(binnerType type);
	void setTallyDomainType(tallyDomainType type);
	void setBinEdges(float* edges, int num_edges);
	void setIsotopes(char* isotopes);

	void generateBinEdges(float start, float end, int num_bins, spacingType type);
	void generateBinCenters();

	void tally(float sample);
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
