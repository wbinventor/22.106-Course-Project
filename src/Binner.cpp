/*
 * binner.cpp
 *
 *  Created on: Mar 13, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "Binner.h"


/**
 * Default Binner constructor
 */
Binner::Binner() {

	/* Default name is blank */
	_name = (char*)"";

	_num_threads = 1;

	 /* Sets the default delta between bins to zero */
	_bin_delta = 0;

	/* Default is to tally all isotopes */
	_isotopes = (char*)"all";

	_statistics_compute = false;

	omp_init_lock(&_lock);
}


/**
 * Binner destructor deletes memory for tallies, number of tallies,
 * bin centers and bin edges if they have been created
 */
Binner::~Binner() {

	if (_num_bins != 0) {
		delete [] _tallies;
		delete [] _tally_acc;
		delete [] _tallies_squared;
		delete [] _bin_mu;
		delete [] _bin_variance;
		delete [] _bin_std_dev;
		delete [] _bin_rel_err;
		delete [] _centers;

		/* If the bin edges were generated rather than defined by the user */
		if (_spacing_type != OTHER)
			delete [] _edges;
	}

	omp_destroy_lock(&_lock);
}


/**
 * Returns the name of this Binner as specified by the user
 * @return the Binner's name
 */
char* Binner::getBinnerName() {
	return _name;
}


/**
 * Returns the type of Binner (FLUX, CAPTURE, ABSORPTION, etc)
 * @return the type of Binner (the subclass type)
 */
binnerType Binner::getBinnerType() {
	return _binner_type;
}


/**
 * Returns the type of bin (LINEAR, LOGARITHMIC, OTHER)
 * @return the bin type
 */
spacingType Binner::getBinType() {
	return _spacing_type;
}


/**
 * Returns the type of domain from which samples are binned using
 * this binner
 * @return the tally domain type (X, Y, Z, ENERGY, or TIME)
 */
tallyDomainType Binner::getTallyDomainType() {
	return _tally_domain_type;
}


/**
 * Returns the number of bins
 * @return the number of bins
 */
int Binner::getNumBins() {
	return _num_bins;
}


/**
 * Returns a double array of bin edge values
 * @return array of bin edge values
 */
float* Binner::getBinEdges() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return bin edges for Binner %s since "
				 "the bins have not yet been created", _name);

	 return _edges;
}


/**
 * Returns a double array of bin center values
 * @return array of bin center values
 */
double* Binner::getBinCenters() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return bin centers for Binner %s since the "
				 "bins have not yet been created", _name);

	 return _centers;
}


/**
 * Returns the delta spacing between bins. NOTE: this value is only non-zero
 * for LINEAR and LOGARITHMIC bin types
 * @return the spacing between bins
 */
float Binner::getBinDelta() {

	if (_spacing_type != LINEAR && _spacing_type != LOGARITHMIC)
		log_printf(ERROR, "Cannot return the bin delta for a set"
				" of bins that do not have LINEAR or LOGARITHMIC spaced"
				" bin edges");

	return _bin_delta;
}


/**
 * Returns the distance between the two bin edges which sandwich a neutron
 * @param neutron the neutron of interest
 * @param the width of the bin that this neutron is inside
 */
float Binner::getBinDelta(neutron* neutron) {

	/* If this Binner uses equally spaced bins in linear or logarithmic
	 * space, return the bin delta */
	if (_spacing_type == LINEAR || _spacing_type == LOGARITHMIC)
		return _bin_delta;

	/* If instead this Binner uses irregularly spaced bin edges defined
	 * by a users, compute bin delta of the bin around the sample */
	else {
		int bin_index = getBinIndex(neutron);
		return (_edges[bin_index] - _edges[bin_index-1]);
	}
}


/**
 * Returns a double array of the tallies within each bin
 * @return an array of
 */
double* Binner::getTallies() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return tallies for Binner %s since the "
				 "bins have not yet been created", _name);

	 return _tallies;
}


/**
 * Returns a specific tally for a specific bin
 * @param bin_index the index for the bin of interest
 * @return the tally within that bin
 */
double Binner::getTally(int bin_index) {

	if (bin_index < 0 || bin_index >= _num_bins)
		log_printf(ERROR, "Tried to get a tally for a bin index for Binner %s"
				"which does not exist: %d, num_bins = %d", _name, bin_index,
				_num_bins);

	return _tallies[bin_index];
}


/**
 * Returns a double array of the sum of the squares of the tallies
 * within each bin
 * @return an array of
 */
double* Binner::getTalliesSquared() {
	 if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return tallies squared for Binner %s since "
				 "the bins have not yet been created", _name);

	 return _tallies_squared;
}


/**
 * Returns a specific sum of the tallies squared for a specific bin
 * @param bin_index the index for the bin of interest
 * @return the sum of the squares of the tallies within that bin
 */
double Binner::getTallySquared(int bin_index) {

	if (bin_index < 0 || bin_index >= _num_bins)
		log_printf(ERROR, "Tried to get a tally for a bin index for Binner %s"
				"which does not exist: %d, num_bins = %d", _name, bin_index,
				_num_bins);

	return _tallies[bin_index];
}


/**
 * Returns a pointer to an array of bin averages if they have been
 * computed
 * @return a double array of bin averages for each bin
 */
double* Binner::getBinMu() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for "
				"Binner %s so bin mu cannot be returned", _name);

	return _bin_mu;
}


/**
 * Returns a pointer to an array of bin variances if they have been
 * computed
 * @return a double array of bin variances for each bin
 */
double* Binner::getBinVariance() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for "
				"Binner %s so bin variance cannot be returned", _name);

	return _bin_variance;
}


/**
 * Returns a pointer to an array of bin standard deviations if they have
 * been computed
 * @return a double array of bin standard deviations for each bin
 */
double* Binner::getBinStdDev() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for "
				"Binner %s so bin std dev cannot be returned", _name);

	return _bin_std_dev;
}


/**
 * Returns a pointer to an array of batch relative errors if they have been
 * computed
 * @return a double array of batch relative errors for each bin
 */
double* Binner::getBinRelativeError() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for "
		"Binner %s so bin relative error cannot be returned", _name);

	return _bin_rel_err;
}


/**
 * Finds the bin index for a sample in a set of bins. If the samples
 * is outside the bounds of all bins, it returns infinity
 * @param sample the sample value of interest
 * @return the bin index for the sample
 */
int Binner::getBinIndex(neutron* neutron) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return a bin index for Binner %s since "
				 "the bins have not yet been created", _name);

	/* Set index to infinity to begin with */
	int index = std::numeric_limits<float>::infinity();

	float sample;

	if (_tally_domain_type == X)
		sample = neutron->_x;
	else if (_tally_domain_type == Y)
		sample = neutron->_y;
	else if (_tally_domain_type == Z)
		sample = neutron->_z;
	else if (_tally_domain_type == TIME)
		sample = neutron->_time;
	else
		sample = neutron->_energy;

	/* if the sample is equal to the last bin edge, return the last bin */
	if (sample == _edges[_num_bins])
		return _num_bins-1;

	/* Equally spaced bins */
	if (_spacing_type == LINEAR)
		index = int((sample - _edges[0]) / _bin_delta);

	/* Logarithmically spaced bins */
	else if (_spacing_type == LOGARITHMIC)
		index = int((log10(sample) - log10(_edges[0])) / _bin_delta);

	/* If the bin_type == OTHER then the bin edges were not generated by
	 * generateEqualBinEdges, so use a brute force search to find the bin */
	else {

		/* Loop over all bin edges to find the correct bin index */
		for (int i=0; i <= _num_bins; i++) {
			if (sample >= _edges[i] && sample < _edges[i+1]) {
				index = i;
				break;
			}
		}
	}

	/* If this sample was not contained within a bin set index to infinity*/
	if (index > _num_bins)
		index = std::numeric_limits<float>::infinity();

	return index;
}

/**
 * Finds the bin index for a sample in a set of bins. If the samples
 * is outside the bounds of all bins, it returns infinity
 * @param sample the sample value of interest
 * @return the bin index for the sample
 */
int Binner::getBinIndex(float sample) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot return a bin index for Binner %s since "
				 "the bins have not yet been created", _name);

	/* Set index to infinity to begin with */
	int index = std::numeric_limits<float>::infinity();

	/* if the sample is equal to the last bin edge, return the last bin */
	if (sample == _edges[_num_bins])
		return _num_bins-1;

	/* Equally spaced bins */
	if (_spacing_type == LINEAR)
		index = int((sample - _edges[0]) / _bin_delta);

	/* Logarithmically spaced bins */
	else if (_spacing_type == LOGARITHMIC)
		index = int((log10(sample) - log10(_edges[0])) / _bin_delta);

	/* If the bin_type == OTHER then the bin edges were not generated by
	 * generateEqualBinEdges, so use a brute force search to find the bin */
	else {

		/* Loop over all bin edges to find the correct bin index */
		for (int i=0; i <= _num_bins; i++) {
			if (sample >= _edges[i] && sample < _edges[i+1]) {
				index = i;
				break;
			}
		}
	}

	/* If this sample was not contained within a bin set index to infinity*/
	if (index > _num_bins)
		index = std::numeric_limits<float>::infinity();

	return index;
}


/**
 * Return the isotopes that this binner is meant to tally
 * @return a character array of the isotope's name or "all" for all isotopes
 */
char* Binner::getIsotopes() {
	return _isotopes;
}



/**
 * Sets this Binner's name
 * @param name the name of the Binner
 */
void Binner::setBinnerName(char* name) {
	_name = name;
}


/**
 * Sets the number of threads used in the simulation. The Binner
 * must know this to ensure that asynchronous tallying is accumulated
 * in S1 and S2 counters appropriately
 * @param _num_threads the number of OpenMP threads
 */
void Binner::setNumThreads(int num_threads) {
	_num_threads = num_threads;
	return;
}


/**
 * Set the type of subclass of Binner
 * @param type the type of Binner (FLUX, CAPTURE, ABSORPTION, etc)
 */
void Binner::setBinnerType(binnerType type) {
	_binner_type = type;
}


/**
 * Set the domain over which neutrons are tallied into this Binner
 * @param type the tally domain type (X, Y, Z, TIME, or ENERGY)
 */
void Binner::setTallyDomainType(tallyDomainType type) {
	_tally_domain_type = type;
}


/**
 * Set a user-defined double array of bin edge values
 * @param edges the array of bin edges
 * @param num_bins the number of bins
 */
void Binner::setBinEdges(float* edges, int num_bins) {

	_num_bins = num_bins;
	_edges = edges;
	_spacing_type = OTHER;

	/* Set all tallies to zero by default */
	_tallies = new double[num_bins];
	_tallies_squared = new double[num_bins];
	_tally_acc = new double[num_bins*_num_threads];

	/* Allocate memory for history-based bin statistics */
	_bin_mu = new double[_num_bins];
	_bin_variance = new double[_num_bins];
	_bin_std_dev = new double[_num_bins];
	_bin_rel_err = new double[_num_bins];

	/* Loop over tallies and set to zero */
	for (int i=0; i < _num_bins; i++) {
		_tallies[i] = 0.0;
		_tallies_squared[i] = 0.0;
		for (int j=0; j < _num_threads; j++)
			_tally_acc[i*_num_threads+j] = 0.0;
	}

	/* Create an array of the center values between bins */
	generateBinCenters();
}



/**
 * Set the isotope that this binner is meant to tally
 * @param isotopes a character array of the isotope's name or
 * "all" for all isotopes
 */
void Binner::setIsotopes(char* isotopes) {
	_isotopes = isotopes;
}


/**
 * Generate edges between bins defined by a start and end point
 * @param start first bin edge value
 * @param end last bin edge value
 * @param num_bins the number of bins to be created
 * @param type the type of bins (LINEAR or LOGARITHMIC)
 */
void Binner::generateBinEdges(float start, float end, int num_bins,
												spacingType type) {
	if (start == end)
		log_printf(ERROR, "Unable to create bins for Binner %s between"
				"the same start and end points: %f", _name, start);

	_num_bins = num_bins;
	_spacing_type = type;

	/* Allocate memory for tallies */
	_tallies = new double[num_bins];
	_tallies_squared = new double[num_bins];
	_tally_acc = new double[num_bins*_num_threads];

	/* Allocate memory for history-based bin statistics */
	_bin_mu = new double[_num_bins];
	_bin_variance = new double[_num_bins];
	_bin_std_dev = new double[_num_bins];
	_bin_rel_err = new double[_num_bins];

	/* Set all tallies to zero by default */
	for (int i=0; i < num_bins; i++) {
		_tallies[i] = 0;
		_tallies_squared[i] = 0.0;
		for (int j=0; j < _num_threads; j++)
			_tally_acc[i*_num_threads+j] = 0.0;
	}

	/* Equal spacing between bins */
	if (type == LINEAR) {
		_bin_delta = float(end - start) / float(_num_bins);

		/* Generate points from start to end for each bin edge */
		_edges = linspace<float, float>(start, end, num_bins+1);
	}

	/* Logarithmically equal spacing between bins */
	else if (type == LOGARITHMIC) {
		_bin_delta = float(log10(end) - log10(start)) / float(_num_bins);

		/* Generate points from start to end for each bin edge */
		_edges = logspace<float, float>(start, end, num_bins+1);
	}

	else
		log_printf(ERROR, "Bin type %d is not yet implemented for Binner %s",
															_name, type);

	/* Create an array of the center values between bins */
	generateBinCenters();

	return;
}


/**
 * Compute the center points between bin edges for this Binner's bins
 */
void Binner::generateBinCenters() {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot generate bin centers for Binner %s since "
				 "the bins have not yet been created", _name);

	/* Allocate memory for the bin centers array */
	_centers = new double[_num_bins];

	/* Loop over all bins and find the midpoint between edges */
	for (int i=0; i < _num_bins; i++)
		_centers[i] = (_edges[i] + _edges[i+1]) / 2.0;

	return;
}


/**
 * Tallies unity for a sample
 * @param samples array of samples to tally
 * @param thread_num the thread number making this tally
 */
void Binner::tally(int thread_num, neutron* neutron) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		_tally_acc[bin_index*_num_threads+thread_num]++;
		omp_unset_lock(&_lock);
	}

	return;
}


/**
 * Tallies unity for a sample
 * @param samples array of samples to tally
 * @param thread_num the thread number making this tally
 */
void Binner::tally(int thread_num, float sample) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(sample);

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		_tally_acc[bin_index*_num_threads+thread_num]++;
		omp_unset_lock(&_lock);
	}

	return;
}


/**
 * This method processes all of the tally accumulators for each bin
 * and adds each accumulator, and each accumulator squared, to the
 * tallies and tallies squared. This is akin to the S1 and S2 counters
 * dicussed in 22.106. At the end, the tally accumulators are reset to 0
 */
void Binner::processTallyAccumulators() {

	/* Loop over all bins */
	for (int i=0; i < _num_bins; i++) {

		/* Loop over all threads */
		for (int j=0; j < _num_threads; j++) {
			_tallies[i] += _tally_acc[i*_num_threads+j];
			_tallies_squared[i] += (_tally_acc[i*_num_threads+j] *
									_tally_acc[i*_num_threads+j]);
			_tally_acc[i*_num_threads+j] = 0.0;
		}
	}

	return;
}


void Binner::computeScaledHistoryStatistics(float factor) {

	if (_num_bins == 0)
		log_printf(ERROR, "Cannot compute bin statistics for Binner %s"
				" since the bins have not yet been generated", _name);

	/* Loop over each bin */
	for (int i=0; i < _num_bins; i++) {
		_bin_mu[i] = _tallies[i] / double(factor);
		_bin_variance[i] = (1.0 / (double(factor) - 1.0)) *
						((_tallies_squared[i] / double(factor))
									- (_bin_mu[i]*_bin_mu[i]));
		_bin_std_dev[i] = sqrt(_bin_variance[i]);
		_bin_rel_err[i] = _bin_std_dev[i] / _bin_mu[i];
	}

	_statistics_compute = true;

	return;
}


/**
 * Outputs the bin statistics (if they have been computed) to an
 * ASCII file
 * @param filename the output filename
 */
void Binner::outputHistoryStatistics(const char* filename) {

	if (_num_bins == 0)
		log_printf(ERROR, "Cannot output bin statistics for Binner %s "
				"since the bins have not yet been generated", _name);

	else if (!_statistics_compute)
		log_printf(ERROR, "Cannot output bin statistics for Binner %s "
				"since they have not yet been computed", _name);

	/* Create output file */
	FILE* output_file;
	output_file = fopen(filename, "w");

	/* Print header to output file */
	fprintf(output_file, "Bin center, Mu, Variance, Std Dev, Rel Err\n");

	/* Loop over each bin and print mu, var, std dev and rel err */
	for (int i=0; i < _num_bins; i++) {
		fprintf(output_file, "%1.10f, %1.10f, %1.10f, %1.10f, %1.10f\n",
				_centers[i], _bin_mu[i], _bin_variance[i], _bin_std_dev[i],
														_bin_rel_err[i]);
	}

	fclose(output_file);

	return;
}


/******************************************************************************
 ********************************  FluxBinner  ********************************
 *****************************************************************************/

FluxBinner::FluxBinner() {
	setBinnerType(FLUX);
};

FluxBinner::~FluxBinner() { };


void FluxBinner::weightedTally(neutron* neutron, float sigma_t,
					int energy_index, Material* material, Isotope* isotope) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);
	int thread_num = neutron->_thread_num;

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		_tally_acc[bin_index*_num_threads+thread_num] += neutron->_weight /
																sigma_t;
		omp_unset_lock(&_lock);
	}

	return;
}



/******************************************************************************
 ****************************  CaptureRateBinner  *****************************
 *****************************************************************************/

CaptureRateBinner::CaptureRateBinner() {
	setBinnerType(CAPTURE_RATE);
};


CaptureRateBinner::~CaptureRateBinner() { };

void CaptureRateBinner::weightedTally(neutron* neutron, float sigma_t,
					int energy_index, Material* material, Isotope* isotope) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);
	int thread_num = neutron->_thread_num;

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		if (strcmp(_isotopes, "all"))
			_tally_acc[bin_index*_num_threads+thread_num] +=
					material->getCaptureMacroXS(energy_index) / sigma_t;
		else
			_tally_acc[bin_index*_num_threads+thread_num] +=
			material->getIsotope(_isotopes)->getCaptureXS(energy_index) /
																sigma_t;
		omp_unset_lock(&_lock);
	}

	return;
}



/******************************************************************************
 **************************  AbsorptionRateBinner  ****************************
 *****************************************************************************/

AbsorptionRateBinner::AbsorptionRateBinner() {
	setBinnerType(ABSORPTION_RATE);
};

AbsorptionRateBinner::~AbsorptionRateBinner() { };

void AbsorptionRateBinner::weightedTally(neutron* neutron, float sigma_t,
							int energy_index, Material* material, Isotope* isotope) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);
	int thread_num = neutron->_thread_num;

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		if (strcmp(_isotopes, "all"))
			_tally_acc[bin_index*_num_threads+thread_num] +=
					material->getAbsorbMacroXS(energy_index) / sigma_t;
		else
			_tally_acc[bin_index*_num_threads+thread_num] +=
			material->getIsotope(_isotopes)->getAbsorbXS(energy_index) /
																sigma_t;
		omp_unset_lock(&_lock);
	}

	return;
}



/******************************************************************************
 ****************************  ElasticRateBinner  *****************************
 *****************************************************************************/

ElasticRateBinner::ElasticRateBinner() {
	setBinnerType(ELASTIC_RATE);
};

ElasticRateBinner::~ElasticRateBinner() { };

void ElasticRateBinner::weightedTally(neutron* neutron, float sigma_t,
							int energy_index, Material* material, Isotope* isotope) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);
	int thread_num = neutron->_thread_num;

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		if (strcmp(_isotopes, "all"))
			_tally_acc[bin_index*_num_threads+thread_num] +=
					material->getElasticMacroXS(energy_index) / sigma_t;
		else
			_tally_acc[bin_index*_num_threads+thread_num] +=
			material->getIsotope(_isotopes)->getElasticXS(energy_index) /
																sigma_t;
		omp_unset_lock(&_lock);
	}

	return;
}


/******************************************************************************
 ***************************  InelasticRateBinner  ****************************
 *****************************************************************************/

InelasticRateBinner::InelasticRateBinner() {
	setBinnerType(INELASTIC_RATE);
};

InelasticRateBinner::~InelasticRateBinner() { };

void InelasticRateBinner::weightedTally(neutron* neutron, float sigma_t,
							int energy_index, Material* material, Isotope* isotope) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);
	int thread_num = neutron->_thread_num;

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		if (strcmp(_isotopes, "all"))
			_tally_acc[bin_index*_num_threads+thread_num] +=
					material->getInelasticMacroXS(energy_index) / sigma_t;
		else
			_tally_acc[bin_index*_num_threads+thread_num] +=
			material->getIsotope(_isotopes)->getInelasticXS(energy_index) /
																sigma_t;
		omp_unset_lock(&_lock);
	}

	return;
}


/******************************************************************************
 ****************************  FissionRateBinner  *****************************
 *****************************************************************************/

FissionRateBinner::FissionRateBinner() {
	setBinnerType(FISSION_RATE);
};

FissionRateBinner::~FissionRateBinner() { };

void FissionRateBinner::weightedTally(neutron* neutron, float sigma_t,
							int energy_index, Material* material, Isotope* isotope) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);
	int thread_num = neutron->_thread_num;

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		if (strcmp(_isotopes, "all"))
			_tally_acc[bin_index*_num_threads+thread_num] +=
					material->getFissionMacroXS(energy_index) / sigma_t;
		else
			_tally_acc[bin_index*_num_threads+thread_num] +=
			material->getIsotope(_isotopes)->getFissionXS(energy_index) /
																sigma_t;
		omp_unset_lock(&_lock);
	}

	return;
}


/******************************************************************************
 ***************************  TransportRateBinner  ****************************
 *****************************************************************************/

TransportRateBinner::TransportRateBinner() {
	setBinnerType(TRANSPORT_RATE);
};

TransportRateBinner::~TransportRateBinner() { };

void TransportRateBinner::weightedTally(neutron* neutron, float sigma_t,
							int energy_index, Material* material, Isotope* isotope) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);
	int thread_num = neutron->_thread_num;

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		if (strcmp(_isotopes, "all"))
			_tally_acc[bin_index*_num_threads+thread_num] +=
					material->getTransportMacroXS(energy_index) / sigma_t;
		else
			_tally_acc[bin_index*_num_threads+thread_num] +=
			material->getIsotope(_isotopes)->getTransportXS(energy_index) /
																sigma_t;
		omp_unset_lock(&_lock);
	}

	return;
}


/******************************************************************************
 ***************************  CollisionRateBinner  ****************************
 *****************************************************************************/

CollisionRateBinner::CollisionRateBinner() {
	setBinnerType(COLLISION_RATE);
};

CollisionRateBinner::~CollisionRateBinner() { };

void CollisionRateBinner::weightedTally(neutron* neutron, float sigma_t,
				int energy_index, Material* material, Isotope* isotope){

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);
	int thread_num = neutron->_thread_num;

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
			_tally_acc[bin_index*_num_threads+thread_num] +=
													neutron->_weight;
		omp_unset_lock(&_lock);
	}

	return;
}



/******************************************************************************
 ***************************  DiffusionRateBinner  ****************************
 *****************************************************************************/

DiffusionRateBinner::DiffusionRateBinner() {
	setBinnerType(DIFFUSION_RATE);
};

DiffusionRateBinner::~DiffusionRateBinner() { };

void DiffusionRateBinner::weightedTally(neutron* neutron, float sigma_t,
							int energy_index, Material* material, Isotope* isotope) {

	if (_num_bins == 0)
		 log_printf(ERROR, "Cannot tally weighted sample in Binner %s since "
				 "the bins have not yet been created", _name);

	int bin_index = getBinIndex(neutron);
	int thread_num = neutron->_thread_num;

	if (bin_index >= 0 && bin_index < _num_bins) {
		omp_set_lock(&_lock);
		if (strcmp(_isotopes, "all")) {
			_tally_acc[bin_index*_num_threads+thread_num] +=
				neutron->_weight * 1.0 / (3.0 *
					material->getTransportMacroXS(energy_index) * sigma_t);
		}
		else {
			Isotope* isotope = material->getIsotope(_isotopes);
			float num_density = material->getIsotopeNumDensity(_isotopes);
			float isotope_sigma_tr =
			isotope->getTransportXS(energy_index) * num_density * 1E-24;
			_tally_acc[bin_index*_num_threads+thread_num] +=
			neutron->_weight * 1.0 / (3.0 * isotope_sigma_tr * sigma_t);
		}
		omp_unset_lock(&_lock);
	}

	return;
}
