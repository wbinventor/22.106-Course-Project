/*
 * BatchBinSet.cpp
 *
 *  Created on: Mar 17, 2012
 *      Author: William Boyd
 *				MIT, Course 22
 *              wboyd@mit.edu
 */

#include "BatchBinSet.h"


/**
 * Default constructor for a BatchBinSet
 */
BatchBinSet::BatchBinSet() {
	/* Default for empty batch bin set */
	_num_batches = 0;
	_statistics_compute = false;
}


/**
 * BatchBinSet destructor deletes memory for arrays of binners,
 * and batch statistics if that memory has been allocated
 */
BatchBinSet::~BatchBinSet() {

	if (_num_batches != 0) {
		delete [] _binners;
		delete [] _batch_mu;
		delete [] _batch_variance;
		delete [] _batch_std_dev;
		delete [] _batch_rel_err;
	}
}


/**
 * Returns a Binner class pointer for a particular batch
 * @param batch the batch number
 * @return a pointer to a Binner class
 */
Binner* BatchBinSet::getBinner(int batch) {

	if (_num_batches == 0)
		log_printf(ERROR, "Unable to return binner %d since the binners"
				" for this BatchBinSet have not yet been created", batch);

	else if (batch < 0 || batch >= _num_batches)
		log_printf(ERROR, "Unable to return binner %d since this BatchBinSet"
				" has %d batches", batch, _num_batches);

	return &_binners[batch];
}


/**
 * Returns a pointer to an array of batch averages if they have been
 * computed
 * @return a float array of batch averages for each bin
 */
float* BatchBinSet::getBatchMu() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for this "
				"BatchBinSet so batch mu cannot be returned");

	return _batch_mu;
}


/**
 * Returns a pointer to an array of batch variances if they have been
 * computed
 * @return a float array of batch variances for each bin
 */
float* BatchBinSet::getBatchVariance() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for this "
				"BatchBinSet so batch variance cannot be returned");

	return _batch_variance;
}


/**
 * Returns a pointer to an array of batch standard deviations if they have
 * been computed
 * @return a float array of batch standard deviations for each bin
 */
float* BatchBinSet::getBatchStdDev() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for this "
				"BatchBinSet so batch std dev cannot be returned");

	return _batch_std_dev;
}


/**
 * Returns a pointer to an array of batch relative errors if they have been
 * computed
 * @return a float array of batch relative errors for each bin
 */
float* BatchBinSet::getBatchRelativeError() {

	if (!_statistics_compute)
		log_printf(ERROR, "Statistics have not yet been computed for this "
				"BatchBinSet so batch relative error cannot be returned");

	return _batch_rel_err;
}


/**
 * Creates a set of binners each with the same set of bins between a
 * certain start and end point
 * @param start starting bin value
 * @param end ending bin value
 * @param num_bins the number of bins for each Binner
 * @param num_batches the number of batches (number of Binner classes)
 * @param bin_type the type of bin (EQUAL, LOGARITHMIC, OTHER)
 * @param tally_type the type of tally (FLUX, ABSORB_RATE, COLLISION_RATE)
 */
void BatchBinSet::createBinners(float start, float end, int num_bins,
			int num_batches, binnerType binner_type, spacingType bin_type,
												tallyDomainType tally_type) {

	_num_batches = num_batches;
	_num_bins = num_bins;

	/* Allocate memory for binners and batch statistics arrays */
	switch(binner_type) {
	case FLUX:
		_binners = new FluxBinner[_num_batches];
		break;
	case CAPTURE_RATE:
		_binners = new CaptureRateBinner[_num_batches];
		break;
	case ABSORPTION_RATE:
		_binners = new AbsorptionRateBinner[_num_batches];
		break;
	case ELASTIC_RATE:
		_binners = new ElasticRateBinner[_num_batches];
		break;
	case INELASTIC_RATE:
		_binners = new InelasticRateBinner[_num_batches];
		break;
	case FISSION_RATE:
		_binners = new FissionRateBinner[_num_batches];
		break;
	case TRANSPORT_RATE:
		_binners = new TransportRateBinner[_num_batches];
		break;
	case COLLISION_RATE:
		_binners = new CollisionRateBinner[_num_batches];
		break;
	case DIFFUSION_RATE:
		_binners = new DiffusionRateBinner[_num_batches];
		break;
	}

	_batch_mu = new float[_num_bins];
	_batch_variance = new float[_num_bins];
	_batch_std_dev = new float[_num_bins];
	_batch_rel_err = new float[_num_bins];

	/* Create a Binner for each batch */
	for (int b=0; b < num_batches; b++) {
		_binners[b].generateBinEdges(start, end, _num_bins, bin_type);
		_binners[b].setTallyDomainType(tally_type);
	}

	return;
}


/**
 * Creates a set of binners each with the same set of bins defined by
 * an array of bin edges
 * @param bin_edges array of edges which define the bins
 * @param num_bins the number of bins for each Binner
 * @param num_batches the number of batches (number of Binner classes)
 * @param tally_type the type of tally (FLUX, ABSORB_RATE, COLLISION_RATE)
 */
void BatchBinSet::createBinners(float* bin_edges, int num_bins,
		int num_batches, binnerType binner_type, tallyDomainType tally_type) {

	_num_batches = num_batches;
	_num_bins = num_bins;

	/* Allocate memory for binners and batch statistics arrays */
	switch(binner_type) {
	case FLUX:
		_binners = new FluxBinner[_num_batches];
		break;
	case CAPTURE_RATE:
		_binners = new CaptureRateBinner[_num_batches];
		break;
	case ABSORPTION_RATE:
		_binners = new AbsorptionRateBinner[_num_batches];
		break;
	case ELASTIC_RATE:
		_binners = new ElasticRateBinner[_num_batches];
		break;
	case INELASTIC_RATE:
		_binners = new InelasticRateBinner[_num_batches];
		break;
	case FISSION_RATE:
		_binners = new FissionRateBinner[_num_batches];
		break;
	case TRANSPORT_RATE:
		_binners = new TransportRateBinner[_num_batches];
		break;
	case COLLISION_RATE:
		_binners = new CollisionRateBinner[_num_batches];
		break;
	case DIFFUSION_RATE:
		_binners = new DiffusionRateBinner[_num_batches];
		break;
	}

	_batch_mu = new float[_num_bins];
	_batch_variance = new float[_num_bins];
	_batch_std_dev = new float[_num_bins];
	_batch_rel_err = new float[_num_bins];

	/* Create a Binner for each batch */
	for (int b=0; b < num_batches; b++) {
		_binners[b].setBinEdges(bin_edges, _num_bins);
		_binners[b].setTallyDomainType(tally_type);
	}

	return;
}


/**
 * Computes average, variance, standard deviation and relative error for each
 * bin over the set of batches
 */
void BatchBinSet::computeBatchStatistics() {

	log_printf(NORMAL, "Computing batch statistics...");

	if (_num_batches == 0)
		log_printf(ERROR, "Cannot compute batch statistics since the binners"
				" for this BatchBinSet have not yet been generated");

	double s1 = 0.0;
	double s2 = 0.0;

	/* Loop over each bin */
	for (int i=0; i < _num_bins; i++) {

		/* Initialize statistics to zero */
		_batch_mu[i] = 0.0;
		_batch_variance[i] = 0.0;
		_batch_std_dev[i] = 0.0;
		_batch_rel_err[i] = 0.0;

		/* Accumulate in s1, s2 counters */
		for (int j=0; j < _num_batches; j++) {
			s1 += _binners[j].getTally(i);
			s2 += _binners[j].getTally(i) * _binners[j].getTally(i);
			log_printf(NORMAL, "j = %d, s1 = %f, s2 = %f", j, s1, s2);
		}

		/* Compute batch average */
		_batch_mu[i] = s1 / _num_batches;

		log_printf(NORMAL, "Batch mu = %f", _batch_mu[i]);

		/* Compute batch variance */
		_batch_variance[i] = (1.0 / (float(_num_batches) - 1.0)) *
				(s2 / float(_num_batches) - (_batch_mu[i]*_batch_mu[i]));

		log_printf(NORMAL, "_batch_variance[i] = %f", _batch_variance[i]);

		_batch_std_dev[i] = sqrt(_batch_variance[i]);
		_batch_rel_err[i] = _batch_std_dev[i] / _batch_mu[i];

//		/* Accumulate flux from each batch */
//		for (int j=0; j < _num_batches; j++)
//			_batch_mu[i] += _binners[j].getTally(i);
//
//		/* Compute average flux for this bin */
//		_batch_mu[i] /= float(_num_batches);
//
//		/* Compute the variance for this bin */
//		for (int j=0; j < _num_batches; j++) {
//			_batch_variance[i] += (_binners[j].getTally(i) - _batch_mu[i])
//					* (_binners[j].getTally(i) - _batch_mu[i]);
//		}
//		_batch_variance[i] /= float(_num_batches);
//
//		/* Compute the standard deviation for this bin */
//		_batch_std_dev[i] = sqrt(_batch_variance[i]);
//
//		/* Compute the relative error for this bin */
//		_batch_rel_err[i] = _batch_std_dev[i] / _batch_mu[i];
	}

	_statistics_compute = true;

	return;
}


/**
 * Computes average, variance, standard deviation and relative error for each
 * bin over the set of batches. This method first scales each bin value by
 * a scaling factor
 * @param scale_factor the factor to scale each bin value by
 */
void BatchBinSet::computeScaledBatchStatistics(float scale_factor) {

	if (_num_batches == 0)
		log_printf(ERROR, "Cannot compute batch statistics since the binners"
				" for this BatchBinSet have not yet been generated");

	double s1 = 0.0;
	double s2 = 0.0;

	/* Loop over each bin */
	for (int i=0; i < _num_bins; i++) {

		/* Initialize statistics to zero */
		_batch_mu[i] = 0.0;
		_batch_variance[i] = 0.0;
		_batch_std_dev[i] = 0.0;
		_batch_rel_err[i] = 0.0;

		/* Accumulate in s1, s2 counters */
		for (int j=0; j < _num_batches; j++) {
			s1 += _binners[j].getTally(i) / scale_factor;
			s2 += (_binners[j].getTally(i) / scale_factor) *
						(_binners[j].getTally(i) / scale_factor);
		}

		/* Compute batch average */
		_batch_mu[i] = s1 / _num_batches;

		/* Compute batch variance */
		_batch_variance[i] = (1.0 / (float(_num_batches) - 1.0)) *
				(s2 / float(_num_batches) - (_batch_mu[i]*_batch_mu[i]));

		_batch_std_dev[i] = sqrt(_batch_variance[i]);
		_batch_rel_err[i] = _batch_std_dev[i] / _batch_mu[i];

//		/* Accumulate flux from each batch */
//		for (int j=0; j < _num_batches; j++)
//			_batch_mu[i] += _binners[j].getTally(i) / scale_factor;
//
//		/* Compute average flux for this bin */
//		_batch_mu[i] /= float(_num_batches);
//
//		/* Compute the variance for this bin */
//		for (int j=0; j < _num_batches; j++) {
//			_batch_variance[i] += (_binners[j].getTally(i) / scale_factor
//			- _batch_mu[i]) * (_binners[j].getTally(i) / scale_factor
//												- _batch_mu[i]);
//		}
//		_batch_variance[i] /= float(_num_batches);
//
//		/* Compute the standard deviation for this bin */
//		_batch_std_dev[i] = sqrt(_batch_variance[i]);
//
//		/* Compute the relative error for this bin */
//		_batch_rel_err[i] = _batch_std_dev[i] / _batch_mu[i];
	}

	_statistics_compute = true;

	return;
}


/**
 * Outputs the batch statistics (if they have been computed) to an
 * ASCII file
 * @param filename the output filename
 */
void BatchBinSet::outputBatchStatistics(const char* filename) {

	if (_num_batches == 0)
		log_printf(ERROR, "Cannot output batch statistics since the binners"
				" for this BatchBinSet have not yet been generated");

	else if (!_statistics_compute)
		log_printf(ERROR, "Cannot output batch statistics since they have"
				" not yet been computed for this BatchBinSet");

	/* Create output file */
	FILE* output_file;
	output_file = fopen(filename, "w");

	/* Print header to output file */
	fprintf(output_file, "Bin center, Mu, Variance, Std Dev, Rel Err\n");

	/* Loop over each bin and print mu, var, std dev and rel err */
	for (int i=0; i < _num_bins; i++) {
		fprintf(output_file, "%1.10f, %1.10f, %1.10f, %1.10f, %1.10f\n",
				_binners[0].getBinCenters()[i], _batch_mu[i],
				_batch_variance[i], _batch_std_dev[i], _batch_rel_err[i]);
	}

	fclose(output_file);

	return;
}
