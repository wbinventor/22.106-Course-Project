/*
 * BaysoySubasi.cpp
 *
 *  Created on: Apr 26, 2012
 *      Author: wboyd
 */


#include "BaysoySubasi.h"

void BaysoySubasi(Options* options) {

	int num_neutrons = options->getNumNeutrons();
	int num_batches = options->getNumBatches();

	BatchBinSet* detector_coll_rate = new BatchBinSet();
	BatchBinSet* tot_coll_rate_x = new BatchBinSet();
	BatchBinSet* tot_coll_rate_y = new BatchBinSet();
	BatchBinSet* tot_coll_rate_z = new BatchBinSet();
	BatchBinSet* tot_coll_rate_rz = new BatchBinSet();

	detector_coll_rate->createBinners(1E-5, 1E7, 1,
			num_batches, COLLISION_RATE, LINEAR, ENERGY);
	tot_coll_rate_x->createBinners(-25.0, 25.0, 200, num_batches,
			COLLISION_RATE, LINEAR, X);
	tot_coll_rate_y->createBinners(-25.0, 25.0, 200, num_batches,
			COLLISION_RATE, LINEAR, Y);
	tot_coll_rate_z->createBinners(-25.0, 25.0, 200, num_batches,
			COLLISION_RATE, LINEAR, Z);
	tot_coll_rate_rz->createBinners(0, 25.0, 200, num_batches,
			COLLISION_RATE, LINEAR, R_Z);

	double N_A = 6.023E23;		/* Avogadro's number */

	/* Dry air density at 20 C and 1 atm */
	float rho_air = 1.2041E-3;
	float rho_N14_air = rho_air * 0.781 * N_A / 14.0;
	float rho_O16_air = rho_air * 0.21 * N_A / 16.0;
	float rho_Ar40_air = rho_air * 0.09 * N_A / 40.0;

	/* He3 detector material densities from Baysoy and Subasi */
	float rho_detector_gas = 7.3E-4;
	float rho_He3_detector_gas = rho_detector_gas * 1.0 * N_A / 3.0;

	/* TNT density from Esheikh, Habbani, ElAgib */
	float rho_tnt = 1.8;
	float rho_H1_tnt = rho_tnt * 0.022 * N_A / 1.0;
	float rho_C12_tnt = rho_tnt * 0.370 * N_A / 12.0;
	float rho_N14_tnt = rho_tnt * 0.185 * N_A / 14.0;
	float rho_O16_tnt = rho_tnt * 0.423 * N_A / 16.0;

	/* Soil isotopic from paper by Elsheikh, Habbani, ElAgib,
	 * 2011 (g/cm^3) with natural isotope fractions taken from Wikipedia */

	/* Dry porous soil */
	float rho_soil = 1.1810;
	float rho_H1_soil = rho_soil * 0.015 * N_A / 1.0;
	float rho_O16_soil = rho_soil * 0.529 * N_A / 16.0;
	float rho_Si28_soil = rho_soil * 0.243 * 0.92223 * N_A / 28.0;
	float rho_Si29_soil = rho_soil * 0.243 * 0.04685 * N_A / 29.0;
	float rho_Si30_soil = rho_soil * 0.243 * 0.03092 * N_A / 30.0;
	float rho_Al27_soil = rho_soil * 0.071 * N_A / 27;
	float rho_Fe54_soil = rho_soil * 0.044 * 0.05845 * N_A / 54.0;
	float rho_Fe56_soil = rho_soil * 0.044 * 0.91754 * N_A / 56.0;
	float rho_Fe57_soil = rho_soil * 0.044 * 0.02119 * N_A / 57.0;
	float rho_Fe58_soil = rho_soil * 0.044 * 0.00282 * N_A / 58.0;
	float rho_Ca40_soil = rho_soil * 0.032 * 0.96941 * N_A / 40.0;
	float rho_Ca42_soil = rho_soil * 0.032 * 0.00647 * N_A / 42.0;
	float rho_Ca44_soil = rho_soil * 0.032 * 0.02086 * N_A / 44.0;
	float rho_K39_soil = rho_soil * 0.023 * 0.932581 * N_A / 39.0;
	float rho_K41_soil = rho_soil * 0.023 * 0.067302 * N_A / 41.0;
	float rho_Na23_soil = rho_soil * 0.025 * N_A / 23.0;
	float rho_Mg24_soil = rho_soil * 0.018 * 0.78958 * N_A / 24.0;
	float rho_Mg25_soil = rho_soil * 0.018 * 0.0999986 * N_A / 25.0;
	float rho_Mg26_soil =  rho_soil * 0.018 * 0.10987 * N_A / 26.0;

	/* Dry dense soil */
	//float rho_soil = 1.7714;
	//float rho_H1_soil = rho_soil * 0.015 * N_A / 16.0;
	//float rho_O16_soil = rho_soil * 0.529 * N_A / 1.0;
	//float rho_Si28_soil = rho_soil * 0.243 * 0.92223 * N_A / 28.0;
	//float rho_Si29_soil = rho_soil * 0.243 * 0.04685 * N_A / 29.0;
	//float rho_Si30_soil = rho_soil * 0.243 * 0.03092 * N_A / 30.0;
	//float rho_Al27_soil = rho_soil * 0.071 * N_A / 27;
	//float rho_Fe54_soil = rho_soil * 0.044 * 0.05845 * N_A / 54.0;
	//float rho_Fe56_soil = rho_soil * 0.044 * 0.91754 * N_A / 56.0;
	//float rho_Fe57_soil = rho_soil * 0.044 * 0.02119 * N_A / 57.0;
	//float rho_Fe58_soil = rho_soil * 0.044 * 0.00282 * N_A / 58.0;
	//float rho_Ca40_soil = rho_soil * 0.032 * 0.96941 * N_A / 40.0;
	//float rho_Ca42_soil = rho_soil * 0.032 * 0.00647 * N_A / 42.0;
	//float rho_Ca44_soil = rho_soil * 0.032 * 0.02086 * N_A / 44.0;
	//float rho_K39_soil = rho_soil * 0.023 * 0.932581 * N_A / 39.0;
	//float rho_K41_soil = rho_soil * 0.023 * 0.067302 * N_A / 41.0;
	//float rho_Na23_soil = rho_soil * 0.025 * N_A / 23.0;
	//float rho_Mg24_soil = rho_soil * 0.018 * 0.78958 * N_A / 24.0;
	//float rho_Mg25_soil = rho_soil * 0.018 * 0.0999986 * N_A / 25.0;
	//float rho_Mg26_soil =  rho_soil * 0.018 * 0.10987 * N_A / 26.0;
	//
	/* Wet porous soil */
	//float rho_soil = 1.3957;
	//float rho_H1_soil = rho_soil * 0.030 * N_A / 1.0;
	//float rho_O16_soil = rho_soil * 0.585 * N_A / 16.0;
	//float rho_Si28_soil = rho_soil * 0.205 * 0.92223 * N_A / 28.0;
	//float rho_Si29_soil = rho_soil * 0.205 * 0.04685 * N_A / 29.0;
	//float rho_Si30_soil = rho_soil * 0.205 * 0.03092 * N_A / 30.0;
	//float rho_Al27_soil = rho_soil * 0.060 * N_A / 27;
	//float rho_Fe54_soil = rho_soil * 0.037 * 0.05845 * N_A / 54.0;
	//float rho_Fe56_soil = rho_soil * 0.037 * 0.91754 * N_A / 56.0;
	//float rho_Fe57_soil = rho_soil * 0.037 * 0.02119 * N_A / 57.0;
	//float rho_Fe58_soil = rho_soil * 0.037 * 0.00282 * N_A / 58.0;
	//float rho_Ca40_soil = rho_soil * 0.027 * 0.96941 * N_A / 40.0;
	//float rho_Ca42_soil = rho_soil * 0.027 * 0.00647 * N_A / 42.0;
	//float rho_Ca44_soil = rho_soil * 0.027 * 0.02086 * N_A / 44.0;
	//float rho_K39_soil = rho_soil * 0.019 * 0.932581 * N_A / 39.0;
	//float rho_K41_soil = rho_soil * 0.019 * 0.067302 * N_A / 41.0;
	//float rho_Na23_soil = rho_soil * 0.021 * N_A / 23.0;
	//float rho_Mg24_soil = rho_soil * 0.015 * 0.78958 * N_A / 24.0;
	//float rho_Mg25_soil = rho_soil * 0.015 * 0.0999986 * N_A / 25.0;
	//float rho_Mg26_soil =  rho_soil * 0.015 * 0.10987 * N_A / 26.0;
	//
	/* Wet dense soil */
	//float rho_soil = 2.0935;
	//float rho_H1_soil = rho_soil * 0.030 * N_A / 1.0;
	//float rho_O16_soil = rho_soil * 0.585 * N_A / 16.0;
	//float rho_Si28_soil = rho_soil * 0.205 * 0.92223 * N_A / 28.0;
	//float rho_Si29_soil = rho_soil * 0.205 * 0.04685 * N_A / 29.0;
	//float rho_Si30_soil = rho_soil * 0.205 * 0.03092 * N_A / 30.0;
	//float rho_Al27_soil = rho_soil * 0.060 * N_A / 27;
	//float rho_Fe54_soil = rho_soil * 0.037 * 0.05845 * N_A / 54.0;
	//float rho_Fe56_soil = rho_soil * 0.037 * 0.91754 * N_A / 56.0;
	//float rho_Fe57_soil = rho_soil * 0.037 * 0.02119 * N_A / 57.0;
	//float rho_Fe58_soil = rho_soil * 0.037 * 0.00282 * N_A / 58.0;
	//float rho_Ca40_soil = rho_soil * 0.027 * 0.96941 * N_A / 40.0;
	//float rho_Ca42_soil = rho_soil * 0.027 * 0.00647 * N_A / 42.0;
	//float rho_Ca44_soil = rho_soil * 0.027 * 0.02086 * N_A / 44.0;
	//float rho_K39_soil = rho_soil * 0.019 * 0.932581 * N_A / 39.0;
	//float rho_K41_soil = rho_soil * 0.019 * 0.067302 * N_A / 41.0;
	//float rho_Na23_soil = rho_soil * 0.021 * N_A / 23.0;
	//float rho_Mg24_soil = rho_soil * 0.015 * 0.78958 * N_A / 24.0;
	//float rho_Mg25_soil = rho_soil * 0.015 * 0.0999986 * N_A / 25.0;
	//float rho_Mg26_soil =  rho_soil * 0.015 * 0.10987 * N_A / 26.0;
	//

	/* Declare and initialize all isotopes in problem */
	Isotope Al27, Ar40, B10, B11, C12, Ca40, Ca42, Ca44, Fe54, Fe56, Fe57, Fe58;
	Isotope H1, He3, K39, K41, Mg24, Mg25, Mg26, N14, Na23, O16, Si28, Si29, Si30;
	Material soil, air, detector_gas, tnt;
	ZCylinder mine, atmosphere, dirt, detector;

	/* Aluminum 27 */
	Al27.setA(27);
	Al27.setIsotopeType("Al27");
	Al27.loadXS("../pendf/Al27-capture.txt", CAPTURE, "\t");
	Al27.loadXS("../pendf/Al27-elastic.txt", ELASTIC, "\t");
	Al27.setElasticAngleType(ISOTROPIC_CM);

	/* Argon 40 */
	Ar40.setA(40);
	Ar40.setIsotopeType("Ar40");
	Ar40.loadXS("../pendf/Ar40-capture.txt", CAPTURE, "\t");
	Ar40.loadXS("../pendf/Ar40-elastic.txt", ELASTIC, "\t");
	Ar40.setElasticAngleType(ISOTROPIC_CM);

	/* Boron 10 */
	B10.setA(10);
	B10.setIsotopeType("B10");
	B10.loadXS("../pendf/B10-capture.txt", CAPTURE, "\t");
	B10.loadXS("../pendf/B10-elastic.txt", ELASTIC, "\t");
	B10.setElasticAngleType(ISOTROPIC_CM);

	/* Boron 11 */
	B11.setA(11);
	B11.setIsotopeType("B11");
	B11.loadXS("../pendf/B11-capture.txt", CAPTURE, "\t");
	B11.loadXS("../pendf/B11-elastic.txt", ELASTIC, "\t");
	B11.setElasticAngleType(ISOTROPIC_CM);

	/* Carbon 12 */
	C12.setA(12);
	C12.setIsotopeType("C12");
	C12.loadXS("../pendf/C12-capture.txt", CAPTURE, "\t");
	C12.loadXS("../pendf/C12-elastic.txt", ELASTIC, "\t");
	C12.setElasticAngleType(ISOTROPIC_CM);

	/* Calcium 40 */
	Ca40.setA(40);
	Ca40.setIsotopeType("Ca40");
	Ca40.loadXS("../pendf/Ca40-capture.txt", CAPTURE, "\t");
	Ca40.loadXS("../pendf/Ca40-elastic.txt", ELASTIC, "\t");
	Ca40.setElasticAngleType(ISOTROPIC_CM);

	/* Calcium 42 */
	Ca42.setA(42);
	Ca42.setIsotopeType("Ca42");
	Ca42.loadXS("../pendf/Ca42-capture.txt", CAPTURE, "\t");
	Ca42.loadXS("../pendf/Ca42-elastic.txt", ELASTIC, "\t");
	Ca42.setElasticAngleType(ISOTROPIC_CM);

	/* Calcium 44 */
	Ca44.setA(44);
	Ca44.setIsotopeType("Ca44");
	Ca44.loadXS("../pendf/Ca44-capture.txt", CAPTURE, "\t");
	Ca44.loadXS("../pendf/Ca44-elastic.txt", ELASTIC, "\t");
	Ca44.setElasticAngleType(ISOTROPIC_CM);

	/* Iron 54 */
	Fe54.setA(54);
	Fe54.setIsotopeType("Fe54");
	Fe54.loadXS("../pendf/Fe54-capture.txt", CAPTURE, "\t");
	Fe54.loadXS("../pendf/Fe54-elastic.txt", ELASTIC, "\t");
	Fe54.setElasticAngleType(ISOTROPIC_CM);

	/* Iron 56 */
	Fe56.setA(56);
	Fe56.setIsotopeType("Fe56");
	Fe56.loadXS("../pendf/Fe56-capture.txt", CAPTURE, "\t");
	Fe56.loadXS("../pendf/Fe56-elastic.txt", ELASTIC, "\t");
	Fe56.setElasticAngleType(ISOTROPIC_CM);

	/* Iron 57 */
	Fe57.setA(57);
	Fe57.setIsotopeType("Fe57");
	Fe57.loadXS("../pendf/Fe57-capture.txt", CAPTURE, "\t");
	Fe57.loadXS("../pendf/Fe57-elastic.txt", ELASTIC, "\t");
	Fe57.setElasticAngleType(ISOTROPIC_CM);

	/* Iron 58 */
	Fe58.setA(58);
	Fe58.setIsotopeType("Fe58");
	Fe58.loadXS("../pendf/Fe58-capture.txt", CAPTURE, "\t");
	Fe58.loadXS("../pendf/Fe58-elastic.txt", ELASTIC, "\t");
	Fe58.setElasticAngleType(ISOTROPIC_CM);

	/* Hydrogen 1 */
	H1.setA(1);
	H1.setIsotopeType("H1");
	H1.loadXS("../pendf/H1-capture.txt", CAPTURE, "\t");
	H1.loadXS("../pendf/H1-elastic.txt", ELASTIC, "\t");
	H1.setElasticAngleType(ISOTROPIC_CM);

	/* Helium 3 */
	He3.setA(3);
	He3.setIsotopeType("He3");
	He3.loadXS("../pendf/He3-reprocessed_capture.txt", CAPTURE, "\t");
	He3.loadXS("../pendf/He3-elastic.txt", ELASTIC, "\t");
	He3.setElasticAngleType(ISOTROPIC_CM);
	He3.plotXS(1E-5, 1E7, 1000, CAPTURE, ELASTIC, NULL);

	/* Potassium 39 */
	K39.setA(39);
	K39.setIsotopeType("K39");
	K39.loadXS("../pendf/K39-capture.txt", CAPTURE, "\t");
	K39.loadXS("../pendf/K39-elastic.txt", ELASTIC, "\t");
	K39.setElasticAngleType(ISOTROPIC_CM);

	/* Potassium 41 */
	K41.setA(41);
	K41.setIsotopeType("K41");
	K41.loadXS("../pendf/K41-capture.txt", CAPTURE, "\t");
	K41.loadXS("../pendf/K41-elastic.txt", ELASTIC, "\t");
	K41.setElasticAngleType(ISOTROPIC_CM);

	/* Magnesium 24 */
	Mg24.setA(24);
	Mg24.setIsotopeType("Mg24");
	Mg24.loadXS("../pendf/Mg24-capture.txt", CAPTURE, "\t");
	Mg24.loadXS("../pendf/Mg24-elastic.txt", ELASTIC, "\t");
	Mg24.setElasticAngleType(ISOTROPIC_CM);

	/* Magnesium 25 */
	Mg25.setA(25);
	Mg25.setIsotopeType("Mg25");
	Mg25.loadXS("../pendf/Mg25-capture.txt", CAPTURE, "\t");
	Mg25.loadXS("../pendf/Mg25-elastic.txt", ELASTIC, "\t");
	Mg25.setElasticAngleType(ISOTROPIC_CM);

	/* Magnesium 26 */
	Mg26.setA(26);
	Mg26.setIsotopeType("Mg26");
	Mg26.loadXS("../pendf/Mg26-capture.txt", CAPTURE, "\t");
	Mg26.loadXS("../pendf/Mg26-elastic.txt", ELASTIC, "\t");
	Mg26.setElasticAngleType(ISOTROPIC_CM);

	/* Nitrogen 14 */
	N14.setA(14);
	N14.setIsotopeType("N14");
	N14.loadXS("../pendf/N14-capture.txt", CAPTURE, "\t");
	N14.loadXS("../pendf/N14-elastic.txt", ELASTIC, "\t");
	N14.setElasticAngleType(ISOTROPIC_CM);

	/* Sodium 23 */
	Na23.setA(23);
	Na23.setIsotopeType("Na23");
	Na23.loadXS("../pendf/Na23-capture.txt", CAPTURE, "\t");
	Na23.loadXS("../pendf/Na23-elastic.txt", ELASTIC, "\t");
	Na23.setElasticAngleType(ISOTROPIC_CM);

	/* Oxygen 16 */
	O16.setA(16);
	O16.setIsotopeType("O16");
	O16.loadXS("../pendf/O16-capture.txt", CAPTURE, "\t");
	O16.loadXS("../pendf/O16-elastic.txt", ELASTIC, "\t");
	O16.setElasticAngleType(ISOTROPIC_CM);

	/* Silicon 28 */
	Si28.setA(28);
	Si28.setIsotopeType("Si28");
	Si28.loadXS("../pendf/Si28-capture.txt", CAPTURE, "\t");
	Si28.loadXS("../pendf/Si28-elastic.txt", ELASTIC, "\t");
	Si28.setElasticAngleType(ISOTROPIC_CM);

	/* Silicon 29 */
	Si29.setA(29);
	Si29.setIsotopeType("Si29");
	Si29.loadXS("../pendf/Si29-capture.txt", CAPTURE, "\t");
	Si29.loadXS("../pendf/Si29-elastic.txt", ELASTIC, "\t");
	Si29.setElasticAngleType(ISOTROPIC_CM);

	/* Silicon 30 */
	Si30.setA(30);
	Si30.setIsotopeType("Si30");
	Si30.loadXS("../pendf/Si30-capture.txt", CAPTURE, "\t");
	Si30.loadXS("../pendf/Si30-elastic.txt", ELASTIC, "\t");
	Si30.setElasticAngleType(ISOTROPIC_CM);

	/* Load materials */
	soil.setMaterialName("soil");
	soil.addIsotope(&H1, rho_H1_soil);
	soil.addIsotope(&O16, rho_O16_soil);
	soil.addIsotope(&Si28, rho_Si28_soil);
	soil.addIsotope(&Si29, rho_Si29_soil);
	soil.addIsotope(&Si30, rho_Si30_soil);
	soil.addIsotope(&Al27, rho_Al27_soil);
	soil.addIsotope(&Fe54, rho_Fe54_soil);
	soil.addIsotope(&Fe56, rho_Fe56_soil);
	soil.addIsotope(&Fe57, rho_Fe57_soil);
	soil.addIsotope(&Fe58, rho_Fe58_soil);
	soil.addIsotope(&Ca40, rho_Ca40_soil);
	soil.addIsotope(&Ca42, rho_Ca42_soil);
	soil.addIsotope(&Ca44, rho_Ca44_soil);
	soil.addIsotope(&K39, rho_K39_soil);
	soil.addIsotope(&K41, rho_K41_soil);
	soil.addIsotope(&Na23, rho_Na23_soil);
	soil.addIsotope(&Mg24, rho_Mg24_soil);
	soil.addIsotope(&Mg25, rho_Mg25_soil);
	soil.addIsotope(&Mg26, rho_Mg26_soil);
	soil.rescaleCrossSections(1E-5, 2E6, 50000);

	air.setMaterialName("air");
	air.addIsotope(&N14, rho_N14_air);
	air.addIsotope(&O16, rho_O16_air);
	air.addIsotope(&Ar40, rho_Ar40_air);
	air.rescaleCrossSections(1E-5, 2E6, 50000);

	detector_gas.setMaterialName("detector gas");
	detector_gas.addIsotope(&He3, rho_He3_detector_gas);
	detector_gas.rescaleCrossSections(1E-5, 2E6, 50000);

	tnt.setMaterialName("tnt");
	tnt.addIsotope(&H1, rho_H1_tnt);
	tnt.addIsotope(&C12, rho_C12_tnt);
	tnt.addIsotope(&N14, rho_N14_tnt);
	tnt.addIsotope(&O16, rho_O16_tnt);
	tnt.rescaleCrossSections(1E-5, 2E6, 50000);

	/* Setup atmosphere cylinder */
	ZCircle* atmosphere_top = new ZCircle();
	atmosphere_top->setBoundaryType(VACUUM);
	atmosphere_top->setRadius(25.0);
	atmosphere_top->setX0(0.0);
	atmosphere_top->setY0(0.0);
	atmosphere_top->setZ0(75.0);

	ZCircle* atmosphere_dirt_interface = new ZCircle();
	atmosphere_dirt_interface->setBoundaryType(INTERFACE);
	atmosphere_dirt_interface->setRadius(25.0);
	atmosphere_dirt_interface->setX0(0.0);
	atmosphere_dirt_interface->setY0(0.0);
	atmosphere_dirt_interface->setZ0(0.0);

	OpenZCylinder* atmosphere_cylinder = new OpenZCylinder();
	atmosphere_cylinder->setBoundaryType(VACUUM);
	atmosphere_cylinder->setX0(0.0);
	atmosphere_cylinder->setY0(0.0);
	atmosphere_cylinder->setZ0(0.0);
	atmosphere_cylinder->setRadius(25.0);
	atmosphere_cylinder->setZLeft(0.0);
	atmosphere_cylinder->setZRight(75.0);

	atmosphere.setRegionName((char*)"atmosphere");
	atmosphere.setMaterial(&air);
	atmosphere.setOpenZCylinder(atmosphere_cylinder);
	atmosphere.setZLeftCircle(atmosphere_dirt_interface);
	atmosphere.setZRightCircle(atmosphere_top);

	/* Set up dirt cylinder */
	ZCircle* dirt_bottom = new ZCircle();
	dirt_bottom->setBoundaryType(VACUUM);
	dirt_bottom->setRadius(25.0);
	dirt_bottom->setX0(0.0);
	dirt_bottom->setY0(0.0);
	dirt_bottom->setZ0(-50.0);

	OpenZCylinder* dirt_cylinder = new OpenZCylinder();
	dirt_cylinder->setBoundaryType(VACUUM);
	dirt_cylinder->setX0(0.0);
	dirt_cylinder->setY0(0.0);
	dirt_cylinder->setZ0(0.0);
	dirt_cylinder->setRadius(25.0);
	dirt_cylinder->setZLeft(-50.0);
	dirt_cylinder->setZRight(0.0);

	dirt.setRegionName((char*)"dirt");
	dirt.setMaterial(&soil);
	dirt.setOpenZCylinder(dirt_cylinder);
	dirt.setZLeftCircle(dirt_bottom);
	dirt.setZRightCircle(atmosphere_dirt_interface);

	atmosphere_top->setLeftRegion(&atmosphere);
	atmosphere_dirt_interface->setLeftRegion(&dirt);
	atmosphere_dirt_interface->setRightRegion(&atmosphere);
	dirt_bottom->setRightRegion(&dirt);
	atmosphere_cylinder->setLeftRegion(&atmosphere);
	dirt_cylinder->setLeftRegion(&dirt);

	/* Set up mine cylinder */
	ZCircle* mine_top = new ZCircle();
	mine_top->setBoundaryType(INTERFACE);
	mine_top->setRadius(3.5);
	mine_top->setX0(0.0);
	mine_top->setY0(0.0);
	mine_top->setZ0(-10.0);

	ZCircle* mine_bottom = new ZCircle();
	mine_bottom->setBoundaryType(INTERFACE);
	mine_bottom->setRadius(3.5);
	mine_bottom->setX0(0.0);
	mine_bottom->setY0(0.0);
	mine_bottom->setZ0(-12.0);

	OpenZCylinder* mine_cylinder = new OpenZCylinder();
	mine_cylinder->setBoundaryType(INTERFACE);
	mine_cylinder->setRadius(3.5);
	mine_cylinder->setX0(0.0);
	mine_cylinder->setY0(0.0);
	mine_cylinder->setZ0(-5.0);
	mine_cylinder->setZLeft(-12.0);
	mine_cylinder->setZRight(-10.0);

	mine.setRegionName((char*)"mine");
	mine.setMaterial(&tnt);
	mine.setOpenZCylinder(mine_cylinder);
	mine.setZLeftCircle(mine_bottom);
	mine.setZRightCircle(mine_top);

	mine_top->setLeftRegion(&mine);
	mine_top->setRightRegion(&dirt);
	mine_bottom->setLeftRegion(&dirt);
	mine_bottom->setRightRegion(&mine);
	mine_cylinder->setLeftRegion(&mine);
	mine_cylinder->setLeftRegion(&mine);
	mine_cylinder->setRightRegion(&dirt);

	dirt.addInteriorRegion(&mine);

	/* Set up simple detector region */
	ZCircle* detector_top = new ZCircle();
	detector_top->setBoundaryType(INTERFACE);
	detector_top->setRadius(5.0);
	detector_top->setX0(0.0);
	detector_top->setY0(0.0);
	detector_top->setZ0(25.0);

	ZCircle* detector_bottom = new ZCircle();
	detector_bottom->setBoundaryType(INTERFACE);
	detector_bottom->setRadius(5.0);
	detector_bottom->setX0(0.0);
	detector_bottom->setY0(0.0);
	detector_bottom->setZ0(10.0);

	OpenZCylinder* detector_cylinder = new OpenZCylinder();
	detector_cylinder->setBoundaryType(INTERFACE);
	detector_cylinder->setRadius(5.0);
	detector_cylinder->setX0(0.0);
	detector_cylinder->setY0(0.0);
	detector_cylinder->setZ0(10.0);
	detector_cylinder->setZLeft(10.0);
	detector_cylinder->setZRight(25.0);

	detector.setRegionName((char*)"detector");
	detector.setMaterial(&detector_gas);
	detector.setZLeftCircle(detector_bottom);
	detector.setZRightCircle(detector_top);
	detector.setOpenZCylinder(detector_cylinder);

	detector_top->setLeftRegion(&detector);
	detector_top->setRightRegion(&atmosphere);
	detector_bottom->setLeftRegion(&atmosphere);
	detector_bottom->setRightRegion(&detector);
	detector_cylinder->setLeftRegion(&detector);
	detector_cylinder->setRightRegion(&atmosphere);

	atmosphere.addInteriorRegion(&detector);

	/* Set up implicit capture regions */
	if (options->useImplicitCapture()) {
		float w_low = options->getWeightLow();
		float w_avg = options->getWeightAvg();
		atmosphere.useImplicitCapture(w_low, w_avg);
		dirt.useImplicitCapture(w_low, w_avg);
		mine.useImplicitCapture(w_low, w_avg);
		detector.useImplicitCapture(w_low, w_avg);
	}

	/* Set up forced collision in detector */
	if (options->useForcedCollision()) {
		float w_low = options->getWeightLow();
		float w_avg = options->getWeightAvg();
		detector.useForcedCollision(w_low, w_avg);
	}

	/* Loop over batches */
	for (int i=0; i < num_batches; i++) {

		/* Load atmosphere with neutrons */
		for (int j=0; j < num_neutrons; j++) {
			neutron* new_neutron = initializeNewNeutron();
			new_neutron->_energy = 1E6;
			new_neutron->_mu = -1.0;
			new_neutron->_phi = 0.0;
			new_neutron->_thread_num = 1;
			new_neutron->_time = 0.0;
			new_neutron->_weight = 1.0;
			new_neutron->_x = 0.0;
			new_neutron->_y = 0.0;
			new_neutron->_z = 5.0;

			atmosphere.addNeutron(new_neutron);
		}

		/* Clear the Binners for each Region */
		atmosphere.clearBinners();
		dirt.clearBinners();
		mine.clearBinners();
		detector.clearBinners();

		/* Add the detector collision rate binner over energy */
		detector.addBinner(detector_coll_rate->getBinner(i));

		/* Add collision rate binners over space to each region */
		atmosphere.addBinner(tot_coll_rate_x->getBinner(i));
		dirt.addBinner(tot_coll_rate_x->getBinner(i));
		mine.addBinner(tot_coll_rate_x->getBinner(i));
		detector.addBinner(tot_coll_rate_x->getBinner(i));

		atmosphere.addBinner(tot_coll_rate_y->getBinner(i));
		dirt.addBinner(tot_coll_rate_y->getBinner(i));
		mine.addBinner(tot_coll_rate_y->getBinner(i));
		detector.addBinner(tot_coll_rate_y->getBinner(i));

		atmosphere.addBinner(tot_coll_rate_z->getBinner(i));
		dirt.addBinner(tot_coll_rate_z->getBinner(i));
		mine.addBinner(tot_coll_rate_z->getBinner(i));
		detector.addBinner(tot_coll_rate_z->getBinner(i));

		atmosphere.addBinner(tot_coll_rate_rz->getBinner(i));
		dirt.addBinner(tot_coll_rate_rz->getBinner(i));
		mine.addBinner(tot_coll_rate_rz->getBinner(i));
		detector.addBinner(tot_coll_rate_rz->getBinner(i));


		int num_alive = num_neutrons;
		log_printf(NORMAL, "Testing neutrons inside of Baysoy Subasi mine/detector geometry");

		while (num_alive != 0) {

			log_printf(NORMAL, "batch # = %d\tnum alive = %d\tin air = %d\tin dirt = %d"
					"\tin mine = %d\tin detector = %d", i, num_alive,
					atmosphere.getNumNeutrons(), dirt.getNumNeutrons(),
					mine.getNumNeutrons(), detector.getNumNeutrons());

			atmosphere.moveNeutrons();
			dirt.moveNeutrons();
			mine.moveNeutrons();
			detector.moveNeutrons();

			atmosphere_cylinder->moveNeutrons();
			dirt_cylinder->moveNeutrons();
			mine_cylinder->moveNeutrons();
			detector_cylinder->moveNeutrons();
			atmosphere_top->moveNeutrons();
			atmosphere_dirt_interface->moveNeutrons();
			dirt_bottom->moveNeutrons();
			mine_top->moveNeutrons();
			mine_bottom->moveNeutrons();
			mine_cylinder->moveNeutrons();
			detector_top->moveNeutrons();
			detector_bottom->moveNeutrons();

			num_alive = atmosphere.getNumNeutrons() + dirt.getNumNeutrons() +
					mine.getNumNeutrons() + detector.getNumNeutrons();
		}

		for (int j=0; j < 50; j++) {
			log_printf(NORMAL, "j = %d, x center = %f, tally = %f", j,
					tot_coll_rate_x->getBinner(i)->getBinCenters()[j],
					tot_coll_rate_x->getBinner(i)->getTally(j));
		}
	}


	detector_coll_rate->computeScaledBatchStatistics(num_neutrons);
	tot_coll_rate_x->computeScaledBatchStatistics(num_neutrons);
	tot_coll_rate_y->computeScaledBatchStatistics(num_neutrons);
	tot_coll_rate_z->computeScaledBatchStatistics(num_neutrons);
	tot_coll_rate_rz->computeScaledBatchStatistics(num_neutrons);

	detector_coll_rate->outputBatchStatistics((const char*)"detector_coll_rate.txt");
	tot_coll_rate_x->outputBatchStatistics((const char*)"tot_coll_rate_x.txt");
	tot_coll_rate_y->outputBatchStatistics((const char*)"tot_coll_rate_z.txt");
	tot_coll_rate_z->outputBatchStatistics((const char*)"tot_coll_rate_y.txt");
	tot_coll_rate_rz->outputBatchStatistics((const char*)"tot_coll_rate_rz.txt");

//	tot_coll_rate_x->plotBatchMu((const char*)"tot_coll_rate_x",
//					(const char*)"x position",
//					(const char*)"Collision rate per source neutron");
//	tot_coll_rate_y->plotBatchMu((const char*)"tot_coll_rate_y",
//					(const char*)"y position",
//					(const char*)"Collision rate per source neutron");
//	tot_coll_rate_z->plotBatchMu((const char*)"tot_coll_rate_z",
//					(const char*)"y position",
//					(const char*)"Collision rate per source neutron");
//	tot_coll_rate_rz->plotBatchMu((const char*)"tot_coll_rate_r",
//					(const char*)"r position",
//					(const char*)"Collision rate per source neutron");
//
	log_printf(NORMAL, "Detector coll rate: mu = %1.10f\t var = %1.10f",
			detector_coll_rate->getBatchMu()[0],
			detector_coll_rate->getBatchVariance()[0]);

	delete detector_coll_rate;
	delete tot_coll_rate_x;
	delete tot_coll_rate_y;
	delete tot_coll_rate_z;
	delete tot_coll_rate_rz;
	delete atmosphere_top;
	delete atmosphere_dirt_interface;
	delete dirt_bottom;
	delete atmosphere_cylinder;
	delete dirt_cylinder;
	delete mine_top;
	delete mine_bottom;
	delete mine_cylinder;
	delete detector_top;
	delete detector_bottom;
	delete detector_cylinder;

	return;
}
