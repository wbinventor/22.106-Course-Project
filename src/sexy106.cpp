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
#include "Isotope.h"
#include "Material.h"
#include "Neutron.h"
#include "Region.h"


//TODO: Test different region types
//TODO: Create neutron source
//TODO: Tie different types of regions together

int main(int argc, const char **argv) {

	Options options(argc, argv);
	Timer timer;
	log_setlevel(options.getVerbosity());

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
	float rho_soil;
	float rho_H1_soil, rho_O16_soil, rho_Si28_soil, rho_Si29_soil;
	float rho_Si30_soil, rho_Al27_soil, rho_Fe54_soil, rho_Fe56_soil;
	float rho_Fe57_soil, rho_Fe58_soil, rho_Ca40_soil, rho_Ca42_soil;
	float rho_Ca44_soil, rho_K39_soil, rho_K41_soil, rho_Mg24_soil;
	float rho_Na23_soil, rho_Mg25_soil, rho_Mg26_soil;

	if (options.getSoilType() == DRY_POROUS) {
		rho_soil = 1.1810;
		rho_H1_soil = rho_soil * 0.015 * N_A / 1.0;
		rho_O16_soil = rho_soil * 0.529 * N_A / 16.0;
		rho_Si28_soil = rho_soil * 0.243 * 0.92223 * N_A / 28.0;
		rho_Si29_soil = rho_soil * 0.243 * 0.04685 * N_A / 29.0;
		rho_Si30_soil = rho_soil * 0.243 * 0.03092 * N_A / 30.0;
		rho_Al27_soil = rho_soil * 0.071 * N_A / 27;
		rho_Fe54_soil = rho_soil * 0.044 * 0.05845 * N_A / 54.0;
		rho_Fe56_soil = rho_soil * 0.044 * 0.91754 * N_A / 56.0;
		rho_Fe57_soil = rho_soil * 0.044 * 0.02119 * N_A / 57.0;
		rho_Fe58_soil = rho_soil * 0.044 * 0.00282 * N_A / 58.0;
		rho_Ca40_soil = rho_soil * 0.032 * 0.96941 * N_A / 40.0;
		rho_Ca42_soil = rho_soil * 0.032 * 0.00647 * N_A / 42.0;
		rho_Ca44_soil = rho_soil * 0.032 * 0.02086 * N_A / 44.0;
		rho_K39_soil = rho_soil * 0.023 * 0.932581 * N_A / 39.0;
		rho_K41_soil = rho_soil * 0.023 * 0.067302 * N_A / 41.0;
		rho_Na23_soil = rho_soil * 0.025 * N_A / 23.0;
		rho_Mg24_soil = rho_soil * 0.018 * 0.78958 * N_A / 24.0;
		rho_Mg25_soil = rho_soil * 0.018 * 0.0999986 * N_A / 25.0;
		rho_Mg26_soil =  rho_soil * 0.018 * 0.10987 * N_A / 26.0;
	}
	else if (options.getSoilType() == DRY_DENSE) {
		rho_soil = 1.7714;
		rho_H1_soil = rho_soil * 0.015 * N_A / 16.0;
		rho_O16_soil = rho_soil * 0.529 * N_A / 1.0;
		rho_Si28_soil = rho_soil * 0.243 * 0.92223 * N_A / 28.0;
		rho_Si29_soil = rho_soil * 0.243 * 0.04685 * N_A / 29.0;
		rho_Si30_soil = rho_soil * 0.243 * 0.03092 * N_A / 30.0;
		rho_Al27_soil = rho_soil * 0.071 * N_A / 27;
		rho_Fe54_soil = rho_soil * 0.044 * 0.05845 * N_A / 54.0;
		rho_Fe56_soil = rho_soil * 0.044 * 0.91754 * N_A / 56.0;
		rho_Fe57_soil = rho_soil * 0.044 * 0.02119 * N_A / 57.0;
		rho_Fe58_soil = rho_soil * 0.044 * 0.00282 * N_A / 58.0;
		rho_Ca40_soil = rho_soil * 0.032 * 0.96941 * N_A / 40.0;
		rho_Ca42_soil = rho_soil * 0.032 * 0.00647 * N_A / 42.0;
		rho_Ca44_soil = rho_soil * 0.032 * 0.02086 * N_A / 44.0;
		rho_K39_soil = rho_soil * 0.023 * 0.932581 * N_A / 39.0;
		rho_K41_soil = rho_soil * 0.023 * 0.067302 * N_A / 41.0;
		rho_Na23_soil = rho_soil * 0.025 * N_A / 23.0;
		rho_Mg24_soil = rho_soil * 0.018 * 0.78958 * N_A / 24.0;
		rho_Mg25_soil = rho_soil * 0.018 * 0.0999986 * N_A / 25.0;
		rho_Mg26_soil =  rho_soil * 0.018 * 0.10987 * N_A / 26.0;
	}
	else if (options.getSoilType() == WET_POROUS) {
		rho_soil = 1.3957;
		rho_H1_soil = rho_soil * 0.030 * N_A / 1.0;
		rho_O16_soil = rho_soil * 0.585 * N_A / 16.0;
		rho_Si28_soil = rho_soil * 0.205 * 0.92223 * N_A / 28.0;
		rho_Si29_soil = rho_soil * 0.205 * 0.04685 * N_A / 29.0;
		rho_Si30_soil = rho_soil * 0.205 * 0.03092 * N_A / 30.0;
		rho_Al27_soil = rho_soil * 0.060 * N_A / 27;
		rho_Fe54_soil = rho_soil * 0.037 * 0.05845 * N_A / 54.0;
		rho_Fe56_soil = rho_soil * 0.037 * 0.91754 * N_A / 56.0;
		rho_Fe57_soil = rho_soil * 0.037 * 0.02119 * N_A / 57.0;
		rho_Fe58_soil = rho_soil * 0.037 * 0.00282 * N_A / 58.0;
		rho_Ca40_soil = rho_soil * 0.027 * 0.96941 * N_A / 40.0;
		rho_Ca42_soil = rho_soil * 0.027 * 0.00647 * N_A / 42.0;
		rho_Ca44_soil = rho_soil * 0.027 * 0.02086 * N_A / 44.0;
		rho_K39_soil = rho_soil * 0.019 * 0.932581 * N_A / 39.0;
		rho_K41_soil = rho_soil * 0.019 * 0.067302 * N_A / 41.0;
		rho_Na23_soil = rho_soil * 0.021 * N_A / 23.0;
		rho_Mg24_soil = rho_soil * 0.015 * 0.78958 * N_A / 24.0;
		rho_Mg25_soil = rho_soil * 0.015 * 0.0999986 * N_A / 25.0;
		rho_Mg26_soil =  rho_soil * 0.015 * 0.10987 * N_A / 26.0;
	}
	else {
		rho_soil = 2.0935;
		rho_H1_soil = rho_soil * 0.030 * N_A / 1.0;
		rho_O16_soil = rho_soil * 0.585 * N_A / 16.0;
		rho_Si28_soil = rho_soil * 0.205 * 0.92223 * N_A / 28.0;
		rho_Si29_soil = rho_soil * 0.205 * 0.04685 * N_A / 29.0;
		rho_Si30_soil = rho_soil * 0.205 * 0.03092 * N_A / 30.0;
		rho_Al27_soil = rho_soil * 0.060 * N_A / 27;
		rho_Fe54_soil = rho_soil * 0.037 * 0.05845 * N_A / 54.0;
		rho_Fe56_soil = rho_soil * 0.037 * 0.91754 * N_A / 56.0;
		rho_Fe57_soil = rho_soil * 0.037 * 0.02119 * N_A / 57.0;
		rho_Fe58_soil = rho_soil * 0.037 * 0.00282 * N_A / 58.0;
		rho_Ca40_soil = rho_soil * 0.027 * 0.96941 * N_A / 40.0;
		rho_Ca42_soil = rho_soil * 0.027 * 0.00647 * N_A / 42.0;
		rho_Ca44_soil = rho_soil * 0.027 * 0.02086 * N_A / 44.0;
		rho_K39_soil = rho_soil * 0.019 * 0.932581 * N_A / 39.0;
		rho_K41_soil = rho_soil * 0.019 * 0.067302 * N_A / 41.0;
		rho_Na23_soil = rho_soil * 0.021 * N_A / 23.0;
		rho_Mg24_soil = rho_soil * 0.015 * 0.78958 * N_A / 24.0;
		rho_Mg25_soil = rho_soil * 0.015 * 0.0999986 * N_A / 25.0;
		rho_Mg26_soil =  rho_soil * 0.015 * 0.10987 * N_A / 26.0;
	}


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
	He3.loadXS("../pendf/He3-capture.txt", CAPTURE, "\t");
	He3.loadXS("../pendf/He3-elastic.txt", ELASTIC, "\t");
	He3.setElasticAngleType(ISOTROPIC_CM);

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


	/* Setup Surfaces and Regions */
	XPlane* x_left = new XPlane();
	XPlane* x_right = new XPlane;
	YPlane* y_left = new YPlane();
	YPlane* y_right = new YPlane;
	ZPlane* z_left = new ZPlane();
	ZPlane* z_right = new ZPlane;

	XCircle* x_circle_left = new XCircle();
	XCircle* x_circle_right = new XCircle();
	YCircle* y_circle_left = new YCircle();
	YCircle* y_circle_right = new YCircle();
	ZCircle* z_circle_left = new ZCircle();
	ZCircle* z_circle_right = new ZCircle();

	Sphere* sphere = new Sphere();

	OpenXCylinder* x_cylinder = new OpenXCylinder();
	OpenYCylinder* y_cylinder = new OpenYCylinder();
	OpenZCylinder* z_cylinder = new OpenZCylinder();

	x_left->setX(0.0);
	x_right->setX(5.0);
	y_left->setY(0.0);
	y_right->setY(5.0);
	z_left->setZ(0.0);
	z_right->setZ(5.0);

	x_left->setBoundaryType(VACUUM);
	x_right->setBoundaryType(VACUUM);
	y_left->setBoundaryType(VACUUM);
	y_right->setBoundaryType(VACUUM);
	z_left->setBoundaryType(VACUUM);
	z_right->setBoundaryType(VACUUM);

	x_circle_left->setX0(0.0);
	x_circle_left->setY0(0.0);
	x_circle_left->setZ0(0.0);
	x_circle_left->setRadius(2.5);

	x_circle_right->setX0(0.0);
	x_circle_right->setY0(0.0);
	x_circle_right->setZ0(0.0);
	x_circle_right->setRadius(2.5);

	y_circle_left->setX0(0.0);
	y_circle_left->setY0(0.0);
	y_circle_left->setZ0(0.0);
	y_circle_left->setRadius(2.5);

	y_circle_right->setX0(0.0);
	y_circle_right->setY0(0.0);
	y_circle_right->setZ0(0.0);
	y_circle_right->setRadius(2.5);

	z_circle_left->setX0(0.0);
	z_circle_left->setY0(0.0);
	z_circle_left->setZ0(0.0);
	z_circle_left->setRadius(2.5);

	z_circle_right->setX0(0.0);
	z_circle_right->setY0(0.0);
	z_circle_right->setZ0(0.0);
	z_circle_right->setRadius(2.5);

	x_circle_left->setBoundaryType(VACUUM);
	x_circle_right->setBoundaryType(VACUUM);
	y_circle_left->setBoundaryType(VACUUM);
	y_circle_right->setBoundaryType(VACUUM);
	z_circle_left->setBoundaryType(VACUUM);
	z_circle_right->setBoundaryType(VACUUM);

	sphere->setBoundaryType(VACUUM);
	sphere->setX0(2.0);
	sphere->setY0(2.0);
	sphere->setZ0(0.0);
	sphere->setRadius(1.5);

	x_cylinder->setBoundaryType(VACUUM);
	x_cylinder->setX0(0.0);
	x_cylinder->setY0(2.0);
	x_cylinder->setZ0(2.0);
	x_cylinder->setRadius(1.5);
	x_cylinder->setXLeft(0.0);
	x_cylinder->setXLeft(5.0);

	y_cylinder->setBoundaryType(VACUUM);
	y_cylinder->setX0(0.0);
	y_cylinder->setY0(2.0);
	y_cylinder->setZ0(2.0);
	y_cylinder->setRadius(1.5);
	y_cylinder->setYLeft(0.0);
	y_cylinder->setYLeft(5.0);

	z_cylinder->setBoundaryType(VACUUM);
	z_cylinder->setX0(0.0);
	z_cylinder->setY0(2.0);
	z_cylinder->setZ0(2.0);
	z_cylinder->setRadius(1.5);
	z_cylinder->setZLeft(0.0);
	z_cylinder->setZLeft(5.0);


	/* First test region */
	Parallelepiped test;
	test.setRegionName((char*)"test cube");
	test.setXLeftSurface(x_left);
	test.setXRightSurface(x_right);
	test.setYLeftSurface(y_left);
	test.setYRightSurface(y_right);
	test.setZLeftSurface(z_left);
	test.setZRightSurface(z_right);
	test.setMaterial(&soil);
//	test.setMaterial(&air);

	x_left->setRightRegion(&test);
	x_right->setLeftRegion(&test);
	y_left->setRightRegion(&test);
	y_right->setLeftRegion(&test);
	z_left->setRightRegion(&test);
	z_right->setLeftRegion(&test);

	/* Load cube with */
	for (int i=0; i < 1000; i++) {
		neutron* new_neutron = initializeNewNeutron();
		new_neutron->_energy = 1E6;
		new_neutron->_mu = -1.0;
		new_neutron->_phi = 0.0;
		new_neutron->_thread_num = 1;
		new_neutron->_time = 0.0;
		new_neutron->_weight = 1.0;
		new_neutron->_x = 2.5;
		new_neutron->_y = 2.5;
		new_neutron->_z = 2.5;

		test.addNeutron(new_neutron);
	}

	int num_alive = 1000;
	while (num_alive != 0) {
		log_printf(NORMAL, "num alive = %d", num_alive);
		test.moveNeutrons();
		num_alive = test.getNumNeutrons();

	}



	timer.printSplits();

	return 0;
}
