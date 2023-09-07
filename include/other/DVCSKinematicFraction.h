#ifndef DVCS_KINEMATIC_FRACTION_H
#define DVCS_KINEMATIC_FRACTION_H

#include <random>
#include <utility>

#include "../../include/other/BaseObject.h"

/*
 * Get '\Delta_{allowed} / \Delta_{all}' fraction, where '\Delta_{all}' is a (\Delta xB, \Delta Q2, \Delta t, \Delta y) volume, 
 * and '\Delta_{allowed}' is a subdomain of this volume allowed by DVCS kinematics. 
 */
class DVCSKinematicFraction : public BaseObject{

public:

	//constructor taking beam energies
	DVCSKinematicFraction(double Elepton, double Eproton, size_t nPoints = 1E3);

	//get fraction
	double getFraction(const std::pair<double, double>& xBRange, 
		const std::pair<double, double>& Q2Range, 
		const std::pair<double, double>& tAbsRange, 
		const std::pair<double, double>& yRange);

private:

	/*
	 *	Check if kinematics valid. 
	 */
	bool isKinematicsValid(double xB, double Q2, double t) const; 

	double m_E;	//lepton energy in fix target mode
	size_t m_nPoints;	//number of points probing the phase-space 
	std::default_random_engine m_generator; //random generator
};

#endif