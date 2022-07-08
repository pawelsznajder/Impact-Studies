#include "../../include/other/DVCSKinematicFraction.h"

DVCSKinematicFraction::DVCSKinematicFraction(double E1, double E2, size_t nPoints) : BaseObject("DVCSKinematicFraction"){

	//masses
	const double PROTON_MASS = 0.938272013;
	const double ELECTRON_MASS = 0.510998910e-3;

	//fixed target equivalent
    double p1 = sqrt(pow(E1, 2) - pow(ELECTRON_MASS, 2));
    double p2 = sqrt(pow(E2, 2) - pow(PROTON_MASS, 2));

    m_E = (E1 * E2 + p1 * p2) / PROTON_MASS;

    //number of points
    m_nPoints = nPoints;
}

double DVCSKinematicFraction::getFraction(const std::pair<double, double>& xBRange, 
		const std::pair<double, double>& Q2Range, 
		const std::pair<double, double>& tAbsRange, 
		const std::pair<double, double>& yRange){

	//masses
	const double PROTON_MASS = 0.938272013;

	//distributions
	std::uniform_real_distribution<double> xBDistribution(xBRange.first, xBRange.second);
	std::uniform_real_distribution<double> Q2Distribution(Q2Range.first, Q2Range.second);
	std::uniform_real_distribution<double> tDistribution(tAbsRange.first, tAbsRange.second);

	//points
	size_t nAll = 0;
	size_t nAccepted = 0;

	//loop over points
	for(;;){

		//check limit
		if(nAll == m_nPoints) break;

		//generate numbers
		double xB = xBDistribution(m_generator);
		double Q2 = Q2Distribution(m_generator);
		double t = -1 * tDistribution(m_generator);

		//count 
		nAll++;

		//check y
		double y = Q2 / (2 * xB * PROTON_MASS * m_E);
		if(y < yRange.first || y > yRange.second) continue;

		//check limit
		if(! isKinematicsValid(xB, Q2, t)) continue;

		//count
		nAccepted++;
	}

	//return
	return nAccepted/double(nAll);
}

bool DVCSKinematicFraction::isKinematicsValid(double xB, double Q2, double t) const{

	//masses
	const double PROTON_MASS = 0.938272013;

	//variables	
 	double epsilon = 2 * xB * PROTON_MASS / sqrt(Q2);
 	double y = Q2 / (2 * xB * PROTON_MASS * m_E);
    double eps2 = epsilon * epsilon;
    double epsroot = sqrt(1 + eps2);
    double tfactor = -Q2 / (4 * xB * (1 - xB) + eps2);
    double tmin = tfactor * (2 * (1 - xB) * (1 - epsroot) + eps2);
    double tmax = tfactor * (2 * (1 - xB) * (1 + epsroot) + eps2);
    double xBmin = 2 * Q2 * m_E / PROTON_MASS / (4 * m_E * m_E - Q2);

    if (xB < xBmin || xB > 1.) {
	  return false;
    }

    if (t > tmin || t < tmax) {
        return false;
    }

    if (y < 0. || y > 1.) {
      return false;
    }

    return true;
}