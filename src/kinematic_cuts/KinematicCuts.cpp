#include "../../include/kinematic_cuts/KinematicCuts.h"

namespace KinematicCuts{

bool checkKinematicCuts(const DVCSEvent& event, KinematicsType::Type kinematicsType, int skip){

	double xB = event.getXB(kinematicsType);
	double Q2 = event.getQ2(kinematicsType);
	double mT = fabs(event.getT(kinematicsType));
	double y = event.getY(kinematicsType);

	/*
	generated:
	    <param name="range_xB" value="0.00001|0.7" />
        <param name="range_t" value="-1.6|-0.01" />
        <param name="range_Q2" value="1.0|100.0" />
        <param name="range_phi" value="0.03|6.253185" />
        <param name="range_phiS" value="0.03|6.253185" />
        <param name="range_y" value="0.01|0.9" />
	*/

	if(skip != 0 && (mT < 0.05 || mT > 1.2)) return false;
	if(skip != 1 && (y < 0.05 || y > 0.6)) return false;
	if(skip != 2 && (Q2 < 1. || Q2 > 100.)) return false;
	if(skip != 3 && (xB < 0.00001 || xB > 0.7)) return false;


	return true;
}

}