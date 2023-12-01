#include "../../include/kinematic_cuts/KinematicCuts.h"

namespace KinematicCuts{

bool checkKinematicCuts(const DVCSEvent& event, KinematicsType::Type kinematicsType){

	double xB = event.getXB(kinematicsType);
	double Q2 = event.getQ2(kinematicsType);
	double mT = fabs(event.getT(kinematicsType));
	double y = event.getY(kinematicsType);

	if(xB < 1.E-4 || xB > 0.7) return false;
	if(Q2 < 1. || Q2 > 100.) return false;
	if(mT < 0. || mT > 1.2) return false;
	if(y < 0.05 || y > 0.6) return false;

	return true;
}

}