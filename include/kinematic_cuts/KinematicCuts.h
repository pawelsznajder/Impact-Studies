#ifndef KINEMATIC_CUTS_H
#define KINEMATIC_CUTS_H

#include "../include/event/DVCSEvent.h"

namespace KinematicCuts{

/*
 * Check kinematic cuts.
 */
bool checkKinematicCuts(const DVCSEvent& event, KinematicsType::Type kinematicsType);

}

#endif