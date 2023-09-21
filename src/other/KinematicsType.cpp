#include "../../include/other/KinematicsType.h"

KinematicsType::KinematicsType()
    : BaseObject("KinematicsType"), m_type(KinematicsType::UNDEFINED) {
}

KinematicsType::KinematicsType(Type type) : BaseObject("KinematicsType"), m_type(type) {

}

KinematicsType::KinematicsType(const KinematicsType &other)
    : BaseObject("KinematicsType"), m_type(other.m_type) {
}

KinematicsType::operator KinematicsType::Type() const { return m_type; }

std::string KinematicsType::toString() const {

  switch (m_type) {

  case Observed:
    return "Observed";
    break;
  case True:
    return "True";
    break;
  case Born:
    return "Born";
    break;

  default:
    return "UNDEFINED";
  }
}

KinematicsType::Type KinematicsType::getType() const { 
  return m_type; 
}

void KinematicsType::setType(Type type) { 
  m_type = type; 
}

