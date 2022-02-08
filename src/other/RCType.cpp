#include "../../include/other/RCType.h"

RCType::RCType()
    : BaseObject("RCType"), m_type(RCType::UNDEFINED) {

    }

RCType::RCType(Type type) : BaseObject("RCType"), m_type(type) {

}

RCType::RCType(const RCType &other)
    : BaseObject("RCType"), m_type(other.m_type) {
}

RCType::operator RCType::Type() const { return m_type; }

std::string RCType::toString() const {

  switch (m_type) {

  case Born:
    return "Born";
    break;
  case ISR:
    return "ISR";
    break;
  case FSR:
    return "FSR";
    break;

  default:
    return "UNDEFINED";
  }
}

RCType::Type RCType::getType() const { 
  return m_type; 
}

void RCType::setType(Type type) { 
  m_type = type; 
}

