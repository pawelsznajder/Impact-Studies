#include "../../include/other/SubProcessType.h"

SubProcessType::SubProcessType()
    : BaseObject("SubProcessType"), m_type(SubProcessType::UNDEFINED) {
  }

SubProcessType::SubProcessType(Type type) : BaseObject("SubProcessType"), m_type(type) {
}

SubProcessType::SubProcessType(const SubProcessType &other)
    : BaseObject("SubProcessType"), m_type(other.m_type) {
}

SubProcessType::operator SubProcessType::Type() const { return m_type; }

std::string SubProcessType::toString() const {

  switch (m_type) {

  case BH:
    return "BH";
    break;
  case INT:
    return "INT";
    break;
  case DVCS:
    return "DVCS";
    break;

  default:
    return "UNDEFINED";
  }
}

SubProcessType::Type SubProcessType::getType() const { 
  return m_type; 
}

void SubProcessType::setType(Type type) { 
  m_type = type; 
}

int SubProcessType::getSubProcessTypeMaskFromStdString(const std::string& str){

  int result = 0;

  if(str.find("BH") != std::string::npos) result |= SubProcessType::BH;
  if(str.find("INT") != std::string::npos) result |= SubProcessType::INT;
  if(str.find("DVCS") != std::string::npos) result |= SubProcessType::DVCS;

  return result;
}

