#include "../../include/other/BaseObject.h"

BaseObject::BaseObject(const std::string& className){
	m_className = className;
}

BaseObject::~BaseObject(){
}

const std::string& BaseObject::getClassName() const{
	return m_className;
}