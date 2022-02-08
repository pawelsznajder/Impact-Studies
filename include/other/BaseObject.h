#ifndef BASE_OBJECT_H
#define BASE_OBJECT_H

#include <string>

/*
 * Base object class.
 */
class BaseObject{

public:

	//constructor
	BaseObject(const std::string& className);

	//destructor
	virtual ~BaseObject();

	//get class name
	const std::string& getClassName() const;

private:

	std::string m_className;	//class name
};

#endif