#ifndef KINEMATICS_TYPE_H
#define KINEMATICS_TYPE_H

#include "BaseObject.h"

#include <string>

/**
 * Definition of enumeration values for kinematics types.
 */
class KinematicsType : public BaseObject {

public:

  /**
   * Definition of enumerate values corresponding to axis types.
   */
  enum Type {

    UNDEFINED = 0, //!< Undefined type.

    Observed = 1, //!<  Observed kinematics.
    True = 2, //!<  True kinematics (corresponding to observed one).
    Born = 3 //!<  True Born kinematics of the process (entering the evaluation of cross-section).
  };

  /**
   * Default constructor.
   */
  KinematicsType();

  /**
   * Assignment constructor.
   */
  KinematicsType(Type type);

  /**
   * Copy constructor.
   */
  KinematicsType(const KinematicsType &other);

  /**
   * Automatic cast to enum.
   */
  operator Type() const;

  /**
   * Get std::string with name.
   */
  std::string toString() const;

  /**
   * Get type being assigned to a declared object of this class.
   */
  KinematicsType::Type getType() const;

  /**
   * Assign type to a declared object of this class.
   */
  void setType(Type type);

private:

  /**
   * Type associated to a declared object of this class.
   */
  KinematicsType::Type m_type;
};

#endif