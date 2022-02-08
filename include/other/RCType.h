#ifndef RC_TYPE_H
#define RC_TYPE_H

#include "BaseObject.h"

#include <string>

/**
 * Definition of enumeration values for RC types.
 */
class RCType : public BaseObject {

public:

  /**
   * Definition of enumerate values corresponding to axis types.
   */
  enum Type {

    UNDEFINED = 0, //!< Undefined type.

    Born = 1, //!<  Born
    ISR = 2, //!<  Initial state radiation
    FSR = 4  //!<  Final state radiation
  };

  /**
   * Default constructor.
   */
  RCType();

  /**
   * Assignment constructor.
   */
  RCType(Type type);

  /**
   * Copy constructor.
   */
  RCType(const RCType &other);

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
  RCType::Type getType() const;

  /**
   * Assign type to a declared object of this class.
   */
  void setType(Type type);

private:

  /**
   * Type associated to a declared object of this class.
   */
  RCType::Type m_type;
};

#endif