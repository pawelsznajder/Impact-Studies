#ifndef SUBPROCESS_TYPE_H
#define SUBPROCESS_TYPE_H

#include "BaseObject.h"

#include <string>

/**
 * Definition of enumeration values for RC types.
 */
class SubProcessType : public BaseObject {

public:

  /**
   * Definition of enumerate values corresponding to axis types.
   */
  enum Type {

    UNDEFINED = 0, //!< Undefined type.

    BH = 1, //!<  BH
    INT = 2, //!< Interference
    DVCS = 4 //! < DVCS
  };


  /**
   * Get subprocess mask from std::string.
   */
  static int getSubProcessTypeMaskFromStdString(const std::string& str);

  /**
   * Default constructor.
   */
  SubProcessType();

  /**
   * Assignment constructor.
   */
  SubProcessType(Type type);

  /**
   * Copy constructor.
   */
  SubProcessType(const SubProcessType &other);

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
  SubProcessType::Type getType() const;

  /**
   * Assign type to a declared object of this class.
   */
  void setType(Type type);

private:

  /**
   * Type associated to a declared object of this class.
   */
  SubProcessType::Type m_type;
};

#endif