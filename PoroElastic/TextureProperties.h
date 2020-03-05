#ifndef TEXTURE_PROPERTIES_H_
#define TEXTURE_PROPERTIES_H_

#include "MatVec.h"
#include <map>

class TiXmlElement;
class Vec3;


class TextureProperties {
public:
  void parse (const TiXmlElement* elem);
  void printLog() const;
  bool getProperty(const std::string& name, const Vec3& X, double& val) const;

protected:
  struct Property {
    double min;
    double max;
    Matrix3D textureData;
    bool prescaled = false;
  };

  std::map<std::string, Property> properties;
};

#endif
