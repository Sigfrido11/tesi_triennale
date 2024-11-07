#ifndef PARTICLE_TYPE_H
#define PARTICLE_TYPE_H

#include <string>

class ParticleType
{
public:
  std::string GetName() const; // forse si può togliere
  double GetMass() const;
  int GetCharge() const;
  virtual double GetWidth() const;
  virtual void Print() const;
  explicit ParticleType(std::string name, double mass, int charge);

private:
  std::string const name_;
  double const mass_;
  int const charge_;
};
#endif