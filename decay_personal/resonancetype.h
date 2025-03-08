#ifndef RESONANCE_TYPE_H
#define RESONANCE_TYPE_H

#include "particle.h"
#include "particletype.h"
#include "string"

class ResonanceType final : public ParticleType {
public:
  virtual void Print() const override;

   double GetWidth() const override;

  ResonanceType(std::string name, double mass, int charge, double width);

private:
  double const width_;
};

#endif