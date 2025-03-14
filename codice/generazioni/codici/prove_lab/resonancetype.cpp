#include "resonancetype.h"
#include "string"

ResonanceType::ResonanceType(std::string name, double mass, int charge, double width=0)
    : ParticleType(name, mass, charge), width_{width} {}

void ResonanceType::Print() const {
  ParticleType::Print();
  std::cout << "width: " << width_ << '\n';
}

double ResonanceType::GetWidth() const { return width_; }