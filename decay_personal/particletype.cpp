#include "particletype.h"

#include <iostream>

ParticleType::ParticleType(std::string name, double mass, int charge)
    : name_{name}, mass_{mass}, charge_{charge} {};

 std::string ParticleType::GetName() const { return name_; }

 double ParticleType::GetMass() const { return mass_; }

 int ParticleType::GetCharge() const { return charge_; }

 double ParticleType::GetWidth() const { return 0.; }

void ParticleType::Print() const {
  std::cout << name_ << " "
            << "mass: " << mass_ << " kg "
            << "charge: " << charge_ << " e" << '\n';
}