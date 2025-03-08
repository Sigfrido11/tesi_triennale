#include "particle.h"
#include "unordered_map"
#include <algorithm>
#include <cmath>
#include <cstdlib> //for RAND_MAX
#include <iostream>

// inizializzazione membri statici
int const Particle::maxNumParticleType_ = 10;
int Particle::NParticleType_ = 0;
std::unordered_map<std::string, ParticleType *> Particle::particleTypes_;

// costruttori
Particle::Particle(){};

Particle::Particle(std::string name, double px = 0, double py = 0,
                   double pz = 0)
    : name_{name}, px_{px}, py_{py}, pz_{pz} {
  auto search = particleTypes_.find(name);
  if (search == particleTypes_.end()) {
    std::cout << "Particle type not found" << '\n';
  }
}

// Getters
 std::string Particle::GetName() const { return name_; }

 double Particle::GetPx() const { return px_; }
 double Particle::GetPy() const { return py_; }
 double Particle::GetPz() const { return pz_; }

 double Particle::GetMass() const {
  auto search = particleTypes_.find(name_);
  return search->second->GetMass();
}

 int Particle::GetCharge() const { 
  auto search = particleTypes_.find(name_);
  return search->second->GetCharge();
}

 double Particle::GetEnergy() const {
  double const m2 = std::pow(GetMass(), 2);
  double const p2 = std::pow(GetModuleP(), 2);
  return std::sqrt(m2 + p2);
}

 double Particle::GetInvariantMass(Particle & other_dau) const {
  double const sum_e2{
  std::pow(other_dau.GetEnergy() + GetEnergy(), 2)};
  double const pxtot = GetPx()+ other_dau.GetPx();
  double const pytot = GetPy()+ other_dau.GetPy();
  double const pztot = GetPz()+ other_dau.GetPz();
  double const sum_p2{std::pow(pxtot,2) + std::pow(pytot,2) + std::pow(pztot,2)};
  if (sum_e2 < sum_p2){
    std::cout << "error in the energy distribution" << "\n";
  }
  return std::sqrt(sum_e2 - sum_p2);
}

 double Particle::GetModuleP() const {
  double const p2 = std::pow(px_, 2) + std::pow(py_, 2) + std::pow(pz_, 2);
  return std::sqrt(p2);
}

//Setters
void Particle::SetName(std::string new_name) {
  auto search = particleTypes_.find(new_name);
  if (search == particleTypes_.end()) {
    return;
  } else {
    name_ = new_name;
  }
}

void Particle::SetP(double px, double py, double pz) {
  px_ = px;
  py_ = py;
  pz_ = pz;
}

void Particle::Print() const {
  std::cout << "particle index:" << FindParticle(name_) << '\n'
            << "particle type name: " << name_ << '\n'
            << "linear momentum (GeV): (" << px_ << ", " << py_ << ", " << pz_
            << ")" << '\n';
}

int Particle::Decay2Body(Particle &dau1, Particle &dau2) const {
  if (GetMass() == 0.0) {
    printf("Decayment cannot be preformed if mass is zero\n");
    return 1;
  }

  double massMot = GetMass();
  double massDau1 = dau1.GetMass();
  double massDau2 = dau2.GetMass();

  auto search = particleTypes_.find(name_);
  if (search != particleTypes_.end()) { // add width effect

    // gaussian random numbers

    float x1, x2, w, y1;

    double invnum = 1. / RAND_MAX;
    do {
      x1 = 2.0 * std::rand() * invnum - 1.0;
      x2 = 2.0 * std::rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = std::sqrt((-2.0 * std::log(w)) / w);
    y1 = x1 * w;

    massMot += particleTypes_[name_]->GetWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    printf("Decayment cannot be performed because mass is too low in this "
           "channel\n");
    return 2;
  }

  double pout =
      std::sqrt(
          (massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) *
          (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) /
      massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = std::rand() * norm;
  double theta = std::rand() * norm * 0.5 - M_PI / 2.;
  dau1.SetP(pout * std::sin(theta) * std::cos(phi),
            pout * std::sin(theta) * std::sin(phi), pout * std::cos(theta));
  dau2.SetP(-pout * std::sin(theta) * std::cos(phi),
            -pout * std::sin(theta) * std::sin(phi), -pout * std::cos(theta));

  double energy =
      std::sqrt(px_ * px_ + py_ * py_ + pz_ * pz_ + massMot * massMot);

  double bx = px_ / energy;
  double by = py_ / energy;
  double bz = pz_ / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}

void Particle::Boost(double bx, double by, double bz) {

  double energy = GetEnergy();

  // Boost this Lorentz vector
  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / std::sqrt(1.0 - b2);
  double bp = bx * px_ + by * py_ + bz * pz_;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  px_ += gamma2 * bp * bx + gamma * bx * energy;
  py_ += gamma2 * bp * by + gamma * by * energy;
  pz_ += gamma2 * bp * bz + gamma * bz * energy;
}

int Particle::Decay3Body(Particle &dau1, Particle &dau2, Particle &dau3){
  
  // Mother and sum daughter masses. Fail if too close.
  double m0  = GetMass();
  double m1 = dau1.GetMass();
  double m2 = dau1.GetMass();
  double m3 =dau3.GetMass();
  double mSum    = m1 + m2 + m3;
  double mDiff   = m0 - mSum;
  if (mDiff < mSafety) return false;

  // Kinematical limits for 2+3 mass. Maximum phase-space weight.
  double m23Min  = m2 + m3;
  double m23Max  = m0 - m1;
  double p1Max   = 0.5 * sqrtpos( (m0 - m1 - m23Min) * (m0 + m1 + m23Min)
    * (m0 + m1 - m23Min) * (m0 - m1 + m23Min) ) / m0;
  double p23Max  = 0.5 * sqrtpos( (m23Max - m2 - m3) * (m23Max + m2 + m3)
    * (m23Max + m2 - m3) * (m23Max - m2 + m3) ) / m23Max;
  double wtPSmax = 0.5 * p1Max * p23Max;

  // Begin loop over matrix-element corrections.
  double wtME, wtMEmax, wtPS, m23, p1Abs, p23Abs;
  do {
    wtME     = 1.;
    wtMEmax  = 1.;

    // Pick an intermediate mass m23 flat in the allowed range.
    do {
      m23    = m23Min + rndmPtr->flat() * mDiff;

      // Translate into relative momenta and find phase-space weight.
      p1Abs  = 0.5 * sqrtpos( (m0 - m1 - m23) * (m0 + m1 + m23)
        * (m0 + m1 - m23) * (m0 - m1 + m23) ) / m0;
      p23Abs = 0.5 * sqrtpos( (m23 - m2 - m3) * (m23 + m2 + m3)
        * (m23 + m2 - m3) * (m23 - m2 + m3) ) / m23;
      wtPS   = p1Abs * p23Abs;

    // If rejected, try again with new invariant masses.
    } while ( wtPS < rndmPtr->flat() * wtPSmax );

    // Set up m23 -> m2 + m3 isotropic in its rest frame.
    pair<Vec4, Vec4> ps23 = rndmPtr->phaseSpace2(m23, m2, m3);
    prod2.p(ps23.first);
    prod3.p(ps23.second);

    // Set up m0 -> m1 + m23 isotropic in its rest frame.
    pair<Vec4, Vec4> ps123 = rndmPtr->phaseSpace2(m0, m1, m23);
    prod1.p(ps123.first);

    // Boost 2 + 3 to the 0 rest frame.
    prod2.bst( ps123.second, m23 );
    prod3.bst( ps123.second, m23 );

  // If rejected, try again with new invariant masses.
  } while ( wtME < rndmPtr->flat() * wtMEmax );

  // Boost 1 + 2 + 3 to the current frame.
  prod1.bst( decayer.p(), decayer.m() );
  prod2.bst( decayer.p(), decayer.m() );
  prod3.bst( decayer.p(), decayer.m() );

  // Done.
  return true;

  
}

void Particle::AddParticleType(std::string name, double mass, int charge,
                               double width) {
  auto search = particleTypes_.find(name);
  if (search != particleTypes_.end()) {
    return;
  } else {
    if (width == 0) {
      ParticleType *particle = new ParticleType(name, mass, charge);
      particleTypes_[name] = particle;
    } else {
      ResonanceType *particle = new ResonanceType(name, mass, charge, width);
      particleTypes_[name] = particle;
    }
  }
}

void Particle::PrintAllTypes() {
  std::for_each(particleTypes_.begin(), particleTypes_.end(),
                [](std::pair<const std::string, ParticleType*> p) { 
                  p.second->Print(); 
                });
}

 int Particle::FindParticle(std::string name) const {
  return std::distance(particleTypes_.begin(), particleTypes_.find(name));
}