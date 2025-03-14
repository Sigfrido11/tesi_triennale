#include "particle.h"
#include "particletype.h"
#include "resonancetype.h"

#include <iostream>

int main() {
  ParticleType const a{"kaon +", 0.4937, 1};
  ResonanceType const b{"resonance", 0.89166, 0, 0.049};

  // test GetName()
  std::cout << "Nome di a: " << a.GetName() << '\n';
  std::cout << "Nome di b: " << b.GetName() << '\n';

  std::cout << '\n';

  // test GetMass()
  std::cout << "Massa di a: " << a.GetMass() << '\n';
  std::cout << "Massa di b: " << b.GetMass() << '\n';

  std::cout << '\n';

  // test GetCharge()
  std::cout << "Carica di a: " << a.GetCharge() << '\n';
  std::cout << "Carica di b: " << b.GetCharge() << '\n';

  std::cout << '\n';

  // test GetWidth()
  std::cout << "Width di a: " << a.GetWidth() << '\n';
  std::cout << "Width di b: " << b.GetWidth() << '\n';

  std::cout << '\n';

  // test Print()
  std::cout << "a: " << '\n';
  a.Print();
  std::cout << "b: " << '\n';
  b.Print();

  Particle::AddParticleType("pion +", 0.13957, 1);
  Particle::AddParticleType("pion -", 0.13957, -1);
  Particle::AddParticleType("kaon +", 0.4937, 1);
  Particle::AddParticleType("kaon -", 0.49367, -1);
  Particle::AddParticleType("proton +", 0.93827, 1);
  Particle::AddParticleType("proton -", 0.93827, -1);
  Particle::AddParticleType("resonance", 0.89166, 0, 0.050);

  std::cout <<'\n' << "I tipi di particelle aggiunti sono: "<< '\n';
  Particle::PrintAllTypes();

  Particle const p1("kaon +", 350, 250, 200);
  Particle p2("kaon -", 200, 456, 317);

  // test GetName()
  std::cout << '\n' << "Nome di p1: " << p1.GetName() << '\n';

  // test GetPx()
  std::cout << '\n' << "Px di p1: " << p1.GetPx() << '\n';

  // test GetPy()
  std::cout << '\n' << "Py di p1: " << p1.GetPy() << '\n';

  // test GetPz()
  std::cout << '\n' << "Pz di p1: " << p1.GetPz() << '\n';

  // test GetMass()
  std::cout << '\n' << "Massa di p1: " << p1.GetMass() << '\n';

  // test GetEnergy()
  std::cout << '\n' << "Energia di p1: " << p1.GetEnergy() << '\n';
  std::cout << '\n' << "Energia di p2: " << p2.GetEnergy() << '\n';

  // test GetInvariantMass()
  std::cout << '\n' << "Massa invariante di p1 e p2: " << p1.GetInvariantMass(p2) << '\n';

  // test SetName()
  p2.SetName("proton +");
  std::cout << '\n' << "Il nuovo nome di p2 e': " << p2.GetName() << '\n';

  // test SetP()
  p2.SetP(14, 56, 9);
  std::cout << '\n' << "Il nuovo P e': " << p2.GetPx() << ' ' << p2.GetPy() << ' '
            << p2.GetPz() << '\n';

  // test Print()
  std::cout << '\n' << "p1: " << '\n';
  p1.Print();
  std::cout << '\n' << "p2: " << '\n';
  p2.Print();

  std::cout << '\n';
}