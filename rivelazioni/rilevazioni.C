// lo faccio compilato e riciclo come una vagabonda il vecchio codice
// DOPO LA SESSIONE CON PYTHIA FAI TUTTI I DECADIMENTI PER BENE

#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include "particle.h"
#include "particletype.h"
#include "resonancetype.h"
#include <TCanvas.h>
#include <algorithm>
#include <cmath>
#include <iterator>

int main() {
  TH1::AddDirectory(kFALSE);
  TFile *file = new TFile("densFreq.root", "RECREATE");
  TGraph *freqCom = (TGraph *)file->Get("hFreqComu");
  gRandom->SetSeed(3032003);
  Particle::AddParticleType("pion +", 0.13957, 1);
  Particle::AddParticleType("deuteron", 0.49367, 0);
  Particle::AddParticleType("c-deuteron", 3.225000, 0,
                            0.050); // cerca la larghezza di risonanza
  TH1F *pModuleDist = new TH1F("pModuleDist", "pModuleDist", 500, 0, 5);
  TH1F *dModule = new TH1F("dModule", "dModule", 500, 0, 5);
  std::array<double, 3> decProb = {0.3, 0.6, 0.9}; // valori di prova
  double detected{0};

  for (int i{0}; i < 1e6; i++) {
    double phi = gRandom->Uniform(0., 2 * M_PI); // azimuth angle
    double theta = gRandom->Uniform(0., M_PI);   // polar angle
    double randDist = gRandom->Rndm();
    double pM;
    int j{0};
    while (true) {
      if (randDist > freqCom->GetPointY(j)) {
        pM = freqCom->GetPointX(j);
        break;
      }
      j++;
    }
    pModuleDist->Fill(pM);
    double px = pM * std::cos(phi) * std::sin(theta);
    double py = pM * std::sin(phi) * std::sin(theta);
    double pz = pM * std::cos(theta);

    Particle p;
    p.SetName("c-deuteron");
    p.SetP(px, py, pz);
    Particle pD;
    Particle p2;
    // Particle p3;
    // double decayMode = gRandom->Rndm();
    // tentativo di fare decadimenti a più corpi

    int allGood = p.Decay2Body(pD, p2);
    dModule->Fill(pD.GetModuleP());

    if (pD.GetModuleP() < 0.2) {
      double var = gRandom->Rndm();
      if (var < 0.25) {
        detected++;
      }
    }

    else if (pD.GetModuleP() < 0.5) {
      double var = gRandom->Rndm();
      if (var < 0.5) {
        detected++;
      }
    } else if (pD.GetModuleP() < 0.8) {
      double var = gRandom->Rndm();
      if (var < 0.9) {
        detected++;
      }
    }
  }
  std::cout << "nuntio vobis gaudium magnum habemus numero" << detected << '\n';
  pModuleDist->Draw("APE");
  dModule->Draw("APE");
}
