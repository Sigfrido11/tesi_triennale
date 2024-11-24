// lo faccio compilato e riciclo come una vagabonda il vecchio codice
// DOPO LA SESSIONE CON PYTHIA FAI TUTTI I DECADIMENTI PER BENE

#include "Pythia.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include "particle.h"
#include "particletype.h"
#include "resonancetype.h"
#include <TCanvas.h>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>

using namespace Pythia8;

int main() {
  TH1::AddDirectory(kFALSE);
  TFile *file = new TFile("densFreq.root", "READ");
  if (!file || file->IsZombie()) {
    std::cerr << "Errore nell'apertura del file: ";
  }
  TGraph *freqCom = (TGraph *)file->Get("FreqComu");
  if (!freqCom) {
    std::cerr << "Oggetto non trovato nel file." << std::endl;
    file->Close();
  }

  // Inizializzazione di Pythia
  Pythia pythia;

  // Definizione della particella c-deuteron e dei suoi decadimenti
  pythia.readString(
      "ParticleData:magicNumber = 10304122"); // ID particella personalizzata
  pythia.readString(
      "ParticleData:m0[10304122] = 2.0"); // Massa della particella -deuteron (2
                                          // GeV)
  pythia.readString("ParticleData:charge[10304122] = 1"); // Carica

  // Attivazione decadimenti per la particella c-deuteron
  pythia.readString("10304122:all = on");

  // Modalità di decadimento:
  // 1. c-deuteron -> pion + kaon
  pythia.readString(
      "10304122:oneChannel = 1 1 211 321"); // c-deuteron -> pion (ID 211)
                                            // + kaon (ID 321)
  // 2. c-deuteron -> due pioni
  pythia.readString(
      "10304122:oneChannel = 1 2 211 211"); // c-deuteron -> due pioni (ID 211)

  // Preparazione per la simulazione
  pythia.init();

  gRandom->SetSeed(3032003);

  double detected{0};
  TCanvas *c = new TCanvas("c", "Canvas", 800, 600);
  freqCom->Draw("APE");

  for (int i{0}; i < 1e6; i++) {
    double phi = gRandom->Uniform(0., 2 * M_PI); // azimuth angle
    double theta = gRandom->Uniform(0., M_PI);   // polar angle
    double randDist = gRandom->Rndm();
    double pM;
    int j{1};
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

    double pt = std::sqrt(px * px + py * py);
    double E = std::sqrt(pM * pM +
                         3.225000 * 3.225000); // dove 3.225000 massa c-deuteron
    double eta = 0.5 * std::log((E + pz) / (E - pz));

    // Crea un oggetto particella di pythia, attento
    Particle cDeuteron(particellaID, pt, eta, phi, E);

    // Aggiungi la particella all'evento di Pythia
    pythia.event.append(cDeuteron);

    // Esegui il decadimento
    pythia.next();
    double DPM; // deuteron P Module

    for (int i = 0; i < pythia.event.size(); ++i) {
      const Particle &p = pythia.event[i];
         
      }

      if (DPM < 0.2) {
        double var = gRandom->Rndm();
        if (var < 0.25) {
          detected++;
        }
      }

      else if (DPM < 0.5) {
        double var = gRandom->Rndm();
        if (var < 0.5) {
          detected++;
        }
      } else if (DPM < 0.8) {
        double var = gRandom->Rndm();
        if (var < 0.9) {
          detected++;
        }
      }
    }
    std::cout << "nuntio vobis gaudium magnum habemus numero" << detected
              << '\n';
    pModuleDist->Draw("APE");
    dModule->Draw("APE");
  }
