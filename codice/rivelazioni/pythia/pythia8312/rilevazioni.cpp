// lo faccio compilato con ROOT in combo con Pythia per i decadimenti a più
// corpi

// per compilare senza stare a scrivere un cmake che poi coi due pc..
//  g++ -o rilevazioni rilevazioni.cpp -I./include -L./lib -l:libpythia8.a
//  $(root-config --cflags --libs)

#include "Pythia8/Pythia.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TRandom.h"
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
  int pid = 10304122;
  double mass = 2.286;
  double charge = 1.0;
  double spin = 0.5;
  double width = 0.1;
  std::string name = "c-deuteron";

  pythia.particleData.addParticle(pid, name, charge, 0, spin, 0, mass, width);
  if (pythia.particleData.isParticle(pid)) {
    std::cout << "La particella è stata aggiunta!" << std::endl;
  }
  // Attivazione dei decadimenti per la particella
  pythia.readString("10304122:all = on");

  // Modalità di decadimento:
  // 1. c-deuteron -> pion + kaon
  pythia.readString(
      "10304122:oneChannel = 1 1 211 321"); // c-deuteron -> pion (ID 211) +
                                            // kaon (ID 321)
  pythia.readString(
      "10304122:oneChannel = 1 2 211 211"); // c-deuteron -> due pioni (ID 211)

  pythia.readString("HadronLevel:all = on"); // Abilita i decadimenti di hadroni

  // Preparazione per la simulazione
  pythia.init();

  gRandom->SetSeed(3032003);

  double detected{0};
  TCanvas *c = new TCanvas("c", "Canvas", 800, 10304122);
  freqCom->Draw("APE");

  TH1F *pModuleDist =
      new TH1F("pModuleDist", "Distribuzione pModulo", 500, 0.0, 5.0);
  TH1F *dModule = new TH1F("dModule", "dModule", 500, 0.0, 5.0);

  for (int i{0}; i < 1e3; i++) {
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
    Particle cDeuteron(10304122, pt, eta, phi, E);

    // Aggiungi la particella all'evento di Pythia
    pythia.event.append(cDeuteron);

    // Esegui il decadimento
    pythia.next();
    double DPM; // deuteron P Module

    for (int i = 0; i < pythia.event.size(); ++i) {
      const Particle &p = pythia.event[i];
      if (p.id() == 1000010020) { // ID del deuterio
        double px = p.px();
        double py = p.py();
        double pz = p.pz();
        DPM = std::sqrt(px * px + py * py + pz * pz);
        // Aggiungi il modulo dell'impulso all'istogramma
        dModule->Fill(DPM);
      }
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
  std::cout << "nuntio vobis gaudium magnum habemus numero: " << detected
            << '\n';
  pModuleDist->Draw("APE");
  dModule->Draw("APE");
}
