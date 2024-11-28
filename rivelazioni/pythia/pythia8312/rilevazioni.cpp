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

void rilevazioni() {
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

  TH1D *hefficiency = nullptr;
  TFile feff("efficiency.root");
  feff.GetObject("eff", hefficiency);

  // Inizializzazione di Pythia
  Pythia pythia;
  const int pid = 10304122;
  double const mass = 3,225000;
  double const charge = 1.0;
  double const spin = 1;
  double const width = 0.1;
  double const lifetime{}:
  double const color{};
  std::string name = "c-deuteron";
if(0){

  pythia.particleData.addParticle(pid, name, charge, mass, width, lifetime, spin, color);
  if (pythia.particleData.isParticle(pid)) {
    std::cout << "La particella è stata aggiunta!" << std::endl;
  }
}
  // Attivazione dei decadimenti per la particella
  //pythia.readString("10304122:all = on");
  pythia.readString("10304122:all = 2CDeuteron 2CDeuteronBar 1 3 0 3.226");

  // Modalità di decadimento:
   pythia.readString(
       "10304122:tau0 = 0.06");
   pythia.readString(
       "10s304122:oneChannel = 1 .1 0 1000010020 -311 211"); // questa l'ha inventata lui credo

  // 1. c-deuteron -> pion + kaon
   pythia.readString(
       "10304122:oneChannel = 1 0.023 1000010020 -311"); // deuterio e k0 bar

   pythia.readString(
       "10304122:oneChannel = 1 0.005 1000010020 -321 211"); // deuterio k- e pi+

   pythia.readString(
       "10304122:oneChannel = 1 0.0016 1000010020 892"); // deuterio e k*bar

   pythia.readString(
     "10304122:oneChannel = 1 0.0028 1000010020 -321 211"); // deuterio k- e pi+


   pythia.readString(
     "10304122:oneChannel = 1 0.0033 1000010020 -311 111"); // deuterio k0 bar- e pi0


   pythia.readString(
       "10304122:oneChannel = 1 0.0012 1000010020 -311 221"); // deuterio k0 bar
  //     \eta

   pythia.readString(
       "10304122:oneChannel = 1 0.0026 1000010020 -311 221 -221"); // deuterio
  //     kbar0
  //                                                            // p+ p-

   pythia.readString(
       "10304122:oneChannel = 1 0.0034 1000010020 -321 211 111"); // deuterio k-
  //     pi+
  //                                                           // pi0

  pythia.readString("HadronLevel:all = on"); // Abilita i decadimenti di hadroni
  pythia.readString(
      "ProcessLevel:all = off"); // Disabilita il processo inelastico iniziale

  // Preparazione per la simulazione
  pythia.init();

  gRandom->SetSeed(3032003);

  double detected{0};
  TCanvas *c = new TCanvas("c", "Canvas", 800, 10304122);
  freqCom->Draw("APE");

  TH1F *pModuleDist =
      new TH1F("pModuleDist", "Distribuzione pModulo", 500, 0.0, 5.0);
  TH1F *dModule = new TH1F("dModule", "dModule", 500, 0.0, 5.0);
  Particle cDeuteron;
  cDeuteron.id(pid);
  cDeuteron.status(11);
  cDeuteron.m(mass);
  cDeuteron.xProd(0);
  cDeuteron.yProd(0);
  cDeuteron.zProd(0);

  for (int i{0}; i < 10; i++) {
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
    double E = std::sqrt(pM * pM +
                         mass * mass); // dove 3.225000 massa c-deuteron
     // Crea un oggetto particella di pythia, attento
    cDeuteron.e(E);
    cDeuteron.px(px);
    cDeuteron.py(py);
    cDeuteron.pz(pz);
    pythia.event.reset();

    // Aggiungi la particella all'evento di Pythia
    pythia.event.append(cDeuteron);
    // Esegui il decadimento
    pythia.next();
    pythia.event.list();
    double DPM; // deuteron P Module

    bool isDetected = false;
    for (int i = 0; i < pythia.event.size(); ++i) {
      const Particle &p = pythia.event[i];
      if (p.id() == 1000010020) { // ID del deuterio
        double px = p.px();
        double py = p.py();
        double pz = p.pz();

        if (gRandom->Rndm() >=
            hefficiency->GetBinContent(
                hefficiency->GetXaxis()->FindBin(std::sqrt(px * px + py * py))))
          isDetected = true;
        DPM = std::sqrt(px * px + py * py + pz * pz);
        // Aggiungi il modulo dell'impulso all'istogramma
        dModule->Fill(DPM);
      }
    }
  }
  std::cout << "nuntio vobis gaudium magnum habemus numero: " << detected
            << '\n';
  pModuleDist->Draw("APE");
  dModule->Draw("APE");
}


int main(){
  rilevazioni();
}