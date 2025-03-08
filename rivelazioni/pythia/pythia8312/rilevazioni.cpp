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

  TCanvas *c = new TCanvas("c", "Canvas", 800, 10304122);
  freqCom->Draw("APE");

  TH1F *pModuleDist =
      new TH1F("pModuleDist", "Distribuzione pModulo", 500, 0.0, 5.0);
  TH1F *dModule = new TH1F("dModule", "dModule", 500, 0.0, 5.0);

  TH1D *hefficiency = nullptr;
  TFile feff("efficiency.root");
  feff.GetObject("eff", hefficiency);
  // Inizializzazione di Pythia
  Pythia pythia;

  const int cDeuteronId{10304122};
  double const mass{3.225000};
  int const charge{1};
  int const spin{1};
  double const width{0.4};
  double const lifetime{7.3e-13};
  double const color{0.0};

  /*
  dal .h di definizione
    void addParticle(int idIn, string nameIn = " ", int spinTypeIn = 0,
      int chargeTypeIn = 0, int colTypeIn = 0, double m0In = 0.,
      double mWidthIn = 0., double mMinIn = 0., double mMaxIn = 0.,
      double tau0In = 0., bool varWidthIn = false)
  */

  std::string name{"2CDeuteron"};
  pythia.particleData.addParticle(
      cDeuteronId, // Codice PDG della particella personalizzata
      "CDeuteron", // nome
      "anti-cd",
      spin,     // Spin
      charge,   // Carica (in unità di carica elementare)
      0,        // colore
      mass,     // massa
      width,    // larghezza
      mass,     // massa minima (0 vuol dire che non voglio massa variabile)
      mass,     // massa massima (0 vuol dire che non voglio massa variabile)
      lifetime, // vita media
      false);   // variazioni sulla larghezza

  // pythia.particleData.addParticle(1000100200, "20Ne", 6, 30, 0, 19.992440);

  // pythia.particleData.addParticle(cDeuteronId, name, charge, mass, width,
  // lifetime, spin, color);
  if (pythia.particleData.isParticle(cDeuteronId)) {
    std::cout << "La particella è stata aggiunta!" << std::endl;
  }

  // processi
  pythia.readString("HadronLevel:Decay = on"); // decadimenti adroni
  // pythia.readString("ParticleDecays:all = on"); //Attiva il decadimento di
  // tutte le particelle pythia.readString("10304122:Decay = on");
  pythia.readString(
      "ProcessLevel:all = off"); // Disabilita il processo inelastico iniziale

  //  id:addChannel = onMode bRatio meMode product1 product2 ....
  //  method  addChannel( branchingRatio, meMode, product1, ...)
  // adds a decay channel with up to 8 products.
  pythia.readString(
      "10304122:oneChannel = 1 0.023 0 -311 1000010020 "); // deuterio e k0 bar

  pythia.readString(
      "10304122:addChannel = 1 .8 0 1000010020 -311 211"); // questa l'ha
                                                           // inventata lui
                                                           // credo

  // 1. c-deuteron -> pion + kaon
  pythia.readString(
      "10304122:addChannel = 1 0.023 0 1000010020 -311"); // deuterio e k0 bar

  pythia.readString(
      "10304122:addChannel = 1 0.005 0 1000010020 -321 211"); // deuterio k- e
                                                              // pi+

  pythia.readString(
      "10304122:addChannel = 1 0.0016 0 1000010020 892"); // deuterio e k*bar

  pythia.readString(
      "10304122:addChannel = 1 0.0028 0 1000010020 -321 211"); // deuterio k- e
                                                               // pi+

  pythia.readString(
      "10304122:addChannel = 1 0.0033 0 1000010020 -311 111"); // deuterio k0
                                                               // bar- e pi0

  pythia.readString(
      "10304122:addChannel = 1 0.0012 0 1000010020 -311 221"); // deuterio k0
                                                               // bar
                                                               //     \eta

  pythia.readString(
      "10304122:addChannel = 1 0.0026 0 1000010020 -311 221 -221"); // deuterio
                                                                    //     kbar0
  //                                                            // p+ p-

  pythia.readString(
      "10304122:addChannel = 1 0.0034 0 1000010020 -321 211 111"); // deuterio
                                                                   // k-
  //     pi+
  //                                                           // pi0

  // pythia.readString("cDeuteronId:Decay = on");
  //  Preparazione per la simulazione
  pythia.init();

  gRandom->SetSeed(3032003);
  int detected{0};

  /*
  Particle cDeuteron;
  cDeuteron.id(cDeuteronId);
  cDeuteron.status(11);
  cDeuteron.m(mass);
  cDeuteron.xProd(0);
  cDeuteron.yProd(0);
  cDeuteron.zProd(0);
*/
  /*
    pythia.event.append(cDeuteron);
    if (!pythia.next()) {
      std::cerr << "Errore nell'elaborazione dell'evento di Pythia" <<
    std::endl; return;
    }

  int k{0};
  while (k<16) {

    freqCom->GetPointY(k);
     std::cout << "punto y " << randDist << '\n';
  }
*/

  int nDeuteri{0};
  for (int i{0}; i < 1e8; i++) {
    double phi = gRandom->Uniform(0., 2 * M_PI); // azimuth angle
    double theta = gRandom->Uniform(0., M_PI);   // polar angle
    double randDist = gRandom->Rndm();
    double pM;
    int j{1};
    while (true) {
      if (randDist < freqCom->GetPointY(j)) {
        pM = freqCom->GetPointX(j);
        break;
      }
      j++;
    }
    pModuleDist->Fill(pM);
    double px = pM * std::cos(phi) * std::sin(theta);
    double py = pM * std::sin(phi) * std::sin(theta);
    double pz = pM * std::cos(theta);
    double E =
        std::sqrt(pM * pM + mass * mass); // dove 3.225000 massa c-deuteron
    // Crea un oggetto particella di pythia, attento

    Particle cDeuteron(cDeuteronId);
    cDeuteron.m(mass);
    cDeuteron.status(11);
    cDeuteron.e(E);
    cDeuteron.px(px);
    cDeuteron.py(py);
    cDeuteron.pz(pz);
    pythia.event.reset();

    // Aggiungi la particella all'evento di Pythia
    pythia.event.append(cDeuteron);
    // Esegui il decadimento
    pythia.next();
    // pythia.event.list();
    for (int i = 0; i < pythia.event.size(); ++i) {
      const Particle &p = pythia.event[i];
      if (p.id() == 1000010020) { // ID del deuterio
        nDeuteri++;
        double px = p.px();
        double py = p.py();
        double pz = p.pz();
        double const p_perp = std::sqrt(px * px + py * py);
        double const random = gRandom->Rndm();
        double const eff = hefficiency->GetBinContent(
            hefficiency->GetXaxis()->FindBin(p_perp));
        bool const passed = random <= eff;
        if (passed)
        detected++;
        double const DPM = std::sqrt(px * px + py * py + pz * pz);
        // Aggiungi il modulo dell'impulso all'istogramma
        dModule->Fill(DPM);
      }
    }
  }

  std::cout << "numero di deuteri " << nDeuteri << '\n';
  std::cout << "nuntio vobis gaudium magnum habemus numero: " << detected
            << '\n';
  pModuleDist->Draw("APE");
  dModule->Draw("APE");

  TFile *output = new TFile("result.root", "READ");
  pModuleDist->Write();
  dModule->Write();
}

int main() { rilevazioni(); }