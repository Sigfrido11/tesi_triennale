//lo faccio compilato e riciclo come una vagabonda il vecchio codice

#include <TCanvas.h>
#include <algorithm>
#include <cmath>
#include <iterator>
#include "TFile.h"
#include "TH1.h"
#include "TRandom.h"
#include "particle.h"
#include "particletype.h"
#include "resonancetype.h"

int main() {
    TH1::AddDirectory(kFALSE);
    TFile *file = new TFile("histo.root", "RECREATE");
    TGraph *densFreq = (TGraph *)file->Get("densFreq");
    gRandom->SetSeed(3032003);
    Particle::AddParticleType("pion +", 0.13957, 1); 
    Particle::AddParticleType("deuteron", 0.49367, 0);
    Particle::AddParticleType("c-deuteron", 3,225000, 0, 0.050); //cerca la larghezza di risonanza

    const TString fileName[2] = {
        "generazioni/156_mev_4_fm/c-deuteron.dN.dy.dat",
        "generazioni/156_mev_4_fm/anti-c-deuteron.dN.dy.dat"};
    int count{0}; // Contatore per i valori non nulli nella seconda colonna
  
    for (int i{0}; i < n; i++) {
      std::ifstream inputFile(fileName[i]); // Apertura del file
      if (!inputFile) {
        std::cerr << "Impossibile aprire il file: " << fileName[i] << std::endl;
      }
      std::string line;
      while (std::getline(inputFile, line)) {
        std::istringstream ss(line);
        double col1, col2, col3;
        if (ss >> col1 >> col2 >> col3 && col2 != 0) {
        count++;
        }
      } 
}