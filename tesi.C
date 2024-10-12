// Giuseppe Luciano software analisi dati tesi sul c-deuteron

#include "TH1.h"
#include <cmath> // Per std::abs
#include <fstream>
#include <iostream>
#include <sstream>

int main() {
  const int n{0};
  const TString fileName[n] = {"a.txt", "b.txt", "c.txt"};
  std::ifstream *file[n];
  std::string line;
  double mean[n] = (0);
  double error[n] = {0};
  for (int i{0}; i < n; ++i) {
    // Apertura file
    file[i] = new std::ifstream(fileName[i].Data());
    // Verifica se il file è aperto correttamente
    if (!file[i]->is_open()) {
      std::cerr << "Errore nell'apertura del file: " << fileName[i]
                << std::endl;
    }
    // Lettura del file riga per riga

    while (std::getline(*(file[i]), line)) {
      std::stringstream ss(line);
      double col1, col2, col3;
      if (ss >> col1 >> col2 >> col3) {
        // Controllo se la prima colonna è minore di 0.5 in valore assoluto
        if (std::abs(col1) <= 0.5) {
          mean[i] += col2 * col3;
          error[i] += col3;
        }
      }
      (*file)[i].close();
    }
  }
  TGraphErrors *s_t_Prediction = new TGraphErrors("STPredictions", mean, error);
  s_t_Prediction->SetLineColor(1);
  s_t_Prediction->GetYaxis()->SetTitleOffset(1.2);
  s_t_Prediction->GetXaxis()->SetTitleSize(0.04);
  s_t_Prediction->GetYaxis()->SetTitleSize(0.04);
  s_t_Prediction->GetXaxis()->SetTitle("Mass (Gev)");
  s_t_Prediction->GetYaxis()->SetTitle("dN/dy");
  s_t_Prediction->GetXaxis()->CenterTitle(true);
  s_t_Prediction->GetXaxis()->CenterTitle(true);

  TF1* func = new TF1("pred", pred);
  s_t_Prediction->Fit("func");
}
