// Giuseppe Luciano software analisi dati tesi sul c-deuteron

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TLegend.h"
#include "TH1.h"
#include <cmath> 
#include <fstream>
#include <iostream>
#include <sstream>

Double_t Pred(double *x, double *par) {
  double res = par[0] * std::pow(M_E, x[0] * par[1]);
  return res;
}

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
  TF1 *func = new TF1("pred", pred);
  s_t_Prediction->Fit(func);
  s_t_Prediction->GetYaxis()->SetMoreLogLabels();
  s_t_Prediction->GetYaxis()->SetNoExponent();
  gPad->SetLogy(1);

  auto canvas = new TCanvas("prediction", "prediction" , 200, 10, 600, 400);
  cavas->cd(0);
  auto leg = new TLegend(.6, .7, .9, .9);
  leg[i]->SetTextSize(0.04);
  leg[i]->SetBorderSize(0); // no border for legend
  leg[i]->SetFillColor(0);  // fill color is white
  leg[i]->AddEntry(h[i], legName[i], "p");
  leg[i]->AddEntry(f[i], "fit", "l");
  leg[i]->Draw();
}
