// Giuseppe Luciano software analisi dati tesi sul c-deuteron

/*
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
*/
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

Double_t Pred(double *x, double *par) {
  double res = par[0] * std::pow(M_E, x[0] * par[1]);
  return res;
}

void grafici() {
  const int n{14};
  TGraphErrors *graph[n];

  const TString fileName[n] = {
      "generazioni/156_mev_5_fm/p.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-p.dN.dy.dat",
      "generazioni/156_mev_5_fm/n.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-n.dN.dy.dat",
      "generazioni/156_mev_5_fm/d.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-d.dN.dy.dat",
      "generazioni/156_mev_5_fm/H3.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-H3.dN.dy.dat",
      "generazioni/156_mev_5_fm/He3.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-He3.dN.dy.dat",
      "generazioni/156_mev_5_fm/He4.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-He4.dN.dy.dat",
      "generazioni/156_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-c-deuteron.dN.dy.dat"};
  //  "generazioni/156_mev_5_fm/HyperTriton.dN.dy.dat","generazioni/156_mev_5_fm/anti-HyperTriton.dN.dy.dat"};

  TH1F *histo[n];
  TH1F *histoSum[n / 2];

  double mass[n / 2] = {0.938272, 0.939565, 1.87561, 2.80892,
                        2.80839,  3.72738,  3.225}; // m HyperTriton=2.99114
  double integral[n / 2];
  double error[n / 2];

  for (int i{0}; i < n; i = i + 2) {
    graph[i] = new TGraphErrors(fileName[i], "%lg %lg %lg");
    graph[i + 1] = new TGraphErrors(fileName[i + 1], "%lg %lg %lg");
    for (int j{0}; j < graph[i]->GetN(); j++) {
      graph[i]->SetPoint(j, graph[i]->GetPointX(j),
                         graph[i]->GetPointY(j) + graph[i + 1]->GetPointY(j));
      // errore calcolato con la somma in quadratura
      double errY =
          TMath::Sqrt(graph[i]->GetErrorY(j) * graph[i]->GetErrorY(j) +
                      graph[i + 1]->GetErrorY(j) * graph[i + 1]->GetErrorY(j));
      graph[i]->SetPointError(j, 0., errY);
    }
    double partialArea{0};
    double integralError{0};
    for (int k = 0; k < graph[i]->GetN() - 1; ++k) {
      // Estrai i punti x e y
      double x1, y1, x2, y2;
      graph[i]->GetPoint(k, x1, y1);
      graph[i]->GetPoint(k + 1, x2, y2);

      // Applica la formula del trapezio per l'area tra i punti (x1, y1) e (x2,
      // y2)
      double area = 0.5 * (y1 + y2) * (x2 - x1);
      partialArea += area;

      // Estrai gli errori sui punti y
      double y1Error = graph[i]->GetErrorY(k);
      double y2Error = graph[i]->GetErrorY(k + 1);

      // Propagazione degli errori (considerando la somma in quadratura)
      double areaError =
          0.5 * std::sqrt(y1Error * y1Error + y2Error * y2Error) * (x2 - x1);
      integralError +=
          areaError * areaError; // Somma in quadratura degli errori
    }

    error[i / 2] = std::sqrt(integralError); // Errore finale
    integral[i / 2] = partialArea;
  }
  auto canvas = new TCanvas("prediction", "prediction", 200, 10, 600, 400);
  canvas->cd(0);
  canvas->SetLogy();
  TGraphErrors *finalGraph = new TGraphErrors(n / 2, mass, integral, error);
  finalGraph->SetLineColor(1);
  finalGraph->GetYaxis()->SetTitleOffset(1.2);
  finalGraph->GetXaxis()->SetTitleSize(0.04);
  finalGraph->GetYaxis()->SetTitleSize(0.04);
  finalGraph->GetXaxis()->SetTitle("Mass (Gev)");
  finalGraph->GetYaxis()->SetTitle("dN/dy");
  finalGraph->GetXaxis()->CenterTitle(true);
  finalGraph->GetXaxis()->CenterTitle(true);
  finalGraph->SetMarkerStyle(21); // Stile 21 è un cerchio
  finalGraph->SetMarkerSize(1.5); // Cambia la dimensione del marker
  // finalGraph->SetMarkerColor(kRed); // Cambia il colore del marker, ad
  // esempio
  //  rosso
  // finalGraph->GetYaxis()->SetLimits(10, 1e-10);
  // finalGraph->GetYaxis()->SetRangeUser(10, 1e-10);

  const char *labels[n / 2] = {"p", "n", "d", "H3", "He3", "He4", "c-deuteron"};
  finalGraph->Draw("AP");
  for (int i = 0; i < finalGraph->GetN(); ++i) {
    double xLabel = finalGraph->GetPointX(i);
    double yLabel =
        finalGraph->GetPointY(i) * 2.1; // Posiziona l'etichetta sopra il punto
    TLatex *label = new TLatex(xLabel, yLabel, labels[i]);
    label->SetTextSize(.02);    // Dimensione del testo
    label->SetNDC(kFALSE);      // Coordinate del testo non normalizzate
    label->SetTextColor(kBlue); // Colore delle etichette
    label->Draw();              // Disegna l'etichetta
  }
  auto leg1 = new TLegend(0.1, 0.7, 0.3, 0.9);
  leg1->AddEntry(finalGraph, "Temperature 156 Mev", "");
  leg1->AddEntry(finalGraph, "Radius 5 fm", "");
  leg1->AddEntry(finalGraph, "|y|<0.5", "");
  leg1->SetTextSize(0.04);
  leg1->Draw();
  // Aggiorna la canvas per visualizzare i grafici
  gPad->SetLogy(1);

  canvas->Update();
}
void cambiamenti() {
  int const n{14};
  const TString fileName[n] = {
      "generazioni/156_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_3,5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_3,5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/154_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/154_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/158_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/158_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/152_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/152_mev_5_fm/anti-c-deuteron.dN.dy.dat",
  };

  double count[n / 2];

  for (int i{0}; i < n; i = i + 2) {
    std::ifstream inputFile1(fileName[i]); // Apertura del file
    std::ifstream inputFile2(fileName[i + 1]);
    if (!inputFile1) {
      std::cerr << "Impossibile aprire il file: " << fileName[i] << std::endl;
    }
    if (!inputFile2) {
      std::cerr << "Impossibile aprire il file: " << fileName[i + 1]
                << std::endl;
    }

    std::string line;
    int parzCount{0}; // Contatore per i valori non nulli nella seconda colonna

    while (std::getline(inputFile1, line)) {
      std::istringstream ss(line);
      double col1, col2, col3;
      if (ss >> col1 >> col2 >> col3 && col2 != 0) {
        parzCount++;
      }
    }
    while (std::getline(inputFile2, line)) {
      std::istringstream ss(line);
      double col1, col2, col3;
      if (ss >> col1 >> col2 >> col3 && col2 != 0) {
        parzCount++;
      }
    }
    if (i == 0) {
      count[i / 2] = parzCount; // / 2.5e7;
      std::cout << count[i / 2] << " " << parzCount << "\n";
    } else {
      count[i / 2] = parzCount; // / 2e7;
      std::cout << count[i / 2] << " " << parzCount << "\n";
    }
  }

  double volume[n / 2] = {8, 6, 5, 3.5, 5, 5, 5};
  double temperature[n / 2] = {156, 156, 156, 156, 154, 158, 152};

  // Crea un oggetto TGraph2D
  TGraph2D *d2Graph = new TGraph2D(n, volume, temperature, count);

  double a[4] = {8, 6, 5, 3.5};
  double b[4] = {count[0], count[1], count[2], count[3]};
  double c[4] = {156, 154, 158, 152};
  double d[4] = {count[2], count[4], count[5], count[6]};
  TGraph *diffVolume = new TGraph(4, a, b);
  TGraph *diffTemp = new TGraph(3, c, d);
  // Crea un canvas
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

  c1->SetTheta(20); // Imposta l'angolo di vista in azimut
  c1->SetPhi(20);   // Imposta l'angolo di vista in elevazione

  // Disegna il grafico 3D
  d2Graph->SetMarkerStyle(21);    // Marker stile 21 (cerchi pieni)
  d2Graph->SetMarkerSize(1.5);    // Aumenta la dimensione del marker
  d2Graph->SetMarkerColor(kBlue); // Colore del marker blu

  // Aggiungi una griglia per facilitare la lettura dei valori

  // Crea una zona grande in alto (80% dell'altezza)
  TPad *pad1 =
      new TPad("pad1", "Pad grande", 0, 0.33, 1, 1); // xlow, ylow, xup, yup
  pad1->SetBottomMargin(0.02); // Riduci il margine inferiore del pad1
  pad1->Draw();
  // Crea due zone piccole in basso (50% ciascuna della larghezza e 33%
  // dell'altezza)
  TPad *pad2 = new TPad("pad2", "Pad piccolo sinistro", 0, 0, 0.5,
                        0.33); // Parte sinistra in basso
  pad2->SetTopMargin(0.02);    // Riduci il margine superiore del pad2
  pad2->Draw();

  TPad *pad3 = new TPad("pad3", "Pad piccolo destro", 0.5, 0, 1,
                        0.33); // Parte destra in basso
  pad3->SetTopMargin(0.02);    // Riduci il margine superiore del pad3
  pad3->Draw();

  // Seleziona ciascun pad e disegna qualcosa
  pad1->cd();
  d2Graph->SetLineColor(1);
  d2Graph->GetYaxis()->SetTitleOffset(1.2);
  d2Graph->GetXaxis()->SetTitleSize(0.04);
  d2Graph->GetYaxis()->SetTitleSize(0.04);
  d2Graph->GetXaxis()->SetTitle("Radius (fm)");
  d2Graph->GetYaxis()->SetTitle("Temperature (Mev)");
  d2Graph->GetZaxis()->SetTitle("Frequency");
  d2Graph->GetZaxis()->CenterTitle(true);
  d2Graph->GetXaxis()->CenterTitle(true);
  d2Graph->GetXaxis()->CenterTitle(true);
  d2Graph->Draw("P0");
  // Aggiungi etichette ai punti per maggiore chiarezza
  /*
    char *labels[n/2] = {"(156 Mev, 8 fm)",   "(156 Mev, 6 fm)", "(156 Mev, 5
    fm)",
                      "(156 Mev, 3.5 fm)", "(154 Mev, 5 fm)", "(158 Mev, 5 fm)",
                      "(152 Mev, 5 fm)"};
    double *xVals = d2Graph->GetX();
    double *yVals = d2Graph->GetY();
    double *zVals = d2Graph->GetZ();

    for (int i = 0; i < d2Graph->GetN(); ++i) {
      TString label = TString::Format("Point %d", i);
      double x = xVals[i];
      double y = yVals[i];
      double z = zVals[i];
      double x2d, y2d;
      TPad *pad = (TPad*)gPad;
      pad->GetCoord(x, y, z, x2d, y2d);
      TLatex *label = new TLatex(xLabel, yLabel, labels[i]);
      label->SetTextSize(.02); // Dimensione del testo
      label->SetNDC(kFALSE);   // Coordinate del testo non normalizzate
      label->Draw();           // Disegna l'etichetta
    }
  */
  pad2->cd();
  diffVolume->SetLineColor(1);
  diffVolume->GetYaxis()->SetTitleOffset(1.2);
  diffVolume->GetXaxis()->SetTitleSize(0.04);
  diffVolume->GetYaxis()->SetTitleSize(0.04);
  diffVolume->GetXaxis()->SetTitle("Radius (fm)");
  diffVolume->GetYaxis()->SetTitle("Frequency");
  diffVolume->GetXaxis()->CenterTitle(true);
  diffVolume->GetXaxis()->CenterTitle(true);
  diffVolume->Draw("P0");
  diffVolume->SetMarkerStyle(20);     // Marker stile 21 (cerchi pieni)
  diffVolume->SetMarkerSize(.8);      // Aumenta la dimensione del marker
  diffVolume->SetMarkerColor(kGreen); // Colore del marker blu
  auto leg1 = new TLegend(0.1, 0.7, 0.3, 0.9);
  leg1->AddEntry(diffVolume, "fixed temperature 156 Mev", "");
  leg1->SetTextSize(0.04);
  diffVolume->Draw("APE");
  TF1 *fitVolume = new TF1("fitVolume", "[0]*x + [1]", 0, 10);
  fitVolume->SetLineColor(kRed); // Colore della linea di fit
  diffVolume->Fit(fitVolume);
  double slopeVolume = fitVolume->GetParameter(0);     // [0] = pendenza
  double interceptVolume = fitVolume->GetParameter(1); // [1] = intercetta
  double slopeErrorVolume = fitVolume->GetParError(0); // Errore sulla pendenza
  double interceptErrorVolume =
      fitVolume->GetParError(1); // Errore sull'intercetta
  TString fitInfoVolume =
      TString::Format("Slope: %.2f ± %.2f", slopeVolume, slopeErrorVolume);
  TString fitInfoVolume1 = TString::Format(
      "Intercept: %.2f ± %.2f", interceptVolume, interceptErrorVolume);
  leg1->AddEntry((TObject *)0, fitInfoVolume, "");
  leg1->AddEntry((TObject *)0, fitInfoVolume1, "");
  leg1->AddEntry(fitVolume, "Fit Line", "l");
  leg1->Draw();

  pad3->cd();
  diffTemp->SetLineColor(1);
  diffTemp->GetYaxis()->SetTitleOffset(1.2);
  diffTemp->GetXaxis()->SetTitleSize(0.04);
  diffTemp->GetYaxis()->SetTitleSize(0.04);
  diffTemp->GetXaxis()->SetTitle("Temperaure (Mev)");
  diffTemp->GetYaxis()->SetTitle("Frequency");
  diffTemp->GetXaxis()->CenterTitle(true);
  diffTemp->GetXaxis()->CenterTitle(true);
  diffTemp->GetXaxis()->SetRangeUser(150, 160);
  diffTemp->GetXaxis()->SetRange(150, 160);
  auto leg2 = new TLegend(0.1, 0.7, 0.3, 0.9);
  leg2->AddEntry(diffTemp, "fixed radius r=5fm", "");
  leg2->SetTextSize(0.04);
  diffTemp->Draw("P0");
  diffTemp->SetMarkerStyle(20);   // Marker stile 21 (cerchi pieni)
  diffTemp->SetMarkerSize(.8);    // Aumenta la dimensione del marker
  diffTemp->SetMarkerColor(kRed); // Colore del marker blu
  diffTemp->Draw("APE");
  TF1 *fitTemp = new TF1("fitTemp", "[0]*x + [1]", 0, 10);
  fitTemp->SetLineColor(kGreen); // Colore della linea di fit
  diffTemp->Fit(fitTemp);
  double slopeTemp = fitTemp->GetParameter(0);         // [0] = pendenza
  double interceptTemp = fitTemp->GetParameter(1);     // [1] = intercetta
  double slopeErrorTemp = fitTemp->GetParError(0);     // Errore sulla pendenza
  double interceptErrorTemp = fitTemp->GetParError(1); // Errore sull'intercetta
  TString fitInfoTemp =
      TString::Format("Slope: %.2f ± %.2f", slopeTemp, slopeErrorTemp);
  TString fitInfoTemp2 = TString::Format("Intercept: %.2f ± %.2f",
                                         interceptTemp, interceptErrorTemp);

  leg2->AddEntry((TObject *)0, fitInfoTemp, "");
  leg2->AddEntry((TObject *)0, fitInfoTemp2, "");
  leg2->AddEntry(fitTemp, "Fit Line", "l");
  leg2->Draw();

  c1->Update();
}
