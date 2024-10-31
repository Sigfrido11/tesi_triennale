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
  TF1 *fitExp = new TF1("fitExp", "[0] * exp(-[1] * x)", 0, 4);
  fitExp->SetParameter(0, 1e3);
  fitExp->SetParLimits(0, 0, 1e5);
  fitExp->SetParLimits(1, 4, 8);
  finalGraph->Fit(fitExp);

  finalGraph->SetLineColor(1);
  finalGraph->GetYaxis()->SetTitleOffset(1.2);
  finalGraph->GetXaxis()->SetTitleSize(0.04);
  finalGraph->GetYaxis()->SetTitleSize(0.04);
  finalGraph->GetXaxis()->SetTitle("Mass (Gev)");
  finalGraph->GetYaxis()->SetTitle("dN/dy");
  finalGraph->GetXaxis()->CenterTitle(true);
  finalGraph->GetXaxis()->CenterTitle(true);
  finalGraph->SetMarkerStyle(21); // Stile 21 è un cerchio
  finalGraph->SetMarkerSize(1.);  // Cambia la dimensione del marker
  finalGraph->GetYaxis()->SetLimits(1e-7, 10);
  finalGraph->GetYaxis()->SetRangeUser(1e-7, 10);

  const char *labels[n / 2] = {"p", "n", "d", "H3", "He3", "He4", "c-deuteron"};
  finalGraph->Draw("AP");
  for (int i = 0; i < finalGraph->GetN(); ++i) {
    double xLabel = finalGraph->GetPointX(i);
    double yLabel;
    TLatex *label;
    if (i == 1 or i == 4) {
      yLabel = finalGraph->GetPointY(i) /
               4.1; // Posiziona l'etichetta sopra il punto
      label = new TLatex(xLabel, yLabel, labels[i]);
    } else {
      yLabel = finalGraph->GetPointY(i) *
               2.1; // Posiziona l'etichetta sopra il punto
      label = new TLatex(xLabel, yLabel, labels[i]);
    }
    label->SetTextSize(.05);    // Dimensione del testo
    label->SetNDC(kFALSE);      // Coordinate del testo non normalizzate
    label->SetTextColor(kBlue); // Colore delle etichette
    label->Draw();              // Disegna l'etichetta
  }
  // fitExp->Draw();
  auto leg1 = new TLegend(0.7, 0.1, 0.9, 0.3);
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
  int const n{22}; // 24};
  // se aggiungi un file le cose da cambiare sono qui
  const TString fileName[n] = {
      "generazioni/156_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6,5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6,5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5,5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5,5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_4,5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_4,5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_3,5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_3,5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/160_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/160_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/158_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/158_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/154_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/154_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/152_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/152_mev_5_fm/anti-c-deuteron.dN.dy.dat",
  };

  /*
  "generazioni/156_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6,5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6,5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_6_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5,5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5,5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_4,5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_4,5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_4_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_4_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_3,5_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_3,5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/160_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/160_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/158_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/158_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/154_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/154_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/152_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/152_mev_5_fm/anti-c-deuteron.dN.dy.dat",
  };
  */

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
    count[i / 2] = parzCount / 2e7;
    std::cout << count[i / 2] << " " << parzCount << "\n";
  }

  // se aggiungi un file le cose da cambiare sono qui

  double volume[n / 2] = {
      8,   6.5, 6, 5.5, 5, 4.5,
      3.5, 5,   5, 5,   5}; //{8, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 5, 5, 5};
  double temperature[n / 2] = {156, 156, 156, 156, 156, 156,
                               156, 160, 158, 154, 152};
  //{156, 156, 156, 156, 156, 156, 156, 156, 160, 158, 154, 152};

  TGraph2D *d2Graph = new TGraph2D(n / 2, volume, temperature, count);
  double a[7] = {8, 6.5, 6,  5.5,
                 5, 4.5, 3.5}; //{8, 6.5, 6, 5.5, 5, 4.5, 4, 3.5};
  double b[7] = {count[0], count[1], count[2], count[3],
                 count[4], count[5], count[6]}; // count[7]
  double c[5] = {160, 158, 156, 154, 152};
  double d[5] = {count[7], count[8], count[4], count[9], count[10]};
  // count[8], count[9], count[4], count[10],count[11]};

  // se aggiungi un file le cose da cambiare sono qui nel numero di punti
  TGraph *diffVolume = new TGraph(7, a, b);
  TGraph *diffTemp = new TGraph(5, c, d);
  // Crea un canvas
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

  c1->SetTheta(15); // Imposta l'angolo di vista in azimut
  c1->SetPhi(20);   // Imposta l'angolo di vista in elevazione

  // Disegna il grafico 3D
  d2Graph->SetMarkerStyle(21);    // Marker stile 21 (cerchi pieni)
  d2Graph->SetMarkerSize(1.5);    // Aumenta la dimensione del marker
  d2Graph->SetMarkerColor(kBlue); // Colore del marker blu

  // Crea il primo canvas per il pad1
  TPad *pad1 = new TPad("pad1", "Pad grande", 0, 0, 1, 1);
  // pad1->SetBottomMargin(0.05); // Riduci il margine inferiore del pad1
  pad1->Draw();
  pad1->cd();

  d2Graph->SetLineColor(1);
  d2Graph->GetYaxis()->SetTitleOffset(1.2);
  d2Graph->GetXaxis()->SetTitleSize(0.04);
  d2Graph->GetYaxis()->SetTitleSize(0.04);
  d2Graph->GetXaxis()->SetTitle("Radius (fm)");
  d2Graph->GetYaxis()->SetTitle("Temperature (Mev)");
  d2Graph->GetXaxis()->SetLimits(1.5, 10);
  d2Graph->GetXaxis()->SetRangeUser(1.5, 10);
  d2Graph->GetYaxis()->SetLimits(150, 165);
  d2Graph->GetYaxis()->SetRangeUser(150, 165);
  d2Graph->Draw("P0");

  c1->Update(); // Aggiorna il canvas

  // Crea il secondo canvas per il pad2
  TCanvas *c2 = new TCanvas("c2", "Canvas 2", 800, 600);
  TPad *pad2 = new TPad("pad2", "Pad piccolo sinistro", 0, 0, 1, 1);
  pad2->Draw();
  pad2->cd();

  // Disegna il grafico diffVolume
  diffVolume->SetLineColor(1);
  diffVolume->GetYaxis()->SetTitleOffset(1.2);
  diffVolume->GetXaxis()->SetTitleSize(0.04);
  diffVolume->GetYaxis()->SetTitleSize(0.04);
  diffVolume->GetXaxis()->SetTitle("Radius (fm)");
  diffVolume->GetYaxis()->SetTitle("Frequency");
  diffVolume->Draw("P0");
  diffVolume->SetMarkerStyle(20);
  diffVolume->SetMarkerSize(.8);
  diffVolume->SetMarkerColor(kGreen);

  auto leg1 = new TLegend(0.7, 0.1, 0.9, 0.3);
  leg1->AddEntry(diffVolume, "T = 156 Mev", "");
  leg1->SetTextSize(0.04);
  leg1->Draw();

  // (Opzionale) Aggiungi la linea di fit
  TF1 *fitVolume = new TF1("fitVolume", "[0] * log([1] * x) + [2]", 0, 10);
  fitVolume->SetLineColor(kRed);
  fitVolume->SetParameter(0, 1e-6);
  // fitVolume->SetParLimits(0, 0, 1e-5);
  fitVolume->SetParameter(1, 1e3);
  // fitVolume->SetParLimits(1, 0, 1e5);
  diffVolume->Fit(fitVolume);
  leg1->AddEntry(fitVolume, "Fit Line", "l");
  leg1->Draw();
  diffVolume->Draw("APE");

  c2->Update(); // Aggiorna il canvas

  // Crea il terzo canvas per il pad3
  TCanvas *c3 = new TCanvas("c3", "Canvas 3", 800, 600);
  TPad *pad3 = new TPad("pad3", "Pad piccolo destro", 0, 0, 1, 1);
  pad3->Draw();
  pad3->cd();

  // Disegna il grafico diffTemp
  diffTemp->SetLineColor(1);
  diffTemp->GetYaxis()->SetTitleOffset(1.2);
  diffTemp->GetXaxis()->SetTitleSize(0.04);
  diffTemp->GetYaxis()->SetTitleSize(0.04);
  diffTemp->GetXaxis()->SetTitle("Temperature (Mev)");
  diffTemp->GetYaxis()->SetTitle("Frequency");
  diffTemp->Draw("P0");
  diffTemp->SetMarkerStyle(20);
  diffTemp->SetMarkerSize(.8);
  diffTemp->SetMarkerColor(kRed);
  auto leg2 = new TLegend(0.7, 0.1, 0.9, 0.3);
  leg2->AddEntry(diffTemp, "R=5fm", "");
  leg2->SetTextSize(0.04);
  leg2->Draw();

  // (Opzionale) Aggiungi la linea di fit
  TF1 *fitTemp = new TF1("fitTemp", "[0]*x + [1]", 0, 10);
  fitTemp->SetLineColor(kGreen);
  fitTemp->SetParameter(0, 1e-6);
  // fitTemp->SetParLimits(0, 0, 1e-5);
  fitTemp->SetParameter(1, 1e3);
  // fitTemp->SetParLimits(1, 0, 1e5);
  diffTemp->Fit(fitTemp);
  leg2->AddEntry(fitTemp, "Fit Line", "l");
  diffTemp->Draw("APE");
  leg2->Draw();

  c3->Update(); // Aggiorna il canvas
}
