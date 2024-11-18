// Giuseppe Luciano software analisi dati tesi sul c-deuteron

#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TLegend.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>

void grafici() {
  const int n{22};
  TGraphErrors *graph[n];

  const TString fileName[n] = {
      "generazioni/chemical_charm_30/156_mev_5_fm//p.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-p.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/n.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-n.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/d.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-d.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/H3.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-H3.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/He3.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-He3.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/He4.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-He4.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/Lambda.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-Lambda.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/Lambda(c)+.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-Lambda(c)+.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/LambdaNeutron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-LambdaNeutron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/LambdaProton.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-LambdaProton.dN.dy.dat"};

  TH1F *histo[n];
  TH1F *histoSum[n / 2];

  double mass[8] = {0.938272, 0.939565, 1.87561,  2.80892,
                    2.80839,  3.72738,  2.057000, 2.056000};
  double massLambda[3] = {3.225, 1.800000,
                          2.286460}; // massa c-deuteron lamda lambda c
  double integral[8];
  double error[8];
  double integralLambda[3];
  double errorLambda[3];

  for (int i{0}; i < n; i = i + 2) {
    graph[i] = new TGraphErrors(fileName[i], "%lg %lg %lg");
    graph[i + 1] = new TGraphErrors(fileName[i + 1], "%lg %lg %lg");
    for (int j{0}; j < graph[i]->GetN(); j++) {
      double x = graph[i]->GetPointX(j);
      double y1 = graph[i]->GetPointY(j);
      double y2 = graph[i + 1]->GetPointY(j);
      double err1 = graph[i]->GetErrorY(j);
      double err2 = graph[i + 1]->GetErrorY(j);
      graph[i]->SetPoint(j, x, y1 + y2);
      graph[i]->SetPointError(j, 0., TMath::Sqrt(err1 * err1 + err2 * err2));
    }
    double partialArea{0};
    double integralError{0};
    for (int k = 0; k < graph[i]->GetN() - 1; ++k) {
      // Estrai i punti x e y
      double x1, y1, x2, y2;
      graph[i]->GetPoint(k, x1, y1);
      if (std::abs(x1) > 0.5)
        continue;
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
    if (i < 16) {                              // Per i primi 16 grafici
      error[i / 2] = std::sqrt(integralError); // Errore finale
      integral[i / 2] = partialArea;
    } else { // Per i grafici successivi
      int lambdaIndex =
          (i - 16) / 2; // Indice corretto per integralLambda e errorLambda
      errorLambda[lambdaIndex] = std::sqrt(integralError); // Errore finale
      integralLambda[lambdaIndex] = partialArea;
    }
  }

  for (int i = 0; i < 8; ++i) {
    std::cout << "integral[" << i << "] = " << integral[i] << std::endl;
  }

  std::cout << "\nArray integralLambda:" << std::endl;
  for (int i = 0; i < 3; ++i) {
    std::cout << "integralLambda[" << i << "] = " << integralLambda[i]
              << std::endl;
  }

  TMultiGraph *finalGraph = new TMultiGraph;

  // Creazione dei grafici con i dati di input
  TGraphErrors *gNormal = new TGraphErrors(8, mass, integral, error);
  TGraphErrors *gLambda =
      new TGraphErrors(3, massLambda, integralLambda, errorLambda);

  // Configurazione della funzione di fit per entrambi i grafici
  TF1 *fitExp = new TF1("fitExp", "[0] * exp(-[1] * x)", 0, 4);
  TF1 *fitExpLambda = new TF1("fitExpLambda", "[0] * exp(-[1] * x)", 0, 4);

  fitExp->SetParameter(0, 1e3);
  fitExp->SetParLimits(0, 0, 1e5);
  fitExp->SetParameter(1, 5);
  fitExp->SetParLimits(1, 4, 8);

  fitExpLambda->SetParameter(0, 1e3);
  fitExpLambda->SetParLimits(0, 0, 1e5);
  fitExpLambda->SetParameter(1, 5);
  fitExpLambda->SetParLimits(1, 4, 8);

  // Configurazione dei marker per i grafici
  gNormal->SetMarkerStyle(21);
  gNormal->SetMarkerSize(1.0);
  gNormal->SetTitle("Particles; Mass (GeV); dN/dy");

  gLambda->SetMarkerStyle(21);
  gLambda->SetMarkerSize(1.0);
  gLambda->SetTitle("Lambda Particles; Mass (GeV); dN/dy");

  // Canvas per gNormal
  auto canvasNorm = new TCanvas("canvasNorm", "gNormal - Particles", 800, 600);
  canvasNorm->cd();
  canvasNorm->SetLogy();
  gNormal->GetYaxis()->SetLimits(1e-7, 10);
  gNormal->GetXaxis()->SetLimits(0, 6);
  gNormal->Draw("APE");

  // Canvas per gLambda
  auto canvasLambda =
      new TCanvas("canvasLambda", "gLambda - Particles", 800, 600);
  canvasLambda->cd();
  canvasLambda->SetLogy();
  gLambda->GetYaxis()->SetLimits(1e-7, 10);
  gLambda->GetXaxis()->SetLimits(0, 6);
  gLambda->Draw("APE");

  // Canvas per il grafico combinato
  auto canvasFinal = new TCanvas("canvasFinal", "Combined Graph", 800, 600);
  canvasFinal->cd();
  canvasFinal->SetLogy();
  finalGraph->Add(gNormal, "P");
  finalGraph->Add(gLambda, "P");
  finalGraph->SetTitle("Combined Graph; Mass (GeV); dN/dy");
  finalGraph->Draw("APE");

  // Aggiunta delle etichette ai punti
  const char *labels[n / 2] = {
      "p",          "n",      "d",         "H3", "He3", "He4", "Lambda-Neutron",
      "c-deuteron", "Lambda", "Lambda(c)+"};

  for (int i = 0; i < gNormal->GetN(); ++i) {
    double x, y;
    gNormal->GetPoint(i, x, y);
    TLatex *label = new TLatex(x, y * 1.5, labels[i]);
    label->SetTextSize(0.03);
    label->Draw();
  }

  for (int i = 0; i < gLambda->GetN(); ++i) {
    double x, y;
    gLambda->GetPoint(i, x, y);
    TLatex *label = new TLatex(x, y * 1.5, labels[i + 8]);
    label->SetTextSize(0.03);
    label->Draw();
  }

  // Aggiunta della legenda
  auto leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry(gNormal, "Normal Particles", "P");
  leg->AddEntry(gLambda, "Lambda Particles", "P");
  leg->Draw();

  // Configurazione del grafico finale
  finalGraph->GetYaxis()->SetTitleOffset(1.2);
  finalGraph->GetXaxis()->SetTitleSize(0.04);
  finalGraph->GetYaxis()->SetTitleSize(0.04);
  finalGraph->GetXaxis()->SetTitle("Mass (GeV)");
  finalGraph->GetYaxis()->SetTitle("dN/dy");
  finalGraph->GetXaxis()->CenterTitle(true);
  finalGraph->GetYaxis()->CenterTitle(true);
  finalGraph->GetYaxis()->SetLimits(1e-7, 10);
  finalGraph->GetYaxis()->SetRangeUser(1e-7, 10);

  // Aggiornamento della canvas per visualizzare i risultati
  canvasFinal->Update();
}

void primo_grafico() {
  const int n{22};
  TGraphErrors *graph[n];

  const TString fileName[n] = {
      "generazioni/chemical_charm_30/156_mev_5_fm//p.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-p.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/n.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-n.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/d.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-d.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/H3.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-H3.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/He3.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-He3.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/He4.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-He4.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/Lambda.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-Lambda.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/Lambda(c)+.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-Lambda(c)+.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/LambdaNeutron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/"
      "anti-LambdaNeutron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/LambdaProton.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/"
      "anti-LambdaProton.dN.dy.dat"};

  TH1F *histo[n];
  TH1F *histoSum[n / 2];

  double mass[n / 2] = {0.938272, 0.939565, 1.87561, 2.80892,  2.80839, 3.72738,
                        2.057000, 2.056000, 3.225,   1.800000, 2.286460};
  double integral[8];
  double error[8];

  for (int i{0}; i < n; i = i + 2) {
    graph[i] = new TGraphErrors(fileName[i], "%lg %lg %lg");
    graph[i + 1] = new TGraphErrors(fileName[i + 1], "%lg %lg %lg");
    for (int j{0}; j < graph[i]->GetN(); j++) {
      double x = graph[i]->GetPointX(j);
      double y1 = graph[i]->GetPointY(j);
      double y2 = graph[i + 1]->GetPointY(j);
      double err1 = graph[i]->GetErrorY(j);
      double err2 = graph[i + 1]->GetErrorY(j);
      graph[i]->SetPoint(j, x, y1 + y2);
      graph[i]->SetPointError(j, 0., TMath::Sqrt(err1 * err1 + err2 * err2));
    }
    double partialArea{0};
    double integralError{0};
    for (int k = 0; k < graph[i]->GetN() - 1; ++k) {
      // Estrai i punti x e y
      double x1, y1, x2, y2;
      graph[i]->GetPoint(k, x1, y1);
      if (std::abs(x1) > 0.5)
        continue;
      graph[i]->GetPoint(k + 1, x2, y2);

      // Applica la formula del trapezio per l'area tra i punti (x1, y1) e
      // (x2, y2)
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
    // Per i primi 16 grafici
    error[i / 2] = std::sqrt(integralError); // Errore finale
    integral[i / 2] = partialArea;
  }

  for (int i = 0; i < n / 2; ++i) {
    std::cout << "integral[" << i << "] = " << integral[i] << std::endl;
  }

  TGraphErrors *finalGraph = new TGraphErrors(n / 2, mass, integral, error);

  TF1 *fitExp = new TF1("fitExp", "[0] * exp(-[1] * x)", 0, 4);

  fitExp->SetParameter(0, 1e3);
  fitExp->SetParLimits(0, 0, 1e5);
  fitExp->SetParameter(1, 5);
  fitExp->SetParLimits(1, 4, 8);

  auto canvasFinal = new TCanvas("canvasFinal", "Combined Graph", 800, 600);
  canvasFinal->cd();
  canvasFinal->SetLogy();
  finalGraph->SetTitle("Combined Graph; Mass (GeV); dN/dy");
  finalGraph->Draw("APE");

  // Aggiunta delle etichette ai punti
  const char *labels[n / 2] = {
      "p",          "n",      "d",         "H3", "He3", "He4", "Lambda-Neutron",
      "c-deuteron", "Lambda", "Lambda(c)+"};

  for (int i = 0; i < finalGraph->GetN(); ++i) {
    double x, y;
    finalGraph->GetPoint(i, x, y);
    TLatex *label = new TLatex(x, y * 1.5, labels[i]);
    label->SetTextSize(0.03);
    label->Draw();
  }

  // Configurazione del grafico finale
  finalGraph->GetYaxis()->SetTitleOffset(1.2);
  finalGraph->GetXaxis()->SetTitleSize(0.04);
  finalGraph->GetYaxis()->SetTitleSize(0.04);
  finalGraph->GetXaxis()->SetTitle("Mass (GeV)");
  finalGraph->GetYaxis()->SetTitle("dN/dy");
  finalGraph->GetXaxis()->CenterTitle(true);
  finalGraph->GetYaxis()->CenterTitle(true);
  finalGraph->GetYaxis()->SetLimits(1e-7, 10);
  finalGraph->GetYaxis()->SetRangeUser(1e-7, 10);

  // Aggiornamento della canvas per visualizzare i risultati
  canvasFinal->Update();
}

void cambiamenti() {
  int const n{24};
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

  double count[n / 2];
  double countError[n / 2];

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
    countError[i / 2] = std::sqrt(count[i / 2]) / 2e7;
    std::cout << count[i / 2] << "+/-" << countError[i / 2] << "\n";
  }

  double volume[n / 2] = {8, 6.5, 6, 5.5, 5, 4.5, 4, 3.5, 5, 5, 5};
  double temperature[n / 2] = {156, 156, 156, 156, 156, 156,
                               156, 156, 160, 158, 154, 152};

  TGraph2D *d2Graph = new TGraph2D(n / 2, volume, temperature, count);
  double radius[8] = {8, 6.5, 6, 5.5, 5, 4.5, 4, 3.5};
  double countRadius[8] = {count[0], count[1], count[2], count[3],
                           count[4], count[5], count[6], count[7]};
  double errorVol[8] = {countError[0], countError[1], countError[2],
                        countError[3], countError[4], countError[5],
                        countError[6]};

  double temp[5] = {160, 158, 156, 154, 152};
  double countTemp[5] = {count[8], count[9], count[4], count[10], count[11]};
  double errorTemp[5] = {countError[8], countError[9], countError[4],
                         countError[10], countError[11]};

  TGraphErrors *diffVolume = new TGraphErrors(8, radius, countRadius, errorVol);
  TGraphErrors *diffTemp = new TGraphErrors(5, temp, countTemp, errorTemp);
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

  c1->Update();

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
  // TF1 *fitVolume = new TF1("fitVolume", "[0] * log([1] * x) + [2]", 0, 10);
  TF1 *fitVolume = new TF1("fitVolume", "[0] * x + [1]", 0, 10);
  fitVolume->SetLineColor(kRed);
  fitVolume->SetParameter(0, 1e-6);
  fitVolume->SetParLimits(0, 0, 3e-6);
  fitVolume->SetParameter(1, 1e-3);
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

void densFreq() {

  TH1::AddDirectory(kFALSE);
  TFile *file = new TFile("rivelazioni/densFreq.root", "RECREATE");
  const TString fileName[4] = {
      "generazioni/chemical_charm_30/156_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/c-deuteron.dN.dp.dat",
      "generazioni/chemical_charm_30/156_mev_5_fm/anti-c-deuteron.dN.dp.dat"};

  int count{0}; // Contatore per i valori non nulli nella seconda colonna
  int bin1{0};  // conto il numero di bin di entrambi i dati
  int bin2{0};
  for (int i{0}; i < 2; i++) {
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
      if (i == 0) {
        bin1++;
      }
      if (i == 1) {
        bin2++;
      }
    }
  }
  std::cout << "conteggi " << count << '\n';
  assert(bin1 == bin2);
  std::ifstream inputP(fileName[2]);     // Apertura del file
  std::ifstream inputAntiP(fileName[3]); // Apertura del file
  if (!inputP) {
    std::cerr << "Impossibile aprire il file: " << fileName[2] << std::endl;
  }
  if (!inputAntiP) {
    std::cerr << "Impossibile aprire il file: " << fileName[3] << std::endl;
  }
  std::string line1;
  std::string line2;
  int i{0};
  double P[bin1];   // impulso
  double val[bin1]; // freq
  double err[bin1]; // errfreq
  while (std::getline(inputAntiP, line1) && std::getline(inputP, line2)) {
    std::istringstream s1(line1);
    std::istringstream s2(line2);
    double col1_1, colVal1, colErr1;
    double col2_1, colVal2, colErr2;
    if (s1 >> col1_1 >> colVal1 >> colErr1 &&
        s2 >> col2_1 >> colVal2 >> colErr2) {
      if (col1_1 != col2_1) {
        std::cout << "colonne disaccoppiate" << '\n';
        assert(false);
      } else {
        P[i] = col1_1;
      }
      val[i] = 0.5 * (colVal1 + colVal2) / count;
      err[i] = std::sqrt(colErr1 * colErr1 + colErr2 * colErr2) / count;
      i++;
    }
  }
  inputP.close();
  inputAntiP.close();

  double comulative[bin1];
  comulative[0] = val[0] * (P[1] - P[0]); // Primo intervallo
  for (int j = 1; j < bin1 - 1; ++j) {
    comulative[j] = comulative[j - 1] + val[j] * (P[j + 1] - P[j]);
  }

  // Normalizzazione
  double normalization = comulative[bin1 - 2];
  for (int j = 0; j < bin1 - 1; ++j) {
    comulative[j] /= normalization;
  }

  TGraph *hFreq = new TGraph(bin1, P, val);
  hFreq->SetName("freq");
  hFreq->GetXaxis()->SetTitle("Impulso");
  hFreq->GetYaxis()->SetTitle("Frequenza");
  TCanvas *c1 = new TCanvas("Frequenza", "Frequenza", 200, 10, 600, 400);
  c1->cd();
  hFreq->Write();
  hFreq->Draw("APE");

  TGraph *hFreqComu = new TGraph(bin1, P, comulative);
  hFreqComu->SetName("hFreqComu");
  hFreqComu->GetXaxis()->SetTitle("Impulso");
  hFreqComu->GetYaxis()->SetTitle("Frequenza comulativa");
  TCanvas *c2 = new TCanvas("hFreqComu", "hFreqComu", 200, 10, 600, 400);
  c2->cd();
  hFreqComu->Write();
  hFreqComu->Draw("APE");

  file->Close();
}
