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
  const int n{14};
  TGraphErrors *graph[n];
  /*
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
    //
    "generazioni/156_mev_5_fm/HyperTriton.dN.dy.dat","generazioni/156_mev_5_fm/anti-HyperTriton.dN.dy.dat"};

  */

  const TString fileName[n] = {
      "generazioni/155_mev_50_yc_8_fm/p.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/anti-p.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/n.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/anti-n.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/d.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/anti-d.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/H3.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/anti-H3.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/He3.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/anti-He3.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/He4.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/anti-He4.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/155_mev_50_yc_8_fm/anti-c-deuteron.dN.dy.dat"};
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

  double a[1] = {mass[n / 2]};
  double b[1] = {integral[n / 2]};
  TGraph *cGraph = new TGraph(1, a, b);
  TF1 *cFit = new TF1("c-Fit", "[0] * exp(-[1] * x) + [2]", 0, 4);
  cFit->FixParameter(0, fitExp->GetParameter(0));
  cFit->FixParameter(1, fitExp->GetParameter(1));
  cGraph->Fit(cFit);

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
  finalGraph->Draw("APE");
  cFit->Draw();
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
  TFile *file = new TFile("densFreq.root");
  const TString fileName[4] = {
      "generazioni/156_mev_4_fm/c-deuteron.dN.dy.dat",
      "generazioni/156_mev_4_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/156_mev_4_fm/c-deuteron.dN.dp.dat",
      "generazioni/156_mev_4_fm/anti-c-deuteron.dN.dp.dat"};

  int count{0}; // Contatore per i valori non nulli nella seconda colonna
  int bin1{0};
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
  double val[bin1]; // densfreq
  double err[bin1]; // errdensfreq
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
  TGraph *hDensFreq = new TGraph(bin1, P, val);
  hDensFreq->SetName("densFreq");
  hDensFreq->GetXaxis()->SetTitle("Impulso");
  hDensFreq->GetYaxis()->SetTitle("Densità di Frequenza");
  double totVal{0.0};
  for (int i = 0; i < bin1 - 1; ++i) {
    totVal += (P[i + 1] - P[i]) * val[i];
}
  double totErr = std::sqrt(std::accumulate(err, err + bin1 - 1, 0.0, [](double acc, double e){
    return acc + e * e;
  }));
  
  if( 1 > totVal-totErr && totVal + totErr >1 ){
    std::cout << "compà t'appost" << '\n';
    std::cout << totVal << " +/- " << totErr << '\n';
  }
  //controllo la normalizzazione
  hDensFreq->Write();
  file->Close();
}
