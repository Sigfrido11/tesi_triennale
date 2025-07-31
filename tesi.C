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

void primo_grafico() {
  const int n{12};
  const int div{8};
  TGraphErrors *graph[n];

  const TString fileName[n] = {
      "generazioni/fugacity_30/156_mev_8_fm/d.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/anti-d.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/H3.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/anti-H3.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/He3.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/anti-He3.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/He4.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/anti-He4.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/Lambda(c)+.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/anti-Lambda(c)+.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/anti-c-deuteron.dN.dy.dat"};
  for (int i{0}; i < n; i++) {
    std::ifstream file(fileName[i].Data());
    std::cout << "file: " << fileName[i] << std::endl;
    if (!file.is_open()) {
      std::cerr << "Impossibile aprire il file: " << fileName[i] << std::endl;
      return;
    }
  }

  TH1F *histo[n];
  TH1F *histoSum[n / 2];
  double mass[n / 2] = { 1.87561, 2.80892, 2.80839, 3.72738, 2.28646, 3.225};
  // massa deuterio trizio he3 he4 lambda c c-deuteron
  double massNorm[div / 2] = { 1.87561, 2.80892, 2.80839, 3.72738};

  double massCharm[(n - div) / 2] = {2.28646,
                                     3.225}; // massa lamda c c-deuteron

  double integral[n / 2];
  double error[n / 2];
  double integralogorm[div / 2];
  double errorNorm[div / 2];
  double integralCharm[(n - div) / 2];
  double errorCharm[(n - div) / 2];

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
    error[i / 2] = std::sqrt(integralError); // Errore finale
    integral[i / 2] = partialArea;
    std::cout << "valori " << i / 2 << " " << integral[i / 2] << "\n";
  }

  for (int i = 0; i < n / 2; ++i) {
    if (i < div / 2) {
      integralogorm[i] = integral[i];
      errorNorm[i] = error[i];
    } else {
      integralCharm[i - div / 2] = integral[i];
      errorCharm[i - div / 2] = error[i];
    }
  }
  // verifica correttezza
  for (int i{0}; i < n / 2; i++) {
    std::cout << "integralogorm[i] " << i << " " << integralogorm[i] << '\n';
  }

  for (int i{0}; i < 2; i++) {
    std::cout << "integralCharm[i] " << i << " " << integralCharm[i] <<" " << errorCharm[i] <<'\n';
  }

  // grafici finali
  TMultiGraph *finalGraph = new TMultiGraph;

  TGraphErrors *gTest = new TGraphErrors(n / 2, mass, integral, nullptr, error);

  TGraphErrors *gNormal =
      new TGraphErrors(div / 2, massNorm, integralogorm, nullptr, errorNorm);
  double cDMass[1]={3.225};
  double cDInt[1]={integralCharm[1]};
  double cDError[1]={errorCharm[1]};
  TGraphErrors *gCharm =
      new TGraphErrors(1,cDMass, cDInt, nullptr, cDError);

  // Configurazione della funzione di fit per entrambi i grafici
  TF1 *fitExp = new TF1("fitExp", "[0] * exp(-[1] * x)", 0, 4);
  TF1 *fitExpCharm = new TF1("fitExpCharm", "[0] * exp(-6.2107 * x)", 0, 4);

  fitExp->SetParameter(0, 1e3);
  //fitExp->SetParLimits(0, 0, 1e5);
  fitExp->SetParameter(1, 5);
  //fitExp->SetParLimits(1, 4, 8);

  fitExpCharm->SetParameter(0, 1e3);
  //fitExpCharm->SetParLimits(0, 0, 1e5);
  fitExpCharm->SetParameter(1, 5);
  //fitExpCharm->SetParLimits(1, 4, 8);

  fitExp->SetLineColor(kRed);
  fitExpCharm->SetLineColor(kGreen);
  gNormal->Fit("fitExp");
  gCharm->Fit("fitExpCharm", "R");
  gCharm->SetPoint(0,2.28646,0.316);
  gCharm->SetPoint(1,3.225,integralCharm[1]);

  // Configurazione dei marker per i grafici
  gNormal->SetMarkerStyle(21);
  gNormal->SetMarkerSize(1.0);
  gNormal->SetTitle("Particles; Mass (GeV); dN/dy");

  gCharm->SetMarkerStyle(21);
  gCharm->SetMarkerSize(1.0);
  gCharm->SetTitle("Charm Particles; Mass (GeV); dN/dy");
  gCharm->GetYaxis()->SetRangeUser(1e-7, 10);

  gTest->SetMarkerStyle(21);
  gTest->SetMarkerSize(1.0);
  gTest->SetTitle("Charm Particles; Mass (GeV); dN/dy");

  // Canvas per gNormal
  auto canvasNorm = new TCanvas("canvasNorm", "gNormal - Particles", 800, 600);
  canvasNorm->cd();
  canvasNorm->SetLogy();
  gNormal->GetYaxis()->SetLimits(1e-7, 10);
  gNormal->GetXaxis()->SetLimits(0, 6);
  gNormal->Draw("APE");

  // Canvas per gCharm
  auto canvasCharm = new TCanvas("canvasCharm", "gCharm - Particles", 800, 600);
  canvasCharm->cd();
  canvasCharm->SetLogy();
  gCharm->GetYaxis()->SetLimits(1e-7, 10);
  gCharm->GetXaxis()->SetLimits(0, 6);
  gCharm->Draw("APE");

  // Canvas per test
  auto canvasTest = new TCanvas("canvastest", "gTest - Particles", 800, 600);
  canvasTest->cd();
  canvasTest->SetLogy();
  gTest->GetXaxis()->SetLimits(0, 6);
  gTest->GetYaxis()->SetLimits(1e-7, 10);
  gTest->Draw("APE");

  // Canvas per il grafico combinato
  auto canvasFinal = new TCanvas("canvasFinal", "Combined Graph", 800, 600);
  canvasFinal->cd();
  canvasFinal->SetLogy();
  finalGraph->Add(gNormal, "P");
  finalGraph->Add(gCharm, "P");
  finalGraph->SetTitle("Combined Graph; Mass (GeV); dN/dy");
  finalGraph->Draw("APE");

  // Aggiunta delle etichette ai punti
  const char *labels[n / 2] = {"d",
                               "H3",
                               "He3",
                               "He4",
                               "Lambda(c)",
                               "c-deuteron"};

  canvasTest->cd();
  for (int i = 0; i < gTest->GetN(); ++i) {
    double x, y;
    gTest->GetPoint(i, x, y);
    TLatex *label = new TLatex(x, y * 1.5, labels[i]);
    label->SetTextSize(0.03);
    canvasTest->cd();
    label->Draw();
    canvasFinal->cd();
    label->Draw();
  }
  
  canvasNorm->cd();
  for (int i = 0; i < gNormal->GetN(); ++i) {
    double x, y;
    gNormal->GetPoint(i, x, y);
    TLatex *label = new TLatex(x, y * 1.5, labels[i]);
    label->SetTextSize(0.03);
    label->Draw();
  }
  
  canvasCharm->cd();
  for (int i = 0; i < 1; ++i) {
    double x, y;
    gCharm->GetPoint(i + div / 2, x, y);
    TLatex *label = new TLatex(x, y * 1.5, labels[i + n / 2]);
    label->SetTextSize(0.03);
    label->Draw();
  }
  
  // Aggiunta della legenda
  auto leg = new TLegend(0.7, 0.7, 0.9, 0.9);
  leg->AddEntry(gNormal, "Normal Particles", "P");
  leg->AddEntry(gCharm, "Charm Particles", "P");
  leg->Draw();
std::cout<< "fin qui tutto bene" << '\n';
  // Configurazione del grafico finale
  gTest->SetMarkerStyle(21);
  gTest->SetMarkerSize(1.0);
  gTest->SetTitle("Particles; Mass (GeV); dN/dy");
  gTest->GetYaxis()->SetLimits(1e-7, 10);
  gTest->GetXaxis()->SetLimits(0, 4);

  // Configurazione del grafico finale
  finalGraph->GetYaxis()->SetTitleOffset(1.2);
  finalGraph->GetXaxis()->SetTitleSize(0.04);
  finalGraph->GetYaxis()->SetTitleSize(0.04);
  finalGraph->GetXaxis()->SetTitle("Mass (GeV)");
  finalGraph->GetYaxis()->SetTitle("dN/dy");
  finalGraph->GetXaxis()->CenterTitle(true);
  finalGraph->GetYaxis()->CenterTitle(true);
  finalGraph->GetYaxis()->SetLimits(1e-7, 1);
  finalGraph->GetYaxis()->SetRangeUser(1e-7, 1);
  finalGraph->GetXaxis()->SetLimits(0, 4);
  finalGraph->GetXaxis()->SetRangeUser(0, 4);
  leg->Draw();

  // Aggiornamento della canvas per visualizzare i risultati
  canvasFinal->Update();
}


void cambiamenti() {
  int const n{62};
  // se aggiungi un file le cose da cambiare sono qui
  const TString fileName[n] = {
      "generazioni/fugacity_30/156_mev_4_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_4_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_6_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_6_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_7_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_7_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_9_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_9_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_10_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_10_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_11_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_11_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_12_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_12_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/150_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/150_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/151_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/151_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/152_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/152_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/153_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/153_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/154_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/154_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/155_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/155_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/157_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/157_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/158_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/158_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/159_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/159_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/160_mev_8_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/160_mev_8_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_24/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_24/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_25/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_25/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_26/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_26/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_27/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_27/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_28/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_28/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_30/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_30/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_31/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_31/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_32/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_32/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_33/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_33/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_34/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_34/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_35/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_35/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_36/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/fugacity_36/anti-c-deuteron.dN.dy.dat",

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
  }
  TGraphErrors *graph[n];
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
          0.5 * std::sqrt(y1Error * y1Error + y2Error * y2Error)*(x2-x1);
      integralError +=
          areaError * areaError; // Somma in quadratura degli errori
            }
    countError[i / 2] = std::sqrt(integralError); // Errore finale
    count[i / 2] = partialArea;
    std::cout << "ingresso " << fileName[i] << " " << i << '\n';
    std::cout << partialArea << " err " << countError[i/2] <<  '\n';
  }

  double volume[n / 2] = {4, 5, 6, 7, 8, 9, 10, 11, 12, 8,
                          8, 8, 8, 8, 8, 8, 8,  8,  8};
  double temperature[n / 2] = {156, 156, 156, 156, 156, 156, 156, 156, 156, 150,
                               151, 152, 153, 154, 155, 157, 158, 159, 160};

  TGraph2D *d2Graph = new TGraph2D(n / 2, volume, temperature, count);
  double radius[9] = {4, 5, 6, 7, 8, 9, 10, 11, 12};
  double countRadius[9] = {count[0], count[1], count[2], count[3], count[4],
                           count[5], count[6], count[7], count[8]};
  double errorVol[9] = {countError[0], countError[1], countError[2],
                        countError[3], countError[4], countError[5],
                        countError[6], countError[7], countError[8]};
  double exV [9] ={0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006};

  double temp[11] = {156, 150, 151, 152, 153, 154, 155, 157, 158, 159, 160};
  double countTemp[11] = {count[4],  count[9],  count[10], count[11],
                          count[12], count[13],   count[14], count[15],
                          count[16], count[17], count[18]};
  double errorTemp[11] = {countError[4],  countError[9],  countError[10],
                          countError[11], countError[12], countError[13],
                          countError[14], countError[15], countError[16],
                          countError[17], countError[18]};
  
  double exT[11]={0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006};
  double logTemp[11] = {std::log(count[4]),std::log(count[9]),std::log(count[10]),std::log(count[11]),std::log(count[12]),std::log(0.000534),std::log(count[14]),std::log(count[15]),std::log(count[16]),std::log(count[17]),std::log(count[18])};

  double fugacity[13] = {24, 25, 26, 27, 28, 30, 31, 32, 33, 34, 35, 36, 29};
  double countFugacity[13] = {count[19], count[20], count[21], count[22],
                              count[23], count[24], count[25], count[26],
                              count[27], count[28], count[29], count[30],
                              count[4]};
  double errorFugacity[13] = {countError[19], countError[20], countError[21],
                              countError[22], countError[23], countError[24],
                              countError[25], countError[26], countError[27],
                              countError[28], countError[29], countError[30],
                              countError[4]};
  double exF[13]={0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006,0.006};
  

  TGraphErrors *diffVolume = new TGraphErrors(9, radius, countRadius, exV,errorVol);
  TGraphErrors *diffTemp = new TGraphErrors(11, temp, countTemp,  exT,errorTemp);
  TGraphErrors *gFugacity =
      new TGraphErrors(13, fugacity, countFugacity, exF,errorFugacity);

  // Crea un canvas
  TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);

  c1->SetTheta(15); // Imposta l'angolo di vista in azimut
  c1->SetPhi(20);   // Imposta l'angolo di vista in elevazione
  c1->SetGridx();  // Griglia solo sull'asse X
  c1->SetGridy();  // Griglia solo sull'asse Y

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
  d2Graph->GetXaxis()->SetTitleOffset(1.3);
  d2Graph->GetYaxis()->SetTitleOffset(1.3);
  d2Graph->GetZaxis()->SetTitleOffset(1.);
  d2Graph->GetXaxis()->SetTitleSize(0.05);
  d2Graph->GetYaxis()->SetTitleSize(0.05);
  d2Graph->GetZaxis()->SetTitleSize(0.05);
  d2Graph->GetXaxis()->SetTitle("Radius (fm)");
  d2Graph->GetYaxis()->SetTitle("Temperature (Mev)");
  d2Graph->SetTitle(" ");
  d2Graph->GetZaxis()->SetTitle("dN/dy");
  d2Graph->GetXaxis()->SetLimits(1.5, 14);
  d2Graph->GetXaxis()->SetRangeUser(1.5, 14);
  d2Graph->GetYaxis()->SetLimits(149, 165);
  d2Graph->GetYaxis()->SetRangeUser(149, 165);
  d2Graph->GetZaxis()->SetLimits(1e-5, 1e-3);
  d2Graph->GetZaxis()->SetRangeUser(1e-5, 1e-5);
  d2Graph->Draw("P0");

  c1->Update();

  TCanvas *c2 = new TCanvas("c2", "Canvas 2", 800, 600);
  TPad *pad2 = new TPad("pad2", "Pad piccolo sinistro", 0, 0, 1, 1);
  pad2->Draw();
  pad2->cd();

  // Disegna il grafico diffVolume
  diffVolume->SetLineColor(1);
  diffVolume->GetYaxis()->SetTitleOffset(1.);
  diffVolume->GetXaxis()->SetTitleSize(0.05);
  diffVolume->GetYaxis()->SetTitleSize(0.05);
  diffVolume->GetXaxis()->SetTitle("Radius (fm)");
  diffVolume->GetYaxis()->SetTitle("dN/dy");
  diffVolume->SetTitle(" ");
  diffVolume->GetXaxis()->SetLabelSize(0.04);  // Dimensione etichette sull'asse X
  diffVolume->GetYaxis()->SetLabelSize(0.04);
  diffVolume->Draw("P0");
  diffVolume->SetMarkerStyle(21);
  diffVolume->SetMarkerSize(.8);
  diffVolume->SetMarkerColor(9);
/*
  auto leg1 = new TLegend(0.7, 0.1, 0.9, 0.3);
  leg1->AddEntry(diffVolume, "T = 156 Mev", "");
  leg1->AddEntry(diffVolume, "|y|<0.5", "");
  leg1->SetTextSize(0.04);
  leg1->Draw();
*/
  // (Opzionale) Aggiungi la linea di fit
  // TF1 *fitVolume = new TF1("fitVolume", "[0] * log([1] * x) + [2]", 0, 10);
  TF1 *fitVolume = new TF1("fitVolume", "[0] * x^3 + [1]", 0, 10);
  fitVolume->SetLineColor(kRed);
  fitVolume->SetParameter(0, 1e-6);
  // fitVolume->SetParLimits(0, 0, 3e-6);
  fitVolume->SetParameter(1, 1e-3);
  // fitVolume->SetParLimits(1, 0, 1e5);
  diffVolume->Fit(fitVolume);
  /*
  leg1->AddEntry(fitVolume, "Fit Line", "l");
  leg1->Draw();
  */
  diffVolume->Draw("APE");

  c2->Update(); // Aggiorna il canvas

  // Crea il terzo canvas per il pad3
  TCanvas *c3 = new TCanvas("c3", "Canvas 3", 800, 600);
  TPad *pad3 = new TPad("pad3", "Pad piccolo destro", 0, 0, 1, 1);;
  pad3->Draw();
  c3->cd();
  
TF1 *fitTemp = new TF1("fitTemp", "[0]*exp(x*[1])+[2]", 140, 165);

// Imposta i parametri iniziali per il fit
fitTemp->SetParameter(0, 4.24449e-14);  // Parametro per l'intercetta (start value)
fitTemp->SetParameter(1, 0.150492); 
fitTemp->SetParameter(2, 0.0000284396);  // Parametro per la pendenza (start value)
fitTemp->SetLineColor(kGreen);   // Imposta il colore della linea di fit
//fitTemp->SetParameter(2,1);
// fitTemp->SetParameter(0,0);
// fitTemp->SetParameter(2,2e-6);

  // Disegna il grafico diffTemp
  diffTemp->SetLineColor(1);
  diffTemp->GetYaxis()->SetTitleOffset(1.);
  diffTemp->GetXaxis()->SetTitleSize(0.05);
  diffTemp->GetYaxis()->SetTitleSize(0.05);
  diffTemp->GetXaxis()->SetTitle("Temperature (Mev)");
  diffTemp->GetYaxis()->SetTitle("dN/dy");
  diffTemp->SetTitle(" ");
  diffTemp->GetXaxis()->SetLabelSize(0.04);  // Dimensione etichette sull'asse X
  diffTemp->GetYaxis()->SetLabelSize(0.04);
  diffTemp->Fit(fitTemp);
  diffTemp->Draw("APE");
  diffTemp->SetMarkerStyle(21);
  diffTemp->SetMarkerSize(.8);
  diffTemp->SetMarkerColor(kRed);
  auto leg2 = new TLegend(0.7, 0.1, 0.9, 0.3);
  leg2->AddEntry(diffTemp, "R=8fm", "");
  leg2->AddEntry(diffTemp, "|y|<0.5", "");
  leg2->SetTextSize(0.04);
  //leg2->Draw();

  // (Opzionale) Aggiungi la linea di fit


  c3->Update(); // Aggiorna il canvas

  TF1 *fitFug = new TF1("fugacity", "[0]*x + [1]", 0, 10);
  fitFug->SetLineColor(kBlue);
  fitFug->SetParameter(0, 1e-3);
  // fitFug->SetParLimits(0, 0, 1e-5);
  fitFug->SetParameter(1, 1e-3);
  //fitFug->SetParLimits(1, 0, 1e-4);
  // fitFug->SetParameter(2, 1e-3);
  // fitFug->SetParLimits(2, 0, 1e-3);

  gFugacity->Fit(fitFug);
  // leg2->AddEntry(fitTemp, "Fit Line", "l");

  TCanvas *c4 = new TCanvas("c4", "Canvas 4", 800, 600);
  // Disegna il grafico diffTemp
  gFugacity->SetLineColor(1);
  gFugacity->GetYaxis()->SetTitleOffset(1.);
  gFugacity->GetXaxis()->SetTitleSize(0.05);
  gFugacity->GetYaxis()->SetTitleSize(0.05);
  gFugacity->GetXaxis()->SetTitle("Charm fugacity");
  gFugacity->GetYaxis()->SetTitle("dN/dy");
  gFugacity->SetTitle(" ");
  gFugacity->GetXaxis()->SetLabelSize(0.04);  // Dimensione etichette sull'asse X
  gFugacity->GetYaxis()->SetLabelSize(0.04);  // Dimensione etichette sull'asse 
  gFugacity->SetMarkerStyle(21);
  gFugacity->SetMarkerSize(.8);
  gFugacity->SetMarkerColor(kRed);
  gFugacity->Draw("APE");
  auto leg3 = new TLegend(0.7, 0.1, 0.9, 0.3);
  leg3->AddEntry(diffTemp, "R=8fm", "");
  leg3->AddEntry(diffTemp, "t=156 MeV", "");
  leg3->AddEntry(diffTemp, "|y|<0.5", "");
  leg3->SetTextSize(0.04);
  //leg3->Draw();
}
/*
void densFreq() {

  TH1::AddDirectory(kFALSE);
  TFile *file =
      new TFile("rivelazioni/pythia/pythia8312/densFreq.root", "RECREATE");
  const TString fileName[4] = {
      "generazioni/fugacity_30/156_mev_5_fm/c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_5_fm/anti-c-deuteron.dN.dy.dat",
      "generazioni/fugacity_30/156_mev_5_fm/c-deuteron.dN.dp.dat",
      "generazioni/fugacity_30/156_mev_5_fm/anti-c-deuteron.dN.dp.dat"};

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

  TGraph *Freq = new TGraph(bin1, P, val);
  Freq->SetName("freq");
  Freq->GetXaxis()->SetTitle("Impulso");
  Freq->GetYaxis()->SetTitle("Frequenza");
  TCanvas *c1 = new TCanvas("Frequenza", "Frequenza", 200, 10, 600, 400);
  c1->cd();
  Freq->Write();
  Freq->Draw("APE");

  TGraph *FreqComu = new TGraph(bin1, P, comulative);
  FreqComu->SetName("FreqComu");
  FreqComu->GetXaxis()->SetTitle("Impulso");
  FreqComu->GetYaxis()->SetTitle("Frequenza comulativa");
  TCanvas *c2 = new TCanvas("hFreqComu", "hFreqComu", 200, 10, 600, 400);
  c2->cd();
  FreqComu->Write();
  FreqComu->Draw("APE");

  file->Close();
}
*/  
