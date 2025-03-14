#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include <TCanvas.h>
#include <iostream>

void setStyle()
{
  gROOT->SetStyle("Plain");
  gStyle->SetPalette(57);
  gStyle->SetOptTitle(0);
}

void analyse()
{
  TH1::AddDirectory(kFALSE);
  TFile *file = new TFile("histo.root");
  TH1F *gen_particles = (TH1F *)file->Get("gen_particles");
  TH1F *azimuth = (TH1F *)file->Get("azimuth");
  TH1F *polar_angle = (TH1F *)file->Get("polar_angle");
  TH1F *p_module = (TH1F *)file->Get("p_module");
  TH1F *p_trans = (TH1F *)file->Get("p_trans");
  TH1F *energy = (TH1F *)file->Get("energy");
  TH1F *all_inv_mass = (TH1F *)file->Get("all_inv_mass");
  TH1F *same_charge_inv_mass = (TH1F *)file->Get("same_charge_inv_mass");
  TH1F *opposite_charge_inv_mass =
      (TH1F *)file->Get("opposite_charge_inv_mass");
  TH1F *pi_k_same = (TH1F *)file->Get("pi_k_same_charge_inv_mass");
  TH1F *pi_k_opposite = (TH1F *)file->Get("pi_k_opposite_charge_inv_mass");
  TH1F *dec_inv_mass = (TH1F *)file->Get("dec_inv_mass");
  file->Close();

  if (gen_particles->GetEntries() == 1e7)
  {
    std::cout << "entries in gen particles histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in gen particles histo: unexpected value" << '\n';
  }
  if (azimuth->GetEntries() == 1e7)
  {
    std::cout << "entries in azimuth histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in azimuth histo: unexpected value" << '\n';
  }
  if (polar_angle->GetEntries() == 1e7)
  {
    std::cout << "entries in polar angle histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in polar angle histo: unexpected value" << '\n';
  }
  if (p_module->GetEntries() == 1e7)
  {
    std::cout << "entries in p module histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in p module histo: unexpected value" << '\n';
  }
  if (p_trans->GetEntries() == 1e7)
  {
    std::cout << "entries in p trans histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in p trans histo: unexpected value" << '\n';
  }
  if (energy->GetEntries() == 1e7)
  {
    std::cout << "entries in energy histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in energy histo: unexpected value" << '\n';
  }
  if (all_inv_mass->GetEntries() < 1e5 * 60 * 119 && //number of combinations of n (maximum)
      all_inv_mass->GetEntries() > 1e5 * 40 * 79)  //number of combinations of n (minimum)
  {
    std::cout << "entries in gen particles histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in gen particles histo: unexpected value" << '\n';
  }
  if (same_charge_inv_mass->GetEntries() < 1e5 * 70 * 71 && //two times sum of n numbers (maximum)
      same_charge_inv_mass->GetEntries() > 1e5 * 45 * 44) //two times sum of n numbers (minimum)
  {
    std::cout << "entries in same charge inv mass histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in same charge inv mass histo: unexpected value" << same_charge_inv_mass->GetEntries()
              << '\n';
  }
  if (opposite_charge_inv_mass->GetEntries() < 1e5 * 70 * 50 && //n positives times n negatives (maximum)
      opposite_charge_inv_mass->GetEntries() > 1e5 * 45 * 45) //n positives times n negatives (minimum)
  {
    std::cout << "entries in opposite charge inv mass histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in opposite charge inv mass histo: unexpected value" << opposite_charge_inv_mass->GetEntries()
              << '\n';
  }
  if (pi_k_same->GetEntries() < 1e5 * 60 * 25 * 2 && //two times n k same times n pi same (maximum)
      pi_k_same->GetEntries() > 1e5 * 30 * 3 * 2) //two times n k same times n pi same (minimum)
  {
    std::cout << "entries in pi k same histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in pi k same histo: unexpected value" << '\n';
  }
  if (pi_k_opposite->GetEntries() < 1e5 * 60 * 25 * 2 && //two times n k opposite times n pi opposite (maximum)
      pi_k_opposite->GetEntries() > 1e5 * 30 * 3 * 2) //two times n k opposite times n pi opposite (minimum)
  {
    std::cout << "entries in pi k opposite histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in pi k opposite histo: unexpected value" << '\n';
  }
  if (dec_inv_mass->GetEntries() < 1e6 && dec_inv_mass->GetEntries() > 1e4)
  {
    std::cout << "entries in dec inv mass histo: ok" << '\n';
  }
  else
  {
    std::cout << "entries in dec inv mass histo: unexpected value" << '\n';
  }

  if (gen_particles->GetBinContent(1) - 3 * gen_particles->GetBinError(1) <=
          gen_particles->GetEntries() * 0.4 &&
      gen_particles->GetBinContent(1) + 3 * gen_particles->GetBinError(1) >=
          gen_particles->GetEntries() * 0.4)
  {
    std::cout << "n pi+: ok" << '\n';
  }
  else
  {
    std::cout << "n pi+: unexpected" << '\n';
  }
  if (gen_particles->GetBinContent(2) - 3 * gen_particles->GetBinError(2) <=
          gen_particles->GetEntries() * 0.4 &&
      gen_particles->GetBinContent(2) + 3 * gen_particles->GetBinError(2) >=
          gen_particles->GetEntries() * 0.4)
  {
    std::cout << "n pi-: ok" << '\n';
  }
  else
  {
    std::cout << "n pi-: unexpected" << '\n';
  }
  if (gen_particles->GetBinContent(3) - 3 * gen_particles->GetBinError(3) <=
          gen_particles->GetEntries() * 0.05 &&
      gen_particles->GetBinContent(3) + 3 * gen_particles->GetBinError(3) >=
          gen_particles->GetEntries() * 0.05)
  {
    std::cout << "n k+: ok" << '\n';
  }
  else
  {
    std::cout << "n k+: unexpected" << '\n';
  }
  if (gen_particles->GetBinContent(4) - 3 * gen_particles->GetBinError(4) <=
          gen_particles->GetEntries() * 0.05 &&
      gen_particles->GetBinContent(4) + 3 * gen_particles->GetBinError(4) >=
          gen_particles->GetEntries() * 0.05)
  {
    std::cout << "n k-: ok" << '\n';
  }
  else
  {
    std::cout << "n k-: unexpected" << '\n';
  }
  if (gen_particles->GetBinContent(5) - 3 * gen_particles->GetBinError(5) <=
          gen_particles->GetEntries() * 0.045 &&
      gen_particles->GetBinContent(5) + 3 * gen_particles->GetBinError(5) >=
          gen_particles->GetEntries() * 0.045)
  {
    std::cout << "n p+: ok" << '\n';
  }
  else
  {
    std::cout << "n p+: unexpected" << '\n';
  }
  if (gen_particles->GetBinContent(6) - 3 * gen_particles->GetBinError(6) <=
          gen_particles->GetEntries() * 0.045 &&
      gen_particles->GetBinContent(6) + 3 * gen_particles->GetBinError(6) >=
          gen_particles->GetEntries() * 0.045)
  {
    std::cout << "n p-: ok" << '\n';
  }
  else
  {
    std::cout << "n p-: unexpected" << '\n';
  }
  if (gen_particles->GetBinContent(7) - 3 * gen_particles->GetBinError(7) <=
          gen_particles->GetEntries() * 0.01 &&
      gen_particles->GetBinContent(7) + 3 * gen_particles->GetBinError(7) >=
          gen_particles->GetEntries() * 0.01)
  {
    std::cout << "n k*: ok" << '\n';
  }
  else
  {
    std::cout << "n k*: unexpected" << '\n';
  }

  gStyle->SetOptFit(1111);

  TF1 *f1 = new TF1("f1", "[0]", 0, 10);
  azimuth->Fit(f1);
  std::cout << "azimuth: uniform fit parameter: " << f1->GetParameter(0)
            << "+/-" << f1->GetParError(0) << '\n';
  std::cout << "ChiSquare/NDF; " << f1->GetChisquare() / f1->GetNDF() << '\n';
  std::cout << "probability: " << f1->GetProb() << '\n';
  polar_angle->Fit(f1);
  std::cout << "polar angle: uniform fit parameter: " << f1->GetParameter(0)
            << "+/-" << f1->GetParError(0) << '\n';
  std::cout << "ChiSquare/NDF " << f1->GetChisquare() / f1->GetNDF() << '\n';
  std::cout << "probability: " << f1->GetProb() << '\n';

  f1 = new TF1("f2", "expo(0)", 0, 1);
  p_module->Fit(f1);
  if ((-f1->GetParameter(1) - 1) / f1->GetParError(1) < 2 &&
      (-f1->GetParameter(1) - 1) / f1->GetParError(1) > -2)
  {
    std::cout << "p module distrution is as expexted" << '\n';
  }
  else
  {
    std::cout << "unexpexted results for p module distriution" << '\n';
  }
  std::cout << "p module: expo fit parameter: " << f1->GetParameter(1) << "+/-"
            << -f1->GetParError(1) << '\n';
  std::cout << "ChiSquare/NDF " << f1->GetChisquare() / f1->GetNDF() << '\n';
  std::cout << "probability: " << f1->GetProb() << '\n';

  TH1F *diff = new TH1F(
      "difference_inv_mass", "difference_inv_mass", 1000, 0, 6.);
  diff->Sumw2();
  diff->Add(opposite_charge_inv_mass, same_charge_inv_mass, 1., -1.);
  TH1F *diff_pi_k = new TH1F(
      "difference_inv_mass_pi_k", "difference_inv_mass_pi_k", 1000, 0, 6.);
  diff_pi_k->Sumw2();
  diff_pi_k->Add(pi_k_opposite, pi_k_same, 1., -1.);

  std::cout << "maximum: "
            << '\n'
            << "from diff opposite same all: "
            << diff->GetBinCenter(diff->GetMaximumBin())
            << '\n'
            << "from diff opposite same pi-k: "
            << diff_pi_k->GetBinCenter(diff_pi_k->GetMaximumBin())
            << '\n'
            << "from decay product: "
            << dec_inv_mass->GetBinCenter(dec_inv_mass->GetMaximumBin())
            << '\n';

  f1 = new TF1("f3", "gaus(0)", 0, 10);
  dec_inv_mass->Fit(f1);
  std::cout << "from decay products: "
            << '\n'
            << "k* mass: " << f1->GetParameter(1) << "+/-" << f1->GetParError(1)
            << '\n'
            << "k* width: " << f1->GetParameter(2) << "+/-" << f1->GetParError(2)
            << '\n'
            << "gauss width: " << f1->GetParameter(0) << "+/-" << f1->GetParError(0)
            << '\n';
  std::cout << "ChiSquare/NDF " << f1->GetChisquare() / f1->GetNDF() << '\n';
  std::cout << "probability: " << f1->GetProb() << '\n';
  diff->Fit(f1);
  std::cout << "from difference opposite and same charge inv mass: "
            << '\n'
            << "k* mass: " << f1->GetParameter(1) << "+/-" << f1->GetParError(1)
            << '\n'
            << "k* width: " << f1->GetParameter(2) << "+/-" << f1->GetParError(2)
            << '\n'
            << "gauss width: " << f1->GetParameter(0) << "+/-" << f1->GetParError(0)
            << '\n';
  std::cout << "ChiSquare/NDF " << f1->GetChisquare() / f1->GetNDF() << '\n';
  std::cout << "probability: " << f1->GetProb() << '\n';
  diff_pi_k->Fit(f1);
  std::cout << "from difference opposite and same charge pions and kaons inv mass: "
            << '\n'
            << "k* mass: " << f1->GetParameter(1) << "+/-" << f1->GetParError(1)
            << '\n'
            << "k* width: " << f1->GetParameter(2) << "+/-" << f1->GetParError(2)
            << '\n'
            << "gauss width: " << f1->GetParameter(0) << "+/-" << f1->GetParError(0)
            << '\n';
  std::cout << "ChiSquare/NDF " << f1->GetChisquare() / f1->GetNDF() << '\n';
  std::cout << "probability: " << f1->GetProb() << '\n';

  diff->GetXaxis()->SetTitle("invariant mass (GeV)");
  diff->GetXaxis()->SetTitleSize(0.045);
  diff->GetYaxis()->SetTitle("events");
  diff->GetYaxis()->SetTitleSize(0.045);
  diff_pi_k->GetXaxis()->SetTitle("invariant mass (GeV)");
  diff_pi_k->GetXaxis()->SetTitleSize(0.045);
  diff_pi_k->GetYaxis()->SetTitle("events");
  diff_pi_k->GetYaxis()->SetTitleSize(0.045);

  TCanvas *c1 = new TCanvas("c1", "n gen, p module and angles", 200, 10, 600, 400);
  c1->Divide(2, 2);
  c1->cd(1);
  gen_particles->Draw("HEP");
  c1->cd(2);
  p_module->Draw("HEP");
  c1->cd(3);
  azimuth->Draw("HEP");
  c1->cd(4);
  polar_angle->Draw("HEP");
  c1->Print("N_P_Angles.pdf");
  c1->Print("N_P_Angles.C");
  c1->Print("N_P_Angles.root");

  TCanvas *c2 = new TCanvas("c2", "p trans", 200, 10, 600, 400);
  p_trans->Draw("HEP");
  c2->Print("PTrans.pdf");
  c2->Print("PTrans.C");
  c2->Print("PTrans.root");

  TCanvas *c3 = new TCanvas("c3", "energy and all inv mass", 200, 10, 600, 400);
  c3->Divide(1, 2);
  c3->cd(1);
  energy->Draw("HEP");
  c3->cd(2);
  all_inv_mass->Draw("HEP");
  c3->Print("Energy_AllInvMass.pdf");
  c3->Print("Energy_AllInvMass.C");
  c3->Print("Energy_AllInvMass.root");

  TCanvas *c4 = new TCanvas("c4", "inv mass same and opposite (all and pi-k)", 200, 10, 600, 400);
  c4->Divide(2, 2);
  c4->cd(1);
  same_charge_inv_mass->Draw("HEP");
  c4->cd(2);
  opposite_charge_inv_mass->Draw("HEP");
  c4->cd(3);
  pi_k_same->Draw("HEP");
  c4->cd(4);
  pi_k_opposite->Draw("HEP");
  c4->Print("InvMassSameOpposite.pdf");
  c4->Print("InvMassSameOpposite.C");
  c4->Print("InvMassSameOpposite.root");

  TCanvas *c5 = new TCanvas("c5", "invariant mass from decay and differences", 200, 10, 600, 400);
  c5->Divide(1, 3);
  c5->cd(1);
  dec_inv_mass->Draw("HEP");
  c5->cd(2);
  diff->Draw("HEP");
  c5->cd(3);
  diff_pi_k->Draw("HEP");
  c5->Print("Decay_Diffs.pdf");
  c5->Print("Decay_Diffs.C");
  c5->Print("Decay_Diffs.root");
}