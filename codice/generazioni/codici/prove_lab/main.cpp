#include <TCanvas.h>
#include <algorithm>
#include <cmath>
#include <iterator>
#include <TFile.h>
#include <TH1.h>
#include <TRandom.h>
#include "particle.h"
#include "particletype.h"
#include "resonancetype.h"

int main()
{
  TH1::AddDirectory(kFALSE);
  TFile *file = new TFile("histo.root", "RECREATE");
  gRandom->SetSeed(200769);
  Particle::AddParticleType("pion +", 0.13957, 1);
  Particle::AddParticleType("pion -", 0.13957, -1);
  Particle::AddParticleType("kaon +", 0.49367, 1);
  Particle::AddParticleType("kaon -", 0.49367, -1);
  Particle::AddParticleType("proton +", 0.93827, 1);
  Particle::AddParticleType("proton -", 0.93827, -1);
  Particle::AddParticleType("resonance", 0.89166, 0, 0.050);
  std::array<Particle, 120> EventParticles{};

  TH1F *gen_particles = new TH1F("gen_particles", "gen_particles", 7, 0, 7);

  TH1F *azimuth = new TH1F("azimuth", "azimuth", 1000, 0., 2 * M_PI);

  TH1F *polar_angle = new TH1F("polar_angle", "polar_angle", 1000, 0., M_PI);

  TH1F *p_module = new TH1F("p_module", "p_module", 1000, 0, 5);

  TH1F *p_trans = new TH1F("p_trans", "p_trans", 1000, 0, 5);

  TH1F *energy = new TH1F("energy", "energy", 1000, 0, 6);

  TH1F *all_inv_mass = new TH1F("all_inv_mass", "all_inv_mass", 1000, 0, 6);

  all_inv_mass->Sumw2();

  TH1F *same_charge_inv_mass =
      new TH1F("same_charge_inv_mass", "same_charge_inv_mass", 1000, 0, 6);

  same_charge_inv_mass->Sumw2();

  TH1F *opposite_charge_inv_mass = new TH1F(
      "opposite_charge_inv_mass", "opposite_charge_inv_mass", 1000, 0, 6);
  opposite_charge_inv_mass->Sumw2();
  TH1F *pi_k_same = new TH1F("pi_k_same_charge_inv_mass",
                             "pi_k_same_charge_inv_mass", 1000, 0, 6);
  pi_k_same->Sumw2();
  TH1F *pi_k_opposite = new TH1F("pi_k_opposite_charge_inv_mass",
                                 "pi_k_opposite_charge_inv_mass", 1000, 0, 6);
  pi_k_opposite->Sumw2();
  TH1F *dec_inv_mass = new TH1F("dec_inv_mass", "dec_inv_mass", 1000, 0, 4.5);
  dec_inv_mass->Sumw2();

  int star_num{0};

  for (int i{0}; i < 1e5; i++)
  {
    if (star_num > 20)
    {
      std::cout << "too much k* has been generated" << '\n';
    }
    star_num = 0; // k* counter

    for (int event{0}; event < 1e2; event++)
    {
      bool is_star{false};                         // tells if a k* is found
      double type = gRandom->Rndm();               // type selector
      double phi = gRandom->Uniform(0., 2 * M_PI); // azimuth angle
      double theta = gRandom->Uniform(0., M_PI);   // polar angle
      double p_m = gRandom->Exp(1.);               // P module
      double px = p_m * std::cos(phi) * std::sin(theta);
      double py = p_m * std::sin(phi) * std::sin(theta);
      double pz = p_m * std::cos(theta);

      azimuth->Fill(phi); // filling histo
      polar_angle->Fill(theta);
      p_module->Fill(p_m);
      p_trans->Fill(std::sqrt(std::pow(px, 2) + std::pow(py, 2)));

      Particle p;
      Particle p1;
      Particle p2;

      if (type <= 0.4)
      {
        p.SetName("pion +");
        gen_particles->Fill(0.5);
      }

      else if (type <= 0.8)
      {
        p.SetName("pion -");
        gen_particles->Fill(1.5);
      }

      else if (type <= 0.85)
      {
        p.SetName("kaon +");
        gen_particles->Fill(2.5);
      }

      else if (type <= 0.90)
      {
        p.SetName("kaon -");
        gen_particles->Fill(3.5);
      }

      else if (type <= 0.945)
      {
        p.SetName("proton +");
        gen_particles->Fill(4.5);
      }

      else if (type <= 0.99)
      {
        p.SetName("proton -");
        gen_particles->Fill(5.5);
      }

      else
      {
        p.SetName("resonance");
        gen_particles->Fill(6.5);
        is_star = true;
        double decad{gRandom->Rndm()};
        p.SetP(px, py, pz);
        energy->Fill(p.GetEnergy());
        star_num += 1;
        if (decad <= 0.5)
        {
          p1.SetName("pion +");
          p2.SetName("kaon -");
        }
        else
        {
          p1.SetName("pion -");
          p2.SetName("kaon +");
        }
        int all_good{p.Decay2Body(p1, p2)};
        if (all_good != 0)
        {
          std::cout << "something went wrong during decay" << '\n'
                    << "\n";
        }
        else
        {
          dec_inv_mass->Fill(p1.GetInvariantMass(p2));
        }
      }
      if (!is_star)
      {
        p.SetP(px, py, pz);
        energy->Fill(p.GetEnergy());
        EventParticles[event + star_num] = p;
      }
      else
      {
        EventParticles[event + star_num] = p1;
        EventParticles[event - 1 + star_num] = p2;
      }

      for (int compare{0}; compare < event + star_num; compare++)
      {
        Particle old_particle{EventParticles[compare]};
        if (is_star)
        {
          p = p1;
        }
        all_inv_mass->Fill(p.GetInvariantMass(old_particle));
        if (p.GetCharge() == old_particle.GetCharge())
        { // same charge
          same_charge_inv_mass->Fill(p.GetInvariantMass(old_particle));
          const bool first_cond{p.GetName() == "pion +" &&
                                old_particle.GetName() == "kaon +"};
          const bool second_cond{p.GetName() == "pion -" &&
                                 old_particle.GetName() == "kaon -"};
          if (first_cond || second_cond)
          {
            pi_k_same->Fill(p.GetInvariantMass(old_particle));
          }
        }
        else
        { // different charge
          opposite_charge_inv_mass->Fill(p.GetInvariantMass(old_particle));
          const bool first_cond{p.GetName() == "pion +" &&
                                old_particle.GetName() == "kaon -"};
          const bool second_cond{p.GetName() == "pion -" &&
                                 old_particle.GetName() == "kaon +"};
          if (first_cond || second_cond)
          {
            pi_k_opposite->Fill(p.GetInvariantMass(old_particle));
          }
        }
      }
    }
  }

  gen_particles->GetXaxis()->SetBinLabel(1, "pi+");
  gen_particles->GetXaxis()->SetBinLabel(2, "pi-");
  gen_particles->GetXaxis()->SetBinLabel(3, "k+");
  gen_particles->GetXaxis()->SetBinLabel(4, "k-");
  gen_particles->GetXaxis()->SetBinLabel(5, "p+");
  gen_particles->GetXaxis()->SetBinLabel(6, "p-");
  gen_particles->GetXaxis()->SetBinLabel(7, "k*");

  azimuth->GetXaxis()->SetTitle("angle (rad)");
  azimuth->GetXaxis()->SetTitleSize(0.045);
  azimuth->GetYaxis()->SetTitle("events");
  azimuth->GetYaxis()->SetTitleSize(0.045);
  azimuth->GetYaxis()->SetRangeUser(0.,12000.);

  polar_angle->GetXaxis()->SetTitle("angle (rad)");
  polar_angle->GetXaxis()->SetTitleSize(0.045);
  polar_angle->GetYaxis()->SetTitle("events");
  polar_angle->GetYaxis()->SetTitleSize(0.045);
  polar_angle->GetYaxis()->SetRangeUser(0.,12000.);

  gen_particles->GetXaxis()->SetTitle("particles");
  gen_particles->GetXaxis()->SetTitleSize(0.045);
  gen_particles->GetYaxis()->SetTitle("events");
  gen_particles->GetYaxis()->SetTitleSize(0.045);

  dec_inv_mass->GetXaxis()->SetTitle("invariant mass (GeV)");
  dec_inv_mass->GetXaxis()->SetTitleSize(0.045);
  dec_inv_mass->GetYaxis()->SetTitle("events");
  dec_inv_mass->GetYaxis()->SetTitleSize(0.045);

  pi_k_opposite->GetXaxis()->SetTitle("invariant mass (GeV)");
  pi_k_opposite->GetXaxis()->SetTitleSize(0.045);
  pi_k_opposite->GetYaxis()->SetTitle("events");
  pi_k_opposite->GetYaxis()->SetTitleSize(0.045);

  pi_k_same->GetXaxis()->SetTitle("invariant mass (GeV)");
  pi_k_same->GetXaxis()->SetTitleSize(0.045);
  pi_k_same->GetYaxis()->SetTitle("events");
  pi_k_same->GetYaxis()->SetTitleSize(0.045);

  opposite_charge_inv_mass->GetXaxis()->SetTitle("invariant mass (GeV)");
  opposite_charge_inv_mass->GetXaxis()->SetTitleSize(0.045);
  opposite_charge_inv_mass->GetYaxis()->SetTitle("events");
  opposite_charge_inv_mass->GetYaxis()->SetTitleSize(0.045);

  same_charge_inv_mass->GetXaxis()->SetTitle("invariant mass (GeV)");
  same_charge_inv_mass->GetXaxis()->SetTitleSize(0.045);
  same_charge_inv_mass->GetYaxis()->SetTitle("events");
  same_charge_inv_mass->GetYaxis()->SetTitleSize(0.045);

  all_inv_mass->GetXaxis()->SetTitle("invariant mass (GeV)");
  all_inv_mass->GetXaxis()->SetTitleSize(0.045);
  all_inv_mass->GetYaxis()->SetTitle("events");
  all_inv_mass->GetYaxis()->SetTitleSize(0.045);

  energy->GetXaxis()->SetTitle("energy (GeV)");
  energy->GetXaxis()->SetTitleSize(0.045);
  energy->GetYaxis()->SetTitle("events");
  energy->GetYaxis()->SetTitleSize(0.045);

  p_trans->GetXaxis()->SetTitle("linear momentum (GeV)");
  p_trans->GetXaxis()->SetTitleSize(0.045);
  p_trans->GetYaxis()->SetTitle("events");
  p_trans->GetYaxis()->SetTitleSize(0.045);

  p_module->GetXaxis()->SetTitle("linear momentum (GeV)");
  p_module->GetXaxis()->SetTitleSize(0.045); 
  p_module->GetYaxis()->SetTitle("events");
  p_module->GetYaxis()->SetTitleSize(0.045);

  azimuth->Write();
  polar_angle->Write();
  p_module->Write();
  p_trans->Write();
  gen_particles->Write();
  energy->Write();
  all_inv_mass->Write();
  same_charge_inv_mass->Write();
  opposite_charge_inv_mass->Write();
  pi_k_same->Write();
  pi_k_opposite->Write();
  dec_inv_mass->Write();
  file->Close();
}
