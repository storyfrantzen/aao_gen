// ROOT macro to plot the aao_norad MC cross section vs phi at fixed xB, Q2, t,
// and beam energy.
//
// Run on the server that has ROOT and can execute ./aao_xsec:
//   cd /Users/wangyijie/remote_clas12/Simulation/aao_print/aao_gen/aao_norad
//   root -l -q 'plot_xsec_phi.C(0.3,2.0,-0.2,10.6)'
//
// Optional scan control:
//   root -l -q 'plot_xsec_phi.C(0.3,2.0,-0.2,10.6,0,360,5,"my_phi_plot")'

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TString.h>
#include <TLegend.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

bool RunAAOXsec(double xB, double Q2, double t, double phi, double beamEnergy,
                double &sigma)
{
  TString cmd;
  cmd.Form("./aao_xsec -xB %.12g -Q2 %.12g -t %.12g -phi %.12g -BeamEnergy %.12g",
           xB, Q2, t, phi, beamEnergy);

  FILE *pipe = gSystem->OpenPipe(cmd.Data(), "r");
  if (!pipe) {
    std::cerr << "ERROR: cannot run command: " << cmd << std::endl;
    return false;
  }

  char buffer[512];
  std::string output;
  while (fgets(buffer, sizeof(buffer), pipe)) {
    output += buffer;
  }
  const int status = gSystem->ClosePipe(pipe);
  if (status != 0) {
    std::cerr << "ERROR: aao_xsec failed at phi = " << phi << " deg\n"
              << "Command: " << cmd << "\n"
              << "Output:\n" << output << std::endl;
    return false;
  }

  std::istringstream in(output);
  if (!(in >> sigma)) {
    std::cerr << "ERROR: cannot parse aao_xsec output at phi = " << phi
              << " deg. Output:\n" << output << std::endl;
    return false;
  }

  return true;
}

TString SafeNumber(double value)
{
  TString s;
  s.Form("%.4g", value);
  s.ReplaceAll("-", "m");
  s.ReplaceAll("+", "p");
  s.ReplaceAll(".", "p");
  return s;
}

} // namespace

void plot_xsec_phi(double xB, double Q2, double t, double beamEnergy,
                   double phiMin = 0.0, double phiMax = 360.0,
                   double phiStep = 5.0,
                   const char *outPrefix = "")
{
  if (phiStep <= 0.0) {
    std::cerr << "ERROR: phiStep must be positive." << std::endl;
    return;
  }
  if (phiMax < phiMin) {
    std::cerr << "ERROR: phiMax must be >= phiMin." << std::endl;
    return;
  }

  TString prefix(outPrefix);
  if (prefix.Length() == 0) {
    prefix.Form("xsec_phi_xB%s_Q2%s_t%s_BE%s",
                SafeNumber(xB).Data(), SafeNumber(Q2).Data(),
                SafeNumber(t).Data(), SafeNumber(beamEnergy).Data());
  }

  std::vector<double> phis;
  std::vector<double> sigmas;

  for (double phi = phiMin; phi <= phiMax + 0.5 * phiStep; phi += phiStep) {
    double sigma = 0.0;
    if (!RunAAOXsec(xB, Q2, t, phi, beamEnergy, sigma)) {
      std::cerr << "Stopped before plotting because one point failed." << std::endl;
      return;
    }
    phis.push_back(phi);
    sigmas.push_back(sigma);
    std::cout << "phi = " << phi << " deg, sigma = " << sigma << std::endl;
  }

  const TString datName = prefix + ".dat";
  std::ofstream dat(datName.Data());
  dat << "# xB " << xB << "\n"
      << "# Q2 " << Q2 << "\n"
      << "# t " << t << "\n"
      << "# BeamEnergy " << beamEnergy << "\n"
      << "# phi_deg sigma\n";
  for (size_t i = 0; i < phis.size(); ++i) {
    dat << phis[i] << " " << sigmas[i] << "\n";
  }
  dat.close();

  gStyle->SetOptStat(0);
  auto *canvas = new TCanvas("c_xsec_phi", "AAO MC cross section vs phi", 900, 650);
  canvas->SetGrid();

  auto *graph = new TGraph(static_cast<int>(phis.size()), phis.data(), sigmas.data());
  graph->SetName("g_xsec_phi");
  graph->SetTitle("");
  graph->SetMarkerStyle(20);
  graph->SetMarkerSize(0.9);
  graph->SetMarkerColor(kBlue + 1);
  graph->SetLineColor(kBlue + 1);
  graph->SetLineWidth(2);
  graph->Draw("ALP");

  graph->GetXaxis()->SetTitle("#phi [deg]");
  graph->GetYaxis()->SetTitle("#sigma_{0} used by MC");
  graph->GetXaxis()->CenterTitle();
  graph->GetYaxis()->CenterTitle();

  auto *legend = new TLegend(0.16, 0.78, 0.72, 0.90);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(graph,
                   Form("x_{B}=%.4g, Q^{2}=%.4g GeV^{2}, t=%.4g GeV^{2}, E=%.4g GeV",
                        xB, Q2, t, beamEnergy),
                   "lp");
  legend->Draw();

  canvas->SaveAs(prefix + ".png");
  canvas->SaveAs(prefix + ".pdf");

  std::cout << "\nSaved:\n"
            << "  " << datName << "\n"
            << "  " << prefix << ".png\n"
            << "  " << prefix << ".pdf\n";
}
