#include "drawresonances.h"

#include <cstdlib>
#include <iostream>
#include <sstream>

//#include "TFile.h"
#include "TChain.h"
#include "TCanvas.h"
#include "TH2.h"
#include "TROOT.h"
#include <TPad.h>

DrawResonances::DrawResonances(const double &mmo, const double &mcha, const double &mchb, const double &mchc):
  DalitzPhaseSpace(mmo,mcha,mchb,mchc),
  prefix("/home/vitaly/B0toD0pipi/B0toD0pipiFeas/DPTree_d0pipi"),
  exten(".root"),
  treename("tree"),
  sum_mass_sq(mmo*mmo+mcha*mcha+mchb*mchb+mchc*mchc)
{
  names.push_back("_d2");      labels.push_back("D*2(2460)");
  names.push_back("_d0");      labels.push_back("D*0(2400)");
  names.push_back("_dv1");     labels.push_back("Dv*(2010)");
  names.push_back("_dJ");      labels.push_back("D*J(2760)");
  names.push_back("_rho770");  labels.push_back("rho(770)");
  names.push_back("_rho1450"); labels.push_back("rho(1450)");
  names.push_back("_rho1700"); labels.push_back("rho(1700)");
  names.push_back("_f2");      labels.push_back("f2(1270)");
  names.push_back("_f0bugg");  labels.push_back("f0(500)");
  names.push_back("_f0flatte");labels.push_back("f0(980)");
  names.push_back("_f02020");  labels.push_back("f0(2020)");
  names.push_back("_nr");      labels.push_back("NR");
  names.push_back("_rho770_dv1"); labels.push_back("rho(770) + Dv*(2010)");
  names.push_back("");         labels.push_back("B^{0}#rightarrowD^{0}#pi^{+}#pi^{-} Dalitz plot");
}

DrawResonances::DrawResonances(DalitzModel* model):
  DrawResonances(model->mM(),model->mA(),model->mB(),model->mC())
{
  abaxis = std::string(model->ABaxis());
  acaxis = std::string(model->ACaxis());
  bcaxis = std::string(model->BCaxis());
}

void DrawResonances::Draw(void){
  for(int i=0; i<(int)names.size(); i++) Draw(i);
  return;
}

void DrawResonances::Draw(const int i){
  const std::string fname = prefix + names[i] + exten;
  TChain* tree = new TChain(treename.c_str());
  tree->Add(fname.c_str());
  std::stringstream out;

  const std::string cname = "c" + names[i];
  TCanvas* c1 = new TCanvas(cname.c_str(),cname.c_str(),800,400);
  c1->Draw();
  TPad* pad1 = new TPad("pad1","pad1",0.0,0.0,0.5,0.99);
  TPad* pad2 = new TPad("pad2","pad2",0.5,0.0,0.99,0.99);
  pad1->Draw();
  pad2->Draw();
  pad1->cd();
  tree->Draw("mm:mp");
  TH2F *histo1 = (TH2F*)gROOT->FindObject("htemp");
  histo1->SetName("histo1");
  out.str(""); out << "m_{AC}^{2}:m_{AC}^{2} for " << labels[i];
  histo1->SetTitle(out.str().c_str());
  histo1->GetXaxis()->SetRangeUser(mACsq_min(),mACsq_max());
  histo1->GetYaxis()->SetRangeUser(mABsq_min(),mABsq_max());
  histo1->GetXaxis()->SetTitle(acaxis.c_str());
  histo1->GetYaxis()->SetTitle(abaxis.c_str());
  histo1->GetXaxis()->SetLabelSize(0.05);
  histo1->GetYaxis()->SetLabelSize(0.05);
  histo1->GetXaxis()->SetTitleSize(0.05);
  histo1->GetYaxis()->SetTitleSize(0.05);
  histo1->GetXaxis()->SetTitleOffset(0.85);
  histo1->GetYaxis()->SetTitleOffset(0.85);
  gPad->Update();

  pad2->cd();
  out.str(""); out << sum_mass_sq << "-mm-mp:mp";
  tree->Draw(out.str().c_str());
  TH2F *histo2 = (TH2F*)gROOT->FindObject("htemp");
  histo2->SetName("histo2");
  out.str(""); out << "m_{AC}^{2}:m_{BC}^{2} for " << labels[i];
  histo2->SetTitle(out.str().c_str());
  histo2->GetXaxis()->SetRangeUser(mACsq_min(),mACsq_max());
  histo2->GetYaxis()->SetRangeUser(mBCsq_min(),mBCsq_max());
//  std::cout << mBC_min() << " " << mBCsq_max() << " " << mAC_min() << " " << mACsq_max() << std::endl;
  histo2->GetXaxis()->SetTitle(acaxis.c_str());
  histo2->GetYaxis()->SetTitle(bcaxis.c_str());
  histo2->GetXaxis()->SetLabelSize(0.05);
  histo2->GetYaxis()->SetLabelSize(0.05);
  histo2->GetXaxis()->SetTitleSize(0.05);
  histo2->GetYaxis()->SetTitleSize(0.05);
  histo2->GetXaxis()->SetTitleOffset(0.85);
  histo2->GetYaxis()->SetTitleOffset(0.85);
  gPad->Update();

  out.str(""); out << "pics/DP" << names[i] << ".eps";
  const std::string epsname = out.str();
  c1->Print(epsname.c_str());
  out.str(""); out << "evince " << epsname << " &";
  system(out.str().c_str());
  return;
}
