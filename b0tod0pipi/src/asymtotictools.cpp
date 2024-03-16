#include "asymtotictools.h"

#include <iostream>
#include <cmath>

using namespace std;

AsymtoticTools::AsymtoticTools(SymDalitzModel* model,const int gsize, const int phsize):
  m_dpbin(new SmallDPBin(model)),
  m_model(model),
  m_norm(1),
  m_phsize(phsize)
{
  m_dpbin->SetGridSize(gsize);
  m_h1d = new TH1I("h1","h1",phsize,-M_PI,M_PI);
  const double mmin = model->mABsq_min();
  const double mmax = model->mABsq_max();
  m_h2ph = new TH2D("hph","hph",gsize,mmin,mmax,gsize,mmin,mmax);
  m_h2wg = new TH2D("hwg","hwg",gsize,mmin,mmax,gsize,mmin,mmax);
  m_tree = new TTree("tree","tree");

//  initN();
}

//int AsymtoticTools::initN(void){
//  for(int i=0; i<m_nbins; i++){
//    m_N0.push_back(0);
//    m_N1.push_back(0);
//  }
//  return m_N0.size();
//}

//int AsymtoticTools::SetGridSize(const int i);


TTree* AsymtoticTools::GetTree(void) const{
    return (TTree*)m_tree->Clone();
}

int AsymtoticTools::Fill(void){
  const int nbins = m_model->GetNBins();
  double Kp[nbins];
  double Km[nbins];
  double sqKK[nbins];
  for(int i=0; i<nbins; i++){Kp[i] = 0; Km[i] = 0; sqKK[i] = 0;}

  double mp, mm;
  double deld;
  double P;
  double Pb;
  double wght;
  int dpbin;

  m_tree->Branch("mp",&mp,"mp/D");
  m_tree->Branch("mm",&mm,"mm/D");
  m_tree->Branch("deld",&deld,"deld/D");
  m_tree->Branch("p",&P,"p/D");
  m_tree->Branch("pb",&Pb,"pb/D");
  m_tree->Branch("wght",&wght,"wght/D");
  m_tree->Branch("dpbin",&dpbin,"dpbin/I");

  int ncounts = 0;
//  const double Norm = 4.*1.78943e+07;
  double norm = 0;
  const int gsize = m_dpbin->GSize();
  for(int i=0; i<gsize; i++){
    for(int j=0; j<i; j++){
      if(m_dpbin->SetBin(i,j)) continue;
      m_dpbin->GetCurrentPoint(mp,mm); ncounts++;
      deld  = m_dpbin->Phase();
      P     = m_dpbin->P();// /Norm;
      Pb    = m_dpbin->Pb();// /Norm;
      wght  = m_dpbin->Weight();
//      wght /= Norm; wght /= Norm;
      dpbin = abs(m_model->GetBin(mp,mm));
      m_tree->Fill();

      Kp[dpbin-1] += P;
      Km[dpbin-1] += Pb;
      sqKK[dpbin-1] += sqrt(P*Pb);
      norm += P + Pb;

      cout << i << " " << j << " " << deld << " " << wght << endl;
      m_h1d->Fill(deld,wght);
      const int bin = m_h2ph->GetBin(i+1,j+1);
      m_h2ph->SetBinContent(bin,deld);
      m_h2wg->SetBinContent(bin,wght);
      const int binb = m_h2ph->GetBin(j+1,i+1);
      m_h2ph->SetBinContent(binb,-deld);
      m_h2wg->SetBinContent(binb,wght);
    }
  }
  m_h1d->Print();
  m_h2ph->Print();
  m_h2wg->Print();

  cout << "Dalitz plot scanning complete. Norm = " << norm << endl;
  for(int i=0; i<nbins; i++){
    Kp[i] /= norm; Km[i] /= norm; sqKK[i] /= norm;
    cout <<"bin " << i+1 << ": ";
    cout << "K+ = " << Kp[i] << ", ";
    cout << "K- = " << Km[i] << ", ";
    cout << "SummKK = " << sqKK[i] << ", ";
    cout << "sqrtKK = " << sqrt(Kp[i]*Km[i]) << endl;
  }

  return ncounts;
}

int AsymtoticTools::GetHist(   TH1I* h) const { h = (TH1I*)m_h1d->Clone();  return 0;}
int AsymtoticTools::GetMap(    TH2D* h) const { h = (TH2D*)m_h2ph->Clone(); return 0;}
int AsymtoticTools::GetWeights(TH2D* h) const { h = (TH2D*)m_h2wg->Clone(); return 0;}
TH1I* AsymtoticTools::GetHist(void)    const {return (TH1I*)m_h1d->Clone();}
TH2D* AsymtoticTools::GetMap(void)     const {return (TH2D*)m_h2ph->Clone();}
TH2D* AsymtoticTools::GetWeights(void) const {return (TH2D*)m_h2wg->Clone();}

