#include "toytools.h"

#include "mylibs/libDalitz/src/kspipimodel.h"
#include "mylibs/libDalitz/src/btodk_dtokspipi_model.h"
#include "mylibs/libDalitz/src/dalitzmcintegral.h"
#include "mylibs/libDalitz/src/kuzbtodpipiamp.h"

#include "drawbdparams.h"
#include "drawresonances.h"
#include "btodkfit.h"

//#include "toydbldalitzgen.h"
//#include "toydatamanager.h"
//#include "toydbldalitzfit.h"
//#include "toyvars.h"
//#include "toypainter.h"

#include <fstream>

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TRandom3.h"

#include "RooArgSet.h"
#include "RooConstVar.h"
#include "RooTruthModel.h"
#include "RooGenericPdf.h"
#include "Roo1DTable.h"
#include "RooMCStudy.h"

#include <string>
#include <sstream>

using namespace std;

//typedef vector<double> vectd;
typedef KuzBtoDpipiAmp KAmp;

void ToyTools::InitTools(void){
//  ToyCat::InitCat();
//  ToyVars::InitVars();
//  ToyResModel::InitResolution();
}

void ToyTools::DrawBins(AbsDalitzModel& model,cstr& label){
  ModelIntegral integrator(&model);
  integrator.SetGridSize(1000);
  DrawBDParams drawpar;
  vectd C,S,K,Kb;
  integrator.Calculate(label,C,S,K,Kb);
  cstr fname  = label + string("_binning.txt");
  cstr dpname = label + string("_binning");
  drawpar.DrawBinsmABmAC(fname,dpname);
  drawpar.DrawBinsmABmBC(fname,dpname);
  const string csname = label + string("_CS");
  drawpar.DrawCS(C,S,csname);
  const string kkbname = label + string("_K");
  drawpar.DrawK(K,Kb,kkbname);
  return;
}

void ToyTools::DrawKspipiBins(void){
  KspipiModel model;
  cstr label("kspipi");
  ToyTools::DrawBins(model,label);
  return;
}

void ToyTools::DrawB0toD0pipiBins(void){
  B0toD0pipiModel model;
  cstr label("d0pipi_D2_rho");
  ToyTools::DrawBins(model,label);
  return;
}

void ToyTools::DrawKuzBins(void){
  KAmp model;
  cstr label("d0pipi_kuz");
  ToyTools::DrawBins(model,label);
  return;
}

void ToyTools::GenerateDalitzPlot(const AbsDalitzModel& model, const int NEvents,vectd& mABsqv,vectd& mACsqv){
  DalitzGenerator generator(&model);
  generator.Generate(NEvents,mABsqv,mACsqv);
  return;
}

void ToyTools::GenerateKuzminPlot(const unsigned NEvents,vectd& mABsqv,vectd& mACsqv){
  KAmp model;
  DalitzGenerator generator(&model);
  generator.Generate(NEvents,mABsqv,mACsqv);
  return;
}

void ToyTools::GetAKuzBr(void){
  KAmp model;
  DalitzMCIntegral mcint(&model);
  vectd vals;
  vectd errs;
  mcint.CalcBranchings(vals,errs);
}

void ToyTools::GenerateKspipiDP(const int NEvents,vectd& mABsqv,vectd& mACsqv){
  KspipiModel model;
  ToyTools::GenerateDalitzPlot(model,NEvents,mABsqv,mACsqv);
  return;
}

void ToyTools::addVector(vectd& v1,const vectd& v2){
  v1.insert(v1.end(),v2.begin(),v2.end());
}

void ToyTools::Generate_BDK_DKspipiDP(const int NEvents,vectd& mABsqv,vectd& mACsqv, vector<int>& flvv){
  mABsqv.clear(); mACsqv.clear(); flvv.clear();
  BtoDK_DtoKspipi_Model model(70,120,0.1);
  DalitzGenerator generator(&model);
  TRandom3 rndm;
  rndm.SetSeed(0);
  int flv = 1;
  double mAB,mAC;
  for(int i=0; i<NEvents; i++){
    flv = rndm.Rndm()-0.5 > 0 ? 1 : -1;
    model.SetFlv(flv);
    generator.Generate(mAB,mAC);
    mABsqv.push_back(mAB);
    mACsqv.push_back(mAC);
    flvv.push_back(flv);
  }
  return;
}

void ToyTools::SaveDPTree(const string& label,const vector<double>& mABsqv,const vector<double>& mACsqv){
  stringstream out;
  out.str(""); out << "DPTree_" << label << ".root";
  TFile* file = new TFile(out.str().c_str(),"RECREATE");
  TTree* tree = new TTree("tree","tree");
  double mABsq, mACsq;
  tree->Branch("mp",&mABsq,"mp/d");
  tree->Branch("mm",&mACsq,"mm/d");
  for(unsigned i=0; i<mABsqv.size(); i++){
    mABsq = mABsqv[i]; mACsq = mACsqv[i];
    tree->Fill();
  }
  tree->Write();
  file->Close();
  return;
}

void ToyTools::SaveDPTree(const std::string& label,const std::vector<double>& mABsqv,const std::vector<double>& mACsqv, std::vector<int> &flvv){
  stringstream out;
  out.str(""); out << "DPTree_" << label << ".root";
  TFile* file = new TFile(out.str().c_str(),"RECREATE");
  TTree* tree = new TTree("tree","tree");
  double mABsq, mACsq;
  int flv;
  tree->Branch("mp",&mABsq,"mp/D");
  tree->Branch("mm",&mACsq,"mm/D");
  tree->Branch("flv",&flv,"flv/I");
  for(unsigned i=0; i<mABsqv.size(); i++){
    mABsq = mABsqv[i]; mACsq = mACsqv[i]; flv = flvv[i];
    tree->Fill();
  }
  tree->Write();
  file->Close();
  return;
}

void ToyTools::GenerateEachResonance(DalitzModel& model, const int NEvents){
  DalitzGenerator generator(&model);
  const int NRes = model.ResNum();
  stringstream out;
  vector<int> subres;
  vectd mABsqv;
  vectd mACsqv;
  for(int i=0; i<NRes; i++){
    subres.clear(); subres.push_back(i);
    if(model.ResName(i) == string("rho(770)")){// rho and omega together
      i++; subres.push_back(i);
    }
    mABsqv.clear(); mACsqv.clear();
    model.SetRVec(subres);
    generator.Generate(NEvents,mABsqv,mACsqv);
    out.str(""); out << "DPTree" << i;
    SaveDPTree(out.str(),mABsqv,mACsqv);
  }
  subres.clear(); model.SetRVec(subres);
  return;
}

//void ToyTools::SaveFlvGenTree(const std::string& label,const RooDataSet* ds){
//  stringstream out;
//  out.str(""); out << "data/GenTreeFlv" << label << ".root";
//  TFile* file = new TFile(out.str().c_str(),"RECREATE");
//  TTree* tree= new TTree("tree","tree");
//  double dt;
//  int binb;
//  int bind = 0;
//  int flv;
//  int cp = 0;
//  tree->Branch("dt",&dt,"dt/D");
//  tree->Branch("binb",&binb,"binb/I");
//  tree->Branch("bind",&bind,"bind/I");
//  tree->Branch("flv",&flv,"flv/I");
//  tree->Branch("cp",&cp,"cp/I");
//
//  const int NTot = ds->sumEntries();
//  for(int i=0; i<NTot; i++){
//    const RooArgSet* aset = ds->get(i);
//    dt  = aset->getRealValue("dt");
//    binb= aset->getCatIndex("bin");
//    flv = aset->getCatIndex("flv");
//    tree->Fill();
//  }
//  tree->Write();
//  file->Close();
//}

//void ToyTools::SaveCPGenTree(const std::string& label,const RooDataSet* ds){
//  stringstream out;
//  out.str(""); out << "data/GenTreeCP" << label << ".root";
//  TFile* file = new TFile(out.str().c_str(),"RECREATE");
//  TTree* tree= new TTree("tree","tree");
//  double dt;
//  int binb;
//  int bind = 0;
//  int flv;
//  int cp;
//  tree->Branch("dt",  &dt,  "dt/D");
//  tree->Branch("binb",&binb,"binb/I");
//  tree->Branch("bind",&bind,"bind/I");
//  tree->Branch("flv", &flv, "flv/I");
//  tree->Branch("cp",  &cp,  "cp/I");
//
//  const int NTot = ds->sumEntries();
//  for(int i=0; i<NTot; i++){
//    const RooArgSet* aset = ds->get(i);
//    dt  = aset->getRealValue("dt");
//    binb= aset->getCatIndex("binb");
//    flv = aset->getCatIndex("flv");
//    cp  = aset->getCatIndex("cp");
//    tree->Fill();
//  }
//  tree->Write();
//  file->Close();
//}

//void ToyTools::SaveDblDlzGenTree(const std::string& label,const RooDataSet* ds){
//  stringstream out;
//  out.str(""); out << "data/GenTreeDblDlz" << label << ".root";
//  TFile* file = new TFile(out.str().c_str(),"RECREATE");
//  TTree* tree= new TTree("tree","tree");
//  double dt;
//  int binb;
//  int bind;
//  int flv;
//  int cp = 0;
//  tree->Branch("dt",&dt,"dt/D");
//  tree->Branch("binb",&binb,"binb/I");
//  tree->Branch("bind",&bind,"bind/I");
//  tree->Branch("flv",&flv,"flv/I");
//  tree->Branch("cp",&cp,"cp/I");

//  const int NTot = ds->sumEntries();
//  for(int i=0; i<NTot; i++){
//    const RooArgSet* aset = ds->get(i);
//    dt   = aset->getRealValue("dt");
//    binb = aset->getCatIndex("binb");
//    bind = aset->getCatIndex("bind");
//    flv  = aset->getCatIndex("flv");
//    tree->Fill();
//  }
//  tree->Write();
//  file->Close();
//}

void ToyTools::GammaUBFit(void){
  const string fname("tree");
  vector<string> flist;
  flist.push_back("/home/vitaly/B0toD0pipi/B0toD0pipiFeas/DPTree_bdk.root");
  BtoDKfit fitter(fname,flist);
//  const mnt::MnUserParameterState& pstate = ;
  fitter.Fit();
}

//void ToyTools::LifetimeFit(void){
//  ToyDblDalitzFit fitter;
//  ToyDataManager dmgr;
//  RooDataSet* ds = dmgr.GetFlvDS(0);
//  fitter.FitLifetime(ds);
//  ToyPainter::DrawLifetime(ds,fitter.GetLifetimePdf());
//}

//void ToyTools::FlvFit(void){
//  ToyDblDalitzFit fitter;
//  ToyDataManager dmgr;
//  ToyVars::Fix_dm(false);
//  RooDataSet* ds = dmgr.GetFlvDS(0);
//  fitter.FitFlv(ds);
//  ToyPainter::DrawFlvFit(ds,fitter.GetFlvPdf());
//}

//void ToyTools::CPFit(const bool cpvfit){
//  ToyDblDalitzFit fitter;
//  ToyDataManager dmgr;
//  if(cpvfit){
//    fitter.FixBParams(true);
//    ToyVars::Fix_dm(true);
//    ToyVars::Fix_tau(false);
//    ToyVars::Fix_sin(false);
//    ToyVars::Fix_cos(false);
//  } else{
//    fitter.FixBParams(false);
//    fitter.FixDParams(true);
//    ToyVars::Fix_dm(true);
//    ToyVars::Fix_tau(true);
//    ToyVars::Fix_sin(true);
//    ToyVars::Fix_cos(true);
//  }

//  RooDataSet* ds = dmgr.GetCPDS(0);
//  ds->Print();
//  fitter.FitCP(ds);
//  ToyPainter::DrawCPFit(ds,fitter.GetSimPdf());
//}

//void ToyTools::DblDlzFit(const bool cpvfit){
//  ToyDblDalitzFit fitter;
//  ToyDataManager dmgr;
//  if(cpvfit){
//    fitter.FixBParams(false);
//    fitter.FixDParams(true);
//    ToyVars::Fix_dm(true);
//    ToyVars::Fix_tau(true);
//    ToyVars::Fix_sin(false);
//    ToyVars::Fix_cos(false);
//  } else{
//    fitter.FixBParams(false);
//    fitter.FixDParams(true);
//    ToyVars::Fix_dm(true);
//    ToyVars::Fix_tau(true);
//    ToyVars::Fix_sin(true);
//    ToyVars::Fix_cos(true);
//  }
//  RooDataSet* ds = dmgr.GetDblDlzDS(0);
//  ds->Print();
//  fitter.FitDblDlz(ds);
//  ToyPainter::DrawDblDlzFit(ds,fitter.GetSimPdf());
//}

//void ToyTools::SimFit(const bool fixphases){
//  ToyDblDalitzFit fitter;
//  ToyDataManager dmgr;
//  fitter.FixBParams(fixphases);
////  fitter.ReleaseBBin(7);
//  fitter.FixDParams(true);
//  ToyVars::Fix_dm(true);
//  ToyVars::Fix_tau(true);
//  ToyVars::Fix_sin(false);
//  ToyVars::Fix_cos(false);

//  RooDataSet* ds   = dmgr.GetDblDlzDS(0);
//  RooDataSet* dscp = dmgr.GetCPDS(0);
//  ds->append(*dscp);

//  fitter.FitSimul(ds);
//  ToyPainter::DrawDblDlzFit(ds,fitter.GetSimPdf());
//}

////void ToyTools::DblDlzToyGenAndFit(void){
////  ToyDblDalitzFit fitter;
////  fitter.GetDDPdf();

////}
