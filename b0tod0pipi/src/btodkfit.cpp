#include "btodkfit.h"

using namespace std;

BtoDKfit::BtoDKfit(const string& tname,const vector<string>& flist){
  m_tree = new TChain(tname.c_str());
  for(unsigned i=0; i<flist.size(); i++) m_tree->Add(flist[i].c_str());

  m_model = new BtoDK_DtoKspipi_Model();
  m_evt   = new DPEvt(m_tree,TEvt::In);
  m_fcn   = new gfcn(m_model,m_evt);
  m_fit   = new gfit(m_fcn);
}

const mnt::MnUserParameterState& BtoDKfit::Fit(void){
  m_fit->Fit();
  return m_fit->PState();
}

