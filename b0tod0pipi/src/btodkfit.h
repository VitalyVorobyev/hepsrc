#ifndef BTODKFIT_H
#define BTODKFIT_H

#include "dpevt.h"

#include "mylibs/libDalitzFit/src/ubinfcn.h"
#include "mylibs/libDalitzFit/src/ubingammafit.h"
#include "mylibs/libDalitz/src/btodk_dtokspipi_model.h"

#include <string>
#include <vector>

#include "TChain.h"

typedef UBinFCN<BtoDK_DtoKspipi_Model,DPEvt> gfcn;
typedef UBinGammaFit<gfcn> gfit;

///
/// \brief The BtoDKfit class. Unbinned fit of the angle gamma in B- -> D0 K-, D0 -> Ks0 pi+ pi- decays
///
class BtoDKfit{
public:
  BtoDKfit(const std::string& tname,const std::vector<std::string>& flist);
  const mnt::MnUserParameterState& Fit(void);

private:
  BtoDK_DtoKspipi_Model* m_model;
  gfcn*   m_fcn;
  gfit*   m_fit;
  TChain* m_tree;
  DPEvt*  m_evt;
};

#endif // BTODKFIT_H
