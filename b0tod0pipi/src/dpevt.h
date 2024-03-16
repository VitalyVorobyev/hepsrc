#ifndef DPEVT_H
#define DPEVT_H

#include "mylibs/libTEvt/src/tevt.h"

///
/// \brief The DPEvt class.
///
class DPEvt : public TEvt{
public:
  DPEvt(const bool iniflag = true);
  DPEvt(TTree* tree, const int ioflag);
  DPEvt(TTree* itree, TTree* otree);
  DPEvt(const DPEvt& x);

  DPEvt& operator=(const DPEvt& x);

  double mp(void) const;
  double mm(void) const;
  int    flv(void) const;

  void Set_mp(const double& x);
  void Set_mm(const double& x);
  void Set_flv(const int& x);

private:
  void init(void);

  std::vector<int> m_ind;
  int mp_ind( void) const {return m_ind[0];}
  int mm_ind( void) const {return m_ind[1];}
  int flv_ind(void) const {return m_ind[2];}
};

#endif // DPEVT_H

