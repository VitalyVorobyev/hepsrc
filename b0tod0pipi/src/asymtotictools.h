#ifndef ASYMTOTICTOOLS_H
#define ASYMTOTICTOOLS_H

#include "mylibs/libDalitz/src/smalldpbin.h"
#include "mylibs/libDalitz/src/symdalitzmodel.h"

#include "TH1.h"
#include "TH2.h"
#include "TTree.h"

class AsymtoticTools{
public:
  AsymtoticTools(SymDalitzModel* model,const int gsize = 1000, const int phsize = 100);

//  int SetGridSize(const int i);
  int Fill(void);

  int GetHist(   TH1I* h) const;
  int GetMap(    TH2D* h) const;
  int GetWeights(TH2D* h) const;

  TTree* GetTree(void) const;
  TH1I* GetHist(void)    const;
  TH2D* GetMap(void)     const;
  TH2D* GetWeights(void) const;

private:
//  int initN(void);
  SmallDPBin* m_dpbin;
  SymDalitzModel* m_model;

  TH1I* m_h1d;
  TH2D* m_h2ph;
  TH2D* m_h2wg;
  TTree* m_tree;
  double m_norm;
  int m_phsize;
  double x,y;
//  int m_nbins;
  std::vector<double> m_N0;
  std::vector<double> m_N1;
};

#endif // ASYMTOTICTOOLS_H
