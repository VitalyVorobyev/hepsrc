#ifndef TOYTOOLS_H
#define TOYTOOLS_H

#include "mylibs/libDalitz/src/libdalitz.h" 

//#include "RooDataSet.h"

#include <string>
#include <vector>

/// \brief Class for keeping various routines for double Dalitz toy study
///
/// InitTools() should be called before using

class ToyTools{
public:
  ToyTools() {}
  static void InitTools(void);

  static void DrawBins(AbsDalitzModel& model,const std::string& label);
  static void DrawKspipiBins(void);
  static void DrawB0toD0pipiBins(void);
  static void DrawKuzBins(void);

  static void GenerateKuzminPlot(const unsigned NEvents,std::vector<double>& mABsqv,std::vector<double>& mACsqv);
  static void GenerateDalitzPlot(const AbsDalitzModel& model, const int NEvents,std::vector<double>& mABsqv,std::vector<double>& mACsqv);
  static void GenerateKspipiDP(const int NEvents,std::vector<double>& mABsqv,std::vector<double>& mACsqv);
  static void Generate_BDK_DKspipiDP(const int NEvents, std::vector<double>& mABsqv, std::vector<double>& mACsqv, std::vector<int> &flvv);
  static void GenerateEachResonance(DalitzModel& model, const int NEvents);

  static void GetAKuzBr(void);

  static void SaveDPTree(const std::string& label,const std::vector<double>& mABsqv,const std::vector<double>& mACsqv);
  static void SaveDPTree(const std::string& label,const std::vector<double>& mABsqv,const std::vector<double>& mACsqv, std::vector<int> &flvv);

//  static void SaveFlvGenTree(const std::string& label,const RooDataSet* ds);
//  static void SaveCPGenTree(const std::string& label,const RooDataSet* ds);
//  static void SaveDblDlzGenTree(const std::string& label,const RooDataSet* ds);

  static void addVector(std::vector<double>& v1,const std::vector<double>& v2);

  static void GammaUBFit(void);

//  /// untagged lifetime fit
//  static void LifetimeFit(void);
//  /// tagget lifetime or/and dm fit
//  static void FlvFit(void);
//  /// tagget CPV fit with D decays into CP eigenstates
//  static void CPFit(const bool cpvfit);
//  /// fit of double Dalitz distribution
//  static void DblDlzFit(const bool cpvfit);
//  ///
//  static void SimFit(const bool fixphases);

//  static void DblDlzToyGenAndFit(void);

private:
};

#endif // TOYTOOLS_H

