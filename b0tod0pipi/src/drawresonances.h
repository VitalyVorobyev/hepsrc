#ifndef DRAWRESONANCES_H
#define DRAWRESONANCES_H

#include "../../libDalitz/src/libdalitz.h"

#include <vector>
#include <string>

class DrawResonances : public DalitzPhaseSpace{
public:
  DrawResonances(const double& mmo, const double& mcha, const double& mchb, const double& mchc);
  DrawResonances(DalitzModel* model);
  void Draw(void);
  void Draw(const int i);

private:
  std::vector< std::string > names;
  std::vector< std::string > labels;
  std::string prefix;
  std::string exten;
  std::string treename;

  std::string abaxis;
  std::string acaxis;
  std::string bcaxis;

  double sum_mass_sq;
};

#endif // DRAWRESONANCES_H
