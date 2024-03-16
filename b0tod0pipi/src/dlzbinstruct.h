#ifndef DLZBINSTRUCT_H
#define DLZBINSTRUCT_H

#include <vector>
#include <string>

///
/// The DlzBinStruct class. Class for keeping binned Dalitz plot parameters
///

class DlzBinStruct{
public:
  DlzBinStruct(const int nbins);
  DlzBinStruct(const std::string& fname);

  void SetKp(const int bin,const double& v);
  void SetKn(const int bin,const double& v);
  void SetC(const int bin,const double& v);
  void SetS(const int bin,const double& v);

  double Kp(const int bin) const;
  double Kn(const int bin) const;
  double C(const int bin)  const;
  double S(const int bin)  const;

  int GetParams(const std::string& fname);
  void Clear(void);

private:
  bool CheckBinNumber(const int bin) const;
  int m_nbins;
  std::vector<double> m_C;
  std::vector<double> m_S;
  std::vector<double> m_Kp;
  std::vector<double> m_Kn;
};

#endif // DLZBINSTRUCT_H
