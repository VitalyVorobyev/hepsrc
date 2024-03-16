#include "dlzbinstruct.h"

#include <fstream>
#include <iostream>

using namespace std;

DlzBinStruct::DlzBinStruct(const int nbins):
  m_nbins(nbins)
{
  for(int i=0; i<m_nbins; i++){ m_C.push_back(0); m_S.push_back(0); m_Kp.push_back(0); m_Kn.push_back(0);}
}

DlzBinStruct::DlzBinStruct(const string& fname){
  m_nbins = GetParams(fname);
}

void DlzBinStruct::Clear(void){
  m_C.clear(); m_S.clear(); m_Kp.clear(); m_Kn.clear();
  m_nbins = 0;
}

int DlzBinStruct::GetParams(const string& fname){
  Clear();
  cout << "Getting parameters from " << fname << endl;
  ifstream ifile(fname.c_str(),ifstream::in);
  if(!ifile.is_open()) return -1;
  string line;
  int bin;
  double c,s,kp,kn,q;
  int flag = 6;
  while(flag == 6){
    int i = m_S.size();
    getline(ifile,line);
    flag = sscanf(line.c_str(),"%d: C = %lf, S = %lf, K = %lf, Kb = %lf, Q = %lf",&bin,&c,&s,&kp,&kn,&q);
    if(i != (bin-1) || flag != 6) break;
    m_C.push_back(c); m_S.push_back(s); m_Kp.push_back(kp); m_Kn.push_back(kn);
    cout << i+1 << " C: " << m_C[i] << " S: " << m_S[i] << " Kp: " << m_Kp[i] << " Kn: " << m_Kn[i] << endl;
  }
  cout << "Done!" << endl;
  return m_nbins = m_C.size();
}

bool DlzBinStruct::CheckBinNumber(const int bin) const{
  if(bin>0 && bin<=m_nbins) return true;
  else cout << "Wrong bin number detected: " << bin << endl;
  return false;
}

void DlzBinStruct::SetKp(const int bin,const double& v){
  if(!CheckBinNumber(bin)) return;
  m_Kp[bin-1] = v;
}

void DlzBinStruct::SetKn(const int bin,const double& v){
  if(!CheckBinNumber(bin)) return;
  m_Kn[bin-1] = v;
}

void DlzBinStruct::SetC(const int bin,const double& v){
  if(!CheckBinNumber(bin)) return;
  m_C[bin-1] = v;
}

void DlzBinStruct::SetS(const int bin,const double& v){
  if(!CheckBinNumber(bin)) return;
  m_S[bin-1] = v;
}

double DlzBinStruct::Kp(const int bin) const{
  if(!CheckBinNumber(abs(bin))) return 0;
  return bin>0 ? m_Kp[bin-1] : m_Kn[-bin-1];
}

double DlzBinStruct::Kn(const int bin) const{
  if(!CheckBinNumber(abs(bin))) return 0;
  return bin>0 ? m_Kn[bin-1] : m_Kp[-bin-1];
}

double DlzBinStruct::C(const int bin)  const{
  if(!CheckBinNumber(abs(bin))) return 0;
  return m_C[abs(bin)-1];
}

double DlzBinStruct::S(const int bin)  const{
  if(!CheckBinNumber(abs(bin))) return 0;
  return bin>0 ? m_S[bin-1] : -m_S[-bin-1];
}

