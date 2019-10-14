#ifndef CovPPhiTheta2PxPyPz_def
#define CovPPhiTheta2PxPyPz_def 1
#include "TMatrixDSym.h"
#include "TMath.h"
#include "TMatrixD.h"

TMatrixDSym CovPPhiTheta2PxPyPz(const double* pars,const TMatrixDSym& cov){
  //from P,phi,theta
  double sphi=TMath::Sin(pars[1]),cphi=TMath::Cos(pars[1]);
  double sth=TMath::Sin(pars[2]),cth=TMath::Cos(pars[2]);
  double p=pars[0];

  static TMatrixD tran(3,3);

 
  double dpndpo[9]=
    {cphi*sth,-p*sphi*sth,p*cphi*cth,
     sphi*sth,p*cphi*sth,p*sphi*cth,
     cth,0,p*(-sth)};

  for(int i=0;i<9;i++)
    tran(i/3,i%3)=dpndpo[i];
  
  TMatrixDSym newcov=cov;

  return newcov.Similarity(tran);

}

TMatrixDSym CovPPhiTheta2PxPyPz_incorrect(const double* pars,const TMatrixDSym& cov,int charge){
  //from P,phi,theta
  double sphi=TMath::Sin(pars[1]),cphi=TMath::Cos(pars[1]);
  double sth=TMath::Sin(pars[2]),cth=TMath::Cos(pars[2]);
  double p=pars[0];

  static TMatrixD tran(3,3);

 
  double dpndpo[9]=
    {cphi*sth,-p*sphi*sth,p*cphi*cth,
     sphi*sth,p*cphi*sth,p*sphi*cth,
     cth,0,p*(-sth)};

  for(int i=0;i<9;i++)
    tran(i/3,i%3)=dpndpo[i];
  
  TMatrixDSym newcov=cov;
  newcov(1,0)=newcov(0,1)=TMath::Sign(1,charge)*cov(0,1);
  newcov(2,0)=newcov(0,2)=TMath::Sign(1,charge)*cov(0,2);

  return newcov.Similarity(tran);

}

TMatrixDSym CovPPhiTheta2PxPyPz(const double* pars,const float cov[3][3]){
  static TMatrixDSym mcov(3);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      mcov(i,j)=cov[i][j];
  return CovPPhiTheta2PxPyPz(pars,mcov);
}

TMatrixDSym CovPPhiTheta2PxPyPz_incorrect(const double* pars,const float cov[3][3],int charge){
  static TMatrixDSym mcov(3);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      mcov(i,j)=cov[i][j];
  return CovPPhiTheta2PxPyPz_incorrect(pars,mcov,charge);
}

#endif

