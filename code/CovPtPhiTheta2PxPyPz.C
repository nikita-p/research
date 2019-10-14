#ifndef CovPtPhiTheta2PxPyPz_H
#define CovPtPhiTheta2PxPyPz_H 

#include "TMatrixDSym.h"
#include "TMath.h"
#include "TMatrixD.h"

TMatrixDSym CovPtPhiTheta2PxPyPz(const double* pars,const TMatrixDSym& cov){
  //from Pt,phi,theta
  double ptr=pars[0];
  double sphi=TMath::Sin(pars[1]),cphi=TMath::Cos(pars[1]);
  double ctg=TMath::Cos(pars[2])/TMath::Sin(pars[2]);

  static TMatrixD tran(3,3);

 
  double dpndpo[9]=
    {cphi,-ptr*sphi,0,
     sphi,ptr*cphi,0,
     ctg,0,ptr*(-1.-ctg*ctg)};

  for(int i=0;i<9;i++)
    tran(i/3,i%3)=dpndpo[i];
  
  TMatrixDSym newcov=cov;

  return newcov.Similarity(tran);

}

TMatrixDSym CovPtPhiTheta2PxPyPz_incorrect(const double* pars,const TMatrixDSym& cov,int charge){
  //from Pt,phi,theta
  double ptr=pars[0];
  double sphi=TMath::Sin(pars[1]),cphi=TMath::Cos(pars[1]);
  double ctg=TMath::Cos(pars[2])/TMath::Sin(pars[2]);

  static TMatrixD tran(3,3);

 
  double dpndpo[9]=
    {cphi,-ptr*sphi,0,
     sphi,ptr*cphi,0,
     ctg,0,ptr*(-1.-ctg*ctg)};

  for(int i=0;i<9;i++)
    tran(i/3,i%3)=dpndpo[i];
  
  TMatrixDSym newcov=cov;
  newcov(1,0)=newcov(0,1)=TMath::Sign(1,charge)*cov(0,1);
  newcov(2,0)=newcov(0,2)=TMath::Sign(1,charge)*cov(0,2);

  return newcov.Similarity(tran);

}

TMatrixDSym CovPtPhiTheta2PxPyPz(const double* pars,const float cov[3][3]){
  static TMatrixDSym mcov(3);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      mcov(i,j)=cov[i][j];
  return CovPtPhiTheta2PxPyPz(pars,mcov);
}

TMatrixDSym CovPtPhiTheta2PxPyPz_incorrect(const double* pars,const float cov[3][3],int charge){
  static TMatrixDSym mcov(3);
  for(int i=0;i<3;i++)
    for(int j=0;j<3;j++)
      mcov(i,j)=cov[i][j];
  return CovPtPhiTheta2PxPyPz_incorrect(pars,mcov,charge);
}
#endif

