#include "Cmd3KF.h"
#include <iostream>

double Cmd3KF(double EBeam, vector<KFParticle> &MeasParticle, vector<KFParticle> &RecoParticle) {

  TVector3 CurParticle;

  TFitParticlePxPyPz FitParticle[10];

  TKinFitter fitter;

  int Status = -1000;
  double ChiSq = 100000;

  TFitConstraintEp pX("Px","Px conservation", 0, TFitConstraintEp::pX, 0.0);
  TFitConstraintEp pY("Py","Py conservation", 0, TFitConstraintEp::pY, 0.0);
  TFitConstraintEp pZ("Pz","Pz conservation", 0, TFitConstraintEp::pZ, 0.0);
  TFitConstraintEp  E("E" ,"E conservation" , 0, TFitConstraintEp::E, EBeam/500.0);

  for ( size_t i = 0; i < MeasParticle.size(); i++ ) {

    CurParticle = TVector3(0.001*MeasParticle[i].P.X(),
    0.001*MeasParticle[i].P.Y(),
    0.001*MeasParticle[i].P.Z());

    FitParticle[i] = TFitParticlePxPyPz(&CurParticle,0.001*MeasParticle[i].P.M(),&MeasParticle[i].Cov);
    pX.addParticle(&FitParticle[i]);
    pY.addParticle(&FitParticle[i]);
    pZ.addParticle(&FitParticle[i]);
    E.addParticle(&FitParticle[i]);

    fitter.addMeasParticle(&FitParticle[i]);
  }

  fitter.addConstraint( &pX );
  fitter.addConstraint( &pY );
  fitter.addConstraint( &pZ );
  fitter.addConstraint( &E  );

  fitter.setMaxNbIter( 50 );
  fitter.setMaxDeltaS( 5e-5 );
  fitter.setMaxF( 1e-4 );
  fitter.setVerbosity(0);

  Status = fitter.fit();

  RecoParticle.clear();
  for ( size_t i = 0; i < MeasParticle.size(); i++ ) {
    KFParticle OutPart;
    TLorentzVector POut   = (*FitParticle[i].getCurr4Vec())*1000.0;
    OutPart.P = POut;
    //    OutPart.Cov = (*FitParticle[i].getCovMatrixFit());
    RecoParticle.push_back(OutPart);
  }

  if ( Status == 0 ) ChiSq =  fitter.getS();

  return ChiSq;
};

TMatrixD GetTrErrorMatrix(TLorentzVector P, float ErrMat[3][3]) {

  TMatrixDSym Sym(3);

  TMatrixD Cov(3,3);

  double pars[3];

  pars[0] = P.P()*sin(P.Theta());
  pars[1] = P.Phi();
  pars[2] = P.Theta();

  Sym = CovPtPhiTheta2PxPyPz(pars,ErrMat);

  Cov(0,0) = Sym(0,0)/1E6;
  Cov(0,1) = Sym(0,1)/1E6;
  Cov(0,2) = Sym(0,2)/1E6;
  Cov(1,0) = Cov(0,1);
  Cov(1,1) = Sym(1,1)/1E6;
  Cov(1,2) = Sym(1,2)/1E6;
  Cov(2,0) = Cov(0,2);
  Cov(2,1) = Cov(1,2);
  Cov(2,2) = Sym(2,2)/1E6;

  return Cov;

};

TMatrixD GetPhErrorMatrix(TLorentzVector P, double &SigmaPhE, double &SigmaPhTh, double &SigmaPhPhi) {
  /*
  TMatrixD A(3,3);
  TMatrixD AT(3,3);
  TMatrixD Err(3,3);
  */

  double pars[3];
  float err[3][3];

  TMatrixDSym Sym(3);

  TMatrixD Cov(3,3);

  pars[0] = P.P()/1E3;
  pars[1] = P.Phi();
  pars[2] = P.Theta();

  err[0][0] = SigmaPhE*SigmaPhE/1E6;
  err[0][1] = 0.0;
  err[0][2] = 0.0;
  err[1][0] = 0.0;
  err[1][1] = SigmaPhPhi*SigmaPhPhi;
  err[1][2] = 0.0;
  err[2][0] = 0.0;
  err[2][1] = 0.0;
  err[2][2] = SigmaPhTh*SigmaPhTh;

  Sym = CovPPhiTheta2PxPyPz(pars,err);

  Cov(0,0) = Sym(0,0);
  Cov(0,1) = Sym(0,1);
  Cov(0,2) = Sym(0,2);
  Cov(1,0) = Cov(0,1);
  Cov(1,1) = Sym(1,1);
  Cov(1,2) = Sym(1,2);
  Cov(2,0) = Cov(0,2);
  Cov(2,1) = Cov(1,2);
  Cov(2,2) = Sym(2,2);
  /*
  A(0,0) =  P.X()/P.P();
  A(0,1) =  P.X()/tan(P.Theta())/1E3;
  A(0,2) = -P.X()*tan(P.Phi())/1E3;

  A(1,0) =  P.Y()/P.P();
  A(1,1) =  P.Y()/tan(P.Theta())/1E3;
  A(1,2) =  P.Y()/tan(P.Phi())/1E3;

  A(2,0) =  P.Z()/P.P();
  A(2,1) = -P.Z()/tan(P.Theta())/1E3;
  A(2,2) =  0.0;

  AT.Transpose(A);

  Err(0,0) = SigmaPhE*SigmaPhE/1E6;
  Err(0,1) = 0.0;
  Err(0,2) = 0.0;
  Err(1,0) = 0.0;
  Err(1,1) = SigmaPhTh*SigmaPhTh;
  Err(1,2) = 0.0;
  Err(2,0) = 0.0;
  Err(2,1) = 0.0;
  Err(2,2) = SigmaPhPhi*SigmaPhPhi;

  Cov = A*Err*AT;
  */
  return Cov;
}
