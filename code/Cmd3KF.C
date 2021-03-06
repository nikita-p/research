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
  //TFitConstraintMGaus KMass("KMass", "KMass", 0, 0, 0.4976, 0.01);//

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
  //KMass.addParticles1(&FitParticle[0], &FitParticle[1]);//

  fitter.addConstraint( &pX );
  fitter.addConstraint( &pY );
  fitter.addConstraint( &pZ );
  fitter.addConstraint( &E  );
  //fitter.addConstraint( &KMass );//

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
    RecoParticle.push_back(OutPart);
  }

  if ( Status == 0 ) ChiSq =  fitter.getS();

  return ChiSq;
};

TMatrixD GetTrErrorMatrix(TLorentzVector P, float ErrMat[3][3]) {

  TMatrixDSym Sym(3);

  TMatrixD Cov(3,3);

  double pars[3];

  pars[0] = P.P()*sin(P.Theta())/1e3;
  pars[1] = P.Phi();
  pars[2] = P.Theta();

  ErrMat[0][0] *= 1e-6;
  ErrMat[1][0] *= 1e-3;
  ErrMat[2][0] *= 1e-3;
  ErrMat[0][1] *= 1e-3;
  ErrMat[0][2] *= 1e-3;

  Sym = CovPtPhiTheta2PxPyPz(pars,ErrMat);
/*
  Cov(0,0) = Sym(0,0);//1E6;
  Cov(0,1) = Sym(0,1);//1E6;
  Cov(0,2) = Sym(0,2);//1E6;
  Cov(1,0) = Cov(0,1);
  Cov(1,1) = Sym(1,1);//1E6;
  Cov(1,2) = Sym(1,2);//1E6;
  Cov(2,0) = Cov(0,2);
  Cov(2,1) = Cov(1,2);
  Cov(2,2) = Sym(2,2);//1E6;
*/
  return Sym;//Cov;

};

TMatrixD GetPhErrorMatrix(TLorentzVector P, double SigmaPhE, double SigmaPhTh, double SigmaPhPhi) {

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
/*
  Cov(0,0) = Sym(0,0);
  Cov(0,1) = Sym(0,1);
  Cov(0,2) = Sym(0,2);
  Cov(1,0) = Cov(0,1);
  Cov(1,1) = Sym(1,1);
  Cov(1,2) = Sym(1,2);
  Cov(2,0) = Cov(0,2);
  Cov(2,1) = Cov(1,2);
  Cov(2,2) = Sym(2,2);
*/
  return Sym;//Cov;
}
