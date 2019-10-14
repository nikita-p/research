#ifndef Cmd3KF_H
#define Cmd3KF_H

#include "KinFitter/TKinFitter.h"
#include "KinFitter/TFitParticlePThetaPhi.h"
#include "KinFitter/TFitParticlePxPyPz.h"
#include "KinFitter/TFitConstraintEp.h"

#include "TVector3.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"

#include <vector>

#include "CovPtPhiTheta2PxPyPz.C"
#include "CovPPhiTheta2PxPyPz.C"

class KFParticle {
 public:
  TLorentzVector P;
  TMatrixD Cov;
  KFParticle():Cov(3,3){};
};

/*
typedef struct KFParticle {

  TLorentzVector P;
  float ErrorMatrix[3][3];
};
*/
#endif
