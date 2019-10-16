#define MC_cxx
#include "MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <iostream>

#define mKs 497.614
#define mPi 139.570

using namespace std;

TLorentzVector MC::VectorCreator(double P, double Theta, double Phi, double Mass){
  TLorentzVector V;
  V.SetXYZM(0, 0, P, Mass);
  V.SetTheta(Theta);
  V.SetPhi(Phi);/*
  double eps = 0.0001;
  if( (fabs(V.Theta()-Theta)>eps) || ( (fabs(V.Phi()-Phi)>eps) && (fabs(V.Phi()-Phi-2*TMath::Pi())>eps) && (fabs(V.Phi()-Phi+2*TMath::Pi())>eps) ) || (fabs(V.M()-Mass)>eps) || (fabs(V.P()-P)>eps) ){
    cout << "P: " << V.P() << '\t' << P << endl;
    cout << "Theta: " << V.Theta() << '\t' << Theta << endl;
    cout << "Phi: " << V.Phi() << '\t' << Phi << endl;
    cout << "Mass: " << V.M() << '\t' << Mass << endl;
    cout << "Error" << endl;
    //throw 1; //выбросить исключение
  }*/
  return V;
}

TMatrixD GetJacobianMatix(TLorentzVector& P1, TLorentzVector& P2){
  TLorentzVector PK = P1 + P2;

  double P[3], Th[3], Ph[3];
  Th[0] = PK.Theta();  Th[1] = P1.Theta(); Th[2] = P2.Theta();
  Ph[0] = PK.Phi();  Ph[1] = P1.Phi(); Ph[2] = P2.Phi();
  P[0] = PK.P()*sin(Th[0]);  P[1] = P1.P()*sin(Th[1]); P[2] = P2.P()*sin(Th[2]);

  TMatrixD InvT(3,3); //Inverted transform matrix (P, phi, theta)
  double InvT_arr[9] = {cos(Ph[0]), sin(Ph[0]), 0,
                        -sin(Ph[0])/P[0], cos(Ph[0])/P[0], 0,
                        cos(Th[0])*cos(Ph[0])*sin(Th[0])/P[0], cos(Th[0])*sin(Th[0])*sin(Ph[0])/P[0], -sin(Th[0])*sin(Th[0])/P[0]  };
  InvT.SetMatrixArray( &InvT_arr[0] );

  TMatrixD R(3,6); //Right part (P, phi, theta)[1,2]
  int j;
  for(int i=1; i<=2; i++){
    j = 3*(i-1);
    R(0, j) = cos(Ph[i]); // P/P1
    R(1, j) = sin(Ph[i]); // P/phi1
    R(2, j) = 1/tan(Th[i]); // P/theta1
    R(0, j+1) = -P[i]*sin(Ph[i]); // phi/P1
    R(1, j+1) = P[i]*cos(Ph[i]); // phi/phi1
    R(2, j+1) = 0; // phi/theta1
    R(0, j+2) = 0;
    R(1, j+2) = 0;
    R(2, j+2) = -P[i]/pow(sin(Th[i]),2);
  }

  return InvT*R;
}

pair<double, double> psi_angle(double P_KS) //P [MeV]
{
  //Landau: vol.2 p.56 ex.3
  pair<double, double> p;
  double V = P_KS / sqrt(P_KS * P_KS + mKs * mKs);               //KS speed in Lab System
  double v0 = sqrt(pow(mKs / 2., 2) - pow(mPi, 2)) / (mKs / 2.); //Pi speed in center-of-mass

  if (V < v0)
  {
    p = {2 * atan((v0 / V) * sqrt(1 - V * V)), TMath::Pi()};
    return p;
  }
  if (V > v0 && V < (v0 / sqrt(1 - v0 * v0)))
  {
    p = {0, asin(sqrt((1 - V * V) / (1 - v0 * v0)))};
    return p;
  }
  p = {0, 2 * atan((v0 / V) * sqrt(1 - V * V))};

  return p;
}

double MC::pidedx(double P, double dEdX)    //calculate dEdX for pions
{
  double pidedx = 5.58030e+9/pow(P+40.,3)+2.21228e+3-3.77103e-1*P-dEdX;
  return pidedx;
}

void MC::SetOutputPath(string key)
{
  string filepath = "../outputs/" + key + "/";
  if (fChain == 0 ){
    this->path = (filepath + "failed.root");
    return;
  }

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  char label[100];
  for (Long64_t jentry=0; jentry<1;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    snprintf(label, sizeof(label), "%.2f.root", ebeam);
  }
  this->path = filepath + label;
  cout << "Out path: " << this->path << endl;
  return;
}

void MC::GetSoftPhotonsNumber(string file)
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  bool model = (fChain->GetMaximum("nsim")<1) ? false : true;
  if(!model) return;

  int N_SOFT = 0;
  double EMEAS = -1;
  ofstream o(file.c_str());

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry==0)
      EMEAS = emeas;

    if (Cut(ientry) < 0) continue;
    N_SOFT++;

    if (fabs(EMEAS-emeas)>0.01 ){
      o << EMEAS << ',' << N_SOFT << endl;
      EMEAS = emeas;
      N_SOFT = 0;
    }

  }
  o << EMEAS << ',' << N_SOFT << endl;
  return;
}

int MC::StandardProcedure(Long64_t entry, std::vector<int> goods){ //вернуть отрицательное значение, если не проходит отборы

  if(nks_total<=0) return -1; //нет каонов - нет и отбора

  double P_CUT = 2*(0.0869*emeas-36.53);  //для каждой энергии кат по импульсу(энергии) теперь будет различаться
  double minDiv = TMath::Infinity();

  int bestKs = -1; //пора искать лучший каон. Из всех, найденных процедурой, найдём лучший, с инв.массой наиболее близкой к массе каона.
  for(int i=0; i<nks_total; i++){
    if(fabs(ksminv[i]-mKs)>minDiv) continue;
    minDiv = fabs(ksminv[i]-mKs);
    bestKs = i;
  }

  if(ksvind[bestKs][0]>ksvind[bestKs][1]) // проверить, что первый записанный трек может быть больше второго (добавил на всякий случай)
    cout << "WARNING THERE: " << ksvind[bestKs][0] << '\t' << ksvind[bestKs][1] << endl;
  if( ksvind[bestKs][0]!=goods[0] || ksvind[bestKs][1]!=goods[1] ) return -1; //если хоть один трек процедурного каона не совпадает с нашим хорошим, то к чёрту его

  if( ksalign[bestKs] < 0.8) return -1; //отбор по косинусу между направлением импульса и направлением на пучок в р-фи плоскости

  if( fabs( ksptot[bestKs] - sqrt(emeas*emeas - mKs*mKs) ) > P_CUT ) return -1; //отбор по импульсу каона

  return bestKs;
}

int MC::Kinfit(Long64_t entry, std::vector<int> goods, double& mass_rec, double& chi){ //Кин.фит
  if(nph==0) return 0; //если нет фотонов, то выкинуться

  TLorentzVector PKS[2], KS, KL;
  PKS[0] = VectorCreator(tptot[goods[0]], tth[goods[0]], tphi[goods[0]], mPi);
  PKS[1] = VectorCreator(tptot[goods[1]], tth[goods[1]], tphi[goods[1]], mPi);

  KS = PKS[0] + PKS[1];
  KL = VectorCreator( KS.P(), TMath::Pi() - KS.Theta(), KS.Phi()+TMath::Pi(), mKs );

  //Проверить пространственный угол
  pair<double, double> ang = psi_angle(KS.P());
  if(PKS[0].Angle(PKS[1].Vect())<ang.first*0.95 || PKS[0].Angle(PKS[1].Vect())>ang.second*1.05) return 0;

  TLorentzVector Photon;
  double MIN_ANG = TMath::Infinity();
  int bestPh = -1;
  for(int i=0; i<nph; i++){
    Photon = VectorCreator(phen[i], phth[i], phphi[i], 0);
    if(Photon.Angle(KL.Vect()) < MIN_ANG){
      MIN_ANG = Photon.Angle(KL.Vect());
      bestPh = i;
    }
  }
  KL = VectorCreator( KS.P(), phth[bestPh], phphi[bestPh], mKs);

  vector<KFParticle> InParticles;
  vector<KFParticle> OutParticles;
  KFParticle Par[4];

  InParticles.clear();
  OutParticles.clear();

  Par[0].P = PKS[0];
  Par[0].Cov = GetTrErrorMatrix(PKS[0], terr[goods[0]]);
  InParticles.push_back(Par[0]);

  Par[1].P = PKS[1];
  Par[1].Cov = GetTrErrorMatrix(PKS[1], terr[goods[1]]);
  InParticles.push_back(Par[1]);

  float klerr[3][3];
  for(int i=0; i<3; i++)
  for(int j=0; j<3; j++)
  klerr[i][j] = 0;
  klerr[0][0] = pow( 20, 2);
  klerr[1][1] = pherr[bestPh][2];
  klerr[2][2] = pherr[bestPh][1];

  Par[2].P = KL;
  Par[2].Cov = GetTrErrorMatrix(KL, klerr);
  InParticles.push_back(Par[2]);

  chi = Cmd3KF(emeas, InParticles, OutParticles);
  //cout << KL.P() << '\t' << klerr[0][0] << '\t' << klerr[1][1] << '\t' << klerr[2][2] << '\t' << chi << endl;
  if(chi>100) return 0;
  mass_rec = ( OutParticles[0].P + OutParticles[1].P ).M();
  return 1;
}

void MC::Loop()
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  TFile* f = TFile::Open(path.c_str(), "recreate");
  double BEAM_ENERGY, LABEL;
  double MASS, MASS_REC, IMPULSE, ALIGN, THETA, LEN, CHI; //масса, реконструированная масса, импульс, угол, качество события, хи-квадрат
  int TRIGGER; //триггеры
  double DEDX[2]; // dE/dX
  int PROCEDURE; //процедура, которая отобрала событие (1 - kinfit, 2 - standard, 3 - both)


  //SimpleVars
  std::vector<int> goods;
  int bestKs;

  TTree *t = new TTree("t", "Tree for invariant mass without energy cut");
  t->Branch("label", &LABEL, "label/D");
  t->Branch("be", &BEAM_ENERGY, "be/D");
  t->Branch("m", &MASS, "m/D");
  t->Branch("m_rec", &MASS_REC, "m/D");
  t->Branch("p", &IMPULSE, "p/D");
  t->Branch("align", &ALIGN, "align/D");
  t->Branch("t", &TRIGGER, "t/I");
  t->Branch("proc", &PROCEDURE, "proc/I");
  t->Branch("dedx", &DEDX, "dedx[2]/D");
  t->Branch("theta", &THETA, "theta/D");
  t->Branch("len", &LEN, "len/D");
  t->Branch("chi", &CHI, "chi/D");

  bool model = (fChain->GetMaximum("nsim")<1) ? false : true;
  cout << "Is this model? " << (model ? "Yes" : "No") << endl;

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;


    //Общие условия
    goods = Good_tracks(ientry);
    if(goods.size() != 2 ) continue; //2 хороших трека

    if( tcharge[goods[0]] + tcharge[goods[1]] != 0 ) continue; //если суммарный заряд ненулевой, то выкинуться


    LABEL = ebeam;
    BEAM_ENERGY = emeas;
    DEDX[0] = tdedx[goods[0]];
    DEDX[1] = tdedx[goods[1]];
    TRIGGER = trigbits - 1;

        //инициализировать стандартными значениями переменные для дерева
    PROCEDURE = 0;
    MASS = -1;
    MASS_REC = -1;
    IMPULSE = -1;
    ALIGN = -1;
    THETA = -1;
    LEN = -1;
    CHI = -1;

    //Kinfit
    PROCEDURE += Kinfit(ientry, goods, MASS_REC, CHI);

    //Стандартная процедура
    bestKs = model ? ( (Cut(ientry) < 0) ? -1 : StandardProcedure(ientry, goods) ) : StandardProcedure(ientry, goods); //Специальный отбор на мягкие фотоны для моделирования
    if(bestKs>=0){
      PROCEDURE += 2;
      MASS = ksminv[bestKs];
      IMPULSE = ksptot[bestKs];
      ALIGN = ksalign[bestKs];
      THETA = ksth[bestKs];
      LEN = kslen[bestKs];
    }

    if(PROCEDURE>0)
      t->Fill();

  }
  f->Write();
  return;
}
