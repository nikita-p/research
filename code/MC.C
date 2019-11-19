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

TLorentzVector MC::VectorCreator(double P, double Theta, double Phi, double Mass)
{
  TLorentzVector V;
  V.SetXYZM(0, 0, P, Mass);
  V.SetTheta(Theta);
  V.SetPhi(Phi);
  return V;
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

double MC::pidedx(double P, double dEdX) //calculate dEdX for pions
{
  double pidedx = 5.58030e+9 / pow(P + 40., 3) + 2.21228e+3 - 3.77103e-1 * P - dEdX;
  return pidedx;
}

void MC::SetOutputPath(string key)
{
  string filepath = "../outputs/" + key + "/";
  if (fChain == 0)
  {
    this->path = (filepath + "failed.root");
    return;
  }

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  char label[100];
  for (Long64_t jentry = 0; jentry < 1; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    snprintf(label, sizeof(label), "%.2f.root", ebeam);
  }
  this->path = filepath + label;
  cout << "Out path: " << this->path << endl;
  return;
}

void MC::GetSoftPhotonsNumber(string file){
  if (fChain == 0)
    return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  bool model = (fChain->GetMaximum("nsim") < 1) ? false : true;
  if (!model)
    return;

  int N_SOFT = 0;
  double EMEAS = -1;
  ofstream o(file.c_str());

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    if (jentry == 0)
      EMEAS = emeas;

    if (Cut(ientry) < 0)
      continue;
    N_SOFT++;

    if (fabs(EMEAS - emeas) > 0.01)
    {
      o << EMEAS << ',' << N_SOFT << endl;
      EMEAS = emeas;
      N_SOFT = 0;
    }
  }
  o << EMEAS << ',' << N_SOFT << endl;
  return;
}

int MC::StandardProcedure(Long64_t entry, std::vector<int> goods, bool& passed_align, bool& passed_mom)
{ //вернуть отрицательное значение, если не проходит отборы
  passed_align = false;
  passed_mom = false;

  if (nks_total <= 0)
    return -1; //нет каонов - нет и отбора

  double P_CUT = 2 * (0.0869 * emeas - 36.53); //для каждой энергии кат по импульсу(энергии) теперь будет различаться
  double minDiv = TMath::Infinity();

  int bestKs = -1; //пора искать лучший каон. Из всех, найденных процедурой, найдём лучший, с инв.массой наиболее близкой к массе каона.
  for (int i = 0; i < nks_total; i++)
  {
    if (fabs(ksminv[i] - mKs) > minDiv)
      continue;
    minDiv = fabs(ksminv[i] - mKs);
    bestKs = i;
  }

  if (ksvind[bestKs][0] > ksvind[bestKs][1]) // проверить, что первый записанный трек может быть больше второго (добавил на всякий случай)
    cout << "WARNING THERE: " << ksvind[bestKs][0] << '\t' << ksvind[bestKs][1] << endl;
  if (ksvind[bestKs][0] != goods[0] || ksvind[bestKs][1] != goods[1])
    return -1; //если хоть один трек процедурного каона не совпадает с нашим хорошим, то к чёрту его

  if (ksalign[bestKs] < 0.8)
    return -1; //отбор по косинусу между направлением импульса и направлением на пучок в р-фи плоскости
  passed_align = true;

  if (fabs(ksptot[bestKs] - sqrt(emeas * emeas - mKs * mKs)) > P_CUT)
    return -1; //отбор по импульсу каона
  passed_mom = true;

  return bestKs;
}

int MC::Kinfit(Long64_t entry, std::vector<int> goods, double &mass_rec, double &chi2, double& en_ph, bool& pass_chi2, bool& pass_en_ph)
{ //Кин.фит
  pass_chi2 = false;  chi2 = -1;
  pass_en_ph = false; en_ph = -1;

  if (nph == 0)
    return 0; //если нет фотонов, то выкинуться

  TLorentzVector PKS[2], KS, KL;
  PKS[0] = VectorCreator(tptot[goods[0]], tth[goods[0]], tphi[goods[0]], mPi);
  PKS[1] = VectorCreator(tptot[goods[1]], tth[goods[1]], tphi[goods[1]], mPi);

  KS = PKS[0] + PKS[1];
  KL = VectorCreator(KS.P(), TMath::Pi() - KS.Theta(), KS.Phi() + TMath::Pi(), mKs);

  TLorentzVector Photon;
  double MIN_ANG = 0.1;//TMath::Infinity();
  int bestPh = -1;
  double cluster_energy;
  double max_cluster_energy = -1;

  for (int i = 0; i < nph; i++)
  {
    Photon = VectorCreator(phen0[i], phth0[i], phphi0[i], 0);
    cluster_energy = phen0[i];
    if(cluster_energy > max_cluster_energy)
      max_cluster_energy = cluster_energy;
    if ((Photon.Angle(KL.Vect()) < MIN_ANG)&&(cluster_energy>1.1*(emeas-550)+100 ) )
    {
      MIN_ANG = Photon.Angle(KL.Vect());
      bestPh = i;
      en_ph = cluster_energy;
    }
  }

  if(bestPh<0){
    en_ph = max_cluster_energy;
    return 0;
  }
  
  pass_en_ph = true;

  KL = VectorCreator(KS.P(), phth0[bestPh], phphi0[bestPh], mKs);

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
/*
  float klerr[3][3];
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      klerr[i][j] = 0;
  klerr[0][0] = pow(35, 2);
  klerr[1][1] = pherr[bestPh][2];
  klerr[2][2] = pherr[bestPh][1];
*/
  Par[2].P = KL;
  //Par[2].Cov = GetTrErrorMatrix(KL, klerr);
  Par[2].Cov = GetPhErrorMatrix(KL, 2*KL.P(), pherr[bestPh][1], pherr[bestPh][2]);
  InParticles.push_back(Par[2]);

  chi2 = Cmd3KF(emeas, InParticles, OutParticles);

  KS = (OutParticles[0].P + OutParticles[1].P); //кинфитированный KS
  KL = OutParticles[2].P; //кинфитированный KL
  PKS[0] = OutParticles[0].P;
  PKS[1] = OutParticles[1].P;
  mass_rec = KS.M();
  
  if (chi2 > 10 || chi2 < 0)
    return 0;
  pass_chi2 = true;

  //Проверить пространственный угол!! не факт, что нужно его проверять, поскольку у пионов как раз определены параметры неправильно 
  /*  
  pair<double, double> ang = psi_angle(KS.P());
  if (PKS[0].Angle(PKS[1].Vect()) < ang.first * 0.8 || PKS[0].Angle(PKS[1].Vect()) > ang.second * 1.2 ) //слабый кат
    return 0;
  */ 
  return 1;
}

void MC::Loop()
{
  if (fChain == 0)
    return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  TFile *f = TFile::Open(path.c_str(), "recreate");
  
  //SimpleVars
  std::vector<int> goods;
  int bestKs;

  TTree *t = new TTree("t", "Tree for invariant mass without energy cut");
  double LABEL;       t->Branch("label", &LABEL, "label/D"); //метка дерева (номинальная энергия)
  double BEAM_ENERGY; t->Branch("beam_energy", &BEAM_ENERGY, "beam_energy/D"); //реальная энергия (по комптону)
  int PROCEDURE;      t->Branch("procedure", &PROCEDURE, "procedure/I"); //метка процедуры, через которую прошло событие (1-standard, 2-kinfit, 3-both)
  int EXIT_CODE;      t->Branch("exit_code", &EXIT_CODE, "exit_code/I"); //метка, указывающая на непройденный отбор

  int TRIGGER;        t->Branch("trigger", &TRIGGER, "t/I"); //номер сработавшего триггера
  double MASS;        t->Branch("mass", &MASS, "mass/D"); //масса из стандартной процедуры
  double MASS_REC;    t->Branch("mass_reco", &MASS_REC, "mass_reco/D"); //масса из кинфита
  double ANGLE_KS;    t->Branch("angle_ks", &ANGLE_KS, "angle_ks/D"); //пространственный угол между KS и суммарным импульсом двух пионов

  //Организовать способ извлекать картинки
  TTree *pic_align = new TTree("pic_align", "Tree as a picture of align selection"); //отбор по косинусу
  double ALIGN;   pic_align->Branch("align", &ALIGN, "align/D");
                  pic_align->Branch("mass", &MASS, "mass/D");
  bool PASSED_A;  pic_align->Branch("passed", &PASSED_A, "passed/O");

  TTree *pic_mom = new TTree("pic_mom", "Tree as a picture of momentum selection"); //отбор по импульсу
  double MOMENTUM;   pic_mom->Branch("momentum", &MOMENTUM, "momentum/D");
                     pic_mom->Branch("mass", &MASS, "mass/D");
  bool PASSED_M;     pic_mom->Branch("passed", &PASSED_M, "passed/O");

  TTree *pic_chi = new TTree("pic_chi", "Tree as a picture of kinfit chi-square selection"); //отбор по хи-2 в кинфите
  double CHI2;   pic_chi->Branch("chi2", &CHI2, "chi2/D");
                 pic_chi->Branch("mass_reco", &MASS_REC, "mass_reco/D");
  bool PASSED_C; pic_chi->Branch("passed", &PASSED_C, "passed/O");

  TTree *pic_kl = new TTree("pic_kl", "Tree as a picture of KL cluster energy selection (kinfit)"); //отбор по энергии кластера KL в кинфите
  double KL_EN;   pic_kl->Branch("kl_en", &KL_EN, "kl_en/D");
                  pic_kl->Branch("mass_reco", &MASS_REC, "mass_reco/D");
  bool PASSED_KL; pic_kl->Branch("passed", &PASSED_KL, "passed/O");

  bool model = (fChain->GetMaximum("nsim") < 1) ? false : true;
  cout << "Is this model? " << (model ? "Yes" : "No") << endl;

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    //Общие условия
    goods = Good_tracks(ientry);
    if (goods.size() != 2)
      continue; //2 хороших трека


    if (tcharge[goods[0]] + tcharge[goods[1]] != 0)
      continue; //если суммарный заряд ненулевой, то выкинуться

    LABEL = ebeam;
    BEAM_ENERGY = emeas;
    TRIGGER = trigbits - 1;

    //инициализировать стандартными значениями переменные для дерева
    PROCEDURE = 0;
    ANGLE_KS = -100;
    MASS = -1;
    MASS_REC = -1;
    double pks = -1;
    double mks = -1;

    //Kinfit
    PROCEDURE += Kinfit(ientry, goods, MASS_REC, CHI2, KL_EN, PASSED_C, PASSED_KL);
    if(CHI2>0)    pic_chi->Fill();
    if(KL_EN>0)   pic_kl->Fill();

    //Стандартная процедура
    bestKs = model ? ((Cut(ientry) < 0) ? -1 : StandardProcedure(ientry, goods, PASSED_A, PASSED_M)) : StandardProcedure(ientry, goods, PASSED_A, PASSED_M); //Специальный отбор на мягкие фотоны для моделирования
    if (bestKs >= 0)
    {
      PROCEDURE += 2;
      MASS = ksminv[bestKs];
      MOMENTUM = ksptot[bestKs];
      ALIGN = ksalign[bestKs];
      TLorentzVector KS = VectorCreator(ksptot[bestKs], ksth[bestKs], ksphi[bestKs], mKs);
      int i1 = ksvind[bestKs][0];
      int i2 = ksvind[bestKs][1];
      TLorentzVector Pi1 = VectorCreator(tptot[i1], tth[i1], tphi[i1], mPi);
      TLorentzVector Pi2 = VectorCreator(tptot[i2], tth[i2], tphi[i2], mPi);
      ANGLE_KS = (Pi1+Pi2).Angle(KS.Vect());
      pic_mom->Fill();
      pic_align->Fill();
    }

    if (PROCEDURE > 0)
      t->Fill();
  }
  f->Write();
  return;
}
