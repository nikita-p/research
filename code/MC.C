#define MC_cxx
#include "MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <TMath.h>
#include <TLorentzVector.h>
#include <iostream>
#include <boost/algorithm/string.hpp>

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

double MC::pidedx(double P, double dEdX) //calculate dEdX for pions
{
  double pidedx = 5.58030e+9 / pow(P + 40., 3) + 2.21228e+3 - 3.77103e-1 * P - dEdX;
  return pidedx;
}

std::vector<double> MC::Pcut(double Ebeam) //Mev
{
    double x = Ebeam*2e-3;
    double angle =  0.205/(x-0.732) + 0.14;
    double s1 = 100; // Mev old: 29.2*(x-0.5) + 0.955;
    double s2 = 10; //MeV
    return {angle, s2, 2*s1};
//   return 2 * (0.0869 * Ebeam - 36.53);
}

void MC::SetOutputPath(string key, string input_name)
{
  string filepath = (string)"../outputs/" + key + "/trees" + (SYS ? "_sys" : "") + "/";
  if (fChain == 0)
  {
    this->path = (filepath + "failed.root");
    return;
  }

  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t ent = fChain->GetEntries();
  Long64_t nbytes = 0, nb = 0;

  if(input_name==""){
  char label[100];
      for (Long64_t jentry = 0; jentry < 1; jentry++)
      {
        Long64_t ientry = LoadTree(jentry);
        if (ientry < 0)
          break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;

        snprintf(label, sizeof(label), "%.2f_%d.root", ebeam, runnum);
      }
  this->path = filepath + label;
  }
  else{
      std::vector<std::string> results;
      boost::algorithm::split(results, input_name, boost::is_any_of("/"));
      this->path = filepath + results.back();
  }
  cout << "Out path: " << this->path << endl;
  return;
}

void MC::GetLums(string file) //doesn't work
{
  if (fChain == 0)
    return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;
  
  bool model = (fChain->GetMaximum("nsim") < 1) ? false : true;
  if (!model)
    return;

  double LUM = 0;
  double EMEAS = -1;
  int RUN = -1;
  ofstream o(file.c_str());
  o << "label,lum\n";

  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;

    if (jentry == 0){
      EMEAS = emeas;
      RUN = runnum;
      LUM = 0;//getrunlum(RUN);
    }
      
    if(fabs(EMEAS - emeas) > 0.01){
      o << EMEAS << ',' << LUM << endl;
      EMEAS = emeas;
      RUN = runnum;
      LUM = 0;//getrunlum(RUN);
    }
      
    if(RUN!=runnum){
      RUN = runnum;
      LUM += 0;//getrunlum(RUN);
    }
  }
  o << EMEAS << ',' << LUM << endl;
  return;
}

int MC::StandardProcedure(Long64_t entry, std::vector<int> goods)
{ //вернуть отрицательное значение, если не проходит отборы
  MASS = -1;
  ANGLE_KS = -1;
  M1 = -1;
  M2 = -1;
  PASSED_A = false;
  PASSED_M = false;

  if (nks_total <= 0)
    return 0; //нет каонов - нет и отбора

  std::vector<double> P_CUT = Pcut(emeas); //для каждой энергии кат по импульсу(энергии) теперь будет различаться
  double minDiv = TMath::Infinity();

  int bestKs = -1; //пора искать лучший каон. Из всех, найденных процедурой, найдём лучший, с инв.массой наиболее близкой к массе каона.
  for (int i = 0; i < nks_total; i++)
  {
    if (fabs(ksminv[i] - mKs) > minDiv)
      continue;
    minDiv = fabs(ksminv[i] - mKs);
    bestKs = i;
  }

  if ((ksvind[bestKs][0] >= nt) || (ksvind[bestKs][1] >= nt)) //проверить, что в ksvind записана не дичь
    cout << "ksvind WARNING: nt=" << nt << "\tksvind0=" << ksvind[bestKs][0] << "\tksvind1=" << ksvind[bestKs][1] << endl;
  if (ksvind[bestKs][0] > ksvind[bestKs][1]) // проверить, что первый записанный трек может быть больше второго (добавил на всякий случай)
    cout << "WARNING THERE: " << ksvind[bestKs][0] << '\t' << ksvind[bestKs][1] << endl;
  if (ksvind[bestKs][0] != goods[0] || ksvind[bestKs][1] != goods[1])
    return -1; //если хоть один трек процедурного каона не совпадает с нашим хорошим, то к чёрту его

  MASS = ksminv[bestKs];
  ALIGN = ksalign[bestKs];
  MOMENTUM = ksptot[bestKs];
  THETA_KS = ksth[bestKs];
  PHI_KS = ksphi[bestKs];
  LEN_KS = kslen[bestKs];

  PASSED_A = (ksalign[bestKs] > (SYS ? 0.77 : 0.8) ) ? true : false;    //origin: 0.8                                      //отбор по косинусу между направлением импульса и направлением на пучок в р-фи плоскости
  if(emeas < mKs )
      PASSED_M = false;
  else{
      double p0 = sqrt(emeas * emeas - mKs * mKs);
      double th = P_CUT[0]; //корреляция
      M1 = (MASS - mKs)*cos(th) - (MOMENTUM - p0)*sin(th);
      M2 = (MASS - mKs)*sin(th) + (MOMENTUM - p0)*cos(th);
      bool cut_MOM1 = ( fabs( M1 ) < (SYS ? P_CUT[1]*1.05 : P_CUT[1]) );
//       bool cut_MOM2 = ( fabs( M2 ) < (SYS ? P_CUT[2]*1.05 : P_CUT[2]) );
      PASSED_M = cut_MOM1;//&&cut_MOM2; //отбор по импульсу каона
  }

  TLorentzVector KS = VectorCreator(ksptot[bestKs], ksth[bestKs], ksphi[bestKs], mKs);
  int i1 = ksvind[bestKs][0];
  int i2 = ksvind[bestKs][1];
  TLorentzVector Pi1 = VectorCreator(tptot[i1], tth[i1], tphi[i1], mPi);
  TLorentzVector Pi2 = VectorCreator(tptot[i2], tth[i2], tphi[i2], mPi);
  ANGLE_KS = (Pi1 + Pi2).Angle(KS.Vect());
  
  pic_align->Fill();
  pic_mom->Fill();
      
  if (PASSED_A && PASSED_M)
      return 1;
  return 0;
}

int MC::Kinfit(Long64_t entry, std::vector<int> goods)
{ //Кин.фит
  PASSED_CHI2 = false;
  PASSED_KL = false;
  PASSED_ANGLE = false;
  CHI2 = -1;
  KL_EN = -1;
  ANGLE_DIFF = -1;
  MASS_REC = -1;

  if (nph == 0)
    return 0; //если нет фотонов, то выкинуться

  TLorentzVector pion[2], ks, kl, cluster;
  double angle, max_angle = TMath::Infinity();
  int best_cluster = -1;

  pion[0] = VectorCreator(tptot[goods[0]], tth[goods[0]], tphi[goods[0]], mPi);
  pion[1] = VectorCreator(tptot[goods[1]], tth[goods[1]], tphi[goods[1]], mPi);
  ks = pion[0] + pion[1];

  for (int i = 0; i < nph; i++)
  {
    cluster = VectorCreator(phen0[i], phth0[i], phphi0[i], 0);
    angle = TMath::Pi() - cluster.Angle(ks.Vect());
    if (angle < max_angle)
    {
      best_cluster = i;
      max_angle = angle;
    }
  }

  if (best_cluster < 0)
    return 0;

  kl = VectorCreator(sqrt(emeas * emeas - mKs * mKs), phth0[best_cluster], phphi0[best_cluster], mKs);
  ANGLE_DIFF = TMath::Pi() - kl.Angle(ks.Vect());
  vector<KFParticle> InParticles;
  vector<KFParticle> OutParticles;
  KFParticle Par[4];

  InParticles.clear();
  OutParticles.clear();

  Par[0].P = pion[0];
  Par[0].Cov = GetTrErrorMatrix(pion[0], terr[goods[0]]);
  InParticles.push_back(Par[0]);

  Par[1].P = pion[1];
  Par[1].Cov = GetTrErrorMatrix(pion[1], terr[goods[1]]);
  InParticles.push_back(Par[1]);

  Par[2].P = kl;
  Par[2].Cov = GetPhErrorMatrix(kl, 1e5, pherr[best_cluster][1], pherr[best_cluster][2]);
  InParticles.push_back(Par[2]);

  CHI2 = Cmd3KF(emeas, InParticles, OutParticles);

  MOM_KS = ks.P();
  MOM_SUM = (pion[0].P() + pion[1].P()) / 2.;
  ks = (OutParticles[0].P + OutParticles[1].P); //кинфитированный KS
  kl = OutParticles[2].P;                       //кинфитированный KL
  pion[0] = OutParticles[0].P;
  pion[1] = OutParticles[1].P;
  KL_EN = phen0[best_cluster];
  MASS_REC = ks.M();

  PASSED_KL = (KL_EN > 100) ? true : false; // 1.1 * (emeas - 550) + 100
  PASSED_CHI2 = ((CHI2 > 0) && (CHI2 < 100)) ? true : false;
  PASSED_ANGLE = (ANGLE_DIFF < 0.5) ? true : false;
  PASSED_MOM = (fabs(MOM_KS - sqrt(emeas * emeas - mKs * mKs)) < 2 * (0.0869 * emeas - 36.53)) ? true : false;
  PASSED_MOM_SUM = (fabs(MOM_SUM - sqrt(pow(emeas / 2., 2.) - mPi * mPi)) < 50) ? true : false;

  pic_kinfit->Fill();

  if (PASSED_KL && PASSED_CHI2 && PASSED_ANGLE && PASSED_MOM && PASSED_MOM_SUM)
    return 1;

  return 0;
}

void MC::FillSimParticles(Long64_t entry, std::vector<double> *simparticles)
{
  simparticles->clear();
  for(int i=0; i<nsim; i++)
    if( (simorig[i]==0)&&(simtype[i]!=22) )
      simparticles->push_back( simtype[i] );
  return;
}

double MC::RadiativePhotonsEnergy(Long64_t entry)
{
    double ph_energy = 0;
    for(int i=0; i<nsim; i++)
        if((simtype[i]==22)&&(simorig[i]==0))
        {
            ph_energy += simmom[i];
        }
    return ph_energy;
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

  t = new TTree("t", "Tree for invariant mass without energy cut");
  t->Branch("label", &LABEL, "label/D");                   //метка дерева (номинальная энергия)
  t->Branch("beam_energy", &BEAM_ENERGY, "beam_energy/D"); //реальная энергия (по комптону)
  t->Branch("procedure", &PROCEDURE, "procedure/I");       //метка процедуры, через которую прошло событие (1-standard, 2-kinfit, 3-both)
  t->Branch("trigger", &TRIGGER, "t/I");                   //номер сработавшего триггера
  t->Branch("mass", &MASS, "mass/D");                      //масса из стандартной процедуры
  t->Branch("m1", &M1, "m1/D");                            //независимые оси
  t->Branch("m2", &M2, "m2/D");
  t->Branch("mass_reco", &MASS_REC, "mass_reco/D");        //масса из кинфита
  t->Branch("angle_ks", &ANGLE_KS, "angle_ks/D");          //пространственный угол между KS и суммарным импульсом двух пионов
  t->Branch("theta_ks", &THETA_KS, "theta_ks/D");          //полярный угол KS из стандартной процедуры
  t->Branch("phi_ks", &PHI_KS, "phi_ks/D");          //азимутальный угол KS из стандартной процедуры
  t->Branch("len_ks", &LEN_KS, "len_ks/D");          //длина отлёта KS
  

  //Способ извлекать картинки
  pic_align = new TTree("pic_align", "Tree as a picture of align selection"); //отбор по косинусу
  pic_align->Branch("align", &ALIGN, "align/D");
  pic_align->Branch("mass", &MASS, "mass/D");
  pic_align->Branch("passed", &PASSED_A, "passed/O");

  pic_mom = new TTree("pic_mom", "Tree as a picture of momentum selection"); //отбор по импульсу
  pic_mom->Branch("align", &ALIGN, "align/D");
  pic_mom->Branch("momentum", &MOMENTUM, "momentum/D");
  pic_mom->Branch("mass", &MASS, "mass/D");
  pic_mom->Branch("passed", &PASSED_M, "passed/O");
  pic_mom->Branch("m1", &M1, "m1/D");
  pic_mom->Branch("m2", &M2, "m2/D");

  pic_kinfit = new TTree("pic_kinfit", "Tree as a picture of kinfit selection"); //отборы в кинфите
  pic_kinfit->Branch("kl_en", &KL_EN, "kl_en/D");
  pic_kinfit->Branch("chi2", &CHI2, "chi2/D");
  pic_kinfit->Branch("angle_diff", &ANGLE_DIFF, "angle_diff/D"); //пространственный угол между KL и ближайшим кластером
  pic_kinfit->Branch("mass_reco", &MASS_REC, "mass_reco/D");
  pic_kinfit->Branch("mom_ks", &MOM_KS, "mom_ks/D");
  pic_kinfit->Branch("mom_sum", &MOM_SUM, "mom_sum/D");
  pic_kinfit->Branch("passed_kl", &PASSED_KL, "passed_kl/O");
  pic_kinfit->Branch("passed_chi2", &PASSED_CHI2, "passed_chi2/O");          //bool по chi2
  pic_kinfit->Branch("passed_angle", &PASSED_ANGLE, "passed_angle/O");       //bool по angle_diff
  pic_kinfit->Branch("passed_mom", &PASSED_MOM, "passed_mom/O");             //bool по mom_ks
  pic_kinfit->Branch("passed_mom_sum", &PASSED_MOM_SUM, "passed_mom_sum/O"); //bool по mom_sum

  bool model = (fChain->GetMaximum("nsim") < 1) ? false : true;
  cout << "Is this model? " << (model ? "Yes" : "No") << endl;
  
  std::vector<double> *st = new std::vector<double>(); 
  if(model){
    pic_mom->Branch("simtypes","vector<double>",&st);
    pic_mom->Branch("gamma_energy", &E_GAMMA, "gamma_energy/D");
    
    mc_passed = new TTree("mc_passed", "Passed photons"); //дерево для эффективности регистрации от энергии фотона
    mc_passed->Branch("ph_energy", &PH_ENERGY, "ph_energy/D");
    mc_passed->Branch("passed_cuts", &PASSED_CUTS, "passed_cuts/O");
  }

  double SYM_TYPES[10];
  if(SYS) 
    cout << "SYSTEMATICS MODE\n";
  for (Long64_t jentry = 0; jentry < nentries; jentry++)
  {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0)
      break;
    nb = fChain->GetEntry(jentry);
    nbytes += nb;
    if(!model && (badrun(runnum) == 1)) //check luminosity
        continue;
        
    
    if(model){
        FillSimParticles(ientry, st);
        E_GAMMA = this->RadiativePhotonsEnergy(ientry);
    }
        
    BEAM_ENERGY = emeas;
    LABEL = ebeam;
    //Общие условия
    goods = Good_tracks(ientry);
    if (goods.size() != 2){
        if(model){
            PASSED_CUTS = false;
            PH_ENERGY = this->RadiativePhotonsEnergy(ientry);
            mc_passed->Fill();
        }
        continue; //2 хороших трека
    }

    if (tcharge[goods[0]] + tcharge[goods[1]] != 0)
      continue; //если суммарный заряд ненулевой, то выкинуться

    TRIGGER = trigbits - 1;

    //инициализировать стандартными значениями переменные для дерева
    PROCEDURE = 0;

    //Kinfit
//     PROCEDURE += Kinfit(ientry, goods); //вернуть kinfit на место (пока я с ним не работаю, пусть отдыхает)

    //Стандартная процедура
    PROCEDURE += 2 * StandardProcedure(ientry, goods);
      
    if( model )
    {
        PASSED_CUTS = (PROCEDURE > 1) ? true : false;
        PH_ENERGY = this->RadiativePhotonsEnergy(ientry);
        mc_passed->Fill();
    }

    if (PROCEDURE > 0)
      t->Fill();
  }
  f->Write();
  delete st;
  return;
}
