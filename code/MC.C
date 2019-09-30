#define MC_cxx
#include "MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

#define mKs 497.614
#define mPi 139.570

double pidedx(double P, double dEdX)    //calculate dEdX for pions
{
  double pidedx = 5.58030e+9/pow(P+40.,3)+2.21228e+3-3.77103e-1*P-dEdX;
  return pidedx;
}

void MC::Loop(string file = "train.root")
{
  if (fChain == 0) return;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nbytes = 0, nb = 0;

  TFile* f = TFile::Open(file.c_str(), "recreate");
  double BEAM_ENERGY; //энергия пучка, измеренная лазером, и ошибка
  double MASS, ENERGY, IMPULSE, ALIGN, QUALITY, THETA, LEN; //масса, энергия, импульс, угол, качество события
  int TRIGGER; //триггеры
  double RADIUS[2]; //отлёт от пучков
  double DEDX[2]; // dE/dX
  bool WIN; //отобранные процедурой события
  bool W1[3]; // отдельные отборы
  double P_CUT; //для каждой энергии кат по импульсу(энергии) теперь будет различаться
  double TH_CUT = 0.6; //отбор по углу тета


  //SimpleVars
  pair<int,int> index(0,0);
  int good;

  TTree *t = new TTree("InvMass", "Tree for invariant mass without energy cut");
  t->Branch("be", &emeas, "be/D");
  t->Branch("m", &MASS, "m/D");
  t->Branch("e", &ENERGY, "e/D");
  t->Branch("p", &IMPULSE, "p/D");
  t->Branch("align", &ALIGN, "align/D");
  t->Branch("r", &RADIUS, "r[2]/D");
  t->Branch("win", &WIN, "win/O");
  t->Branch("w1", &W1, "w1[3]/O");
  t->Branch("t", &TRIGGER, "t/I");
  t->Branch("quality", &QUALITY, "quality/D");
  t->Branch("dedx", &DEDX, "dedx[2]/D");
  t->Branch("theta", &THETA, "theta/D");
  t->Branch("len", &LEN, "len/D");

  bool model = (fChain->GetMaximum("nsim")<1) ? false : true;
  cout << "Is this model? " << (model ? "Yes" : "No") << endl;
  //для моделирования отдельный отбор на мягкие фотоны
  int SOFT_PHOTONS = 0; //количество событий с мягкими фотонами (важно для моделирования, можно на него нормировать)
  double SOFT_PHOTONS_MOMENTUM;
  if(model){
    t->Branch("soft_ph", &SOFT_PHOTONS, "soft_ph/I");
    for(Long64_t jentry=0; jentry<nentries; jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      P_CUT = 2*(0.0869*emeas-36.53);
      SOFT_PHOTONS_MOMENTUM = 0;
      int j = 0;
      for(int i=0; i<nsim; i++){
        if(simtype[i]==310){
          SOFT_PHOTONS_MOMENTUM += simmom[i];
          j++;
        }
        if( j==2 )
        std::cout << "Warning\n";
      }
      if( TMath::Abs(SOFT_PHOTONS_MOMENTUM-sqrt(emeas*emeas - mKs*mKs))>P_CUT ) continue; //если импульс KS в событии отличается от импульсе, при энергии KS равной энергии пучка больше чем на P_CUT, то не работать с ним
      SOFT_PHOTONS += 1;
    }
    std::cout << SOFT_PHOTONS << std::endl;
  }

  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    // if (Cut(ientry) < 0) continue;
    //Init vars
    index = {-1, -1}; //индексы хороших треков (отрицательные индексы говорят о том, что треков нет)
    good = 0; //число хороших треков (важно не забыть, что в начале обработки каждого события здесь должен быть 0)
    WIN = false; //попавшее к нам событие пометится как true
    W1[0] = false; W1[1] = false; W1[2] = false; //флаг, обозначающий прохождение некоторого конкретного отбора
    QUALITY = 0; //субъективно-объективный параметр, отражающий качество данного события

    //Специальный отбор на мягкие фотоны для моделирования
    if(model){
      SOFT_PHOTONS_MOMENTUM = 0;
      for(int i=0; i<nsim; i++){
        if(simtype[i]==310)
        SOFT_PHOTONS_MOMENTUM += simmom[i];
      }
      if( TMath::Abs(SOFT_PHOTONS_MOMENTUM-sqrt(emeas*emeas - mKs*mKs))>P_CUT) continue; //если энергия KS в событии отличается от пучка больше чем на 20, то не работать с ним
    }

    //Conditions
    if(nt<2) continue; //нет двух треков, нет и дел с таким событием
    if(nks_total<=0) continue; //нет каонов - нет и отбора

    for(int i=0; i<nt; i++){ //пробегаем по всем трекам из события, яхууу

      if( fabs(tz[i])>10.0 ) continue; //вылетел из пучка
      if( tchi2r[i]>30.0 ) continue; // хи2 хороший
      if( tchi2z[i]>25.0 ) continue;
      if((tth[i]>(TMath::Pi() - TH_CUT))||(tth[i] < TH_CUT)) continue; //летит в детектор

      if( tptotv[i]<40. ) continue; //меньшие импульсы непригодны, т.к. треки закрутятся в дк
      if( tptotv[i]>2*ebeam ) continue; //куда ж ещё больше
      if( tnhit[i]<6 ) continue; //5 уравнений - 5 неизвестных: phi, theta, P, ...  -->>-- я добавил по сравн. с пред. версией 1 хит (стало 6)
      if( fabs(pidedx(tptotv[i], tdedx[i]))>2000 ) continue; //ионизационные потери

      good++;
      if(good>2) break;
      if(good==1) index.first = i;
      if(good==2) index.second = i; //запись в пару, если всё забито, то брэйкнуться отсюда

    }

    if(good!=2) continue; //как результат "for'а", если треков меньше (больше) чем достаточно, делаем кик аут в поисках лучшей жизни
    if( tcharge[index.first] + tcharge[index.second] != 0 ) continue; //если суммарный заряд ненулевой, то выкинуться

    P_CUT = 2*(0.0869*emeas-36.53);
    double alignMin = 0.8;
    double rhoCut = 0.1;
    double minDiv = TMath::Infinity();

    int bestKs = 0; //пора искать лучший каон. Из всех, найденных процедурой, найдём лучший, с инв.массой наиболее близкой к массе каона.
    for(int i=0; i<nks_total; i++){
      if(fabs(ksminv[i]-mKs)>minDiv) continue;
      minDiv = fabs(ksminv[i]-mKs);
      bestKs = i;
    }

    if(ksvind[bestKs][0]<ksvind[bestKs][1]){ //здесь изменение, раньше не было варианта, что первый записанный трек может быть больше второго (добавил на всякий случай)
      if( ksvind[bestKs][0]!=index.first || ksvind[bestKs][1]!=index.second ) continue; //если хоть один трек процедурного каона не совпадает с нашим хорошим, то к чёрту его
    }
    else{
      if( ksvind[bestKs][0]!=index.second || ksvind[bestKs][1]!=index.first ) continue;
    }

    if( ksalign[bestKs] > alignMin) W1[0] = true; //отбор по косинусу между направлением импульса и направлением на пучок в р-фи плоскости

    TLorentzVector K;
    K.SetXYZM(0,0,ksptot[bestKs],mKs); //после фита с общей вершиной
    K.SetTheta(ksth[bestKs]);
    K.SetPhi(ksphi[bestKs]);
    if( fabs( K.P() - sqrt(emeas*emeas - mKs*mKs) ) < P_CUT ) W1[1] = true; //отбор по энергии каона

    if( fabs(trho[index.first])>rhoCut && fabs(trho[index.second])>rhoCut ) W1[2] = true; //отбор по прицельному параметру в р-фи

    WIN = W1[0]&&W1[1]&&W1[2];

    BEAM_ENERGY = emeas;
    MASS = ksminv[bestKs];
    ENERGY = K.E();
    IMPULSE = K.P();
    ALIGN = ksalign[bestKs];
    RADIUS[0] = trho[index.first];
    RADIUS[1] = trho[index.second];
    TRIGGER = trigbits - 1;
    DEDX[0] = tdedx[index.first];
    DEDX[1] = tdedx[index.second];
    THETA = ksth[bestKs];
    LEN = kslen[bestKs];
    QUALITY = ( (fabs(sqrt(emeas*emeas-mKs*mKs) - IMPULSE)/P_CUT) + (fabs(1-ALIGN)/(1-alignMin)) )/2.; // лежит в [0; 1]
  }
  t->Write();
  return;
}
