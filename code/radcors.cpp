#include <TMath.h>
#include <TRandom.h>
#include <TComplex.h>

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include <ctime>
#include <algorithm>

using namespace std;

class RadCor
{
  double Alpha = 0.0072973525376;
  double Threshold;
  vector<double> e_cs;     //энергии, в которых измерены сечения
  vector<double> cs;       //сами сечения
  vector<double> e_out;    //энергии, которые нам нужны
  vector<double> rad_cors; //рад. поправки, которые считаются

  double Delta = 1e-5;
  int NSim = 1e5;
  double X1max = 0.99; //
  double X2max = 0.99; //

  void ReadCrossSec(string CSFile)
  {
    ifstream o(CSFile.c_str());
    int index;
    double ex, csx;
    while (!o.eof())
    {
      o >> ex >> csx; //>> index >> ex >> csx;
      if (o.eof())
        break;
      e_cs.push_back(ex * 1e-3); //перевод энергии в ГэВ
      cs.push_back(csx);         //сечения в нб
    }
    cout << "Successful reading... Cross section File contains " << cs.size() << " records" << endl;
    // for(auto it = cs.begin(); it!=cs.end(); it++)
      // cout << *it << endl;
    return;
  }

  void ReadEnergies(string energies)
  {
    ifstream o(energies.c_str());
    double e;
    while (!o.eof())
    {
      o >> e;
      if (o.eof())
        break;
      e_out.push_back(e * 2e-3); //перевод энергии в ГэВ
    }
    cout << "Successful reading... Energy File contains " << e_out.size() << " records" << endl;
    return;
  }

  double CrossSection(double s)
  {
    double E = sqrt(s);

    int i = int((e_cs.size() + 1) / 2.);
    int step = int(i / 2.);

    i = (E < e_cs[i]) ? (i - step) : (i + step);

    while (step != 1)
    {
      step = int((step + 1) / 2.);
      i = (E < e_cs[i]) ? (i - step) : (i + step);
    }
    i = (E < e_cs[i]) ? i : (i + 1);

    if (i >= e_cs.size())
      i = e_cs.size() - 1;
    if (i < 1)
      i = 1;
    double cross = cs[i - 1] + (cs[i] - cs[i - 1]) * (E - e_cs[i - 1]) / (e_cs[i] - e_cs[i - 1]);
    cross = (cross > 0) ? cross : 0;
    return cross;
  }

  double L_Func(double s)
  {
    return TMath::Log(s / 0.0005 / 0.0005);
  }

  double B_Func(double s)
  {
    return 2.0 * Alpha * (L_Func(s) - 1) / TMath::Pi();
  }

  double D0_Func(double s)
  {
    double Beta = B_Func(s);
    double L = L_Func(s);
    return 1.0 + 3.0 * Beta / 8.0 - TMath::Power(Beta, 2.0) * (L / 3.0 + TMath::Power(TMath::Pi(), 2.0) - 47.0 / 8.0) / 48.0;
  }

  double D_Func(double Z, double s)
  {
    double Factor;
    double Part1, Part2;
    double Beta = B_Func(s);
    double D0 = D0_Func(s);

    Factor = 2.0 * TMath::Power((1 - Z), (1.0 - 0.5 * Beta)) / Beta;

    Part1 = Factor * Beta * (1 + Z) / 4.0;
    Part2 = Factor * (TMath::Power(Beta, 2.0) * (4.0 * (1 + Z) * TMath::Log(1.0 / (1.0 - Z)) + (1 + 3.0 * TMath::Power(Z, 2.0)) * TMath::Log(1.0 / Z) / (1.0 - Z) - 5.0 - Z) / 32.0);

    double fval = D0 - Part1 + Part2;

    return fval;
  }

  void ComputeOne(double s) //GeV^2
  {
    double emeas = (sqrt(s)/2.)*1e3;
    double pb = sqrt( s/4. - pow(0.4976, 2) );
    double dp = (2 * (0.0869 * emeas - 36.53))*1e-3;
    X1max = 2*( 1 - sqrt(1-2*pb*dp/s) ); //40e-3/sqrt(s); 
    cout << "\nX1max: " << (X1max*sqrt(s))*1e3 << " MeV\n";
    X2max = 2*( 1 - sqrt(1-2*pb*dp/s) ); //40e-3/sqrt(s);
    //cout << emeas << '\t' << X1max << '\n';

    vector<int> Ns(4, 0);                      //вместо Nsim1, Nsim2, Nsim3, Nsim4;
    vector<double> Max(4, -TMath::Infinity()); //вместо Max1, Max2, Max3, Max4
    vector<double> Part(4);
    vector<double> Sigma(4, 0);

    double X1, Z1, X2, Z2, X2m;

    double L = L_Func(s);
    double Beta = B_Func(s);
    double D0 = D0_Func(s);

    double K_Factor = 1.0 + Alpha * (TMath::Power(TMath::Pi(), 2.0) / 3.0 - 0.5) / TMath::Pi();

    double Temp1 = TMath::Log(Delta) * Beta / 2.0;
    double Temp2 = TMath::Exp(Temp1);

    //1 area (это константа, поэтому вне цикла)
    Part[0] = TMath::Power(Delta, Beta) * TMath::Power(D0, 2.0) * CrossSection(s * (1.0 - Delta)) * K_Factor;
    Max[0] = Part[0];

    for (int i = 0; i < NSim; i++)
    {
      //2 area
      X1 = TMath::Power(gRandom->Rndm() * (TMath::Power(X1max, Beta / 2.0) - Temp2) + Temp2, 2.0 / Beta);
      Z1 = 1.0 - X1;

      if (s * Z1 < Threshold * Threshold)
        continue;

      Part[1] = Temp2 * D0 * (TMath::Power(X1max, Beta / 2.0) - Temp2) * D_Func(Z1, s) * CrossSection(s * Z1) * K_Factor;
      Max[1] = (Part[1] > Max[1]) ? Part[1] : Max[1];

      //3 area
      X2 = TMath::Power(gRandom->Rndm() * (TMath::Power(X2max, Beta / 2.0) - Temp2) + Temp2, 2.0 / Beta);
      Z2 = 1.0 - X2;

      if (s * Z2 < Threshold * Threshold)
        continue;

      Part[2] = Temp2 * D0 * (TMath::Power(X2max, Beta / 2.0) - Temp2) * D_Func(Z2, s) * CrossSection(s * Z2) * K_Factor;
      Max[2] = (Part[2] > Max[2]) ? Part[2] : Max[2];

      //4 area
      X1 = TMath::Power(gRandom->Rndm() * (TMath::Power(X1max, Beta / 2.0) - Temp2) + Temp2, 2.0 / Beta);
      Z1 = 1.0 - X1;

      X2m = min(X2max, 1.0 - Threshold * Threshold / s / (1 - X1));
      X2 = TMath::Power(gRandom->Rndm() * (TMath::Power(X2m, Beta / 2.0) - Temp2) + Temp2, 2.0 / Beta);
      Z2 = 1.0 - X2;

      if (s * Z1 * Z2 < Threshold * Threshold)
        continue;

      Part[3] = (TMath::Power(X1max, Beta / 2.0) - Temp2) * (TMath::Power(X2m, Beta / 2.0) - Temp2) * D_Func(Z1, s) * D_Func(Z2, s) * CrossSection(s * Z1 * Z2) * K_Factor;
      Max[3] = (Part[3] > Max[3]) ? Part[3] : Max[3];
    }

    //1 area
    Ns[0] = NSim;
    //2 area
    for (int i = 0; i < NSim; i++)
    {
      X1 = TMath::Power(gRandom->Rndm() * (TMath::Power(X1max, Beta / 2.0) - Temp2) + Temp2, 2.0 / Beta);
      Z1 = 1.0 - X1;
      if (s * Z1 < Threshold * Threshold)
        continue;
      Part[1] = Temp2 * D0 * (TMath::Power(X1max, Beta / 2.0) - Temp2) * D_Func(Z1, s) * CrossSection(s * Z1);
      if (gRandom->Rndm() * Max[1] < Part[1])
        Ns[1] += 1;
    }
    //3 area
    for (int i = 0; i < NSim; i++)
    {
      X2 = TMath::Power(gRandom->Rndm() * (TMath::Power(X2max, Beta / 2.0) - Temp2) + Temp2, 2.0 / Beta);
      Z2 = 1.0 - X2;
      if (s * Z2 < Threshold * Threshold)
        continue;
      Part[2] = Temp2 * D0 * (TMath::Power(X2max, Beta / 2.0) - Temp2) * D_Func(Z2, s) * CrossSection(s * Z2);
      if (gRandom->Rndm() * Max[2] < Part[2])
        Ns[2] += 1;
    }
    //4 area
    for (int i = 0; i < NSim; i++)
    {
      X1 = TMath::Power(gRandom->Rndm() * (TMath::Power(X1max, Beta / 2.0) - Temp2) + Temp2, 2.0 / Beta);
      Z1 = 1.0 - X1;
      X2m = min(X2max, 1.0 - Threshold * Threshold / s / (1 - X1));
      X2 = TMath::Power(gRandom->Rndm() * (TMath::Power(X2m, Beta / 2.0) - Temp2) + Temp2, 2.0 / Beta);
      Z2 = 1.0 - X2;
      if (s * Z1 * Z2 < Threshold * Threshold)
        continue;
      Part[3] = (TMath::Power(X1max, Beta / 2.0) - Temp2) * (TMath::Power(X2m, Beta / 2.0) - Temp2) * D_Func(Z1, s) * D_Func(Z2, s) * CrossSection(s * Z1 * Z2);
      if (gRandom->Rndm() * Max[3] < Part[3])
        Ns[3] += 1;
    }

    for (int i = 0; i < 4; i++)
    {
      Sigma[i] = Max[i] * Ns[i] / NSim;
      //cout << "Max" << i << ": " <<  Max[i] << endl;
    }

    double Sigma0 = Sigma[0] + Sigma[1] + Sigma[2] + Sigma[3];

    double RadiativeCorrection;

    if (CrossSection(s) != 0.0)
      RadiativeCorrection = Sigma0 / CrossSection(s);
    else
      RadiativeCorrection = 1.0;

    rad_cors.push_back(RadiativeCorrection);
    cout << "Radiative correction at E = " << sqrt(s) << " GeV is equal to " << RadiativeCorrection << endl;
    return;
  }

  void Compute()
  {
    double e;
    time_t seconds;
    for (auto it = e_out.begin(); it != e_out.end(); it++)
    {
      e = *it;
      seconds = time(NULL);
      ComputeOne(e * e);
      seconds = time(NULL) - seconds;
      cout << "\tTime of this calculation: " << seconds << " seconds\n";
    }
    return;
  }

  void Write(string filename)
  {
    cout << "Write rad corrections to file.. ";
    ofstream o(filename.c_str());
    for (int i = 0; i < e_out.size(); i++)
    {
      o << e_out[i] << "," << rad_cors[i] << endl;
    }
    cout << "Well done\n";
    return;
  }

public:
  RadCor(string energies, string CSFile, double EThreshold)
  {
    Threshold = EThreshold; //GeV
    cout << "EnergyFile name is " << energies << endl;
    cout << "Reaction threshold is " << Threshold << ", GeV" << endl;
    cout << "Cross Section File Name is " << CSFile << endl;
    ReadCrossSec(CSFile);
    ReadEnergies(energies);
  }

  void PrintEnergies()
  {
    cout << "Out energies\n";
    for (auto it = e_out.begin(); it != e_out.end(); it++)
    {
      cout << *it << endl;
    }
    return;
  }

  void Calc(string outputFile)
  {
    cout << "\nRadiative corrections calculation started..\n";
    Compute();
    Write(outputFile);
    return;
  }

  void SetNSim(int NSim)
  {
    this->NSim = NSim;
    return;
  }
  void SetDelta(double Delta)
  {
    this->Delta = Delta;
    return;
  }
  void SetX1max(double X1max)
  {
    this->X1max = X1max;
    return;
  }
  void SetX2max(double X2max)
  {
    this->X2max = X2max;
    return;
  }
};

int radcors()
{
  RadCor rc("../inputs/radcors/energies.dat", "../inputs/radcors/cs.dat", 0.9952);
  rc.SetNSim(1e5);
  rc.Calc("../outputs/radcors.dat");
  return 0;
}