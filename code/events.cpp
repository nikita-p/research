#include "MC.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <iostream>
#include <string>
#include <ctime>
#include <fstream>
#include <TTree.h>
#include <TFile.h>
#include <TMath.h>

#define mKs 497.614
#define mPi 139.570

using namespace std;

class TreeReader
{
public:
  TreeReader()
  {
  }

  void Washing(std::vector<string> inputPath, string key, bool sys=false, bool time = false)
  {
    double time0 = clock() / (CLOCKS_PER_SEC + 0.);
    cout << "The door is closing" << endl;
    MC cl(inputPath, key, sys);
    cout << "Go loop" << endl;
    cl.Loop();
    cout << "Wow, I feel fresh" << endl;
    if (time)
      cout << "Washing time: " << clock() / (CLOCKS_PER_SEC + 0.) - time0 << "seconds." << endl;
    return;
  }

  void WashingFromFile(string filePath, string key, bool sys=false)
  {
    ifstream f(filePath.c_str());
    string input;
    while (f >> input)
    {
      if(input[0]=='#'){
        cout << "Passed " << input << endl;
        continue;
      }
      cout << "Working at file: " << input << endl;
      Washing({input}, key, sys, true);
      cout << "Well done" << endl
           << endl;
    }
    return;
  }

  void WashingModelOldv6(bool sys = false)
  {
    WashingFromFile("../inputs/model_old_v6/trees", "model_old_v6", sys);
    return;
  }
  
  void WashingModelOldv7(bool sys = false)
  {
    WashingFromFile("../inputs/model_old_v7/trees", "model_old_v7", sys);
    return;
  }
  
  void WashingModelNewv6(bool sys = false)
  {
    WashingFromFile("../inputs/model_new_v6/trees", "model_new_v6", sys);
    return;
  }

  void Washing11(bool sys = false)
  {
    WashingFromFile("../inputs/11/trees", "11", sys);
  }
  void Washing12(bool sys = false)
  {
    WashingFromFile("../inputs/12/trees", "12", sys);
  }
  void Washing17(bool sys = false)
  {
    WashingFromFile("../inputs/17/trees", "17", sys);
  }
  void Washing19(bool sys = false)
  {
    WashingFromFile("../inputs/19/trees", "19", sys);
  }
  void WashingOth(std::vector<string> other_files, string folder = "others", bool sys = false)
  {
    for(auto it=other_files.begin(); it!=other_files.end(); it++){
        cout << "File\n";
        Washing({*it}, folder.c_str(), sys, true);
    }
    return;
  }
};

void events_cores(string file, string key)
{
  TreeReader t;
  t.Washing({file}, key, false, true);
  return;
}

void events()
{
  TreeReader t;
  t.WashingOth({
      "root://cmd//sim/tr_ph_run045875_v5.root",
      "root://cmd//sim/tr_ph_run045672_v5.root",
      "root://cmd//sim/tr_ph_run045670_v5.root",
      "root://cmd//sim/tr_ph_run045669_v5.root",
      "root://cmd//sim/tr_ph_run045668_v5.root",
      "root://cmd//sim/tr_ph_run045667_v5.root",
      "root://cmd//sim/tr_ph_run045665_v5.root",
      "root://cmd//sim/tr_ph_run045664_v5.root",
      "root://cmd//sim/tr_ph_run045660_v5.root",
      "root://cmd//sim/tr_ph_run045659_v5.root",
      "root://cmd//sim/tr_ph_run045658_v5.root",
      "root://cmd//sim/tr_ph_run045657_v5.root",
      "root://cmd//sim/tr_ph_run045656_v5.root",
      "root://cmd//sim/tr_ph_run045655_v5.root",
      "root://cmd//sim/tr_ph_run045654_v5.root",
      "root://cmd//sim/tr_ph_run045653_v5.root",
      "root://cmd//sim/tr_ph_run045652_v5.root",
      "root://cmd//sim/tr_ph_run045651_v5.root",
      "root://cmd//sim/tr_ph_run045650_v5.root",
      "root://cmd//sim/tr_ph_run045649_v5.root",
      "root://cmd//sim/tr_ph_run045648_v5.root",
      "root://cmd//sim/tr_ph_run045647_v5.root",
      "root://cmd//sim/tr_ph_run045643_v5.root",
      "root://cmd//sim/tr_ph_run045641_v5.root"      
  }, "others/multihadrons");
  return;
}
