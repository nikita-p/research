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
    MC cl("../inputs/model_old_v6/trees", "model_old_v6", sys);
    string write_file = (string)"../outputs/model_old_v6/soft_ph" + (sys ? "_sys" : "") + ".csv";
//     cl.GetSoftPhotonsNumber(write_file);
    return;
  }
  
  void WashingModelOldv7(bool sys = false)
  {
    WashingFromFile("../inputs/model_old_v7/trees", "model_old_v7", sys);
    MC cl("../inputs/model_old_v7/trees", "model_old_v7", sys);
    string write_file = (string)"../outputs/model_old_v7/soft_ph" + (sys ? "_sys" : "") + ".csv";
//     cl.GetSoftPhotonsNumber(write_file);
    return;
  }
  
  void WashingModelNewv6(bool sys = false)
  {
    WashingFromFile("../inputs/model_new_v6/trees", "model_new_v6", sys);
    MC cl("../inputs/model_new_v6/trees", "model_new_v6", sys);
    string write_file = (string)"../outputs/model_new_v6/soft_ph" + (sys ? (string)"_sys" : (string)"") + (string)".csv";
//     cl.GetSoftPhotonsNumber(write_file);
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
  void WashingOth(std::vector<string> other_files, bool sys = false)
  {
    for(auto it=other_files.begin(); it!=other_files.end(); it++){
        cout << "File\n";
        Washing({*it}, "others", sys, true);
    }
    return;
  }
};

void events()
{
  TreeReader t;
//   t.WashingModelOldv7(false);
//   t.WashingModelNewv6();
   t.Washing19();
//    t.Washing11();
//   t.Washing11(true);
//   t.WashingOth({
//   "root://cmd//sim/tr_ph_run045069_v5.root",
//   "root://cmd//sim/tr_ph_run045078_v5.root",
//   "root://cmd//sim/tr_ph_run045076_v5.root",
//   "root://cmd//sim/tr_ph_run045067_v5.root"});
  return;
}
