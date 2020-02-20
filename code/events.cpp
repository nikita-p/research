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

  void Washing(std::vector<string> inputPath, string key, bool time = false)
  {
    double time0 = clock() / (CLOCKS_PER_SEC + 0.);

    MC cl(inputPath, key);
    cout << "Go loop" << endl;
    cl.Loop();

    if (time)
      cout << "Washing time: " << clock() / (CLOCKS_PER_SEC + 0.) - time0 << "seconds." << endl;
    return;
  }

  void WashingFromFile(string filePath, string key)
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
      Washing({input}, key, true);
      cout << "Well done" << endl
           << endl;
    }
    return;
  }

  void WashingModelOldv6(string allModels)
  {
    WashingFromFile(allModels, "model_old_v6");
    MC cl(allModels, "model_old_v6");
    cl.GetSoftPhotonsNumber("../outputs/model_old_v6/soft_ph.csv");
    return;
  }
  
  void WashingModelOldv7(string allModels)
  {
    WashingFromFile(allModels, "model_old_v7");
    MC cl(allModels, "model_old_v7");
    cl.GetSoftPhotonsNumber("../outputs/model_old_v7/soft_ph.csv");
    return;
  }
  
  void WashingModelNewv6(string allModels)
  {
    WashingFromFile(allModels, "model_new_v6");
    MC cl(allModels, "model_new_v6");
    cl.GetSoftPhotonsNumber("../outputs/model_new_v6/soft_ph.csv");
    return;
  }

  void Washing11(string all11)
  {
    WashingFromFile(all11, "11");
  }
  void Washing12(string all12)
  {
    WashingFromFile(all12, "12");
  }
  void Washing17(string all17)
  {
    WashingFromFile(all17, "17");
  }
  void Washing19(string all19)
  {
    WashingFromFile(all19, "19");
  }
  void WashingOth(string allOth)
  {
    WashingFromFile(allOth, "others");
  }
};

void events()
{
  TreeReader t;
  t.WashingModelOldv7("../inputs/model_old_v7/trees");
  //t.WashingModelNewv6("../inputs/model_new_v6/trees");
  // t.Washing17("../inputs/17/trees");
  // t.Washing11("../inputs/11/trees");
  return;
}
