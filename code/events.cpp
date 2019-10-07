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

class TreeReader{
public:

    TreeReader(){
    }

    void Washing(std::vector<string> inputPath, string key, bool time=false){
        double time0 = clock()/(CLOCKS_PER_SEC+0.);

        MC cl(inputPath, key);
        cl.Loop();

        if(time) cout << "Washing time: " << clock()/(CLOCKS_PER_SEC+0.) - time0 << "seconds." << endl;
        return;
    }

    void WashingFromFile(string filePath, string key){
        ifstream f(filePath.c_str());
        string input;
        while( f >> input ){
            Washing({input}, key, true);
            cout << "Well done" << endl << input << endl << endl;
        }
        return;
    }

    void WashingModel(string allModels){
      WashingFromFile(allModels, "model");
      MC cl(allModels, "model");
      cl.GetSoftPhotonsNumber("../outputs/model/soft_ph.csv");
      return;
    }

    void Washing11(string all11){
      WashingFromFile(all11, "11");
    }
    void Washing12(string all12){
      WashingFromFile(all12, "12");
    }
    void Washing17(string all17){
      WashingFromFile(all17, "17");
    }
    void Washing19(string all19){
      WashingFromFile(all19, "19");
    }

    static Double_t func(Double_t* x, Double_t* par){
      double X = x[0];
      double mean = par[0];
      double sigma = par[1];
      double N = par[2];
      double C = par[3];

      return (N*TMath::Gaus(X, mean, sigma, kTRUE) + C);
    }

    void FitterOne(TH1D* h, string directory, pair<double, double>& N){
      h->Sumw2();
      h->Scale(1, "width");

      TCanvas c("can", "Canvas", 900, 600);
      TF1 f("f", this->func, 450, 550, 4);
      f.SetParameters(mKs, 10, h->GetEntries(), 0.1);
      h->Fit(&f, "ML");
      c.Print((directory+".svg").c_str());
      N.first = f.GetParameter(2);
      N.second = f.GetParError(2);
      return;
    }

    void Fitter(string directory, string lumfile){
      FILE* o; o = fopen(lumfile.c_str(), "r");
      ofstream out((directory + "/number.csv").c_str());
      double e, lum;
      char e_str[100];
      string filename;
      TChain *ch = new TChain("t");
      TH1D* h = new TH1D("h", "Mass distribution", 50, 450, 550);
      pair<double, double> N;


      while( fscanf(o, "%lf,%lf\n", &e, &lum) == 2 ){
        cout << "E: " << e << "\tLum: " << lum << endl;
        snprintf(e_str, sizeof(e_str), "/%.2f", e);
        filename = directory + e_str + ".root";
        h->Reset();
        ch->SetName("t");
        ch->Add(filename.c_str());
        if(ch->GetEntries()<20)
          continue;

        gROOT->SetBatch(kTRUE);
        ch->Draw("m>>h", "m>450&m<550");
        FitterOne(h, directory + "/img" + e_str, N);
        cout << e_str << '\t' << N.first << '\t' << N.second << endl;
        out << e << "," << N.first << "," << N.second << endl;
        ch->Reset();
        gROOT->SetBatch(kFALSE);
      }
    }
};


void events(){
    TreeReader t;
    //t.Washing11("../inputs/11/trees");
    t.Fitter("../outputs/11", "../inputs/11/lum");
    return;
}
