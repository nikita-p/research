#define MC_cxx
#include "MC.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void MC::Loop()
{

  TFile* f = TFile::Open("train.root", "recreate");
  TTree* t = new TTree("t", "Tree");
  t->Branch("nt", &nt, "nt/I");

  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  Long64_t nbytes = 0, nb = 0;
  for (Long64_t jentry=0; jentry<nentries;jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    t->Fill();
    // if (Cut(ientry) < 0) continue;
  }
  t->Write();
  return;
}
