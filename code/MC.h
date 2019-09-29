#ifndef MC_h
#define MC_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class MC {
public :
   TChain          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Float_t         emeas;
   Float_t         demeas;
   Int_t           nks;
   Int_t           ksvind[4][20];   //[nks]
   Int_t           kstype[4];   //[nks]

   // List of branches
   TBranch        *b_emeas;   //!
   TBranch        *b_demeas;   //!
   TBranch        *b_nks;   //!
   TBranch        *b_ksvind;   //!
   TBranch        *b_kstype;   //!

   MC(TChain *tree=0);
   virtual ~MC();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TChain *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef MC_cxx
MC::MC(TChain *tree) : fChain(0)
{
   Init(tree);
}

MC::~MC()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t MC::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t MC::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void MC::Init(TChain *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("emeas", &emeas, &b_emeas);
   fChain->SetBranchAddress("demeas", &demeas, &b_demeas);
   fChain->SetBranchAddress("nks", &nks, &b_nks);
   fChain->SetBranchAddress("ksvind", ksvind, &b_ksvind);
   fChain->SetBranchAddress("kstype", kstype, &b_kstype);
   Notify();
}

Bool_t MC::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void MC::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t MC::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef MC_cxx
