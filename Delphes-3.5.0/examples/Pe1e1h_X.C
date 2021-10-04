/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, store histograms in a root file and print them as image files.

root -l examples/Example2.C'("delphes_output.root")'
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include "TH1.h"
#include "TSystem.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#endif

//------------------------------------------------------------------------------

struct MyPlots
{
  TH1 *fMuonPT[500000];
};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, MyPlots *plots)
{
  TPaveText *comment;

  // book 2 histograms for PT of 1st and 2nd leading jets

  plots->fMuonPT[0] = result->AddHist1D(
    "Muon_PT", "Muon P_{T}",
    "Muon P_{T}, GeV", "number of events",
    80, 0.0, 80.0);

  // book general comment

  comment = result->AddComment(0.64, 0.86, 0.98, 0.98);
  comment->AddText("demonstration plot");
  comment->AddText("produced by Example2.C");

  // attach comment to single histograms

  result->Attach(plots->fMuonPT[0], comment);
}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
{
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  Muon *muon[1];

  Long64_t entry;

  Int_t i;

  int num =0;
  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    num += branchMuon->GetEntriesFast();

    for(i = 0; i < branchMuon->GetEntriesFast(); ++i)
    {
       muon[0] = (Muon*) branchMuon->At(i);
       plots->fMuonPT[0]->Fill(muon[0]->PT);
     }

  }

  cout << "** Muon contains " << num << " entries" << endl;
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, MyPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

//void Pe1e1h_X(const char *inputFile)
void Pe1e1h_X()
{

  char buf[10000];
  ifstream fp;
  fp.open("/afs/ihep.ac.cn/users/m/myliu/scratchfs/Delphes-3.5.0/all_txt/Pqq.txt");
  std::cout<<"while start\t"<<std::endl;

  int N= 13663;
  int filenum=1;
 // TChain *chain;
  ExRootResult *result = new ExRootResult();
  MyPlots *plots = new MyPlots;

  BookHistograms(result, plots);

//  for(int filenum=10;filenum<N;filenum++) {
  while (fp.getline(buf,200)) {

//      char *inputfile = Form("/cefs/higgs/myliu/Data/higgs/E240.Pe1e1h_X/delphes_output_higgs_Pe1e1h_000%d.root",filenum);
      cout<<" inputfile=\t"<<buf<<endl;

      gSystem->Load("libDelphes");

      TChain *chain = new TChain("Delphes");
      //  chain->Add(inputFile);
      chain->Add(buf);

      ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

cout << treeReader->GetEntries() << endl;

      AnalyseEvents(treeReader, plots);

      //  PrintHistograms(result, plots);

      if(filenum==13662) {
          result->Write("Data/result/results_Pqq.root");
      }
      filenum++;
      delete treeReader;
      delete chain;
  }
//  result->Write("Data/result/results_Pe2e2h_X.root");

  std::cout << "** Exiting..." << std::endl;


  delete plots;
std::cout << __LINE__ << std::endl;
  delete result;
std::cout << __LINE__ << std::endl;
//  delete treeReader;
//  delete chain;
}

//------------------------------------------------------------------------------
