/*
Simple macro showing how to access branches from the delphes output root file,
loop over events, and plot simple quantities such as the jet pt and the di-electron invariant
mass.

root -l examples/Example1.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#endif

//------------------------------------------------------------------------------

void Example1(const char *inputFile)
{
  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain chain("Delphes");
  chain.Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
  Long64_t numberOfEntries = treeReader->GetEntries();

cout << "numberOfEntries= " << numberOfEntries << endl;

  // Get pointers to branches used in this analysis
  TClonesArray *branchJet = treeReader->UseBranch("Jet");
  TClonesArray *branchMuon = treeReader->UseBranch("Muon");
  TClonesArray *branchEvent = treeReader->UseBranch("Event");

  // Book histograms
  TH1 *histJetPT = new TH1F("jet_pt", "jet P_{T}", 100, 0.0, 140.0);
  TH1 *histMass = new TH1F("mass", "M_{inv}(mu_{1}, mu_{2})", 100, 0.0, 140.0);

  // Loop over all events
  for(Int_t entry = 0; entry < numberOfEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    
    //HepMCEvent *event = (HepMCEvent*) branchEvent -> At(0);
    //LHEFEvent *event = (LHEFEvent*) branchEvent -> At(0);
    //Float_t weight = event->Weight;

    // If event contains at least 1 jet

    if(branchJet->GetEntries() > 0)
    {
      // Take first jet
      Jet *jet = (Jet*) branchJet->At(0);

      // Plot jet transverse momentum
      histJetPT->Fill(jet->PT);

      // Print jet transverse momentum
//      cout << "Jet pt: "<<jet->PT << endl;
    }

    Muon *muon1, *muon2;

    // If event contains at least 2 electrons
    if(branchMuon->GetEntries() > 1)
    {
      // Take first two electrons
      muon1 = (Muon *) branchMuon->At(0);
      muon2 = (Muon *) branchMuon->At(1);

      // Plot their invariant mass
      histMass->Fill(((muon1->P4()) + (muon2->P4())).M());
    cout << "mass = " << ((muon1->P4()) + (muon2->P4())).T() << endl;
    }
  }

  // Show resulting histograms
  //histJetPT->Draw();
  histMass->Draw();
}

