void Mass()
{

   TFile *f1 = new TFile("result.root");
   TTree *t1 = (TTree *)f1->Get("treeMuon");
   TTree *t2 = (TTree *)f1->Get("treeJet");

   TLorentzVector *Muon_P4 = new TLorentzVector();
   TLorentzVector *Jet_P4  = new TLorentzVector();

   TLorentzVector vec;
   vec.SetPxPyPzE(0.0,0.0,0.0,240.0);

   t1->SetBranchAddress("Muon_P4",&Muon_P4);
   t2->SetBranchAddress("Jet_P4",&Jet_P4);

   TFile *f2 = new TFile("Mass/Z_H_Mass.root","RECREATE");
   TTree *Z_mass = new TTree("Z_mass","Z_mass");
   TTree *H_mass = new TTree("H_mass","H_mass");
   float z_mass;
   float h_mass;
   Z_mass->Branch("z_mass",&z_mass);
   H_mass->Branch("h_mass",&h_mass);


  int entries_Muon = t1->GetEntries();
  int entries_Jet = t2->GetEntries();

  cout << " entries_Muon= " << entries_Muon
       << " entries_Jet = " << entries_Jet
       << endl;

  // Z_mass
  for(int i=0;i<entries_Muon;i++) {
      t1->GetEntry(i);

      if(i%2==0){
          TLorentzVector Recoil1;
          TLorentzVector Recoil2;

          TLorentzVector muon1;// = *Muon_P4;
          muon1.SetPtEtaPhiM(Muon_P4->Pt(),Muon_P4->Eta(),Muon_P4->Phi(),0.0);

          TLorentzVector muon2;
          t1->GetEntry(i+1);

          muon2.SetPtEtaPhiM(Muon_P4->Pt(),Muon_P4->Eta(),Muon_P4->Phi(),0.0);// = *Muon_P4;
          z_mass = (muon1 + muon2).M();
          cout << "Z_mass= " << z_mass << endl;

          Z_mass->Fill();

      }
      else continue;
  }

  // H_mass
  for(int i=0;i<entries_Jet;i++) {
      t2->GetEntry(i);
      cout << " jet_p4= " << Jet_P4->P() << endl;
      if(i%2==0){
          TLorentzVector jet1 = *Jet_P4;
          cout << " jet1_p4= " << jet1.P() << endl;
          TLorentzVector jet2;
          t2->GetEntry(i+1);
          jet2 = *Jet_P4;
          cout << " jet2_p4= " << jet2.P() << endl;
          h_mass = (jet1 + jet2).M();
          H_mass->Fill();
          cout << "H_mass= " << h_mass << endl;
      }
  }


  f2->cd();
  Z_mass->Write();
  H_mass->Write();
  f2->Close();
}



