/*
c++ -o read_01 `root-config --glibs --cflags` -lm read_01.cpp 
./read_01 LHEfile.lhe
*/

#include <cmath>
#include <vector>
#include "LHEF.h"
#include "TLorentzVector.h"
#include <algorithm>
#include <functional>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

using namespace std ;

struct ptSort: public std::binary_function<TLorentzVector, TLorentzVector, bool>
{
  bool operator() (TLorentzVector x, TLorentzVector y)
    {
      return x.Pt () < y.Pt () ;
    }
} ;

double helicityAngle(TLorentzVector p4mother, TLorentzVector p4, TLorentzVector p4daughter){
      TVector3 boost = -p4.BoostVector();
      p4mother.Boost(boost);
      p4daughter.Boost(boost);
      TVector3 p3daughter= p4daughter.Vect();
      TVector3 p3mother= p4mother.Vect();

      return p3daughter.Dot(p3mother)/((p3daughter.Mag())*(p3mother.Mag()));
}

void computeAngles(TLorentzVector p4M11, 
    TLorentzVector p4M12, 
    TLorentzVector p4M21, 
    TLorentzVector p4M22, 
    TLorentzVector all, 
    float& costhetastar, 
    float& costheta1, 
    float& costheta2, 
    float& Phi, 
    float& Phi1){

      //build Z 4-vectors
      TLorentzVector p4Z1 = p4M11 + p4M12;
      TLorentzVector p4Z2 = p4M21 + p4M22;

      //build H 4-vectors
      TLorentzVector p4H = p4Z1 + p4Z2; 

      //build mother of H
      TLorentzVector p4m(0,0,0,250);

      // -----------------------------------

      costhetastar=helicityAngle(p4m,p4H,p4Z1);
      costheta1=helicityAngle(p4H,p4Z1,p4M11);
      costheta2=helicityAngle(p4H,p4Z2,p4M21);

      /*
      // --------------------------- Phi and Phi1 (old phistar1 - azimuthal production angle)
      //    TVector3 boostX = -(thep4H.BoostVector());
      TLorentzVector p4M11_BX( p4M11 );
      TLorentzVector p4M12_BX( p4M12 );
      TLorentzVector p4M21_BX( p4M21 );
      TLorentzVector p4M22_BX( p4M22 );

      p4M11_BX.Boost( boostX );
      p4M12_BX.Boost( boostX );
      p4M21_BX.Boost( boostX );
      p4M22_BX.Boost( boostX );

      TVector3 tmp1 = p4M11_BX.Vect().Cross( p4M12_BX.Vect() );
      TVector3 tmp2 = p4M21_BX.Vect().Cross( p4M22_BX.Vect() );    

      TVector3 normal1_BX( tmp1.X()/tmp1.Mag(), tmp1.Y()/tmp1.Mag(), tmp1.Z()/tmp1.Mag() ); 
      TVector3 normal2_BX( tmp2.X()/tmp2.Mag(), tmp2.Y()/tmp2.Mag(), tmp2.Z()/tmp2.Mag() ); 

      // Phi
      TLorentzVector p4Z1_BX = p4M11_BX + p4M12_BX;    
      float tmpSgnPhi = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normal2_BX) );
      float sgnPhi = tmpSgnPhi/fabs(tmpSgnPhi);
      Phi = sgnPhi * acos( -1.*normal1_BX.Dot( normal2_BX) );

      //////////////

      //TVector3 beamAxis(0,0,1);
      TVector3 beamAxis(all.Px()/all.Vect().Mag(),all.Py()/all.Vect().Mag(),all.Pz()/all.Vect().Mag());
      TVector3 tmp3 = (p4M11_BX + p4M12_BX).Vect();

      TVector3 p3V1_BX( tmp3.X()/tmp3.Mag(), tmp3.Y()/tmp3.Mag(), tmp3.Z()/tmp3.Mag() );
      TVector3 tmp4 = beamAxis.Cross( p3V1_BX );
      TVector3 normalSC_BX( tmp4.X()/tmp4.Mag(), tmp4.Y()/tmp4.Mag(), tmp4.Z()/tmp4.Mag() );

      // Phi1
      float tmpSgnPhi1 = p4Z1_BX.Vect().Dot( normal1_BX.Cross( normalSC_BX) );
      float sgnPhi1 = tmpSgnPhi1/fabs(tmpSgnPhi1);    
      Phi1 = sgnPhi1 * acos( normal1_BX.Dot( normalSC_BX) );    
      */
}


//int lhe_conversion_HZZ(TString argv="outfile_HZZ_SM.lhe") {
int main (int argc, char **argv) {

  // Open a stream connected to an event file:
  if (argc < 2) exit (1) ;
  //std::ifstream ifs(argv);
  std::ifstream ifs(argv[1]);

  // Create the Reader object:
  LHEF::Reader reader(ifs);

  // Copy header and init blocks and write them out.
  if ( reader.outsideBlock.length() ) std::cout << reader.outsideBlock;

  // Print out the header information:
  std::cerr << reader.headerBlock;

  // Print out the addinional comments in the init block:
  std::cerr << reader.initComments;

  // Print out the beam energies:
  std::cerr << "Beam A: " << reader.heprup.EBMUP.first << " GeV, Beam B: "
            << reader.heprup.EBMUP.second << " GeV." << std::endl;

  double weight;
  int nparticles;

  std::vector<TLorentzVector> P4_ElectronPlus;
  std::vector<TLorentzVector> P4_ElectronMinus;
  std::vector<TLorentzVector> P4_Z;
  std::vector<TLorentzVector> P4_Higgs;
  std::vector<TLorentzVector> P4_MuonZPlus;
  std::vector<TLorentzVector> P4_MuonZMinus;
  std::vector<TLorentzVector> P4_Muon1Plus;
  std::vector<TLorentzVector> P4_Muon1Minus;
  std::vector<TLorentzVector> P4_Muon2Plus;
  std::vector<TLorentzVector> P4_Muon2Minus;

  double ElectronPlusPx, ElectronPlusPy, ElectronPlusPz, ElectronPlusE;
  double ElectronMinusPx, ElectronMinusPy, ElectronMinusPz, ElectronMinusE;
  double ZPx, ZPy, ZPz, ZE;
  double HiggsPx, HiggsPy, HiggsPz, HiggsE;
  double MuonZPlusPx, MuonZPlusPy, MuonZPlusPz, MuonZPlusE;
  double MuonZMinusPx, MuonZMinusPy, MuonZMinusPz, MuonZMinusE;
  double Muon1PlusPx, Muon1PlusPy, Muon1PlusPz, Muon1PlusE;
  double Muon1MinusPx, Muon1MinusPy, Muon1MinusPz, Muon1MinusE;
  double Muon2PlusPx, Muon2PlusPy, Muon2PlusPz, Muon2PlusE;
  double Muon2MinusPx, Muon2MinusPy, Muon2MinusPz, Muon2MinusE;

  TFile *outfile = new TFile("trial_HZZ.root","RECREATE");
  TTree *out_tree = new TTree("trialTree","trialTree");
  TH1D *m_higgs = new TH1D("m_higgs","m_higgs",40,123,127);
  TH1D *m_Zboson = new TH1D("m_Zboson","m_Zboson",20,80,100);
  TH1D *cos_theta1 = new TH1D("cos_theta1","cos_theta1",25,-1,1);
  TH1D *cos_theta2 = new TH1D("cos_theta2","cos_theta2",25,-1,1);
  TH1D *cos_thetastar = new TH1D("cos_thetastar","cos_thetastar",25,-1,1);
  TH1D *plot_phi = new TH1D("phi","phi",50,-3.2,3.2);
  TLorentzVector Muon1Plus, Muon1Minus, Muon2Plus, Muon2Minus, MuonZPlus, MuonZMinus;
  TLorentzVector Z, Higgs, ElectronPlus, ElectronMinus;

  out_tree->Branch("weight",&weight,"weight/D");
  out_tree->Branch("nparticles",&nparticles,"nparticles/I");

  out_tree->Branch("ElectronPlusPx",&ElectronPlusPx,"ElectronPlusPx/D");
  out_tree->Branch("ElectronPlusPy",&ElectronPlusPy,"ElectronPlusPy/D");
  out_tree->Branch("ElectronPlusPz",&ElectronPlusPz,"ElectronPlusPz/D");
  out_tree->Branch("ElectronPlusE",&ElectronPlusE,"ElectronPlusE/D");

  out_tree->Branch("ElectronMinusPx",&ElectronMinusPx,"ElectronMinusPx/D");
  out_tree->Branch("ElectronMinusPy",&ElectronMinusPy,"ElectronMinusPy/D");
  out_tree->Branch("ElectronMinusPz",&ElectronMinusPz,"ElectronMinusPz/D");
  out_tree->Branch("ElectronMinusE",&ElectronMinusE,"ElectronMinusE/D");

  out_tree->Branch("ZPx",&ZPx,"ZPx/D");
  out_tree->Branch("ZPy",&ZPy,"ZPy/D");
  out_tree->Branch("ZPz",&ZPz,"ZPz/D");
  out_tree->Branch("ZE",&ZE,"ZE/D");

  out_tree->Branch("HiggsPx",&HiggsPx,"HiggsPx/D");
  out_tree->Branch("HiggsPy",&HiggsPy,"HiggsPy/D");
  out_tree->Branch("HiggsPz",&HiggsPz,"HiggsPz/D");
  out_tree->Branch("HiggsE",&HiggsE,"HiggsE/D");

  out_tree->Branch("MuonZPlusPx",&MuonZPlusPx,"MuonZPlusPx/D");
  out_tree->Branch("MuonZPlusPy",&MuonZPlusPy,"MuonZPlusPy/D");
  out_tree->Branch("MuonZPlusPz",&MuonZPlusPz,"MuonZPlusPz/D");
  out_tree->Branch("MuonZPlusE",&MuonZPlusE,"MuonZPlusE/D");

  out_tree->Branch("MuonZMinusPx",&MuonZMinusPx,"MuonZMinusPx/D");
  out_tree->Branch("MuonZMinusPy",&MuonZMinusPy,"MuonZMinusPy/D");
  out_tree->Branch("MuonZMinusPz",&MuonZMinusPz,"MuonZMinusPz/D");
  out_tree->Branch("MuonZMinusE",&MuonZMinusE,"MuonZMinusE/D");

  out_tree->Branch("Muon1PlusPx",&Muon1PlusPx,"Muon1PlusPx/D");
  out_tree->Branch("Muon1PlusPy",&Muon1PlusPy,"Muon1PlusPy/D");
  out_tree->Branch("Muon1PlusPz",&Muon1PlusPz,"Muon1PlusPz/D");
  out_tree->Branch("Muon1PlusE",&Muon1PlusE,"Muon1PlusE/D");

  out_tree->Branch("Muon1MinusPx",&Muon1MinusPx,"Muon1MinusPx/D");
  out_tree->Branch("Muon1MinusPy",&Muon1MinusPy,"Muon1MinusPy/D");
  out_tree->Branch("Muon1MinusPz",&Muon1MinusPz,"Muon1MinusPz/D");
  out_tree->Branch("Muon1MinusE",&Muon1MinusE,"Muon1MinusE/D");

  out_tree->Branch("Muon2PlusPx",&Muon2PlusPx,"Muon2PlusPx/D");
  out_tree->Branch("Muon2PlusPy",&Muon2PlusPy,"Muon2PlusPy/D");
  out_tree->Branch("Muon2PlusPz",&Muon2PlusPz,"Muon2PlusPz/D");
  out_tree->Branch("Muon2PlusE",&Muon2PlusE,"Muon2PlusE/D");

  out_tree->Branch("Muon2MinusPx",&Muon2MinusPx,"Muon2MinusPx/D");
  out_tree->Branch("Muon2MinusPy",&Muon2MinusPy,"Muon2MinusPy/D");
  out_tree->Branch("Muon2MinusPz",&Muon2MinusPz,"Muon2MinusPz/D");
  out_tree->Branch("Muon2MinusE",&Muon2MinusE,"Muon2MinusE/D");

  P4_ElectronPlus.clear();
  P4_ElectronMinus.clear();
  P4_Z.clear();
  P4_Higgs.clear();
  P4_MuonZPlus.clear();
  P4_MuonZMinus.clear();
  P4_Muon1Plus.clear();
  P4_Muon1Minus.clear();
  P4_Muon2Plus.clear();
  P4_Muon2Minus.clear();

  // Now loop over all events:
  long ieve = 0;
  while ( reader.readEvent() ) 
    {
      ++ieve;
      if (ieve % 100 == 0) std::cerr << "event " << ieve << "\n" ;
  
      int MuonPlus_Counter=0, MuonMinus_Counter=0;
      //PG loop over particles in the event
      for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); iPart++)
        {  

          //if (reader.hepeup.ISTUP.at (iPart) < 0) continue ;  // for removing initial particle beams

          int type = reader.hepeup.IDUP.at (iPart) ;
          TLorentzVector particle 
             (
               reader.hepeup.PUP.at (iPart).at (0), //PG px
               reader.hepeup.PUP.at (iPart).at (1), //PG py
               reader.hepeup.PUP.at (iPart).at (2), //PG pz
               reader.hepeup.PUP.at (iPart).at (3)  //PG E
             ) ;

          weight = reader.hepeup.XWGTUP;
          nparticles = reader.hepeup.NUP;
 
          if (type == -11){ 
            P4_ElectronPlus.push_back(particle) ;
            ElectronPlus = particle;
            ElectronPlusPx = particle.Px();
            ElectronPlusPy = particle.Py();
            ElectronPlusPz = particle.Pz();
            ElectronPlusE = particle.E();
            continue;
            }//positron

          if (type == 11){ 
            P4_ElectronMinus.push_back(particle);
            ElectronMinus = particle;
            ElectronMinusPx = particle.Px();
            ElectronMinusPy = particle.Py();
            ElectronMinusPz = particle.Pz();
            ElectronMinusE = particle.E();
            continue;
          }//electron

          if (type == 23){ 
            P4_Z.push_back(particle);
            Z = particle;
            ZPx = particle.Px();
            ZPy = particle.Py();
            ZPz = particle.Pz();
            ZE = particle.E();
            m_Zboson->Fill(particle.M(),weight);
            continue;
          }

          if (type == 25){ 
            P4_Higgs.push_back(particle) ;
            Higgs = particle;
            HiggsPx = particle.Px();
            HiggsPy = particle.Py();
            HiggsPz = particle.Pz();
            HiggsE = particle.E();
            m_higgs->Fill(particle.M(),weight);
            continue;
          }           

          if (type == -13 && MuonPlus_Counter==0){ 
            P4_MuonZPlus.push_back(particle) ;
            MuonZPlus = particle;
            MuonZPlusPx = particle.Px();
            MuonZPlusPy = particle.Py();
            MuonZPlusPz = particle.Pz();
            MuonZPlusE = particle.E();
            MuonPlus_Counter=1;
            continue;
          }

          if (type == 13 && MuonMinus_Counter==0){ 
            P4_MuonZMinus.push_back(particle) ;
            MuonZMinus = particle;
            MuonZMinusPx = particle.Px();
            MuonZMinusPy = particle.Py();
            MuonZMinusPz = particle.Pz();
            MuonZMinusE = particle.E();
            MuonMinus_Counter=1;
            continue;
          }

          if (type == -13 && MuonPlus_Counter==1){    
            P4_Muon1Plus.push_back(particle) ;
            Muon1Plus = particle;
            Muon1PlusPx = particle.Px();
            Muon1PlusPy = particle.Py();
            Muon1PlusPz = particle.Pz();
            Muon1PlusE = particle.E();
            MuonPlus_Counter=2;
            continue;
          }

          if (type == 13 && MuonMinus_Counter==1){
            P4_Muon1Minus.push_back(particle) ;
            Muon1Minus = particle;
            Muon1MinusPx = particle.Px();
            Muon1MinusPy = particle.Py();
            Muon1MinusPz = particle.Pz();
            Muon1MinusE = particle.E();
            MuonMinus_Counter=2;
            continue;
          }

          if (type == -13 && MuonPlus_Counter==2){    
            P4_Muon2Plus.push_back(particle) ;
            Muon2Plus = particle;
            Muon2PlusPx = particle.Px();
            Muon2PlusPy = particle.Py();
            Muon2PlusPz = particle.Pz();
            Muon2PlusE = particle.E();
            MuonPlus_Counter=3;
            continue;
          }

          if (type == 13 && MuonMinus_Counter==2){
            P4_Muon2Minus.push_back(particle) ;
            Muon2Minus = particle;
            Muon2MinusPx = particle.Px();
            Muon2MinusPy = particle.Py();
            Muon2MinusPz = particle.Pz();
            Muon2MinusE = particle.E();
            MuonMinus_Counter=3;
            continue;
          }

      } //PG loop over particles in the event

    TLorentzVector null;
    float costhetastar;
    float costheta1, costheta2;
    float Phi, Phi1;
    computeAngles(Muon1Plus, Muon1Minus, Muon2Plus, Muon2Minus, null, costhetastar, costheta1, costheta2, Phi, Phi1);
    cos_theta1->Fill(costheta1,weight);
    cos_theta2->Fill(costheta2,weight);
    cos_thetastar->Fill(costhetastar,weight);
    plot_phi->Fill(Phi,weight);

    out_tree->Fill();

    } // Now loop over all events

//  cos_theta1->Draw();
//  cos_theta2->Draw();
//  cos_thetastar->Draw();
//  plot_phi->Draw();
//  m_Zboson->Draw();
//  m_higgs->Draw();
  outfile->Write();
  outfile->Close();
  // Now we are done.
  return 0 ;
}

/*

TCanvas c2
gStyle->SetPalette (1)


.L ReverseCumultive.C 

void ReverseCumulative (TH1* histo) 
{
  double integral = histo->GetBinContent (histo->GetNbinsX ()) ;
  for (int iBin = 0; iBin < histo->GetNbinsX (); iBin++) 
    {
     Êdouble value = histo->GetBinContent (iBin+1) ;
     Êhisto->SetBinContent (iBin+1, integral - value) ; 
    }
}

gStyle->SetPalette (1)
.L TH2FCumulative.C 
*/
