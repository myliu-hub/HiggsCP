#include <iostream>
#include <sstream>
#include <stdexcept>

#include <signal.h>
#include <map>
#include <vector>


#include "TApplication.h"
#include "TROOT.h"

#include "TDatabasePDG.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "TObjArray.h"
#include "TParticlePDG.h"
#include "TStopwatch.h"
#include "TRef.h"
#include "TProcessID.h"

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "classes/DelphesFactory.h"
#include "classes/DelphesLHEFReader.h"
#include "modules/Delphes.h"

#include "ExRootAnalysis/ExRootProgressBar.h"
#include "ExRootAnalysis/ExRootTreeBranch.h"
#include "ExRootAnalysis/ExRootTreeWriter.h"
#endif

using namespace std;

//-----------------------------------------------------

typedef struct JetBtag
{
    int Num;
    int  Btag;
    double Jet_PT;
    JetBtag(int n, int B,double p):Num(n),Btag(B),Jet_PT(p){}
    bool operator < (const JetBtag &j) const
    {
       return Num < j.Num;
    }
} jetBtag,jetBtag2;

std::map<int,jetBtag> jetbtag;
std::map<int,jetBtag2> jetbtag2;

void ReadDelRoot()
{


  char buf[10000];
  std::ifstream fp;
  fp.open("/afs/ihep.ac.cn/users/m/myliu/scratchfs/Delphes-3.5.0/all_txt/Pzz_sl.txt");
  std::cout<<"while start\t"<<std::endl;


  TFile *outputFile = new TFile("/cefs/higgs/myliu/Data/sample/Sample_cut/2muon_2bjet_Pzz_sl.root","RECREATE");
  //TTree *t1 = new TTree("t1","t1");

  TTree *Z_mass = new TTree("Z_mass","Z_mass");
  float z_mass;
  Z_mass->Branch("z_mass",&z_mass);

  // Event 
  TTree *treeEvent = new TTree("treeEvent","treeEvent");
  Long64_t Event_Number;
  float Event_ReadTime,Event_ProcTime;

  treeEvent->Branch("Event_Number",&Event_Number);
  treeEvent->Branch("Event_ReadTime",&Event_ReadTime);
  treeEvent->Branch("Event_ProcTime",&Event_ProcTime);

  // Particle
  TTree *treeParticle = new TTree("treeParticle","treeParticle");

  int Particle_PID,Particle_Status,Particle_IsPU;
  int Particle_M1,Particle_M2,Particle_D1,Particle_D2,Particle_Charge;
  float Particle_Mass,Particle_E,Particle_Px,Particle_Py,Particle_Pz;
  float Particle_P,Particle_PT,Particle_Eta,Particle_Phi,Particle_Rapidity,Particle_T;
  float Particle_X,Particle_Y,Particle_Z;
  TLorentzVector Particle_P4;

  treeParticle->Branch("Particle_PID",&Particle_PID);
  treeParticle->Branch("Particle_Status",&Particle_Status);
  treeParticle->Branch("Particle_IsPU",&Particle_IsPU);
  treeParticle->Branch("Particle_M1",&Particle_M1);
  treeParticle->Branch("Particle_M2",&Particle_M2);
  treeParticle->Branch("Particle_D1",&Particle_D1);
  treeParticle->Branch("Particle_D2",&Particle_D2);
  treeParticle->Branch("Particle_Charge",&Particle_Charge);
  treeParticle->Branch("Particle_Mass",&Particle_Mass);
  treeParticle->Branch("Particle_E",&Particle_E);
  treeParticle->Branch("Particle_Px",&Particle_Px);
  treeParticle->Branch("Particle_Py",&Particle_Py);
  treeParticle->Branch("Particle_Pz",&Particle_Pz);
  treeParticle->Branch("Particle_P",&Particle_P);
  treeParticle->Branch("Particle_PT",&Particle_PT);
  treeParticle->Branch("Particle_Eta",&Particle_Eta);
  treeParticle->Branch("Particle_Phi",&Particle_Phi);
  treeParticle->Branch("Particle_Rapidity",&Particle_Rapidity);
  treeParticle->Branch("Particle_T",&Particle_T);
  treeParticle->Branch("Particle_X",&Particle_X);
  treeParticle->Branch("Particle_Y",&Particle_Y);
  treeParticle->Branch("Particle_Z",&Particle_Z);
  treeParticle->Branch("Particle_P4",&Particle_P4);

  // GenJet
  TTree *treeGenJet = new TTree("treeGenJet","treeGenJet");

  float GenJet_PT,GenJet_Eta,GenJet_Phi,GenJet_T,GenJet_Mass,GenJet_DeltaEta,GenJet_DeltaPhi;
  int GenJet_Flavor,GenJet_FlavorAlgo,GenJet_FlavorPhys,GenJet_BTag,GenJet_BTagAlgo,GenJet_BTagPhys,GenJet_TauTag,GenJet_Charge,GenJet_NCharged,GenJet_NNeutrals,GenJet_NSubJetsTrimmed,GenJet_NSubJetsPruned,GenJet_NSubJetsSoftDropped;
  double GenJet_ExclYmerge23,GenJet_ExclYmerge34,GenJet_ExclYmerge45,GenJet_ExclYmerge56;
  float GenJet_TauWeight,GenJet_EhadOverEem,GenJet_NeutralEnergyFraction,GenJet_ChargedEnergyFraction,GenJet_Beta,GenJet_BetaStar,GenJet_MeanSqDeltaR,GenJet_PTD;
  float GenJet_FracPt[5],GenJet_Tau[5];
  TLorentzVector GenJet_SoftDroppedJet,GenJet_SoftDroppedSubJet1,GenJet_SoftDroppedSubJet2;
  TLorentzVector GenJet_TrimmedP4[5],GenJet_PrunedP4[5],GenJet_SoftDroppedP4[5];
  TRefArray GenJet_Constituents,GenJet_Particles;
  TLorentzVector GenJet_P4,GenJet_Area;

  treeGenJet->Branch("GenJet_PT",&GenJet_PT);
  treeGenJet->Branch("GenJet_Eta",&GenJet_Eta);
  treeGenJet->Branch("GenJet_Phi",&GenJet_Phi);
  treeGenJet->Branch("GenJet_T",&GenJet_T);
  treeGenJet->Branch("GenJet_Mass",&GenJet_Mass);
  treeGenJet->Branch("GenJet_DeltaEta",&GenJet_DeltaEta);
  treeGenJet->Branch("GenJet_DeltaPhi",&GenJet_DeltaPhi);
  treeGenJet->Branch("GenJet_Flavor",&GenJet_Flavor);
  treeGenJet->Branch("GenJet_FlavorAlgo",&GenJet_FlavorAlgo);
  treeGenJet->Branch("GenJet_FlavorPhys",&GenJet_FlavorPhys);
  treeGenJet->Branch("GenJet_BTag",&GenJet_BTag);
  treeGenJet->Branch("GenJet_BTagAlgo",&GenJet_BTagAlgo);
  treeGenJet->Branch("GenJet_BTagPhys",&GenJet_BTagPhys);
  treeGenJet->Branch("GenJet_TauTag",&GenJet_TauTag);
  treeGenJet->Branch("GenJet_TauWeight",&GenJet_TauWeight);
  treeGenJet->Branch("GenJet_Charge",&GenJet_Charge);
  treeGenJet->Branch("GenJet_EhadOverEem",&GenJet_EhadOverEem);
  treeGenJet->Branch("GenJet_NCharged",&GenJet_NCharged);
  treeGenJet->Branch("GenJet_NNeutrals",&GenJet_NNeutrals);
  treeGenJet->Branch("GenJet_NeutralEnergyFraction",&GenJet_NeutralEnergyFraction);
  treeGenJet->Branch("GenJet_ChargedEnergyFraction",&GenJet_ChargedEnergyFraction);
  treeGenJet->Branch("GenJet_Beta",&GenJet_Beta);
  treeGenJet->Branch("GenJet_BetaStar",&GenJet_BetaStar);
  treeGenJet->Branch("GenJet_MeanSqDeltaR",&GenJet_MeanSqDeltaR);
  treeGenJet->Branch("GenJet_PTD",&GenJet_PTD);
  treeGenJet->Branch("GenJet_FracPt",GenJet_FracPt);
  treeGenJet->Branch("GenJet_Tau",GenJet_Tau);
  treeGenJet->Branch("GenJet_SoftDroppedJet",&GenJet_SoftDroppedJet);
  treeGenJet->Branch("GenJet_SoftDroppedSubJet1",&GenJet_SoftDroppedSubJet1);
  treeGenJet->Branch("GenJet_SoftDroppedSubJet2",&GenJet_SoftDroppedSubJet2);
  treeGenJet->Branch("GenJet_TrimmedP4",GenJet_TrimmedP4);
  treeGenJet->Branch("GenJet_PrunedP4",GenJet_PrunedP4);
  treeGenJet->Branch("GenJet_SoftDroppedP4",GenJet_SoftDroppedP4);
  treeGenJet->Branch("GenJet_NSubJetsTrimmed",&GenJet_NSubJetsTrimmed);
  treeGenJet->Branch("GenJet_NSubJetsPruned",&GenJet_NSubJetsPruned);
  treeGenJet->Branch("GenJet_NSubJetsSoftDropped",&GenJet_NSubJetsSoftDropped);
  treeGenJet->Branch("GenJet_ExclYmerge23",&GenJet_ExclYmerge23);
  treeGenJet->Branch("GenJet_ExclYmerge34",&GenJet_ExclYmerge34);
  treeGenJet->Branch("GenJet_ExclYmerge45",&GenJet_ExclYmerge45);
  treeGenJet->Branch("GenJet_ExclYmerge56",&GenJet_ExclYmerge56);
  treeGenJet->Branch("GenJet_Constituents",&GenJet_Constituents);
  treeGenJet->Branch("GenJet_Particles",&GenJet_Particles);
  treeGenJet->Branch("GenJet_P4",&GenJet_P4);
  treeGenJet->Branch("GenJet_Area",&GenJet_Area);


  // GenMissET
  TTree *treeGenMissingET = new TTree("treeGenMissingET","treeGenMissingET");

  float GenMissingET_MET,GenMissingET_Eta,GenMissingET_Phi;
  TLorentzVector GenMissingET_P4;

  treeGenMissingET->Branch("GenMissingET_MET",&GenMissingET_MET);
  treeGenMissingET->Branch("GenMissingET_Eta",&GenMissingET_Eta);
  treeGenMissingET->Branch("GenMissingET_Phi",&GenMissingET_Phi);
  treeGenMissingET->Branch("GenMissingET_P4",&GenMissingET_P4);

  // Track
  TTree *treeTrack = new TTree("treeTrack","treeTrack");

  int Track_PID,Track_Charge,Track_VertexIndex;
  float Track_P,Track_PT,Track_Eta,Track_Phi,Track_CtgTheta,Track_C,Track_Mass;
  float Track_EtaOuter,Track_PhiOuter,Track_T,Track_X,Track_Y,Track_Z;
  float Track_TOuter,Track_XOuter,Track_YOuter,Track_ZOuter,Track_Xd,Track_Yd,Track_Zd;
  float Track_L,Track_D0,Track_DZ,Track_Nclusters,Track_dNdx,Track_ErrorP,Track_ErrorPT;
  float Track_ErrorPhi,Track_ErrorCtgTheta,Track_ErrorT,Track_ErrorD0,Track_ErrorDZ;
  float Track_ErrorC,Track_ErrorD0Phi,Track_ErrorD0C,Track_ErrorD0DZ,Track_ErrorD0CtgTheta;
  float Track_ErrorPhiC,Track_ErrorPhiDZ,Track_ErrorPhiCtgTheta,Track_ErrorCDZ,Track_ErrorCCtgTheta,Track_ErrorDZCtgTheta;
  TRef Track_Particle;
  TLorentzVector Track_P4;
  TMatrixDSym Track_CovarianceMatrix;

  treeTrack->Branch("Track_PID",&Track_PID);
  treeTrack->Branch("Track_Charge",&Track_Charge);
  treeTrack->Branch("Track_P",&Track_P);
  treeTrack->Branch("Track_PT",&Track_PT);
  treeTrack->Branch("Track_Eta",&Track_Eta);
  treeTrack->Branch("Track_Phi",&Track_Phi);
  treeTrack->Branch("Track_CtgTheta",&Track_CtgTheta);
  treeTrack->Branch("Track_C",&Track_C);
  treeTrack->Branch("Track_Mass",&Track_Mass);
  treeTrack->Branch("Track_EtaOuter",&Track_EtaOuter);
  treeTrack->Branch("Track_PhiOuter",&Track_PhiOuter);
  treeTrack->Branch("Track_T",&Track_T);
  treeTrack->Branch("Track_X",&Track_X);
  treeTrack->Branch("Track_Y",&Track_Y);
  treeTrack->Branch("Track_Z",&Track_Z);
  treeTrack->Branch("Track_TOuter",&Track_TOuter);
  treeTrack->Branch("Track_XOuter",&Track_XOuter);
  treeTrack->Branch("Track_YOuter",&Track_YOuter);
  treeTrack->Branch("Track_ZOuter",&Track_ZOuter);
  treeTrack->Branch("Track_Xd",&Track_Xd);
  treeTrack->Branch("Track_Yd",&Track_Yd);
  treeTrack->Branch("Track_Zd",&Track_Zd);
  treeTrack->Branch("Track_L",&Track_L);
  treeTrack->Branch("Track_D0",&Track_D0);
  treeTrack->Branch("Track_DZ",&Track_DZ);
  treeTrack->Branch("Track_Nclusters",&Track_Nclusters);
  treeTrack->Branch("Track_dNdx",&Track_dNdx);
  treeTrack->Branch("Track_ErrorP",&Track_ErrorP);
  treeTrack->Branch("Track_ErrorPT",&Track_ErrorPT);
  treeTrack->Branch("Track_ErrorPhi",&Track_ErrorPhi);
  treeTrack->Branch("Track_ErrorCtgTheta",&Track_ErrorCtgTheta);
  treeTrack->Branch("Track_ErrorT",&Track_ErrorT);
  treeTrack->Branch("Track_ErrorD0",&Track_ErrorD0);
  treeTrack->Branch("Track_ErrorDZ",&Track_ErrorDZ);
  treeTrack->Branch("Track_ErrorC",&Track_ErrorC);
  treeTrack->Branch("Track_ErrorD0Phi",&Track_ErrorD0Phi);
  treeTrack->Branch("Track_ErrorD0C",&Track_ErrorD0C);
  treeTrack->Branch("Track_ErrorD0DZ",&Track_ErrorD0DZ);
  treeTrack->Branch("Track_ErrorD0CtgTheta",&Track_ErrorD0CtgTheta);
  treeTrack->Branch("Track_ErrorPhiC",&Track_ErrorPhiC);
  treeTrack->Branch("Track_ErrorPhiDZ",&Track_ErrorPhiDZ);
  treeTrack->Branch("Track_ErrorPhiCtgTheta",&Track_ErrorPhiCtgTheta);
  treeTrack->Branch("Track_ErrorCDZ",&Track_ErrorCDZ);
  treeTrack->Branch("Track_ErrorCCtgTheta",&Track_ErrorCCtgTheta);
  treeTrack->Branch("Track_ErrorDZCtgTheta",&Track_ErrorDZCtgTheta);
  treeTrack->Branch("Track_Particle",&Track_Particle);
  treeTrack->Branch("Track_VertexIndex",&Track_VertexIndex);
  treeTrack->Branch("Track_P4",&Track_P4);
  treeTrack->Branch("Track_CovarianceMatrix",&Track_CovarianceMatrix);

  // Tower
  TTree *treeTower = new TTree("treeTower","treeTower");

  float Tower_ET,Tower_Eta,Tower_Phi,Tower_E,Tower_T;
  int Tower_NTimeHits;
  float Tower_Eem,Tower_Ehad,Tower_Edges[4];
  TRefArray Tower_Particles;
  TLorentzVector Tower_P4;

  treeTower->Branch("Tower_ET",&Tower_ET);
  treeTower->Branch("Tower_Eta",&Tower_Eta);
  treeTower->Branch("Tower_Phi",&Tower_Phi);
  treeTower->Branch("Tower_E",&Tower_E);
  treeTower->Branch("Tower_T",&Tower_T);
  treeTower->Branch("Tower_NTimeHits",&Tower_NTimeHits);
  treeTower->Branch("Tower_Eem",&Tower_Eem);
  treeTower->Branch("Tower_Ehad",&Tower_Ehad);
  treeTower->Branch("Tower_Edges",Tower_Edges);
  treeTower->Branch("Tower_Particles",&Tower_Particles);
  treeTower->Branch("Tower_P4",&Tower_P4);

  // EFlowTrack
  TTree *treeEFlowTrack = new TTree("treeEFlowTrack","treeEFlowTrack");

  int EFlowTrack_PID,EFlowTrack_Charge,EFlowTrack_VertexIndex;
  float EFlowTrack_P,EFlowTrack_PT,EFlowTrack_Eta,EFlowTrack_Phi,EFlowTrack_CtgTheta,EFlowTrack_C,EFlowTrack_Mass;
  float EFlowTrack_EtaOuter,EFlowTrack_PhiOuter,EFlowTrack_T,EFlowTrack_X,EFlowTrack_Y,EFlowTrack_Z;
  float EFlowTrack_TOuter,EFlowTrack_XOuter,EFlowTrack_YOuter,EFlowTrack_ZOuter,EFlowTrack_Xd,EFlowTrack_Yd,EFlowTrack_Zd;
  float EFlowTrack_L,EFlowTrack_D0,EFlowTrack_DZ,EFlowTrack_Nclusters,EFlowTrack_dNdx,EFlowTrack_ErrorP,EFlowTrack_ErrorPT;
  float EFlowTrack_ErrorPhi,EFlowTrack_ErrorCtgTheta,EFlowTrack_ErrorT,EFlowTrack_ErrorD0,EFlowTrack_ErrorDZ;
  float EFlowTrack_ErrorC,EFlowTrack_ErrorD0Phi,EFlowTrack_ErrorD0C,EFlowTrack_ErrorD0DZ,EFlowTrack_ErrorD0CtgTheta;
  float EFlowTrack_ErrorPhiC,EFlowTrack_ErrorPhiDZ,EFlowTrack_ErrorPhiCtgTheta,EFlowTrack_ErrorCDZ,EFlowTrack_ErrorCCtgTheta,EFlowTrack_ErrorDZCtgTheta;
  TRef EFlowTrack_Particle;
  TLorentzVector EFlowTrack_P4;
  TMatrixDSym EFlowTrack_CovarianceMatrix;

  treeEFlowTrack->Branch("EFlowTrack_PID",&EFlowTrack_PID);
  treeEFlowTrack->Branch("EFlowTrack_Charge",&EFlowTrack_Charge);
  treeEFlowTrack->Branch("EFlowTrack_P",&EFlowTrack_P);
  treeEFlowTrack->Branch("EFlowTrack_PT",&EFlowTrack_PT);
  treeEFlowTrack->Branch("EFlowTrack_Eta",&EFlowTrack_Eta);
  treeEFlowTrack->Branch("EFlowTrack_Phi",&EFlowTrack_Phi);
  treeEFlowTrack->Branch("EFlowTrack_CtgTheta",&EFlowTrack_CtgTheta);
  treeEFlowTrack->Branch("EFlowTrack_C",&EFlowTrack_C);
  treeEFlowTrack->Branch("EFlowTrack_Mass",&EFlowTrack_Mass);
  treeEFlowTrack->Branch("EFlowTrack_EtaOuter",&EFlowTrack_EtaOuter);
  treeEFlowTrack->Branch("EFlowTrack_PhiOuter",&EFlowTrack_PhiOuter);
  treeEFlowTrack->Branch("EFlowTrack_T",&EFlowTrack_T);
  treeEFlowTrack->Branch("EFlowTrack_X",&EFlowTrack_X);
  treeEFlowTrack->Branch("EFlowTrack_Y",&EFlowTrack_Y);
  treeEFlowTrack->Branch("EFlowTrack_Z",&EFlowTrack_Z);
  treeEFlowTrack->Branch("EFlowTrack_TOuter",&EFlowTrack_TOuter);
  treeEFlowTrack->Branch("EFlowTrack_XOuter",&EFlowTrack_XOuter);
  treeEFlowTrack->Branch("EFlowTrack_YOuter",&EFlowTrack_YOuter);
  treeEFlowTrack->Branch("EFlowTrack_ZOuter",&EFlowTrack_ZOuter);
  treeEFlowTrack->Branch("EFlowTrack_Xd",&EFlowTrack_Xd);
  treeEFlowTrack->Branch("EFlowTrack_Yd",&EFlowTrack_Yd);
  treeEFlowTrack->Branch("EFlowTrack_Zd",&EFlowTrack_Zd);
  treeEFlowTrack->Branch("EFlowTrack_L",&EFlowTrack_L);
  treeEFlowTrack->Branch("EFlowTrack_D0",&EFlowTrack_D0);
  treeEFlowTrack->Branch("EFlowTrack_DZ",&EFlowTrack_DZ);
  treeEFlowTrack->Branch("EFlowTrack_Nclusters",&EFlowTrack_Nclusters);
  treeEFlowTrack->Branch("EFlowTrack_dNdx",&EFlowTrack_dNdx);
  treeEFlowTrack->Branch("EFlowTrack_ErrorP",&EFlowTrack_ErrorP);
  treeEFlowTrack->Branch("EFlowTrack_ErrorPT",&EFlowTrack_ErrorPT);
  treeEFlowTrack->Branch("EFlowTrack_ErrorPhi",&EFlowTrack_ErrorPhi);
  treeEFlowTrack->Branch("EFlowTrack_ErrorCtgTheta",&EFlowTrack_ErrorCtgTheta);
  treeEFlowTrack->Branch("EFlowTrack_ErrorT",&EFlowTrack_ErrorT);
  treeEFlowTrack->Branch("EFlowTrack_ErrorD0",&EFlowTrack_ErrorD0);
  treeEFlowTrack->Branch("EFlowTrack_ErrorDZ",&EFlowTrack_ErrorDZ);
  treeEFlowTrack->Branch("EFlowTrack_ErrorC",&EFlowTrack_ErrorC);
  treeEFlowTrack->Branch("EFlowTrack_ErrorD0Phi",&EFlowTrack_ErrorD0Phi);
  treeEFlowTrack->Branch("EFlowTrack_ErrorD0C",&EFlowTrack_ErrorD0C);
  treeEFlowTrack->Branch("EFlowTrack_ErrorD0DZ",&EFlowTrack_ErrorD0DZ);
  treeEFlowTrack->Branch("EFlowTrack_ErrorD0CtgTheta",&EFlowTrack_ErrorD0CtgTheta);
  treeEFlowTrack->Branch("EFlowTrack_ErrorPhiC",&EFlowTrack_ErrorPhiC);
  treeEFlowTrack->Branch("EFlowTrack_ErrorPhiDZ",&EFlowTrack_ErrorPhiDZ);
  treeEFlowTrack->Branch("EFlowTrack_ErrorPhiCtgTheta",&EFlowTrack_ErrorPhiCtgTheta);
  treeEFlowTrack->Branch("EFlowTrack_ErrorCDZ",&EFlowTrack_ErrorCDZ);
  treeEFlowTrack->Branch("EFlowTrack_ErrorCCtgTheta",&EFlowTrack_ErrorCCtgTheta);
  treeEFlowTrack->Branch("EFlowTrack_ErrorDZCtgTheta",&EFlowTrack_ErrorDZCtgTheta);
  treeEFlowTrack->Branch("EFlowTrack_Particle",&EFlowTrack_Particle);
  treeEFlowTrack->Branch("EFlowTrack_VertexIndex",&EFlowTrack_VertexIndex);
  treeEFlowTrack->Branch("EFlowTrack_P4",&EFlowTrack_P4);
  treeEFlowTrack->Branch("EFlowTrack_CovarianceMatrix",&EFlowTrack_CovarianceMatrix);

  // EFlowPhoton
  TTree *treeEFlowPhoton = new TTree("treeEFlowPhoton","treeEFlowPhoton");

  float EFlowPhoton_PT,EFlowPhoton_Eta,EFlowPhoton_Phi,EFlowPhoton_E,EFlowPhoton_T;
  float EFlowPhoton_EhadOverEem,EFlowPhoton_IsolationVar,EFlowPhoton_IsolationVarRhoCorr;
  float EFlowPhoton_SumPtCharged,EFlowPhoton_SumPtNeutral,EFlowPhoton_SumPtChargedPU,EFlowPhoton_SumPt;
  int EFlowPhoton_Status;
  TLorentzVector EFlowPhoton_P4;
  TRefArray EFlowPhoton_Particles;

  treeEFlowPhoton->Branch("EFlowPhoton_PT",&EFlowPhoton_PT);
  treeEFlowPhoton->Branch("EFlowPhoton_Eta",&EFlowPhoton_Eta);
  treeEFlowPhoton->Branch("EFlowPhoton_Phi",&EFlowPhoton_Phi);
  treeEFlowPhoton->Branch("EFlowPhoton_E",&EFlowPhoton_E);
  treeEFlowPhoton->Branch("EFlowPhoton_T",&EFlowPhoton_T);
  treeEFlowPhoton->Branch("EFlowPhoton_EhadOverEem",&EFlowPhoton_EhadOverEem);
  treeEFlowPhoton->Branch("EFlowPhoton_Particles",&EFlowPhoton_Particles);
  treeEFlowPhoton->Branch("EFlowPhoton_IsolationVar",&EFlowPhoton_IsolationVar);
  treeEFlowPhoton->Branch("EFlowPhoton_IsolationVarRhoCorr",&EFlowPhoton_IsolationVarRhoCorr);
  treeEFlowPhoton->Branch("EFlowPhoton_SumPtCharged",&EFlowPhoton_SumPtCharged);
  treeEFlowPhoton->Branch("EFlowPhoton_SumPtNeutral",&EFlowPhoton_SumPtNeutral);
  treeEFlowPhoton->Branch("EFlowPhoton_SumPtChargedPU",&EFlowPhoton_SumPtChargedPU);
  treeEFlowPhoton->Branch("EFlowPhoton_SumPt",&EFlowPhoton_SumPt);
  treeEFlowPhoton->Branch("EFlowPhoton_Status",&EFlowPhoton_Status);
  treeEFlowPhoton->Branch("EFlowPhoton_P4",&EFlowPhoton_P4);
  
  // Photon
  TTree *treePhoton = new TTree("treePhoton","treePhoton");

  float Photon_PT,Photon_Eta,Photon_Phi,Photon_E,Photon_T;
  float Photon_EhadOverEem,Photon_IsolationVar,Photon_IsolationVarRhoCorr;
  float Photon_SumPtCharged,Photon_SumPtNeutral,Photon_SumPtChargedPU,Photon_SumPt;
  int Photon_Status;
  TLorentzVector Photon_P4;
  TRefArray Photon_Particles;

  treePhoton->Branch("Photon_PT",&Photon_PT);
  treePhoton->Branch("Photon_Eta",&Photon_Eta);
  treePhoton->Branch("Photon_Phi",&Photon_Phi);
  treePhoton->Branch("Photon_E",&Photon_E);
  treePhoton->Branch("Photon_T",&Photon_T);
  treePhoton->Branch("Photon_EhadOverEem",&Photon_EhadOverEem);
  treePhoton->Branch("Photon_Particles",&Photon_Particles);
  treePhoton->Branch("Photon_IsolationVar",&Photon_IsolationVar);
  treePhoton->Branch("Photon_IsolationVarRhoCorr",&Photon_IsolationVarRhoCorr);
  treePhoton->Branch("Photon_SumPtCharged",&Photon_SumPtCharged);
  treePhoton->Branch("Photon_SumPtNeutral",&Photon_SumPtNeutral);
  treePhoton->Branch("Photon_SumPtChargedPU",&Photon_SumPtChargedPU);
  treePhoton->Branch("Photon_SumPt",&Photon_SumPt);
  treePhoton->Branch("Photon_Status",&Photon_Status);
  treePhoton->Branch("Photon_P4",&Photon_P4);

  //Electron
  TTree *treeElectron = new TTree("treeElectron","treeElectron");

  float Electron_PT,Electron_Eta,Electron_Phi,Electron_T,Electron_EhadOverEem;
  float Electron_IsolationVar,Electron_IsolationVarRhoCorr,Electron_SumPtCharged;
  float Electron_SumPtNeutral,Electron_SumPtChargedPU,Electron_SumPt;
  float Electron_D0,Electron_DZ,Electron_ErrorD0,Electron_ErrorDZ;
  int Electron_Charge;
  TLorentzVector Electron_P4;
  TRef Electron_Particle;

  treeElectron->Branch("Electron_PT",&Electron_PT);
  treeElectron->Branch("Electron_Eta",&Electron_Eta);
  treeElectron->Branch("Electron_Phi",&Electron_Phi);
  treeElectron->Branch("Electron_T",&Electron_T);
  treeElectron->Branch("Electron_Charge",&Electron_Charge);
  treeElectron->Branch("Electron_EhadOverEem",&Electron_EhadOverEem);
  treeElectron->Branch("Electron_Particle",&Electron_Particle);
  treeElectron->Branch("Electron_IsolationVar",&Electron_IsolationVar);
  treeElectron->Branch("Electron_IsolationVarRhoCorr",&Electron_IsolationVarRhoCorr);
  treeElectron->Branch("Electron_SumPtCharged",&Electron_SumPtCharged);
  treeElectron->Branch("Electron_SumPtNeutral",&Electron_SumPtNeutral);
  treeElectron->Branch("Electron_SumPtChargedPU",&Electron_SumPtChargedPU);
  treeElectron->Branch("Electron_SumPt",&Electron_SumPt);
  treeElectron->Branch("Electron_D0",&Electron_D0);
  treeElectron->Branch("Electron_DZ",&Electron_DZ);
  treeElectron->Branch("Electron_ErrorD0",&Electron_ErrorD0);
  treeElectron->Branch("Electron_ErrorDZ",&Electron_ErrorDZ);
  treeElectron->Branch("Electron_P4",&Electron_P4);

  // Jet
  TTree *treeJet = new TTree("treeJet","treeJet");

  float Jet_PT,Jet_Eta,Jet_Phi,Jet_T,Jet_Mass,Jet_DeltaEta,Jet_DeltaPhi;
  int Jet_Flavor,Jet_FlavorAlgo,Jet_FlavorPhys,Jet_BTag,Jet_BTagAlgo,Jet_BTagPhys,Jet_TauTag,Jet_Charge,Jet_NCharged,Jet_NNeutrals,Jet_NSubJetsTrimmed,Jet_NSubJetsPruned,Jet_NSubJetsSoftDropped;
  double Jet_ExclYmerge23,Jet_ExclYmerge34,Jet_ExclYmerge45,Jet_ExclYmerge56;
  float Jet_TauWeight,Jet_EhadOverEem,Jet_NeutralEnergyFraction,Jet_ChargedEnergyFraction,Jet_Beta,Jet_BetaStar,Jet_MeanSqDeltaR,Jet_PTD;
  float Jet_FracPt[5],Jet_Tau[5];
  TLorentzVector Jet_SoftDroppedJet,Jet_SoftDroppedSubJet1,Jet_SoftDroppedSubJet2;
  TLorentzVector Jet_TrimmedP4[5],Jet_PrunedP4[5],Jet_SoftDroppedP4[5];
  TRefArray Jet_Constituents,Jet_Particles;
  TLorentzVector Jet_P4,Jet_Area;

  treeJet->Branch("GenJet_PT",&Jet_PT);
  treeJet->Branch("Jet_Eta",&Jet_Eta);
  treeJet->Branch("Jet_Phi",&Jet_Phi);
  treeJet->Branch("Jet_T",&Jet_T);
  treeJet->Branch("Jet_Mass",&Jet_Mass);
  treeJet->Branch("Jet_DeltaEta",&Jet_DeltaEta);
  treeJet->Branch("Jet_DeltaPhi",&Jet_DeltaPhi);
  treeJet->Branch("Jet_Flavor",&Jet_Flavor);
  treeJet->Branch("Jet_FlavorAlgo",&Jet_FlavorAlgo);
  treeJet->Branch("Jet_FlavorPhys",&Jet_FlavorPhys);
  treeJet->Branch("Jet_BTag",&Jet_BTag);
  treeJet->Branch("Jet_BTagAlgo",&Jet_BTagAlgo);
  treeJet->Branch("Jet_BTagPhys",&Jet_BTagPhys);
  treeJet->Branch("Jet_TauTag",&Jet_TauTag);
  treeJet->Branch("Jet_TauWeight",&Jet_TauWeight);
  treeJet->Branch("Jet_Charge",&Jet_Charge);
  treeJet->Branch("Jet_EhadOverEem",&Jet_EhadOverEem);
  treeJet->Branch("Jet_NCharged",&Jet_NCharged);
  treeJet->Branch("Jet_NNeutrals",&Jet_NNeutrals);
  treeJet->Branch("Jet_NeutralEnergyFraction",&Jet_NeutralEnergyFraction);
  treeJet->Branch("Jet_ChargedEnergyFraction",&Jet_ChargedEnergyFraction);
  treeJet->Branch("GenJet_Beta",&Jet_Beta);
  treeJet->Branch("Jet_BetaStar",&Jet_BetaStar);
  treeJet->Branch("Jet_MeanSqDeltaR",&Jet_MeanSqDeltaR);
  treeJet->Branch("Jet_PTD",&Jet_PTD);
  treeJet->Branch("Jet_FracPt",Jet_FracPt);
  treeJet->Branch("Jet_Tau",Jet_Tau);
  treeJet->Branch("Jet_SoftDroppedJet",&Jet_SoftDroppedJet);
  treeJet->Branch("Jet_SoftDroppedSubJet1",&Jet_SoftDroppedSubJet1);
  treeJet->Branch("Jet_SoftDroppedSubJet2",&Jet_SoftDroppedSubJet2);
  treeJet->Branch("Jet_TrimmedP4",Jet_TrimmedP4);
  treeJet->Branch("Jet_PrunedP4",Jet_PrunedP4);
  treeJet->Branch("Jet_SoftDroppedP4",Jet_SoftDroppedP4);
  treeJet->Branch("Jet_NSubJetsTrimmed",&Jet_NSubJetsTrimmed);
  treeJet->Branch("Jet_NSubJetsPruned",&Jet_NSubJetsPruned);
  treeJet->Branch("Jet_NSubJetsSoftDropped",&Jet_NSubJetsSoftDropped);
  treeJet->Branch("Jet_ExclYmerge23",&Jet_ExclYmerge23);
  treeJet->Branch("Jet_ExclYmerge34",&Jet_ExclYmerge34);
  treeJet->Branch("Jet_ExclYmerge45",&Jet_ExclYmerge45);
  treeJet->Branch("Jet_ExclYmerge56",&Jet_ExclYmerge56);
  treeJet->Branch("Jet_Constituents",&Jet_Constituents);
  treeJet->Branch("Jet_Particles",&Jet_Particles);
  treeJet->Branch("Jet_P4",&Jet_P4);
  treeJet->Branch("Jet_Area",&Jet_Area);


  // MissingET
  TTree *treeMissingET = new TTree("treeMissingET","treeMissingET");

  float MissingET_MET,MissingET_Eta,MissingET_Phi;
  TLorentzVector MissingET_P4;

  treeMissingET->Branch("MissingET_MET",&MissingET_MET);
  treeMissingET->Branch("MissingET_Eta",&MissingET_Eta);
  treeMissingET->Branch("MissingET_Phi",&MissingET_Phi);
  treeMissingET->Branch("MissingET_P4",&MissingET_P4);

  // ScalarHT
  TTree *treeScalarHT = new TTree("treeScalarHT","treeScalarHT");

  float ScalarHT_HT;

  treeScalarHT->Branch("ScalarHT_HT",&ScalarHT_HT);

  //  Muon
  TTree *treeMuon = new TTree("treeMuon","treeMuon");

  int Muon_Charge;
  TRef Muon_Particle;
  float Muon_PT,Muon_Eta,Muon_Phi,Muon_T,Muon_IsolationVar,Muon_IsolationVarRhoCorr,Muon_SumPtCharged,Muon_SumPtNeutral,Muon_SumPtChargedPU,Muon_SumPt,Muon_D0,Muon_DZ,Muon_ErrorD0,Muon_ErrorDZ;
  TLorentzVector Muon_P4;

  treeMuon->Branch("Muon_Charge",&Muon_Charge);
  treeMuon->Branch("Muon_Particle",&Muon_Particle);
  treeMuon->Branch("Muon_PT",&Muon_PT);
  treeMuon->Branch("Muon_Eta",&Muon_Eta);
  treeMuon->Branch("Muon_Phi",&Muon_Phi);
  treeMuon->Branch("Muon_T",&Muon_T);
  treeMuon->Branch("Muon_IsolationVar",&Muon_IsolationVar);
  treeMuon->Branch("Muon_IsolationVarRhoCorr",&Muon_IsolationVarRhoCorr);
  treeMuon->Branch("Muon_SumPtCharged",&Muon_SumPtCharged);
  treeMuon->Branch("Muon_SumPtNeutral",&Muon_SumPtNeutral);
  treeMuon->Branch("Muon_SumPtChargedPU",&Muon_SumPtChargedPU);
  treeMuon->Branch("Muon_SumPt",&Muon_SumPt);
  treeMuon->Branch("Muon_D0",&Muon_D0);
  treeMuon->Branch("Muon_DZ",&Muon_DZ);
  treeMuon->Branch("Muon_ErrorD0",&Muon_ErrorD0);
  treeMuon->Branch("Muon_ErrorDZ",&Muon_ErrorD0);
  treeMuon->Branch("Muon_P4",&Muon_P4);

  while (fp.getline(buf,200)) {


      std::cout<<" inputfile=\t"<<buf<<std::endl;

      gSystem->Load("libDelphes");

      // Create chain of root trees
      TChain chain("Delphes");
      chain.Add(buf);

      // Create object of class ExRootTreeReader
      ExRootTreeReader *treeReader = new ExRootTreeReader(&chain);
      long int numberOfEntries = treeReader->GetEntries();

      // Get pointers to branches used in this analysis
      TClonesArray *branchMuon = treeReader->UseBranch("Muon");
      TClonesArray *branchJet = treeReader->UseBranch("Jet");
      TClonesArray *branchEvent = treeReader->UseBranch("Event");
      TClonesArray *branchParticle = treeReader->UseBranch("Particle");
      TClonesArray *branchGenJet = treeReader->UseBranch("GenJet");
      TClonesArray *branchGenMissingET = treeReader->UseBranch("GenMissingET");
      TClonesArray *branchTrack = treeReader->UseBranch("Track");
      TClonesArray *branchTower = treeReader->UseBranch("Tower");
      TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");
      TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");
      //      TClonesArray *branchEFlowNeutralHadron = treeReader->UseBranch("EFlowNeutralHadron");
      TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
      TClonesArray *branchElectron = treeReader->UseBranch("Electron");
      TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
      TClonesArray *branchScalarHT = treeReader->UseBranch("ScalarHT");

      Event *event;
      GenParticle *particle;
      Jet *genjet;
      MissingET *genMissingET;
      Track *track;
      Tower *tower;
      Track *eflowTrack;
      Photon * eflowPhoton;
      Photon *photon;
      Electron *electron;
      MissingET *missingET;
      ScalarHT * scalarHT;
      Muon *muon, *muon0, *muon1, *muon2, *muon3, *muon_mass1, *muon_mass2;
      Jet *jet,*jet0,*jet1,*jet2,*jet3;

      float Muon_mass,BJet_mass;

      // Loop over all events
      for(int entry = 0; entry < numberOfEntries; ++entry)
      {

          // Load selected branches with data from specified event
          treeReader->ReadEntry(entry);

           //cout << "  event = " << branchEvent->GetEntries()
           //     << "  particle = " << branchParticle->GetEntries()
           //     << "  GenJet = " << branchGenJet->GetEntries()
           //     << "  GenMissingET = " << branchGenMissingET->GetEntries()
           //     << "  track = " << branchTrack->GetEntries()
           //     << "  Tower  = " << branchTower->GetEntries()
           //     << "  EFlowTrack= " << branchEFlowTrack->GetEntries()
           //     << "  EFlowPhoton= " << branchEFlowPhoton->GetEntries()
           //     << "  Photon= " << branchPhoton->GetEntries()
           //     << "  Electron = " << branchElectron->GetEntries()
           //     << "  MissingET= " << branchMissingET->GetEntries()
           //     << "  ScalarHT = " << branchScalarHT->GetEntries()
           //     << "  Jet =  " << branchJet->GetEntries()
           //     << "  Muon =  " << branchMuon->GetEntries()
           //     << endl;


          // Muon Analyze
          if(branchMuon->GetEntries()>1) {

//              std::cout << "Muon size= " << branchMuon->GetEntries() << std::endl;

              std::vector<double> muon_pt;
              int muon_max[2] = {0};

              for(int i=0;i<branchMuon->GetEntries();i++){

                  muon3 = (Muon*) branchMuon->At(i);

                  muon_pt.push_back(muon3->PT);
//std::cout << "muon pt[" << i << "]= " << muon->PT << std::endl;
              } // for cycle

              int max1 = distance(begin(muon_pt),max_element(muon_pt.begin(),muon_pt.end()));
              muon1 = (Muon*) branchMuon->At(max1);
              muon_max[0] = max1;

              muon_pt.erase(muon_pt.begin()+max1);
              muon_pt.insert(muon_pt.begin()+max1,-999);

              int max2 = distance(begin(muon_pt),max_element(muon_pt.begin(),muon_pt.end()));
              muon2 = (Muon*) branchMuon->At(max2);
              muon_max[1] = max2;

              Muon_mass = (muon1->P4()+muon2->P4()).M();
//cout << "muon mass = " << Muon_mass << endl;
              muon_pt.clear();
          } // muon if
          else continue;

          // Jet Analyze
          if(branchJet->GetEntries()>1) {

              std::cout << "Jet size= " << branchJet->GetEntries() << std::endl;
              std::vector<double> jet_pt;
              int jet_max[2] = {0};
              int bjet_num = 0;

              for(int i=0;i<branchJet->GetEntries();i++){

                  jet = (Jet*) branchJet->At(i);
                  if(jet->BTag==0 || jet->BTag==1) {
                      bjet_num++;
                  }
              }

              if(bjet_num<2) continue;
              else 
              {
                  int Bnum =0;
                  for(int i=0;i<branchJet->GetEntries();i++)
                  {
                      jet = (Jet*) branchJet->At(i);

//                      cout << " jet pt = " << jet->PT << " Bjet = " << jet->BTag << endl;
                      if(jet->BTag==0 || jet->BTag==1) {
                          jetbtag.insert(std::pair<int,jetBtag>(Bnum,jetBtag(i,jet->BTag,jet->PT)));
                          Bnum++;
                      } // if

                          //jet_pt.push_back(jet->PT);

                  }

                  for(int m=0;m<jetbtag.size();m++)
                  {
                      auto mapit = jetbtag.find(m);
                      if(mapit != jetbtag.end())
                      {
                         jet_pt.push_back(mapit->second.Jet_PT);
                      }
                  }

                  int jet_max1 = distance(begin(jet_pt),max_element(jet_pt.begin(),jet_pt.end()));
                  auto mapit1 = jetbtag.find(jet_max1);
                  if(mapit1 != jetbtag.end())
                  {
                      jet1 = (Jet*) branchJet->At(mapit1->second.Num);
                  }

                  jet_pt.erase(jet_pt.begin()+jet_max1);
                  jet_pt.insert(jet_pt.begin(),-999);

                  int jet_max2 = distance(begin(jet_pt),max_element(jet_pt.begin(),jet_pt.end()));
                  auto mapit2 = jetbtag.find(jet_max2);
                  if(mapit2 != jetbtag.end())
                  {
                      jet2 = (Jet*) branchJet->At(mapit2->second.Num);
                  }

                  BJet_mass = (jet1->P4()+jet2->P4()).M();
cout << " jet mass = " << BJet_mass << endl;
                  jetbtag.clear();
                  jet_pt.clear();
              } // jet else

          } //jet if
          else continue;


          // cut
          if((Muon_mass<105&&Muon_mass>75)&&(BJet_mass<150&&BJet_mass>100))
          {
              cout << " muon mass = " << Muon_mass
                  << " jet mass = " << BJet_mass
                  << endl;

              // Event Fill
              cout << " Event size = " << branchEvent->GetEntries() << endl;

              for(int i=0;i<branchEvent->GetEntries();i++)
              {
                  event = (Event *) branchEvent->At(i);

                  Event_Number = event->Number;
                  Event_ReadTime = event->ReadTime;
                  Event_ProcTime = event->ProcTime;

                  treeEvent->Fill();

              } // event fill

              // GenMissingET Fill
              cout << " GenMissingET size = " << branchGenMissingET->GetEntries() << endl;

              for(int i=0;i<branchGenMissingET->GetEntries();i++)
              {
                  genMissingET = (MissingET *) branchGenMissingET->At(i);

                  GenMissingET_MET = genMissingET->MET;
                  GenMissingET_Eta = genMissingET->Eta;
                  GenMissingET_Phi = genMissingET->Phi;
                  GenMissingET_P4.SetPtEtaPhiM(genMissingET->P4().Pt(),genMissingET->P4().Eta(),genMissingET->P4().Phi(),genMissingET->P4().M());

                  treeGenMissingET->Fill();

              } // GenMissingET Fill

              // MissingET Fill
              cout << " MissingET size = " << branchMissingET->GetEntries() << endl;

              for(int i=0;i<branchMissingET->GetEntries();i++)
              {
                  missingET = (MissingET *) branchMissingET->At(i);

                  MissingET_MET = missingET->MET;
                  MissingET_Eta = missingET->Eta;
                  MissingET_Phi = missingET->Phi;
                  MissingET_P4.SetPtEtaPhiM(missingET->P4().Pt(),missingET->P4().Eta(),missingET->P4().Phi(),missingET->P4().M());

                  treeMissingET->Fill();
              } // MissingET Fill

              // ScalarHT Fill
              cout << " ScalarHT size = " << branchScalarHT->GetEntries() << endl;
              for(int i=0;i<branchScalarHT->GetEntries();i++)
              {
                  scalarHT = (ScalarHT *) branchScalarHT->At(i);

                  ScalarHT_HT = scalarHT->HT;

                  treeScalarHT->Fill();

              } // ScalarHT Fill

              // Particle Fill
              cout << " Particle size  = " << branchParticle->GetEntries() << endl;

              for(int i=0;i<branchParticle->GetEntries();i++) {
                  particle = (GenParticle *) branchParticle->At(i);

                  Particle_PID = particle->PID;
                  Particle_Status = particle->Status;
                  Particle_IsPU = particle->IsPU;
                  Particle_M1 = particle->M1;
                  Particle_M2 = particle->M2;
                  Particle_D1 = particle->D1;
                  Particle_D2 = particle->D2;
                  Particle_Charge = particle->Charge;
                  Particle_Mass = particle->Mass;
                  Particle_E = particle->E;
                  Particle_Px = particle->Px;
                  Particle_Py = particle->Py;
                  Particle_Pz = particle->Pz;
                  Particle_P  = particle->P;
                  Particle_PT = particle->PT;
                  Particle_Eta = particle->Eta;
                  Particle_Phi = particle->Phi;
                  Particle_Rapidity = particle->Rapidity;
                  Particle_T = particle->T;
                  Particle_X = particle->X;
                  Particle_Y = particle->Y;
                  Particle_Z = particle->Z;
                  Particle_P4.SetPxPyPzE(particle->P4().Px(),particle->P4().Py(),particle->P4().Pz(),particle->P4().E());
                  treeParticle->Fill();


              } // for particle

              // GenJet Fill
              cout << " GetJet size = " << branchGenJet->GetEntries() << endl;

              for(int i =0;i<branchGenJet->GetEntries();i++)
              {

                  genjet = (Jet *) branchGenJet->At(i);

                  GenJet_PT          =  genjet->PT;
                  GenJet_Eta         =  genjet->Eta;
                  GenJet_Phi         =  genjet->Phi;
                  GenJet_T           =  genjet->T;
                  GenJet_Mass        =  genjet->Mass;
                  GenJet_DeltaEta    =  genjet->DeltaEta;
                  GenJet_DeltaPhi    =  genjet->DeltaPhi;
                  GenJet_Flavor      =  genjet->Flavor;
                  GenJet_FlavorAlgo  =  genjet->FlavorAlgo;
                  GenJet_FlavorPhys  =  genjet->FlavorPhys;
                  GenJet_BTag        =  genjet->BTag;
                  GenJet_BTagAlgo    =  genjet->BTagAlgo;
                  GenJet_BTagPhys    =  genjet->BTagPhys;
                  GenJet_TauTag      =  genjet->TauTag;
                  GenJet_TauWeight   =  genjet->TauWeight;
                  GenJet_Charge      =  genjet->Charge;
                  GenJet_EhadOverEem =  genjet->EhadOverEem;
                  GenJet_NCharged    =  genjet->NCharged;
                  GenJet_NNeutrals   =  genjet->NNeutrals;
                  GenJet_NeutralEnergyFraction = genjet->NeutralEnergyFraction;
                  GenJet_ChargedEnergyFraction = genjet->ChargedEnergyFraction;
                  GenJet_Beta        =  genjet->Beta;
                  GenJet_BetaStar    =  genjet->BetaStar;
                  GenJet_MeanSqDeltaR=  genjet->MeanSqDeltaR;
                  GenJet_PTD         =  genjet->PTD;
                  GenJet_SoftDroppedJet.SetXYZT(genjet->SoftDroppedJet.X(),genjet->SoftDroppedJet.Y(),genjet->SoftDroppedJet.Z(),genjet->SoftDroppedJet.T());
                  GenJet_SoftDroppedSubJet1.SetXYZT(genjet->SoftDroppedSubJet1.X(),genjet->SoftDroppedSubJet1.Y(),genjet->SoftDroppedSubJet1.Z(),genjet->SoftDroppedSubJet1.T());
                  GenJet_SoftDroppedSubJet2.SetXYZT(genjet->SoftDroppedSubJet2.X(),genjet->SoftDroppedSubJet2.Y(),genjet->SoftDroppedSubJet2.Z(),genjet->SoftDroppedSubJet2.T());
                  GenJet_NSubJetsTrimmed = genjet->NSubJetsTrimmed;
                  GenJet_NSubJetsPruned  = genjet->NSubJetsPruned;
                  GenJet_NSubJetsSoftDropped = genjet->NSubJetsSoftDropped;
                  GenJet_ExclYmerge23        = genjet->ExclYmerge23;
                  GenJet_ExclYmerge34        = genjet->ExclYmerge34;
                  GenJet_ExclYmerge45        = genjet->ExclYmerge45;
                  GenJet_ExclYmerge56        = genjet->ExclYmerge56;
                  GenJet_P4.SetPtEtaPhiM(genjet->P4().Pt(),genjet->P4().Eta(),genjet->P4().Phi(),genjet->P4().M());
                  GenJet_Area.SetXYZT(genjet->Area.X(),genjet->Area.Z(),genjet->Area.Z(),genjet->Area.T());
                  for(int n=0;n<(genjet->Particles).GetEntriesFast();n++)
                  {
                      GenJet_Particles.Add(genjet->Particles.At(n));
                  }
                  for(int n=0;n<(genjet->Constituents).GetEntriesFast();n++)
                  {
                    GenJet_Constituents.Add(genjet->Constituents.At(n));
                  }


                  // for jet ntuple
                  for(int n =0;n<5;n++){
                      GenJet_FracPt[n]   =  genjet->FracPt[n];
                      GenJet_Tau[n]      =  genjet->Tau[n];
                      GenJet_TrimmedP4[n].SetXYZT(genjet->TrimmedP4[n].X(),genjet->TrimmedP4[n].Y(),genjet->TrimmedP4[n].Z(),genjet->TrimmedP4[n].T());
                      GenJet_PrunedP4[n].SetXYZT(genjet->PrunedP4[n].X(),genjet->PrunedP4[n].Y(),genjet->PrunedP4[n].Z(),genjet->PrunedP4[n].T());
                      GenJet_SoftDroppedP4[i+n].SetXYZT(genjet->SoftDroppedP4[n].X(),genjet->SoftDroppedP4[n].Y(),genjet->SoftDroppedP4[n].Z(),genjet->SoftDroppedP4[n].T());
                  } // for jet ntuple

                  treeGenJet->Fill();
                  GenJet_Constituents.Clear();
                  GenJet_Particles.Clear();

              } // for genjet

              // Jet Fill
              cout << "Jet size = " << branchJet->GetEntries() << endl;

              std::vector<double> jet_pt2;
              int jet_max[2] = {0};
              int Bnum2 =0;
              for(int j=0;j<branchJet->GetEntries();j++)
              {
                  jet0 = (Jet*) branchJet->At(j);
                  if(jet0->BTag==0 || jet0->BTag==1) {
                     jetbtag2.insert(std::pair<int,jetBtag2>(Bnum2,jetBtag2(j,jet0->BTag,jet0->PT)));
                     Bnum2++;
                     }
              }
              for(int l=0;l<jetbtag2.size();l++)
              {
                  auto mapit0 = jetbtag2.find(l);
                  if(mapit0 != jetbtag2.end())
                  {
                      jet_pt2.push_back(mapit0->second.Jet_PT);
                  }
              }

              int Jet_max1 = distance(begin(jet_pt2),max_element(jet_pt2.begin(),jet_pt2.end()));
              auto mapJet1 = jetbtag2.find(Jet_max1);
              if(mapJet1 != jetbtag2.end())
              {
                  jet_max[0]  = mapJet1->second.Num;
              }

              jet_pt2.erase(jet_pt2.begin()+Jet_max1);
              jet_pt2.insert(jet_pt2.begin(),-999);

              int Jet_max2 = distance(begin(jet_pt2),max_element(jet_pt2.begin(),jet_pt2.end()));
              auto mapJet2 = jetbtag2.find(Jet_max2);
              if(mapJet2 != jetbtag2.end())
              {
                  jet_max[1]  = mapJet2->second.Num;
              }


              for(int s =0;s<2;s++)//branchJet->GetEntries();i++)
              {

                  jet = (Jet *) branchJet->At(jet_max[s]);

                  Jet_PT          =  jet->PT;
                  Jet_Eta         =  jet->Eta;
                  Jet_Phi         =  jet->Phi;
                  Jet_T           =  jet->T;
                  Jet_Mass        =  jet->Mass;
                  Jet_DeltaEta    =  jet->DeltaEta;
                  Jet_DeltaPhi    =  jet->DeltaPhi;
                  Jet_Flavor      =  jet->Flavor;
                  Jet_FlavorAlgo  =  jet->FlavorAlgo;
                  Jet_FlavorPhys  =  jet->FlavorPhys;
                  Jet_BTag        =  jet->BTag;
                  Jet_BTagAlgo    =  jet->BTagAlgo;
                  Jet_BTagPhys    =  jet->BTagPhys;
                  Jet_TauTag      =  jet->TauTag;
                  Jet_TauWeight   =  jet->TauWeight;
                  Jet_Charge      =  jet->Charge;
                  Jet_EhadOverEem =  jet->EhadOverEem;
                  Jet_NCharged    =  jet->NCharged;
                  Jet_NNeutrals   =  jet->NNeutrals;
                  Jet_NeutralEnergyFraction = jet->NeutralEnergyFraction;
                  Jet_ChargedEnergyFraction = jet->ChargedEnergyFraction;
                  Jet_Beta        =  jet->Beta;
                  Jet_BetaStar    =  jet->BetaStar;
                  Jet_MeanSqDeltaR=  jet->MeanSqDeltaR;
                  Jet_PTD         =  jet->PTD;
                  Jet_SoftDroppedJet.SetXYZT(jet->SoftDroppedJet.X(),jet->SoftDroppedJet.Y(),jet->SoftDroppedJet.Z(),jet->SoftDroppedJet.T());
                  Jet_SoftDroppedSubJet1.SetXYZT(jet->SoftDroppedSubJet1.X(),jet->SoftDroppedSubJet1.Y(),jet->SoftDroppedSubJet1.Z(),jet->SoftDroppedSubJet1.T());
                  Jet_SoftDroppedSubJet2.SetXYZT(jet->SoftDroppedSubJet2.X(),jet->SoftDroppedSubJet2.Y(),jet->SoftDroppedSubJet2.Z(),jet->SoftDroppedSubJet2.T());
                  Jet_NSubJetsTrimmed = jet->NSubJetsTrimmed;
                  Jet_NSubJetsPruned  = jet->NSubJetsPruned;
                  Jet_NSubJetsSoftDropped = jet->NSubJetsSoftDropped;
                  Jet_ExclYmerge23        = jet->ExclYmerge23;
                  Jet_ExclYmerge34        = jet->ExclYmerge34;
                  Jet_ExclYmerge45        = jet->ExclYmerge45;
                  Jet_ExclYmerge56        = jet->ExclYmerge56;
                  Jet_P4.SetPtEtaPhiM(jet->P4().Pt(),jet->P4().Eta(),jet->P4().Phi(),jet->P4().M());
                  Jet_Area.SetXYZT(jet->Area.X(),jet->Area.Z(),jet->Area.Z(),jet->Area.T());
                  for(int n=0;n<(jet->Particles).GetEntriesFast();n++)
                  {
                      Jet_Particles.Add(jet->Particles.At(n));
                  }
                  for(int n=0;n<(jet->Constituents).GetEntriesFast();n++)
                  {
                      Jet_Constituents.Add(jet->Constituents.At(n));
                  }


                  // for jet ntuple
                  for(int n =0;n<5;n++){
                      Jet_FracPt[n]   =  jet->FracPt[n];
                      Jet_Tau[n]      =  jet->Tau[n];
                      Jet_TrimmedP4[n].SetXYZT(jet->TrimmedP4[n].X(),jet->TrimmedP4[n].Y(),jet->TrimmedP4[n].Z(),jet->TrimmedP4[n].T());
                      Jet_PrunedP4[n].SetXYZT(jet->PrunedP4[n].X(),jet->PrunedP4[n].Y(),jet->PrunedP4[n].Z(),jet->PrunedP4[n].T());
                      Jet_SoftDroppedP4[n].SetXYZT(jet->SoftDroppedP4[n].X(),jet->SoftDroppedP4[n].Y(),jet->SoftDroppedP4[n].Z(),jet->SoftDroppedP4[n].T());
                  } // for jet ntuple

                  treeJet->Fill();
                  Jet_Particles.Clear();
                  Jet_Constituents.Clear();

              } // for jet

              jet_pt2.clear();
              jetbtag2.clear();

              // Muon Fill
              cout << "muon size = " << branchMuon->GetEntries() << endl;


              std::vector<double> Muon_pt;
              int Muon_max[2] = {0};

              for(int i=0;i<branchMuon->GetEntries();i++){

                  muon0 = (Muon*) branchMuon->At(i);

                  Muon_pt.push_back(muon0->PT);
                  //std::cout << "muon pt[" << i << "]= " << muon->PT << std::endl;
              } // for cycle

              int Max1 = distance(begin(Muon_pt),max_element(Muon_pt.begin(),Muon_pt.end()));
              Muon_max[0] = Max1;

              Muon_pt.erase(Muon_pt.begin()+Max1);
              Muon_pt.insert(Muon_pt.begin()+Max1,-999);

              int Max2 = distance(begin(Muon_pt),max_element(Muon_pt.begin(),Muon_pt.end()));
              Muon_max[1] = Max2;

              muon_mass1 = (Muon *) branchMuon->At(Muon_max[0]);
              muon_mass2 = (Muon *) branchMuon->At(Muon_max[1]);
//cout << " muon1 p = " << muon_mass1->P4().P() << " muon1 e = " << muon_mass1->P4().E() <<  "muon1 M ="  << muon_mass1->P4().M() << endl;
//cout << " muon1 p = " << muon_mass2->P4().P() << " muon2 E = " << muon_mass2->P4().E() <<  "muon2 M ="  << muon_mass2->P4().M() << endl;
              z_mass = (muon_mass1->P4()+muon_mass2->P4()).M();
              Z_mass->Fill();

              for(int m=0;m<2;m++){
                  muon = (Muon *) branchMuon->At(Muon_max[m]);

                  Muon_PT = muon->PT;
                  Muon_Eta = muon->Eta;
                  Muon_Phi = muon->Phi;
                  Muon_T = muon->T;
                  Muon_Charge = muon->Charge;
                  Muon_Particle = muon->Particle;
                  Muon_IsolationVar = muon->IsolationVar;
                  Muon_IsolationVarRhoCorr = muon->IsolationVarRhoCorr;
                  Muon_SumPtCharged = muon->SumPtCharged;
                  Muon_SumPtNeutral = muon->SumPtNeutral;
                  Muon_SumPtChargedPU = muon->SumPtChargedPU;
                  Muon_SumPt = muon->SumPt;
                  Muon_D0 = muon->D0;
                  Muon_DZ = muon->DZ;
                  Muon_ErrorD0 = muon->ErrorD0;
                  Muon_ErrorDZ = muon->ErrorDZ;
                  //Muon_P4.SetPtEtaPhiM(muon->P4().Pt(),muon->P4().Eta(),muon->P4().Phi(),muon->P4().M());
                  Muon_P4.SetPtEtaPhiM(muon->P4().Pt(),muon->P4().Eta(),muon->P4().Phi(),105.0);

                  treeMuon->Fill();

              } // for muon
              Muon_pt.clear();

              // Track Fill
              cout << " Track Size = " << branchTrack->GetEntries() << endl;

              for(int i=0;i<branchTrack->GetEntries();i++)
              {
                  track = (Track *) branchTrack->At(i);

                  Track_PID       =  track->PID;
                  Track_Charge    =  track->Charge;
                  Track_P         =  track->P;
                  Track_PT        =  track->PT;
                  Track_Eta       =  track->Eta;
                  Track_Phi       =  track->Phi;
                  Track_CtgTheta  =  track->CtgTheta;
                  Track_C         =  track->C;
                  Track_Mass      =  track->Mass;
                  Track_EtaOuter  =  track->EtaOuter;
                  Track_PhiOuter  =  track->PhiOuter;
                  Track_T         =  track->T;
                  Track_X         =  track->X;
                  Track_Y         =  track->Y;
                  Track_Z         =  track->Z;
                  Track_TOuter    =  track->TOuter;
                  Track_XOuter    =  track->XOuter;
                  Track_YOuter    =  track->YOuter;
                  Track_ZOuter    =  track->ZOuter;
                  Track_L         =  track->L;
                  Track_Xd        =  track->Xd;
                  Track_Yd        =  track->Yd;
                  Track_Zd        =  track->Zd;
                  Track_D0        =  track->D0;
                  Track_DZ        =  track->DZ;
                  Track_Nclusters =  track->Nclusters;
                  Track_dNdx      =  track->dNdx;
                  Track_ErrorP    =  track->ErrorP;
                  Track_ErrorPT   =  track->ErrorPT;
                  Track_ErrorPhi  =  track->ErrorPhi;
                  Track_ErrorCtgTheta = track->ErrorCtgTheta;
                  Track_ErrorT    =  track->ErrorT;
                  Track_ErrorD0   =  track->ErrorD0;
                  Track_ErrorDZ   =  track->ErrorDZ;
                  Track_ErrorC    =  track->ErrorC;
                  Track_ErrorD0Phi       =  track->ErrorD0Phi;
                  Track_ErrorD0C         =  track->ErrorD0C;
                  Track_ErrorD0DZ        =  track->ErrorD0DZ;
                  Track_ErrorD0CtgTheta  =  track->ErrorD0CtgTheta;
                  Track_ErrorPhiC        =  track->ErrorPhiC;
                  Track_ErrorPhiDZ       =  track->ErrorPhiDZ;
                  Track_ErrorPhiCtgTheta =  track->ErrorPhiCtgTheta;
                  Track_ErrorCDZ         =  track->ErrorCDZ;
                  Track_ErrorCCtgTheta   =  track->ErrorCCtgTheta;
                  Track_ErrorDZCtgTheta  =  track->ErrorDZCtgTheta;
                  Track_Particle         =  track->Particle;
                  Track_VertexIndex      =  track->VertexIndex;
                  Track_P4.SetPtEtaPhiM(track->P4().Pt(),track->P4().Eta(),track->P4().Phi(),track->P4().M());
                  //const Double_t * mp1 = track->CovarianceMatrix().GetMatrixArray();
                  //int a = track->CovarianceMatrix().GetRowLwb();
                  //int b = track->CovarianceMatrix().GetRowUpb();
                  //int c = track->CovarianceMatrix().GetColLwb();
                  //int d = track->CovarianceMatrix().GetColUpb();
                  //int f = track->CovarianceMatrix().GetNcols();
                  //for(int n=track->CovarianceMatrix().GetRowLwb();n<=track->CovarianceMatrix().GetRowUpb();n++){
                  //    for(int m=track->CovarianceMatrix().GetColLwb();m<=track->CovarianceMatrix().GetColUpb();m++){

                  //         (*Track_CovarianceMatrix)[m][n] = mp1[(n-a)*f+(m-c)];
                  ////mp1[(n-a)*f+(m-c)];
                  //    }
                  //}

                  treeTrack->Fill();

              } //for track Fill

              // EFlowTrack Fill
              cout << " EFlowTrack size = " << branchEFlowTrack->GetEntries() << endl;

              for(int i =0;i<branchEFlowTrack->GetEntries();i++)
              {
                  eflowTrack = (Track *) branchEFlowTrack->At(i);

                  EFlowTrack_PID       =  eflowTrack->PID;
                  EFlowTrack_Charge    =  eflowTrack->Charge;
                  EFlowTrack_P         =  eflowTrack->P;
                  EFlowTrack_PT        =  eflowTrack->PT;
                  EFlowTrack_Eta       =  eflowTrack->Eta;
                  EFlowTrack_Phi       =  eflowTrack->Phi;
                  EFlowTrack_CtgTheta  =  eflowTrack->CtgTheta;
                  EFlowTrack_C         =  eflowTrack->C;
                  EFlowTrack_Mass      =  eflowTrack->Mass;
                  EFlowTrack_EtaOuter  =  eflowTrack->EtaOuter;
                  EFlowTrack_PhiOuter  =  eflowTrack->PhiOuter;
                  EFlowTrack_T         =  eflowTrack->T;
                  EFlowTrack_X         =  eflowTrack->X;
                  EFlowTrack_Y         =  eflowTrack->Y;
                  EFlowTrack_Z         =  eflowTrack->Z;
                  EFlowTrack_TOuter    =  eflowTrack->TOuter;
                  EFlowTrack_XOuter    =  eflowTrack->XOuter;
                  EFlowTrack_YOuter    =  eflowTrack->YOuter;
                  EFlowTrack_ZOuter    =  eflowTrack->ZOuter;
                  EFlowTrack_L         =  eflowTrack->L;
                  EFlowTrack_Xd        =  eflowTrack->Xd;
                  EFlowTrack_Yd        =  eflowTrack->Yd;
                  EFlowTrack_Zd        =  eflowTrack->Zd;
                  EFlowTrack_D0        =  eflowTrack->D0;
                  EFlowTrack_DZ        =  eflowTrack->DZ;
                  EFlowTrack_Nclusters =  eflowTrack->Nclusters;
                  EFlowTrack_dNdx      =  eflowTrack->dNdx;
                  EFlowTrack_ErrorP    =  eflowTrack->ErrorP;
                  EFlowTrack_ErrorPT   =  eflowTrack->ErrorPT;
                  EFlowTrack_ErrorPhi  =  eflowTrack->ErrorPhi;
                  EFlowTrack_ErrorCtgTheta = eflowTrack->ErrorCtgTheta;
                  EFlowTrack_ErrorT    =  eflowTrack->ErrorT;
                  EFlowTrack_ErrorD0   =  eflowTrack->ErrorD0;
                  EFlowTrack_ErrorDZ   =  eflowTrack->ErrorDZ;
                  EFlowTrack_ErrorC    =  eflowTrack->ErrorC;
                  EFlowTrack_ErrorD0Phi       =  eflowTrack->ErrorD0Phi;
                  EFlowTrack_ErrorD0C         =  eflowTrack->ErrorD0C;
                  EFlowTrack_ErrorD0DZ        =  eflowTrack->ErrorD0DZ;
                  EFlowTrack_ErrorD0CtgTheta  =  eflowTrack->ErrorD0CtgTheta;
                  EFlowTrack_ErrorPhiC        =  eflowTrack->ErrorPhiC;
                  EFlowTrack_ErrorPhiDZ       =  eflowTrack->ErrorPhiDZ;
                  EFlowTrack_ErrorPhiCtgTheta =  eflowTrack->ErrorPhiCtgTheta;
                  EFlowTrack_ErrorCDZ         =  eflowTrack->ErrorCDZ;
                  EFlowTrack_ErrorCCtgTheta   =  eflowTrack->ErrorCCtgTheta;
                  EFlowTrack_ErrorDZCtgTheta  =  eflowTrack->ErrorDZCtgTheta;
                  EFlowTrack_Particle         =  eflowTrack->Particle;
                  EFlowTrack_VertexIndex      =  eflowTrack->VertexIndex;
                  EFlowTrack_P4.SetPtEtaPhiM(eflowTrack->P4().Pt(),eflowTrack->P4().Eta(),eflowTrack->P4().Phi(),eflowTrack->P4().M());
                  //EFlowTrack_CovarianceMatrix = eflowTrack->CovarianceMatrix();
                  //const Double_t * mp1 = eflowTrack->CovarianceMatrix().GetMatrixArray();
                  //int a = eflowTrack->CovarianceMatrix().GetRowLwb();
                  //int b = eflowTrack->CovarianceMatrix().GetRowUpb();
                  //int c = eflowTrack->CovarianceMatrix().GetColLwb();
                  //int d = eflowTrack->CovarianceMatrix().GetColUpb();
                  //int f = eflowTrack->CovarianceMatrix().GetNcols();
                  //for(int n=eflowTrack->CovarianceMatrix().GetRowLwb();n<=eflowTrack->CovarianceMatrix().GetRowUpb();n++){
                  //    for(int m=eflowTrack->CovarianceMatrix().GetColLwb();m<=eflowTrack->CovarianceMatrix().GetColUpb();m++){
                  //        EFlowTrack_CovarianceMatrix(m,n) = mp1[(n-a)*f+(m-c)];
                  //    }
                  //}

                  treeEFlowTrack->Fill();

              } // for eflowtrack fill

              // Tower fill
              cout << " tower size= " << branchTower->GetEntries() << endl;

              for(int i=0;i<branchTower->GetEntries();i++)
              {
                  tower = (Tower *) branchTower->At(i);

                  Tower_ET   =  tower->ET;
                  Tower_Eta  =  tower->Eta;
                  Tower_Phi  =  tower->Phi;
                  Tower_E    =  tower->E;
                  Tower_T    =  tower->T;
                  Tower_NTimeHits  =  tower->NTimeHits;
                  Tower_Eem  =  tower->Eem;
                  Tower_Ehad  =  tower->Ehad;
                  for(int n=0;n<(tower->Particles).GetEntriesFast();n++)
                  {
                      Tower_Particles.Add(tower->Particles.At(n));
                  }
                  //Tower_Particles = tower->Particles;
                  Tower_P4.SetPtEtaPhiM(tower->P4().Pt(),tower->P4().Eta(),tower->P4().Phi(),tower->P4().M());
                  for(int n=0;n<4;n++){
                      Tower_Edges[n] = tower->Edges[n];
                  }

                  treeTower->Fill();
                  Tower_Particles.Clear();

              } // for tower fill

              // Photon Fill
              cout << " Photon size = " << branchPhoton->GetEntries() << endl;

              for(int i=0;i<branchPhoton->GetEntries();i++)
              {
                  photon = (Photon *) branchPhoton->At(i);

                  Photon_PT            =  photon->PT;
                  Photon_Eta           =  photon->Eta;
                  Photon_Phi           =  photon->Phi;
                  Photon_E             =  photon->E;
                  Photon_T             =  photon->T;
                  Photon_EhadOverEem   =  photon->EhadOverEem;
                  for(int n=0;n<(photon->Particles).GetEntriesFast();n++)
                  {
                      Photon_Particles.Add(photon->Particles.At(n));
                  }
                  Photon_IsolationVar  =  photon->IsolationVar;
                  Photon_IsolationVarRhoCorr = photon->IsolationVarRhoCorr;
                  Photon_SumPtCharged  =  photon->SumPtCharged;
                  Photon_SumPtNeutral  =  photon->SumPtNeutral;
                  Photon_SumPtChargedPU=  photon->SumPtChargedPU;
                  Photon_SumPt         =  photon->SumPt;
                  Photon_Status        =  photon->Status;
                  Photon_P4.SetPtEtaPhiM(photon->P4().Pt(),photon->P4().Eta(),photon->P4().Phi(),photon->P4().M());

                  treePhoton->Fill();
                  Photon_Particles.Clear();


              } // for Photon Fill

              // EFlowPhoton fill
              cout << " EFlowPhoton size = " << branchEFlowPhoton->GetEntries() << endl;

              for(int i=0;i<branchEFlowPhoton->GetEntries();i++)
              {
                  eflowPhoton = (Photon *) branchEFlowPhoton->At(i);

                  EFlowPhoton_PT            =  eflowPhoton->PT;
                  EFlowPhoton_Eta           =  eflowPhoton->Eta;
                  EFlowPhoton_Phi           =  eflowPhoton->Phi;
                  EFlowPhoton_E             =  eflowPhoton->E;
                  EFlowPhoton_T             =  eflowPhoton->T;
                  EFlowPhoton_EhadOverEem   =  eflowPhoton->EhadOverEem;
                  for(int n=0;n<(eflowPhoton->Particles).GetEntriesFast();n++)
                  {
                      EFlowPhoton_Particles.Add(eflowPhoton->Particles.At(n));
                  }
                  EFlowPhoton_IsolationVar  =  eflowPhoton->IsolationVar;
                  EFlowPhoton_IsolationVarRhoCorr = eflowPhoton->IsolationVarRhoCorr;
                  EFlowPhoton_SumPtCharged  =  eflowPhoton->SumPtCharged;
                  EFlowPhoton_SumPtNeutral  =  eflowPhoton->SumPtNeutral;
                  EFlowPhoton_SumPtChargedPU=  eflowPhoton->SumPtChargedPU;
                  EFlowPhoton_SumPt         =  eflowPhoton->SumPt;
                  EFlowPhoton_Status        =  eflowPhoton->Status;
                  EFlowPhoton_P4.SetPtEtaPhiM(eflowPhoton->P4().Pt(),eflowPhoton->P4().Eta(),eflowPhoton->P4().Phi(),eflowPhoton->P4().M());

                  treeEFlowPhoton->Fill();


              } // for EFlowPhoton fill

              // Electron Fill
              cout << " Electron size= " << branchElectron->GetEntries() << endl;

              for(int i=0;i<branchElectron->GetEntries();i++)
              {
                  electron = (Electron *) branchElectron->At(i);

                  Electron_PT            =  electron->PT;
                  Electron_Eta           =  electron->Eta;
                  Electron_Phi           =  electron->Phi;
                  Electron_T             =  electron->T;
                  Electron_Charge        =  electron->Charge;
                  Electron_EhadOverEem   =  electron->EhadOverEem;
                  Electron_Particle      =  electron->Particle;
                  Electron_IsolationVar  =  electron->IsolationVar;
                  Electron_IsolationVarRhoCorr = electron->IsolationVarRhoCorr;
                  Electron_SumPtCharged  =  electron->SumPtCharged;
                  Electron_SumPtNeutral  =  electron->SumPtNeutral;
                  Electron_SumPtChargedPU=  electron->SumPtChargedPU;
                  Electron_SumPt         =  electron->SumPt;
                  Electron_D0            =  electron->D0;
                  Electron_DZ            =  electron->DZ;
                  Electron_ErrorD0       =  electron->ErrorD0;
                  Electron_ErrorDZ       =  electron->ErrorDZ;
                  //Electron_P4.SetPtEtaPhiM(electron->P4().Pt(),electron->P4().Eta(),electron->P4().Phi(),electron->P4().M());
                  Electron_P4.SetPtEtaPhiM(electron->P4().Pt(),electron->P4().Eta(),electron->P4().Phi(),0.51);

                  treeElectron->Fill();

              } // electron fill

          } // cut if

      } //for
  } // while

  std::cout<<"while end\t"<<std::endl;

  outputFile->cd();
  //  treeMuon->Write();
cout << " open outputFile " << endl;
  treeEvent->Write();
cout << " write event " << endl;
  treeParticle->Write();
cout << " write Particle " << endl;
  treeGenJet->Write();
cout << " write GenJet " << endl;
  treeGenMissingET->Write();
cout << " write GenMissingET " << endl;
  treeTrack->Write();
cout << " write Track " << endl;
  treeTower->Write();
cout << " write Tower " << endl;
  treeEFlowTrack->Write();
cout << " write EFlowTrack " << endl;
  treeEFlowPhoton->Write();
cout << " write EFlowPhoton " << endl;
  treePhoton->Write();
cout << " write Photon " << endl;
  treeElectron->Write();
cout << " write Electron " << endl;
  treeJet->Write();
cout << " write Jet " << endl;
  treeMissingET->Write();
 cout << " write MissingET " << endl;
  treeScalarHT->Write();
cout << " write ScalarHT " << endl;
  treeMuon->Write();
cout << " write Muon " << endl;
  Z_mass->Write();
  outputFile->Close();
cout << " close outputFile " << endl;

}

