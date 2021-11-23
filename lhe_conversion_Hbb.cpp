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

int Accepted(TLorentzVector& mu){
    if(mu.Pt()>3 && mu.CosTheta()<0.85)
      return 1;
    else
      return 0;
}


double helicityAngle(TLorentzVector p4mother, TLorentzVector p4, TLorentzVector p4daughter){
    TVector3 boost = -p4.BoostVector();
    p4mother.Boost(boost);
    p4daughter.Boost(boost);
    TVector3 p3daughter= p4daughter.Vect();
    TVector3 p3mother= p4mother.Vect();

    return p3daughter.Dot(p3mother)/((p3daughter.Mag())*(p3mother.Mag()));
}

float helicityAimuthalAngle(TLorentzVector p4b1, TLorentzVector p4b2, TLorentzVector p4d1, TLorentzVector p4d2){

    TLorentzVector p4m=p4d1+p4d2;
    TVector3 zaxis=(p4m.Vect()).Unit();
    TVector3 beta=(-1./p4m.E())*p4m.Vect();
    p4d1.Boost(beta);
    p4d2.Boost(beta);
    p4b1.Boost(beta);
    p4b2.Boost(beta);
    TVector3 yaxis=((p4b1.Vect()).Cross(p4b2.Vect())).Unit();
    TVector3 xaxis=(yaxis.Cross(zaxis)).Unit();
    float phi= TMath::ATan2((p4d1.Vect()).Dot(yaxis),(p4d1.Vect()).Dot(xaxis));
    return phi;
}

void computeAngles(TLorentzVector p4M11, 
            TLorentzVector p4M12, 
            TLorentzVector p4M21,
            TLorentzVector p4M22,
            float& costhetastar, 
            float& costheta1, 
            float& costheta2, 
            float& Phi, 
            float& Phi1){

    //build Z 4-vectors
    TLorentzVector p4Z = p4M11 + p4M12;

    //build H 4-vectors
    TLorentzVector p4H = p4M21+p4M22; 

    TLorentzVector p4E = p4Z + p4H;
    TLorentzVector p4E_target(0-p4E.Px(),0-p4E.Py(),0-p4E.Pz(),p4E.E()); 

    //build mother of H
    TLorentzVector p4m(0,0,0,240);
    TLorentzVector p4m_target(0-p4m.Px(),0-p4m.Py(),0-p4m.Pz(),p4m.E());

    // -----------------------------------

    TLorentzVector p4h(0-p4Z.Px(),0-p4Z.Py(),0-p4Z.Pz(),p4Z.E());

    costhetastar=helicityAngle(p4m,p4h,p4M21);
    costheta1=p4Z.CosTheta();//helicityAngle(p4m,p4m,p4Z);
    costheta2=helicityAngle(p4m,p4Z,p4M11);
    Phi=helicityAimuthalAngle(p4E,p4E_target,p4M11,p4M12);
    Phi1=helicityAimuthalAngle(p4E,p4E_target,p4M21,p4M22);

}


//int lhe_conversion_HZZ(TString argv="100kevt1000.lhe") {
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

    bool DEBUG=false;
//    bool DEBUG=true;
    double weight;
    int nparticles;

    //momentum for Z, H, Z (on shell), Zs (offshell) and the decays (muon)
    std::vector< std::pair<int,TLorentzVector> >  P4_Particle;
    std::vector<TLorentzVector> P4_Beamp;
    std::vector<TLorentzVector> P4_Beamm;
    std::vector<TLorentzVector> P4_Z;
    std::vector<TLorentzVector> P4_H;
    std::vector<TLorentzVector> P4_H_b;
    std::vector<TLorentzVector> P4_H_bs;
    std::vector<TLorentzVector> P4_Z_Mup;
    std::vector<TLorentzVector> P4_Z_Mum;

    double Px_Beamp, Py_Beamp, Pz_Beamp, E_Beamp; // e+
    double Px_Beamm, Py_Beamm, Pz_Beamm, E_Beamm; // e-
    double Px_Z, Py_Z, Pz_Z, E_Z, P_Z, M_Z, Pt_Z; // Z
    double Px_H, Py_H, Pz_H, E_H, P_H, M_H, Pt_H; // Higgs
    double Px_H_b, Py_H_b, Pz_H_b, E_H_b, P_H_b, M_H_b, Pt_H_b; // Higgs->b1
    double Px_H_bs, Py_H_bs, Pz_H_bs, E_H_bs, P_H_bs, M_H_bs, Pt_H_bs; // Higgs->b2
    double Px_Z_Mup, Py_Z_Mup, Pz_Z_Mup, E_Z_Mup, P_Z_Mup, M_Z_Mup, Pt_Z_Mup; // Z->mu(on)
    double Px_Z_Mum, Py_Z_Mum, Pz_Z_Mum, E_Z_Mum, P_Z_Mum, M_Z_Mum, Pt_Z_Mum; // Z->mu(off)
    float costheta1, costheta2, costhetas;
    float phi, phi1;

    TFile *outfile = new TFile("trial_HZZ.root","RECREATE");
    TTree *out_tree = new TTree("trialTree","trialTree");
    TH1D *m_Higgs = new TH1D("m_Higgs","m_higgs",80,123,127); 
    TH1D *m_Zboson = new TH1D("m_Zboson","m_Zboson",20,80,100);
    TH1D *m_HZboson = new TH1D("m_HZboson","H->ZZstar",20,80,100);
    TH1D *m_HZstar= new TH1D("m_HZstar","H->ZZstar",90,10,100);
    TH1D *cos_theta1 = new TH1D("cos_theta1","cos_theta1",25,-1,1);
    TH1D *cos_theta = new TH1D("cos_theta","cos_theta1 and cos_theta2",25,-1,1);
    TH1D *cos_theta2 = new TH1D("cos_theta2","cos_theta2",25,-1,1);
    TH1D *cos_thetastar = new TH1D("cos_thetastar","cos_thetastar",25,-1,1);
    TH1D *plot_phi = new TH1D("phi","phi",50,-3.2,3.2);
    TH1D *plot_phi1 = new TH1D("phi1","phi1",50,-3.2,3.2);

    out_tree->Branch("weight",&weight,"weight/D");
    out_tree->Branch("nparticles",&nparticles,"nparticles/I");
    out_tree->Branch("costheta1",&costheta1,"costheta1/F");
    out_tree->Branch("costheta2",&costheta2,"costheta2/F");
    out_tree->Branch("costhetas",&costhetas,"costhetas/F");
    out_tree->Branch("phi",&phi,"phi/F");
    out_tree->Branch("phi1",&phi1,"phi1/F");

    // e+
    out_tree->Branch("Px_Beamp",&Px_Beamp,"Px_Beamp/D");
    out_tree->Branch("Py_Beamp",&Py_Beamp,"Py_Beamp/D");
    out_tree->Branch("Pz_Beamp",&Pz_Beamp,"Pz_Beamp/D");
    out_tree->Branch("E_Beamp",&E_Beamp,"E_Beamp/D");

    // e-
    out_tree->Branch("Px_Beamm",&Px_Beamm,"Px_Beamm/D");
    out_tree->Branch("Py_Beamm",&Py_Beamm,"Py_Beamm/D");
    out_tree->Branch("Pz_Beamm",&Pz_Beamm,"Pz_Beamm/D");
    out_tree->Branch("E_Beamm",&E_Beamm,"E_Beamm/D");

    // Z
    out_tree->Branch("Px_Z",&Px_Z,"Px_Z/D");
    out_tree->Branch("Py_Z",&Py_Z,"Py_Z/D");
    out_tree->Branch("Pz_Z",&Pz_Z,"Pz_Z/D");
    out_tree->Branch("E_Z",&E_Z,"E_Z/D");
    out_tree->Branch("P_Z",&P_Z,"P_Z/D");
    out_tree->Branch("M_Z",&M_Z,"M_Z/D");
    out_tree->Branch("Pt_Z",&Pt_Z,"Pt_Z/D");

    // Higgs
    out_tree->Branch("Px_H",&Px_H,"Px_H/D");
    out_tree->Branch("Py_H",&Py_H,"Py_H/D");
    out_tree->Branch("Pz_H",&Pz_H,"Pz_H/D");
    out_tree->Branch("E_H",&E_H,"E_H/D");
    out_tree->Branch("P_H",&P_H,"P_H/D");
    out_tree->Branch("M_H",&M_H,"M_H/D");
    out_tree->Branch("Pt_H",&Pt_H,"Pt_H/D");

    // Higgs->b1
    out_tree->Branch("Px_H_b",&Px_H_b,"Px_H_b/D");
    out_tree->Branch("Py_H_b",&Py_H_b,"Py_H_b/D");
    out_tree->Branch("Pz_H_b",&Pz_H_b,"Pz_H_b/D");
    out_tree->Branch("E_H_b",&E_H_b,"E_H_b/D");
    out_tree->Branch("P_H_b",&P_H_b,"P_H_b/D");
    out_tree->Branch("M_H_b",&M_H_b,"M_H_b/D");
    out_tree->Branch("Pt_H_b",&Pt_H_b,"Pt_H_b/D");

    // Higgs->b2
    out_tree->Branch("Px_H_bs",&Px_H_bs,"Px_H_bs/D");
    out_tree->Branch("Py_H_bs",&Py_H_bs,"Py_H_bs/D");
    out_tree->Branch("Pz_H_bs",&Pz_H_bs,"Pz_H_bs/D");
    out_tree->Branch("E_H_bs",&E_H_bs,"E_H_bs/D");
    out_tree->Branch("P_H_bs",&P_H_bs,"P_H_bs/D");
    out_tree->Branch("M_H_bs",&M_H_bs,"M_H_bs/D");
    out_tree->Branch("Pt_H_bs",&Pt_H_bs,"Pt_H_bs/D");

    // Z->mu(on)
    out_tree->Branch("Px_Z_Mup",&Px_Z_Mup,"Px_Z_Mup/D");
    out_tree->Branch("Py_Z_Mup",&Py_Z_Mup,"Py_Z_Mup/D");
    out_tree->Branch("Pz_Z_Mup",&Pz_Z_Mup,"Pz_Z_Mup/D");
    out_tree->Branch("E_Z_Mup",&E_Z_Mup,"E_Z_Mup/D");
    out_tree->Branch("P_Z_Mup",&P_Z_Mup,"P_Z_Mup/D");
    out_tree->Branch("M_Z_Mup",&M_Z_Mup,"M_Z_Mup/D");
    out_tree->Branch("Pt_Z_Mup",&Pt_Z_Mup,"Pt_Z_Mup/D");

    // Z->mu(off)
    out_tree->Branch("Px_Z_Mum",&Px_Z_Mum,"Px_Z_Mum/D");
    out_tree->Branch("Py_Z_Mum",&Py_Z_Mum,"Py_Z_Mum/D");
    out_tree->Branch("Pz_Z_Mum",&Pz_Z_Mum,"Pz_Z_Mum/D");
    out_tree->Branch("E_Z_Mum",&E_Z_Mum,"E_Z_Mum/D");
    out_tree->Branch("P_Z_Mum",&P_Z_Mum,"P_Z_Mum/D");
    out_tree->Branch("M_Z_Mum",&M_Z_Mum,"M_Z_Mum/D");
    out_tree->Branch("Pt_Z_Mum",&Pt_Z_Mum,"Pt_Z_Mum/D");

    // Now loop over all events:
    long ieve = 0;
    while ( reader.readEvent() ) 
    {
        ++ieve;
        if (ieve % 100 == 0) std::cerr << "event " << ieve << "\n" ;

        P4_Particle.clear();
        P4_Beamp.clear();
        P4_Beamm.clear();
        P4_Z.clear();
        P4_H.clear();
        P4_H_b.clear();
        P4_H_bs.clear();
        P4_Z_Mup.clear();
        P4_Z_Mum.clear();
        //PG loop over particles in the event
        for (int iPart = 0 ; iPart < reader.hepeup.IDUP.size (); iPart++)
        {  

            //if (reader.hepeup.ISTUP.at (iPart) < 0) continue ;  // for removing initial particle beams

            int type = reader.hepeup.IDUP.at (iPart) ;
            int mother=reader.hepeup.MOTHUP.at(iPart).first;
            TLorentzVector particle 
                (
                 reader.hepeup.PUP.at (iPart).at (0), //PG px
                 reader.hepeup.PUP.at (iPart).at (1), //PG py
                 reader.hepeup.PUP.at (iPart).at (2), //PG pz
                 reader.hepeup.PUP.at (iPart).at (3)  //PG E
                ) ;

            weight = reader.hepeup.XWGTUP;
            nparticles = reader.hepeup.NUP;
            if(DEBUG)
                std::cout<<"pdgid "<<type<<" "<<mother<<" "<<particle.E()<<" "<<particle.M()<<std::endl;
            P4_Particle.push_back(std::make_pair(type,particle)) ;
            if (type == -11){
                P4_Beamp.push_back(particle) ;
                Px_Beamp= particle.Px();
                Py_Beamp= particle.Py();
                Pz_Beamp= particle.Pz();
                E_Beamp= particle.E();
                continue;
            }//positron

            if (type == 11){ 
                P4_Beamm.push_back(particle) ;
                Px_Beamm= particle.Px();
                Py_Beamm= particle.Py();
                Pz_Beamm= particle.Pz();
                E_Beamm= particle.E();
                continue;
            }//electron

            if (type == 23 && fabs(P4_Particle.at(mother-1).first) == 11){ 
                P4_Z.push_back(particle);
                Px_Z = particle.Px();
                Py_Z = particle.Py();
                Pz_Z = particle.Pz();
                E_Z = particle.E();
                P_Z = sqrt(Px_Z*Px_Z+Py_Z*Py_Z+Pz_Z*Pz_Z);
                M_Z = sqrt(E_Z*E_Z-P_Z*P_Z);
                Pt_Z = sqrt(Px_Z*Px_Z+Py_Z*Py_Z);
                continue;
	    } // Z

            if (type == 25 && fabs(P4_Particle.at(mother-1).first) == 11){ 
                P4_H.push_back(particle) ;
                Px_H = particle.Px();
                Py_H = particle.Py();
                Pz_H = particle.Pz();
                E_H = particle.E();
                P_H = sqrt(Px_H*Px_H+Py_H*Py_H+Pz_H*Pz_H);
                M_H = sqrt(E_H*E_H-P_H*P_H);
                Pt_H = sqrt(Px_H*Px_H+Py_H*Py_H);
                continue;
            } // Higgs

            if (type == -13 && P4_Particle.at(mother-1).first == 23 &&P4_Z_Mup.size()==0){ 
                P4_Z_Mup.push_back(particle) ;
                Px_Z_Mup = particle.Px();
                Py_Z_Mup = particle.Py();
                Pz_Z_Mup = particle.Pz();
                E_Z_Mup = particle.E();
                P_Z_Mup = sqrt(Px_Z_Mup*Px_Z_Mup+Py_Z_Mup*Py_Z_Mup+Pz_Z_Mup*Pz_Z_Mup);
                M_Z_Mup = sqrt(E_Z_Mup*E_Z_Mup-P_Z_Mup*P_Z_Mup);
                Pt_Z_Mup = sqrt(Px_Z_Mup*Px_Z_Mup+Py_Z_Mup*Py_Z_Mup);
                continue;
            } // Z->mu(on)

            if (type == 13 && P4_Particle.at(mother-1).first == 23 &&P4_Z_Mum.size()==0){ 
                P4_Z_Mum.push_back(particle) ;
                Px_Z_Mum = particle.Px();
                Py_Z_Mum = particle.Py();
                Pz_Z_Mum = particle.Pz();
                E_Z_Mum = particle.E();
                P_Z_Mum = sqrt(Px_Z_Mum*Px_Z_Mum+Py_Z_Mum*Py_Z_Mum+Pz_Z_Mum*Pz_Z_Mum);
                M_Z_Mum = sqrt(E_Z_Mum*E_Z_Mum-P_Z_Mum*P_Z_Mum);
                Pt_Z_Mum = sqrt(Px_Z_Mum*Px_Z_Mum+Py_Z_Mum*Py_Z_Mum);
                continue;
            } // Z->mu(off)

            if (type == 5 && P4_Particle.at(mother-1).first == 25){ 

                P4_H_b.push_back(particle) ;
                Px_H_b = particle.Px();
                Py_H_b = particle.Py();
                Pz_H_b = particle.Pz();
                E_H_b = particle.E();
                P_H_b = sqrt(Px_H_b*Px_H_b+Py_H_b*Py_H_b+Pz_H_b*Pz_H_b);
                M_H_b = sqrt(E_H_b*E_H_b-P_H_b*P_H_b);
                Pt_H_b = sqrt(Px_H_b*Px_H_b+Py_H_b*Py_H_b);
                continue;
            }

            if (type == -5 && P4_Particle.at(mother-1).first == 25){ 

                P4_H_bs.push_back(particle) ;
                Px_H_bs = particle.Px();
                Py_H_bs = particle.Py();
                Pz_H_bs = particle.Pz();
                E_H_bs = particle.E();
                P_H_bs = sqrt(Px_H_bs*Px_H_bs+Py_H_bs*Py_H_bs+Pz_H_bs*Pz_H_bs);
                M_H_bs = sqrt(E_H_bs*E_H_bs-P_H_bs*P_H_bs);
                Pt_H_bs = sqrt(Px_H_bs*Px_H_bs+Py_H_bs*Py_H_bs);
                continue;
            }


        } //PG loop over particles in the event
        if(DEBUG)
          std::cout<<" mass "<<P4_Z.at(0).M()<<" "<<P4_H.at(0).M()<<" "<<P4_Z_Mup.at(0).M()<<" "<<P4_Z_Mum.at(0).M()<<" "<<P4_H_b.at(0).M()<<" "<<P4_H_bs.at(0).M()<<std::endl;

        computeAngles(P4_Z_Mup.at(0),P4_Z_Mum.at(0),P4_H_b.at(0),P4_H_bs.at(0),costhetas, costheta1, costheta2, phi, phi1);
        m_Higgs->Fill((P4_H_b.at(0)+P4_H_bs.at(0)).M(),weight);
        m_Zboson->Fill((P4_Z_Mup.at(0)+P4_Z_Mum.at(0)).M(),weight);
        if(DEBUG)
        std::cout<<P4_H_bs.at(0).M()<<std::endl;
        cos_theta1->Fill(costheta1,weight);
        cos_theta->Fill(costheta1,weight/2.);
        cos_theta->Fill(costheta2,weight/2.);
        cos_theta2->Fill(costheta2,weight);
        cos_thetastar->Fill(costhetas,weight);
        plot_phi->Fill(phi,weight);
        plot_phi1->Fill(phi1,weight);

        out_tree->Fill();

    } // Now loop over all events

    outfile->Write();
    outfile->Close();
    // Now we are done.
    return 0 ;
}

