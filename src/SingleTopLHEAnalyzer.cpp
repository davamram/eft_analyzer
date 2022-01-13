#define SingleTopLHEAnalyzer_cxx
#include "../include/SingleTopLHEAnalyzer.hpp"

#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

using namespace std;

// double weight_SM_size;
// double weight_ctwi_m5_size, weight_ctwi_m2_size, weight_ctwi_m1_size, weight_ctwi_p1_size, weight_ctwi_p2_size, weight_ctwi_p5_size;
// double weight_cbwi_m5_size, weight_cbwi_m2_size, weight_cbwi_m1_size, weight_cbwi_p1_size, weight_cbwi_p2_size, weight_cbwi_p5_size;

void SingleTopLHEAnalyzer::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L SingleTopLHEAnalyzer.C
//      Root > SingleTopLHEAnalyzer t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

   float weight;
   float sinTheta, cosTheta;
   float sinThetaStar, cosThetaStar;
   float sinPhiStar, cosPhiStar, PhiStar;
   float lepton_E_Wframe;
   float top_pt, W_pt, lepton_pt;
   float top_mass, W_mass, W_transverse_mass;
   float nature_lepton;
   float weight_SM;
   float weight_cptb_m10, weight_cptb_m5, weight_cptb_p5, weight_cptb_p10;
   float weight_cptbi_m10, weight_cptbi_m5, weight_cptbi_p5, weight_cptbi_p10;
   float weight_ctw_m2, weight_ctw_m1, weight_ctw_p1, weight_ctw_p2;
   float weight_ctwi_m2, weight_ctwi_m1, weight_ctwi_p1, weight_ctwi_p2;
   float weight_cbw_m2, weight_cbw_m1, weight_cbw_p1, weight_cbw_p2;
   float weight_cbwi_m2, weight_cbwi_m1, weight_cbwi_p1, weight_cbwi_p2;
  //  Ctwi = -2 -> cbwi = {-2,-1,1,2}
   float weight_ctwi_m2_cbwi_m2, weight_ctwi_m2_cbwi_m1, weight_ctwi_m2_cbwi_p1, weight_ctwi_m2_cbwi_p2;
  //  Ctwi = -1 -> cbwi = {-2,-1,1,2}
   float weight_ctwi_m1_cbwi_m2, weight_ctwi_m1_cbwi_m1, weight_ctwi_m1_cbwi_p1, weight_ctwi_m1_cbwi_p2;
  // Ctwi = 1 -> cbwi = {-2,-2,1,2}
   float weight_ctwi_p1_cbwi_m2, weight_ctwi_p1_cbwi_m1, weight_ctwi_p1_cbwi_p1, weight_ctwi_p1_cbwi_p2;
  // Ctwi = 2 -> cbwi = {-2,-1,1,2}
   float weight_ctwi_p2_cbwi_m2, weight_ctwi_p2_cbwi_m1, weight_ctwi_p2_cbwi_p1, weight_ctwi_p2_cbwi_p2;

   double weight_sum[9];

   TFile* fOutput = new TFile("output.root","RECREATE");
   TTree* tOutput = new TTree("Tree","Tree");

   tOutput->Branch("nature_lepton",&nature_lepton,"nature_lepton/F");
   tOutput->Branch("weight",&weight,"weight/F");
   tOutput->Branch("cosTheta",&cosTheta,"cosTheta/F");
   tOutput->Branch("sinTheta",&sinTheta,"sinTheta/F");
   tOutput->Branch("cosThetaStar",&cosThetaStar,"cosThetaStar/F");
   tOutput->Branch("sinThetaStar",&sinThetaStar,"sinThetaStar/F");
   tOutput->Branch("cosPhiStar",&cosPhiStar,"cosPhiStar/F");
   tOutput->Branch("sinPhiStar",&sinPhiStar,"sinPhiStar/F");
   tOutput->Branch("PhiStar",&PhiStar,"PhiStar/F");
   tOutput->Branch("lepton_E_Wframe",&lepton_E_Wframe,"lepton_E_Wframe/F");
   tOutput->Branch("top_pt",&top_pt,"top_pt/F");
   tOutput->Branch("W_pt",&W_pt,"W_pt/F");
   tOutput->Branch("lepton_pt",&lepton_pt,"lepton_pt/F");
   tOutput->Branch("W_mass", &W_mass, "W_mass/F");
   tOutput->Branch("top_mass", &top_mass, "top_mass/F");
   tOutput->Branch("W_transverse_mass", &W_transverse_mass, "W_transverse_mass/F");

   tOutput->Branch("weight_SM",&weight_SM,"weight_SM/F");

   tOutput->Branch("weight_cptb_m10",&weight_cptb_m10,"weight_cptb_m10/F");
   tOutput->Branch("weight_cptb_m5",&weight_cptb_m5,"weight_cptb_m5/F");
   tOutput->Branch("weight_cptb_p5",&weight_cptb_p5,"weight_cptb_p5/F");
   tOutput->Branch("weight_cptb_p10",&weight_cptb_p10,"weight_cptb_p10/F");

   tOutput->Branch("weight_cptbi_m10",&weight_cptbi_m10,"weight_cptbi_m10/F");
   tOutput->Branch("weight_cptbi_m5",&weight_cptbi_m5,"weight_cptbi_m5/F");
   tOutput->Branch("weight_cptbi_p5",&weight_cptbi_p5,"weight_cptbi_p5/F");
   tOutput->Branch("weight_cptbi_p10",&weight_cptbi_p10,"weight_cptbi_p10/F");

   tOutput->Branch("weight_ctw_m2",&weight_ctw_m2,"weight_ctw_m2/F");
   tOutput->Branch("weight_ctw_m1",&weight_ctw_m1,"weight_ctw_m1/F");
   tOutput->Branch("weight_ctw_p1",&weight_ctw_p1,"weight_ctw_p1/F");
   tOutput->Branch("weight_ctw_p2",&weight_ctw_p2,"weight_ctw_p2/F");

   tOutput->Branch("weight_ctwi_m2",&weight_ctwi_m2,"weight_ctwi_m2/F");
   tOutput->Branch("weight_ctwi_m1",&weight_ctwi_m1,"weight_ctwi_m1/F");
   tOutput->Branch("weight_ctwi_p1",&weight_ctwi_p1,"weight_ctwi_p1/F");
   tOutput->Branch("weight_ctwi_p2",&weight_ctwi_p2,"weight_ctwi_p2/F");

   tOutput->Branch("weight_cbw_m2",&weight_cbw_m2,"weight_cbw_m2/F");
   tOutput->Branch("weight_cbw_m1",&weight_cbw_m1,"weight_cbw_m1/F");
   tOutput->Branch("weight_cbw_p1",&weight_cbw_p1,"weight_cbw_p1/F");
   tOutput->Branch("weight_cbw_p2",&weight_cbw_p2,"weight_cbw_p2/F");

   tOutput->Branch("weight_cbwi_m2",&weight_cbwi_m2,"weight_cbwi_m2/F");
   tOutput->Branch("weight_cbwi_m1",&weight_cbwi_m1,"weight_cbwi_m1/F");
   tOutput->Branch("weight_cbwi_p1",&weight_cbwi_p1,"weight_cbwi_p1/F");
   tOutput->Branch("weight_cbwi_p2",&weight_cbwi_p2,"weight_cbwi_p2/F");

   tOutput->Branch("weight_ctwi_m2_cbwi_m2",&weight_ctwi_m2_cbwi_m2,"weight_ctwi_m2_cbwi_m2/F");
   tOutput->Branch("weight_ctwi_m2_cbwi_m1",&weight_ctwi_m2_cbwi_m1,"weight_ctwi_m2_cbwi_m1/F");
   tOutput->Branch("weight_ctwi_m2_cbwi_p1",&weight_ctwi_m2_cbwi_p1,"weight_ctwi_m2_cbwi_p1/F");
   tOutput->Branch("weight_ctwi_m2_cbwi_p2",&weight_ctwi_m2_cbwi_p2,"weight_ctwi_m2_cbwi_p2/F");

   tOutput->Branch("weight_ctwi_m1_cbwi_m2",&weight_ctwi_m1_cbwi_m2,"weight_ctwi_m1_cbwi_m2/F");
   tOutput->Branch("weight_ctwi_m1_cbwi_m1",&weight_ctwi_m1_cbwi_m1,"weight_ctwi_m1_cbwi_m1/F");
   tOutput->Branch("weight_ctwi_m1_cbwi_p1",&weight_ctwi_m1_cbwi_p1,"weight_ctwi_m1_cbwi_p1/F");
   tOutput->Branch("weight_ctwi_m1_cbwi_p2",&weight_ctwi_m1_cbwi_p2,"weight_ctwi_m1_cbwi_p2/F");

   tOutput->Branch("weight_ctwi_p1_cbwi_m2",&weight_ctwi_p1_cbwi_m2,"weight_ctwi_p1_cbwi_m2/F");
   tOutput->Branch("weight_ctwi_p1_cbwi_m1",&weight_ctwi_p1_cbwi_m1,"weight_ctwi_p1_cbwi_m1/F");
   tOutput->Branch("weight_ctwi_p1_cbwi_p1",&weight_ctwi_p1_cbwi_p1,"weight_ctwi_p1_cbwi_p1/F");
   tOutput->Branch("weight_ctwi_p1_cbwi_p2",&weight_ctwi_p1_cbwi_p2,"weight_ctwi_p1_cbwi_p2/F");

   tOutput->Branch("weight_ctwi_p2_cbwi_m2",&weight_ctwi_p2_cbwi_m2,"weight_ctwi_p2_cbwi_m2/F");
   tOutput->Branch("weight_ctwi_p2_cbwi_m1",&weight_ctwi_p2_cbwi_m1,"weight_ctwi_p2_cbwi_m1/F");
   tOutput->Branch("weight_ctwi_p2_cbwi_p1",&weight_ctwi_p2_cbwi_p1,"weight_ctwi_p2_cbwi_p1/F");
   tOutput->Branch("weight_ctwi_p2_cbwi_p2",&weight_ctwi_p2_cbwi_p2,"weight_ctwi_p2_cbwi_p2/F");


   if (fChain == 0) return;

	//Long64_t nentries = 100;
   Long64_t nentries = fChain->GetEntriesFast();
   //cout<<"nentries= "<<nentries<<endl;

   TLorentzVector Ptop;
   TLorentzVector Pb;
   TLorentzVector Pw;
   TLorentzVector Pl;
   TLorentzVector Pnu;
   TLorentzVector Pqspec;
   double Pl_ID;

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) 
   {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

	  //cout << "nParticles="<<Particle_<<endl;
	  for (int i=0; i<Particle_; i++){

      //If on wants to know the sum of the weights use a cout of this variables (sum of weights is linked to the Xsection)
       if (Rwgt_>0){
        // SM 
         weight_SM = Rwgt_Weight[0];
        // Cptb
         weight_cptb_m10 = Rwgt_Weight[1];
         weight_cptb_m5 = Rwgt_Weight[2];
         weight_cptb_p5 = Rwgt_Weight[3];
         weight_cptb_p10 = Rwgt_Weight[4]; 
        // Cptbi
         weight_cptbi_m10 = Rwgt_Weight[5];
         weight_cptbi_m5 = Rwgt_Weight[6];
         weight_cptbi_p5 = Rwgt_Weight[7];
         weight_cptbi_p10 = Rwgt_Weight[8];
        // Ctw
         weight_ctw_m2 = Rwgt_Weight[9];
         weight_ctw_m1 = Rwgt_Weight[10];
         weight_ctw_p1 = Rwgt_Weight[11];
         weight_ctw_p2 = Rwgt_Weight[12];
        // Ctwi
         weight_ctwi_m2 = Rwgt_Weight[13];
         weight_ctwi_m1 = Rwgt_Weight[14];
         weight_ctwi_p1 = Rwgt_Weight[15];
         weight_ctwi_p2 = Rwgt_Weight[16];
        // Cbw
         weight_cbw_m2 = Rwgt_Weight[17];
         weight_cbw_m1 = Rwgt_Weight[18];
         weight_cbw_p1 = Rwgt_Weight[19];
         weight_cbw_p2 = Rwgt_Weight[20];
        // Cbwi
         weight_cbwi_m2 = Rwgt_Weight[21];
         weight_cbwi_m1 = Rwgt_Weight[22];
         weight_cbwi_p1 = Rwgt_Weight[23];
         weight_cbwi_p2 = Rwgt_Weight[24];
        //  OffGrid ctwi/cbwi
        weight_ctwi_m2_cbwi_m2 = Rwgt_Weight[25];
        weight_ctwi_m2_cbwi_m1 = Rwgt_Weight[26];
        weight_ctwi_m2_cbwi_p1 = Rwgt_Weight[27];
        weight_ctwi_m2_cbwi_p2 = Rwgt_Weight[28];

        weight_ctwi_m1_cbwi_m2 = Rwgt_Weight[29];
        weight_ctwi_m1_cbwi_m1 = Rwgt_Weight[30];
        weight_ctwi_m1_cbwi_p1 = Rwgt_Weight[31];
        weight_ctwi_m1_cbwi_p2 = Rwgt_Weight[32];

        weight_ctwi_p1_cbwi_m2 = Rwgt_Weight[33];
        weight_ctwi_p1_cbwi_m1 = Rwgt_Weight[34];
        weight_ctwi_p1_cbwi_p1 = Rwgt_Weight[35];
        weight_ctwi_p1_cbwi_p2 = Rwgt_Weight[36];

        weight_ctwi_p2_cbwi_m2 = Rwgt_Weight[37];
        weight_ctwi_p2_cbwi_m1 = Rwgt_Weight[38];
        weight_ctwi_p2_cbwi_p1 = Rwgt_Weight[39];
        weight_ctwi_p2_cbwi_p2 = Rwgt_Weight[40];

      }

		  if (TMath::Abs(Particle_PID[i])==24 && Particle_Status[i]==2)
			   Pw.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

		  if (TMath::Abs(Particle_PID[i])==5 && Particle_Status[i]==1)
			  Pb.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

		  if (TMath::Abs(Particle_PID[i])<=5 && Particle_Status[i]==1)
			  Pqspec.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

		  if ((TMath::Abs(Particle_PID[i])==11 || TMath::Abs(Particle_PID[i])==13) && Particle_Status[i]==1)
      {
        Pl.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
        Pl_ID = abs(Particle_PID[i]);
        if (Pl_ID == 11)
        {
          nature_lepton = 1;
        }
        if (Pl_ID == 13)
        {
          nature_lepton = 2;
        }
      }


      if ((TMath::Abs(Particle_PID[i])==12 || TMath::Abs(Particle_PID[i])==14) && Particle_Status[i]==1)
  		  Pnu.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);

	  }

    Ptop = Pb + Pl + Pnu;

	  weight = Event_Weight[jentry];

	  /* SELECTION */
    // M2 selection
  	if ((Pl.Pt()<35 || TMath::Abs(Pl.Eta())>1.479) && Pl_ID==11) continue; //Electron
    if ((Pl.Pt()<26 || TMath::Abs(Pl.Eta())>2.4) && Pl_ID==13) continue; //Muon
	  if (Pqspec.Pt()<40 || TMath::Abs(Pqspec.Eta())>4.7) continue; //Jet 1st selection
    if (abs(Pqspec.Eta()) < 3.0 && 2.7 < abs(Pqspec.Eta()) && Pqspec.Pt() < 60) continue; //Jet 2nd selection
	  if (Pb.Pt()<40 || TMath::Abs(Pb.Eta())>2.4) continue;


    // STreco selection
  	// if ((Pl.Pt()<32 || (TMath::Abs(Pl.Eta()) > 1.4442 && TMath::Abs(Pl.Eta()) < 1.5660)) && Pl_ID==11) continue; //Electron Streco Selection
    // if ((Pl.Pt()<30 || TMath::Abs(Pl.Eta())>2.4) && Pl_ID==13) continue; //Muon Streco Selection [Pt(2016 and 2018: 26Gev ; 2017: 30Gev)]
    // if (Pqspec.Pt()<40 || TMath::Abs(Pqspec.Eta())>4.7) continue; //Jet 1st selection
    // if (TMath::Abs(Pqspec.Eta()) >= 2.4 && Pqspec.Pt() < 60) continue; //Jet Streco 2nd selection for |eta|>=2.4
	  // if (Pb.Pt()<40 || TMath::Abs(Pb.Eta())>2.5) continue; //Streco Selection [Eta for 2016 >2.4 and for 2017/2018 >2.5]
    

	  /* ANGLE RECONSTRUCTION */

	  TVector3 InvTopBoost;  InvTopBoost.SetXYZ(-Ptop.Px()/Ptop.E(),-Ptop.Py()/Ptop.E(),-Ptop.Pz()/Ptop.E());
	  TVector3 InvWBoost; InvWBoost.SetXYZ(-Pw.Px()/Pw.E(),-Pw.Py()/Pw.E(),-Pw.Pz()/Pw.E());

	  Pw.Boost(InvTopBoost);
	  Pb.Boost(InvTopBoost);
	  Pqspec.Boost(InvTopBoost);
	  //cout << "Top rest frame, W: Px="<<Pw.Px()<<" Py="<<Pw.Py()<<" Pz="<<Pw.Pz()<<endl;
	  //cout << "Top rest frame, b: Px="<<Pb.Px()<<" Py="<<Pb.Py()<<" Pz="<<Pb.Pz()<<endl;
	  //cout << "Top rest frame, q spec: Px="<<Pqspec.Px()<<" Py="<<Pqspec.Py()<<" Pz="<<Pqspec.Pz()<<endl;

	  TVector3 Zdir = Pw.Vect().Unit();
	  TVector3 PqspecUnit = Pqspec.Vect().Unit();
	  TVector3 Ydir = PqspecUnit.Cross(Zdir).Unit();
	  TVector3 Xdir = Ydir.Cross(Zdir);

	  sinTheta = PqspecUnit.Cross(Zdir).Mag();
	  cosTheta = PqspecUnit.Dot(Zdir);
	  if (sinTheta<0) sinTheta=-sinTheta;
	  //cout << "Angle Theta: cosTheta="<<cosTheta<<" sinTheta="<<sinTheta<<endl;

	  //cout << "Top rest frame, Xdir: X="<<Xdir.X()<<" Y="<<Xdir.Y()<<" Z="<<Xdir.Z()<< " Mag="<< Xdir.Mag() <<endl;
	  //cout << "Top rest frame, Ydir: X="<<Ydir.X()<<" Y="<<Ydir.Y()<<" Z="<<Ydir.Z()<< " Mag="<< Ydir.Mag() <<endl;
	  //cout << "Top rest frame, Zdir: X="<<Zdir.X()<<" Y="<<Zdir.Y()<<" Z="<<Zdir.Z()<< " Mag="<< Zdir.Mag() <<endl;

	  Pl.Boost(InvWBoost);
	  lepton_E_Wframe = Pl.E();

	  TVector3 PlUnit = Pl.Vect().Unit();
	  cosThetaStar = PlUnit.Dot(Zdir);
	  sinThetaStar = PlUnit.Cross(Zdir).Mag();
	  if (sinThetaStar<0) sinThetaStar=-sinThetaStar;
	  //cout << "Angle ThetaStar: cosThetaStar="<<cosThetaStar<<" sinThetaStar="<<sinThetaStar<<endl;

	  TVector3 PlUnit_PlaneXY = (PlUnit - (PlUnit.Dot(Zdir))*Zdir).Unit();
	  cosPhiStar = PlUnit_PlaneXY.Dot(Xdir);
	  sinPhiStar = PlUnit_PlaneXY.Dot(Ydir);
      //PhiStar;// = TMath::ACos(cosPhiStar);
	  if (sinPhiStar>0) PhiStar = TMath::ACos(cosPhiStar);
	  if (sinPhiStar<0) PhiStar = 2*TMath::Pi()-TMath::ACos(cosPhiStar);
	  //cout << "Angle PhiStar: cosPhiStar="<<cosPhiStar<<" sinPhiStar="<<sinPhiStar<<" PhiStar="<<PhiStar<<endl;
   
    top_pt = Ptop.Pt();
    W_pt = Pw.Pt();
    lepton_pt = Pl.Pt();
    top_mass = Ptop.M();
    W_mass = Pw.M();
    W_transverse_mass = Pw.Mt();
    
    tOutput->Fill();
   }

   tOutput->Write();
   fOutput->Close();
}





/*
         weight_SM_size=Rwgt_Weight[0];

         weight_ctwi_p1_size=Rwgt_Weight[1];

         weight_ctwi_m1_size=Rwgt_Weight[2];

         weight_ctwi_p2_size=Rwgt_Weight[3];

         weight_ctwi_m2_size=Rwgt_Weight[4];

         weight_cbwi_p1_size=Rwgt_Weight[5];

         weight_cbwi_m1_size=Rwgt_Weight[6];

         weight_cbwi_p2_size=Rwgt_Weight[7];

         weight_cbwi_m2_size=Rwgt_Weight[8];

    weight_SM_size+=Rwgt_Weight[0];
    weight_ctwi_p1_size+=Rwgt_Weight[1];
    weight_ctwi_m1_size+=Rwgt_Weight[2];
    weight_ctwi_p2_size+=Rwgt_Weight[3];
    weight_ctwi_m2_size+=Rwgt_Weight[4];
    weight_cbwi_p1_size+=Rwgt_Weight[5];
    weight_cbwi_m1_size+=Rwgt_Weight[6];
    weight_cbwi_p2_size+=Rwgt_Weight[7];
    weight_cbwi_m2_size+=Rwgt_Weight[8];


   cout<<"weight_SM_size= "<<weight_SM_size<<endl;
   cout<<"weight_ctwi_p1_size= "<<weight_ctwi_p1_size<<endl;
   cout<<"weight_ctwi_p2_size= "<<weight_ctwi_p2_size<<endl;
   cout<<"weight_ctwi_m1_size= "<<weight_ctwi_m1_size<<endl;
   cout<<"weight_ctwi_m2_size= "<<weight_ctwi_m2_size<<endl;
   cout<<"weight_cbwi_p1_size= "<<weight_cbwi_p1_size<<endl;
   cout<<"weight_cbwi_p2_size= "<<weight_cbwi_p2_size<<endl;
   cout<<"weight_cbwi_m1_size= "<<weight_cbwi_m1_size<<endl;
   cout<<"weight_cbwi_m2_size= "<<weight_cbwi_m2_size<<endl;



  cout<<"weight_SM_size= "<<Rwgt_Weight[0]<<endl;
  cout<<"weight_ctwi_p1_size= "<<Rwgt_Weight[1]<<endl;
  cout<<"weight_ctwi_p2_size= "<<Rwgt_Weight[2]<<endl;
  cout<<"weight_ctwi_m1_size= "<<Rwgt_Weight[3]<<endl;
  cout<<"weight_ctwi_m2_size= "<<Rwgt_Weight[4]<<endl;
  cout<<"weight_cbwi_p1_size= "<<Rwgt_Weight[5]<<endl;
  cout<<"weight_cbwi_p2_size= "<<Rwgt_Weight[6]<<endl;
  cout<<"weight_cbwi_m1_size= "<<Rwgt_Weight[7]<<endl;
  cout<<"weight_cbwi_m2_size= "<<Rwgt_Weight[8]<<endl;

      cout<<"weight_sum0= "<<weight_sum[0]<<endl;
      cout<<"weight_sum1= "<<weight_sum[1]<<endl;
      cout<<"weight_sum2= "<<weight_sum[2]<<endl;
      cout<<"weight_sum3= "<<weight_sum[3]<<endl;
      cout<<"weight_sum4= "<<weight_sum[4]<<endl;
      cout<<"weight_sum5= "<<weight_sum[5]<<endl;
      cout<<"weight_sum6= "<<weight_sum[6]<<endl;
      cout<<"weight_sum7= "<<weight_sum[7]<<endl;
      cout<<"weight_sum8= "<<weight_sum[8]<<endl;




      weight_sum[0]+=Rwgt_Weight[0];
      weight_sum[1]+=Rwgt_Weight[1];
      weight_sum[2]+=Rwgt_Weight[2];
      weight_sum[3]+=Rwgt_Weight[3];
      weight_sum[4]+=Rwgt_Weight[4];
      weight_sum[5]+=Rwgt_Weight[5];
      weight_sum[6]+=Rwgt_Weight[6];
      weight_sum[7]+=Rwgt_Weight[7];
      weight_sum[8]+=Rwgt_Weight[8];

    weight_SM_size=Rwgt_Weight[0];
    weight_ctwi_p1_size=Rwgt_Weight[1];
    weight_ctwi_m1_size=Rwgt_Weight[2];
    weight_ctwi_p2_size=Rwgt_Weight[3];
    weight_ctwi_m2_size=Rwgt_Weight[4];
    weight_cbwi_p1_size=Rwgt_Weight[5];
    weight_cbwi_m1_size=Rwgt_Weight[6];
    weight_cbwi_p2_size=Rwgt_Weight[7];
    weight_cbwi_m2_size=Rwgt_Weight[8];

    weight_SM_size = 0;
    weight_ctwi_p1_size = 0;
    weight_ctwi_m1_size = 0;
    weight_ctwi_p2_size = 0;
    weight_ctwi_m2_size = 0;
    weight_cbwi_p1_size = 0;
    weight_cbwi_m1_size = 0;
    weight_cbwi_p2_size = 0;
    weight_cbwi_m2_size = 0;




*/
