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
#include <string>

#include <chrono>
#include <math.h>

#define nbLepton 4

using namespace std;
string part="";


// double weight_SM_size;
// double weight_ctwi_m5_size, weight_ctwi_m2_size, weight_ctwi_m1_size, weight_ctwi_p1_size, weight_ctwi_p2_size, weight_ctwi_p5_size;
// double weight_cbwi_m5_size, weight_cbwi_m2_size, weight_cbwi_m1_size, weight_cbwi_p1_size, weight_cbwi_p2_size, weight_cbwi_p5_size;

void SingleTopLHEAnalyzer::Loop()
{
  cout<<part<<endl;
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


   float invMassLepton_Z;
   float invMassLepton_W;
   float invMassB;
   float invMass4Lepton2B;
   float invMass4Lepton;
   float invMassDeltaR;
   float invMass2LeptonW2B;

   float ptLepton1;
   float ptLepton2;
   float ptLepton3;
   float ptLepton4;

   float sumPtLeptonVect_Z;
   float sumPtLeptonVect_W;
   float sumPtLeptonScal_Z;
   float sumPtLeptonScal_W;
   float sumPtBVect;
   float sumPtBScal;

   float angleBtwLeptonZ;
   float angleBtwLeptonZ2;
   float thetaStar1;
   float thetaStar2;
   float deltaRLeptonZ, deltaRLeptonZm, deltaRLeptonZp, deltaRLeptonZr;
   float deltaPhiLeptonZ;


   float natureLeptonZ;
   float natureLeptonW;
   float natureLeptonW_1;
   float natureLeptonW_2;

   float leptonID;
   float pLepton3ID;
   float pLepton4ID;
   float pB1ID;
   float pB2ID;

   float mother1Lepton1;
   float mother1Lepton2;
   float mother2Lepton1;
   float mother2Lepton2;

   int countPart=0;
   int countSame=0;

   int test=0;

   TFile* fOutput = new TFile("output.root","UPDATE");
   TTree* tOutput = new TTree("TreePt","TreePt");
   TTree* tOutput_E = new TTree("TreeElec","TreeElec");
   TTree* tOutput_V = new TTree("TreeMuon","TreeMuon");
   TTree* tOutput_W = new TTree("TreeNatureW","TreeNatureW");

   tOutput->Branch("invMassLepton_Z",&invMassLepton_Z,"invMassLepton_Z/F");
   tOutput->Branch("invMassLepton_W",&invMassLepton_W,"invMassLepton_W/F");
   tOutput->Branch("ptLepton1",&ptLepton1,"ptLepton1/F");
   tOutput->Branch("ptLepton2",&ptLepton2,"ptLepton2/F");
   tOutput->Branch("ptLepton3",&ptLepton3,"ptLepton3/F");
   tOutput->Branch("ptLepton4",&ptLepton4,"ptLepton4/F");
   tOutput->Branch("sumPtLeptonVect_Z",&sumPtLeptonVect_Z,"sumPtLeptonVect_Z/F");
   tOutput->Branch("sumPtLeptonVect_W",&sumPtLeptonVect_W,"sumPtLeptonVect_W/F");
   tOutput->Branch("sumPtLeptonScal_Z",&sumPtLeptonScal_Z,"sumPtLeptonScal_Z/F");
   tOutput->Branch("sumPtLeptonScal_W",&sumPtLeptonScal_W,"sumPtLeptonScal_W/F");
   tOutput->Branch("sumPtBVect",&sumPtBVect,"sumPtBVect/F");
   tOutput->Branch("sumPtBScal",&sumPtBScal,"sumPtBScal/F");
   tOutput->Branch("natureLeptonZ",&natureLeptonZ, "natureLeptonZ/F");
   tOutput->Branch("invMassB", &invMassB, "invMassB/F");
   tOutput->Branch("invMass4Lepton2B", &invMass4Lepton2B, "invMassB/F");
   tOutput->Branch("invMass4Lepton", &invMass4Lepton, "invMass4Lepton/F");
   tOutput->Branch("invMassDeltaR", &invMassDeltaR, "invMassDeltaR/F");
   tOutput->Branch("angleBtwLeptonZ", &angleBtwLeptonZ, "angleBtwLeptonZ/F");
   tOutput->Branch("thetaStar1", &thetaStar1, "thetaStar1");
   tOutput->Branch("thetaStar2", &thetaStar2, "thetaStar2");
   tOutput->Branch("deltaRLeptonZ", &deltaRLeptonZ, "deltaRLeptonZ");
   tOutput->Branch("deltaRLeptonZr", &deltaRLeptonZr, "deltaRLeptonZr");
   tOutput->Branch("deltaRLeptonZp", &deltaRLeptonZp, "deltaRLeptonZp");
   tOutput->Branch("deltaRLeptonZm", &deltaRLeptonZm, "deltaRLeptonZm");

   tOutput->Branch("deltaPhiLeptonZ", &deltaPhiLeptonZ, "deltaPhiLeptonZ");

   tOutput->Branch("invMass2LeptonW2B", &invMass2LeptonW2B, "invMass2LeptonW2B/F");



   tOutput_E->Branch("invMassLeptonPt_Elec_Z", &invMassLepton_Z, "invMassLeptonPt_Elec_Z/F");

   tOutput_V->Branch("invMassLepton_Z",&invMassLepton_Z,"invMassLepton_Z/F");
   tOutput_V->Branch("invMassLepton_W",&invMassLepton_W,"invMassLepton_W/F");
   tOutput_V->Branch("ptLepton1",&ptLepton1,"ptLepton1/F");
   tOutput_V->Branch("ptLepton2",&ptLepton2,"ptLepton2/F");
   tOutput_V->Branch("ptLepton3",&ptLepton3,"ptLepton3/F");
   tOutput_V->Branch("ptLepton4",&ptLepton4,"ptLepton4/F");
   tOutput_V->Branch("sumPtLeptonVect_Z",&sumPtLeptonVect_Z,"sumPtLeptonVect_Z/F");
   tOutput_V->Branch("sumPtLeptonVect_W",&sumPtLeptonVect_W,"sumPtLeptonVect_W/F");
   tOutput_V->Branch("sumPtLeptonScal_Z",&sumPtLeptonScal_Z,"sumPtLeptonScal_Z/F");
   tOutput_V->Branch("sumPtLeptonScal_W",&sumPtLeptonScal_W,"sumPtLeptonScal_W/F");
   tOutput_V->Branch("sumPtBVect",&sumPtBVect,"sumPtBVect/F");
   tOutput_V->Branch("sumPtBScal",&sumPtBScal,"sumPtBScal/F");
   tOutput_V->Branch("natureLeptonZ",&natureLeptonZ, "natureLeptonZ/F");
   tOutput_V->Branch("invMassB", &invMassB, "invMassB/F");
   tOutput_V->Branch("invMass4Lepton2B", &invMass4Lepton2B, "invMassB/F");
   tOutput_V->Branch("invMass4Lepton", &invMass4Lepton, "invMass4Lepton/F");
   tOutput_V->Branch("invMassDeltaR", &invMassDeltaR, "invMassDeltaR/F");
   tOutput_V->Branch("angleBtwLeptonZ", &angleBtwLeptonZ, "angleBtwLeptonZ/F");
   tOutput_V->Branch("invMass2LeptonW2B", &invMass2LeptonW2B, "invMass2LeptonW2B/F");
   tOutput_V->Branch("thetaStar1", &thetaStar1, "thetaStar1");
   tOutput_V->Branch("thetaStar2", &thetaStar2, "thetaStar2");
   tOutput_V->Branch("deltaRLeptonZ", &deltaRLeptonZ, "deltaRLeptonZ");
   tOutput_V->Branch("deltaRLeptonZr", &deltaRLeptonZr, "deltaRLeptonZr");
   tOutput_V->Branch("deltaRLeptonZp", &deltaRLeptonZp, "deltaRLeptonZp");
   tOutput_V->Branch("deltaRLeptonZm", &deltaRLeptonZm, "deltaRLeptonZm");
   tOutput_V->Branch("deltaPhiLeptonZ", &deltaPhiLeptonZ, "deltaPhiLeptonZ");




   tOutput_W->Branch("natureLeptonW", &natureLeptonW, "natureLeptonW/F");


   if (fChain == 0) return;

	//Long64_t nentries = 100;
   Long64_t nentries = fChain->GetEntriesFast();

   TLorentzVector *pLepton = new TLorentzVector[nbLepton];
   TLorentzVector pLeptonPair;
   TLorentzVector pLepton1;
   TLorentzVector pLepton2;
   TLorentzVector pLepton3;
   TLorentzVector pLepton4;
   TLorentzVector pLeptonTot;
   TLorentzVector pDeltaR;

   TLorentzVector pB1;
   TLorentzVector pB2;
   TLorentzVector pBPair;
   TLorentzVector pTop1;
   TLorentzVector pTop2;

   TLorentzVector pLeptonPairW;
   TLorentzVector pLeptonB;

   int tabMother1Lepton[nbLepton];
   int tabMother2Lepton[nbLepton];
   int tabLeptonID[nbLepton];
   int tabPosLepton[nbLepton];



   int count=0;
   int countDiffMother=0;


   Long64_t nbytes = 0, nb = 0;
   auto start = std::chrono::steady_clock::now();
   for (Long64_t jentry=0; jentry<nentries;jentry++)
   {
      // Progression bar
      auto end = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed = end - start;
      float progress = (float)jentry / nentries * 100;
      float time_per_event = elapsed.count() / (jentry + 1);
      float remaining_time = time_per_event * (nentries - jentry - 1);
      if((jentry % (nentries/10000)==0)) std::cout << "\rProgress: [" << std::string(progress / 2, '#') << std::string(50 - progress / 2, ' ') << "] " << (float)((int)progress*100)/100 << "% | Estimated remaining time: " << round(remaining_time) <<"s    "<<flush;


      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      //if (jentry > 10) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      int countLepton=0;
      int countB=0;

	  for (int i=0; i<Particle_; i++){

      if (abs(Particle_PID[i])==5){
        if (!countB){
          pB1.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
          pB1ID=Particle_PID[i];
        }
        if (countB){
          pB2.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
          pB2ID=Particle_PID[i];
        }
        countB++;
      }

      // Keep only leptons

		  if ((abs(Particle_PID[i])==11 || abs(Particle_PID[i])==13) && Particle_Status[i]==1)
      {

        // Be warned that it make the assumption that only 4 lepton were detected
        // Wich is the case here specificaly

        // We store the lepton position
        pLepton[countLepton].SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
        tabMother1Lepton[countLepton]=Particle_PID[Particle_Mother1[i]];
        tabMother2Lepton[countLepton]=Particle_PID[Particle_Mother2[i]];
        tabLeptonID[countLepton]=Particle_PID[i];
        tabPosLepton[countLepton]=i;


        // Then move to the next case

        countLepton++;

      }

	  }
    // Here place what you want to do with particle proprties



    // Event selection for Muons (Electrons)
    // Pt > 40 for the first one
    // Pt > 10 and |eta|< 2.4 (2.5)

    float pMax=-1;
    for(int i=0; i<nbLepton;i++){
      if(pLepton[i].Pt()>pMax){
        pMax=pLepton[i].Pt();
      }
    }
    if (pMax<40) continue;

    int cut=0;
    for(int i=0; i<nbLepton;i++){
      if ((pLepton[i].Pt()<10 || fabs(pLepton[i].Eta())>2.4) && fabs(tabLeptonID[i])==13) cut++; //Muon
      if ((pLepton[i].Pt()<10 || fabs(pLepton[i].Eta())>2.5) && fabs(tabLeptonID[i])==11) cut++; //Electron
    }
    if (cut) continue;
    count++;

    //###################################
    //Some variable to test
    //###################################
    // min deltaR between b and lepton
    float deltaR1, deltaR2;
    float deltaR = pB1.DeltaR(pLepton[0]);
    TLorentzVector pBDeltaR = pB1;
    TLorentzVector pLeptonDeltaR = pLepton[0];
    for (int i=0;i<nbLepton;i++){
      deltaR1=pB1.DeltaR(pLepton[i]);
      deltaR2=pB2.DeltaR(pLepton[i]);
      if (deltaR1<deltaR){
        deltaR=deltaR1;
        pBDeltaR=pB1;
        pLeptonDeltaR=pLepton[i];
      }
      if (deltaR2<deltaR){
        deltaR=deltaR2;
        pBDeltaR=pB2;
        pLeptonDeltaR=pLepton[i];
      }
    }
    pDeltaR=pLeptonDeltaR+pBDeltaR;
    invMassDeltaR=pDeltaR.M();



if (part=="mad"){
      //###############################
      //Part One : Use madgraph information on particle history
      //###############################
        int countUsefulPair=0;
        int countUnUsufulPair=0;

    	  for (int i=0; i<Particle_; i++){


          // Keep only leptons

    		  if ((abs(Particle_PID[i])==11 || abs(Particle_PID[i])==13) && Particle_Status[i]==1)
          {


            // Those who come from the same mother that created b and thos who come from W are unuseful
            if (abs(Particle_PID[Particle_Mother1[i]])!=24 &&
               not (abs(Particle_PID[Particle_Mother1[i]])==6 &&
               abs(Particle_PID[i-1])==5))
               {
                 countUsefulPair++;

              // Here place all you want to do with the lepton pair

                if (countUsefulPair==1) {
                pLepton1.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
                natureLeptonZ=abs(Particle_PID[i]);
                }
                if (countUsefulPair==2) pLepton2.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);


            }
            else{
              if (!countUnUsufulPair) pLepton3.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
              if (countUnUsufulPair) pLepton4.SetPxPyPzE(Particle_Px[i], Particle_Py[i], Particle_Pz[i], Particle_E[i]);
              countUnUsufulPair++;
            }
          }

    	  }


      }
    //###############################
    //Part Two : Blind to particle history
    //###############################
if (part=="pt"){
    // We find with a very ugly way the two higher momentum

    // We look for the two highest Pt-lepton-sum of the same nature
    int tabPosPartner[nbLepton]={0};
    int tabPosCounterPartner[nbLepton]={0};
    int nbPartner=0;
    int nbCounterPartner=0;
    float ptMax=0;

    for(int i=0; i<nbLepton; i++){
      // Set the lepton as the one that's interesting
      pLepton1=pLepton[i];
      natureLeptonZ=abs(tabLeptonID[i]);
      mother1Lepton1=tabMother1Lepton[i];
      mother2Lepton1=tabMother2Lepton[i];

      // This lepton has a partner ?
        // If it has one, does it have a twin ?
          // If yes, we take the higher pt
      for(int j=i; j<nbLepton; j++){
        if(tabLeptonID[i]==tabLeptonID[j]){
          if (pLepton1.Pt()<pLepton[j].Pt()) {
            pLepton1=pLepton[j];
            natureLeptonZ=abs(tabLeptonID[j]);
            mother1Lepton1=tabMother1Lepton[j];
            mother2Lepton1=tabMother2Lepton[j];
          }
        }
        if(tabLeptonID[i]!=-tabLeptonID[j]) continue;
        tabPosPartner[nbPartner]=i;
        nbPartner++;
      }
      // If not, we look at another one
      if(nbPartner==0) continue;

      // If exactly one, we look for the other two
      if(nbPartner==1){
        for(int k=0; k<nbLepton; k++){
          if(abs(tabLeptonID[k])==abs(tabLeptonID[i])) continue;
          tabPosCounterPartner[nbCounterPartner]=k;
          nbCounterPartner++;
        }
        // One counter Partner means we had two first-particle and aleady made the selection
          // One partner also mean that the partner is the first one
        if(nbCounterPartner==0){
          pLepton2=pLepton[tabPosPartner[0]];
          mother1Lepton2=tabMother1Lepton[tabPosPartner[0]];
          mother2Lepton2=tabMother2Lepton[tabPosPartner[0]];
        }
        if(nbCounterPartner==1){
          pLepton2=pLepton[tabPosPartner[0]];
          mother1Lepton2=tabMother1Lepton[tabPosPartner[0]];
          mother2Lepton2=tabMother2Lepton[tabPosPartner[0]];
        }
        // If 2 counterpartner, they may be in pair. (We could have look for that earlier but it's clearer this way)
          // If they are, we sum the pt and take the higher
          // Else, it's the first pair
        if(nbCounterPartner==2){
          if(tabLeptonID[tabPosCounterPartner[0]]==-tabLeptonID[tabPosCounterPartner[1]]){
            if(pLepton1.Pt()+pLepton[tabPosPartner[0]].Pt() > pLepton[tabPosCounterPartner[0]].Pt()+pLepton[tabPosCounterPartner[1]].Pt()){
              pLepton2=pLepton[tabPosPartner[0]];
              mother1Lepton2=tabMother1Lepton[tabPosPartner[0]];
              mother2Lepton2=tabMother2Lepton[tabPosPartner[0]];
            }
            else {
              pLepton1=pLepton[tabPosCounterPartner[0]];
              pLepton2=pLepton[tabPosCounterPartner[1]];
              natureLeptonZ=abs(tabLeptonID[tabPosCounterPartner[0]]);
              mother1Lepton1=tabMother1Lepton[tabPosCounterPartner[0]];
              mother2Lepton1=tabMother2Lepton[tabPosCounterPartner[0]];
              mother1Lepton2=tabMother1Lepton[tabPosCounterPartner[1]];
              mother2Lepton2=tabMother2Lepton[tabPosCounterPartner[1]];
            }
          }
        }
      }

      // If 2 or 3 partners, we take the higher-pt one
      else{
        pMax=-1;
        for(int k=0; k<nbPartner; k++){
          if(pLepton[tabPosPartner[k]].Pt()>pMax) {
            pMax=pLepton[tabPosPartner[k]].Pt();
            pLepton2=pLepton[tabPosPartner[k]];
            mother1Lepton2=tabMother1Lepton[tabPosPartner[k]];
            mother2Lepton2=tabMother2Lepton[tabPosPartner[k]];
          }
        }
      }
    }
    if (pLepton1.Pt()==pLepton2.Pt()) countSame++;
    // End of the selection of the higher-Pt couple

    if (mother1Lepton1!=mother1Lepton2 || mother2Lepton1!=mother2Lepton2) {
      countDiffMother++;
    }

    // Take the invMass for the others leptons
    // We take the Pts that aren't alredy used
    bool use=0;
    for(int i=0; i<nbLepton; i++){
      if (pLepton[i]==pLepton1 || pLepton[i]==pLepton2) continue;
      if (use){
        pLepton4=pLepton[i];
        pLepton4ID=tabLeptonID[i];
      }
      if (!use){
        pLepton3=pLepton[i];
        pLepton3ID=tabLeptonID[i];
        use=1;
      }

      natureLeptonW=abs(tabLeptonID[i]);
      tOutput_W->Fill();
    }
  }

    //###################################
    // Part 3 : Analyze
    //###################################



    // Part for Z
    pLeptonPair = pLepton1 + pLepton2;
    invMassLepton_Z=pLeptonPair.M();
    ptLepton1=pLepton1.Pt();
    ptLepton2=pLepton2.Pt();
    sumPtLeptonScal_Z=(ptLepton1+ptLepton2);
    sumPtLeptonVect_Z=pLeptonPair.Pt();

    // Part for W
    pLeptonPairW = pLepton3 + pLepton4;
    invMassLepton_W=pLeptonPairW.M();
    ptLepton3=pLepton3.Pt();
    ptLepton4=pLepton4.Pt();
    sumPtLeptonScal_W=(ptLepton3+ptLepton4);
    sumPtLeptonVect_W=pLeptonPairW.Pt();

    // Part for B
    pBPair = pB1+pB2;
    pLeptonB=pBPair + pLeptonPairW + pLeptonPair;
    pLeptonTot=pLeptonPair+pLeptonPairW;
    invMassB = pBPair.M();
    invMass4Lepton2B=pLeptonB.M();
    invMass4Lepton=pLeptonTot.M();



    // Part for B + lepton W

    invMass2LeptonW2B=(pLeptonPairW+pBPair).M();


    // Angle
    angleBtwLeptonZ=pLepton1.Angle(pLepton2.Vect());
    deltaRLeptonZ=pLepton1.DeltaR(pLepton2);

    deltaRLeptonZm=-99;
    deltaRLeptonZp=-99;
    deltaRLeptonZr=-99;

    if (invMassLepton_Z<76) deltaRLeptonZm=deltaRLeptonZ;
    if (invMassLepton_Z>106) deltaRLeptonZp=deltaRLeptonZ;
    else deltaRLeptonZr=deltaRLeptonZ;

    deltaPhiLeptonZ=pLepton1.DeltaPhi(pLepton2);

    TVector3 invBoost(-pLeptonPair.Px()/pLeptonPair.E(),-pLeptonPair.Py()/pLeptonPair.E(),-pLeptonPair.Pz()/pLeptonPair.E());

    pLepton1.Boost(invBoost);
    pLepton2.Boost(invBoost);

    TVector3 Zdir=pLeptonPair.Vect().Unit();
    TVector3 Ydir(-(pLeptonPair.Py()+pLeptonPair.Pz())/pLeptonPair.Px(), 1, 1);
    Ydir=Ydir.Unit();
    TVector3 Xdir(Ydir.Cross(Zdir));

    thetaStar1=TMath::ACos((pLepton1.Vect().Unit()).Dot(Zdir));
    thetaStar2=TMath::ACos((pLepton2.Vect().Unit()).Dot(Zdir));




    // ##################################
    // Fill
    //###################################

    tOutput->Fill();
    if (natureLeptonZ==11) tOutput_E->Fill();
    if (natureLeptonZ==13) tOutput_V->Fill();
   }
   cout<<"Lepton from different mother : "<<countDiffMother<<endl;
   cout<<"Same lepton : "<<countSame<<endl;



   //####################################
   //Part 3 : Angle Reconstruction
   //####################################
   /*
   */




//#####END#######
   tOutput->Write();
   tOutput_E->Write();
   tOutput_V->Write();
   tOutput_W->Write();
   fOutput->Close();
}
// Should have create an object lepton to keep all the properties
//Or a dictionary with the position of the lepton !

// Fill two times the tree if same lepton
