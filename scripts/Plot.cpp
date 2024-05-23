#include "../include/Plot.hpp"
#define nbVariableMass 12
#define nbVariableAngle 8

int main ()
{
  string Path = "data/madgraph/output/";
  string variableMass[nbVariableMass]={"invMassLepton_Z", "invMassLepton_W", "sumPtLeptonVect_Z", "sumPtLeptonVect_W", "sumPtLeptonScal_Z", "sumPtLeptonScal_W", "sumPtBVect", "sumPtBScal", "invMassB", "invMass4Lepton2B", "invMassDeltaR", "invMass2LeptonW2B"};
  string variableAngle[nbVariableAngle]={"angleBtwLeptonZ", "thetaStar1", "thetaStar2", "deltaRLeptonZ","deltaRLeptonZp","deltaRLeptonZm","deltaRLeptonZr", "deltaPhiLeptonZ"};


  //TTree* SM = FileReader(Path + "output_run01sm.root", "TreePt");
  //TTree* DIM6 = FileReader(Path + "output_EFT_CTL2_-5.root", "TreePt");

  TTree* SM_E = FileReader(Path + "output_run01sm.root", "TreeElec");
  TTree* DIM6_E = FileReader(Path + "output_EFT_CTL2_-5.root", "TreeElec");


  TTree* SM_V = FileReader(Path + "output_run01sm.root", "TreeMuon");
  TTree* DIM6_V = FileReader(Path + "output_EFT_CTL2_-5.root", "TreeMuon");



  // For-loop when I manage to put everything in the same tree

  //Compare_3Histos(SM,DIM6,DIM6_2,variable[0],100,0,500,"","","","M [GeV]","Events","legendUpRight","Models","SM","DIM6@(-2,-2)", "DIM6@(5,0)","eft","results/" + variable[0]);
  //Compare_3Histos(SM_E,DIM6_E,DIM6_E_2,variable[1],100,0,500,"","","","M [GeV]","Events","legendUpRight","Models","SM","DIM6@(-2,-2)", "DIM6@(5,0)","eft","results/" + variable[1]);
  //Compare_3Histos(SM_V,DIM6_V,DIM6_V_2,variable[2],100,0,500,"","","","M [GeV]","Events","legendUpRight","Models","SM","DIM6@(5,0)", "DIM6@(5,0)","eft","results/" + variable[2]);
  //Compare_3Histos(SM,DIM6,DIM6_2,variable[3],4,11,14,"","","","PID","Events","legendUpRight","Models","SM","DIM6@(-2,-2)","DIM6@(5,0)","eft","results/" + variable[3]);
  //Compare_4Histos(SM,DIM6,DIM6_2,DIM6_3,variable[3],4,11,14,"","","","","PID","Events","legendUpRight","Models","SM","DIM6@(5,5)","DIM6@(-2,-2)","DIM6@(5,0)","results/" + variable[3]);
  extern bool log_fity;
  log_fity=true;
  for(int i=0; i<nbVariableMass;i++){
    //Compare_2Histos(SM, DIM6, variableMass[i], 20,0,500,"","","M [GeV]","Events","legendUpRight","Models","SM","DIM6 ctl2=-5", "results/" + variableMass[i] + "duo");
    Compare_2Histos(SM_V, DIM6_V, variableMass[i], 20,0,500,"","","M [GeV]","Events","legendUpRight","Models","SM","DIM6 ctl2=-5", "results/MADlog/" + variableMass[i] + "duoV");


  }

  extern bool flows;
  log_fity=false;
  flows=true;
  for (int i=0;i<nbVariableAngle;i++){
  //Compare_2Histos(SM, DIM6, variableAngle[i], 50,0,4,"","","Angle","Events","legendUpRight","Models","SM","DIM6 ctl2=-5", "results/" + variableAngle[i]+ "duo");
  Compare_2Histos(SM_V, DIM6_V, variableAngle[i], 50,0,6,"","","Angle","Events","legendUpRight","Models","SM","DIM6 ctl2=-5", "results/MADnolog/" + variableAngle[i] + "duoV");

  }
  flows=true;
  for(int i=0; i<nbVariableMass;i++){
    //Compare_2Histos(SM, DIM6, variableMass[i], 10,0,500,"","","M [GeV]","Events","legendUpRight","Models","SM","DIM6 ctl2=-5", "results/" + variableMass[i] + "duo");
    Compare_2Histos(SM_V, DIM6_V, variableMass[i], 10,0,500,"","","M [GeV]","Events","legendUpRight","Models","SM","DIM6 ctl2=-5", "results/MADnolog/" + variableMass[i] + "duoV");


  }


  //Compare_2Histos(SM_V, DIM6_V, variable[2], 100,0,1000,"","","M [GeV]","Events","legendUpRight","Models","SM","DIM6@5&0", "results/" + variable[2] + "duo");
  return 0;
}
