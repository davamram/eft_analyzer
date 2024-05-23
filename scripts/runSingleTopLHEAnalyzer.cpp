
#include "../include/SingleTopLHEAnalyzer.hpp"
#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <sstream>
using namespace std;


int main(int argc, char **argv)
{
  extern string part;
  part=argv[1];
  int nSamples = 1;
  int modes = 1;
  string* suffix = new string[nSamples];
  string* prefix = new string[modes];

prefix[0]="run01";
//prefix[1]="madspin_";
//prefix[0]="EVENTS_t_channel_madspin_Rwt_ctwi_cbwi_";

suffix[0] = "sm";
//suffix[1] = "eft";
// suffix[2] = "cbw_p1";
// suffix[3] = "cbw_p2";
// suffix[4] = "ctw_p1";
// suffix[5] = "ctw_p2";
// suffix[6] = "ctw_m1";
// suffix[7] = "ctw_m2";
// suffix[8] = "cptb_m10";
// suffix[9] = "cptb_m5";
// suffix[10] = "cptb_p5";
// suffix[11] = "cptb_p10";
// suffix[12] = "cptbi_m10";
// suffix[13] = "cptbi_m5";
// suffix[14] = "cptbi_p5";
// suffix[15] = "cptbi_p10";



  TFile** fInput = new TFile*[nSamples];
  TTree** tInput = new TTree*[nSamples];
  SingleTopLHEAnalyzer** singleTopLHEAnalyzer = new SingleTopLHEAnalyzer*[nSamples];

  for (int j=0; j<modes; j++)
  {
    for (int i=0; i<nSamples; i++)
    {
      //cout<<"Cooking "+prefix[j] + suffix[i]+".root"<<endl;

     ///////////////////////////Choose Input Path files//////////////////////////////////////
     string inputPath = "data/madgraph/";
     //string inputName = inputPath + prefix[j] + suffix[i] + ".root";
     //string inputName ="EFT_CTL2_-5.root";
     string inputName ="run01sm.root";

    // string outputName = "output_EFT_CTL2_-5.root";
     string outputName = "output_run01sm.root";
    
     fInput[i] = new TFile((inputPath+inputName).c_str(),"READ");
     tInput[i] = (TTree*) fInput[i]->Get("LHEF");
     singleTopLHEAnalyzer[i] = new SingleTopLHEAnalyzer(tInput[i]);
     singleTopLHEAnalyzer[i]->Loop(part);

     ///////////////////////////Choose Output Path files//////////////////////////////////////
     string outputPath = "data/madgraph/output/";
     string commandline = "mv output.root " + outputPath + outputName;
     cout<<"Creating "<<outputPath + outputName<<endl;

     system(commandline.c_str());
    }
  }



return 0;
}
