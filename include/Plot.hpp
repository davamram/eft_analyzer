#include <TFile.h>
#include <TTree.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <sstream>
#include "TStyle.h"
#include <iostream>
#include <TStyle.h>
#include <string>
#include <TGraph.h>
#include <TF1.h>
#include <fstream>
#include <TSystem.h>
using namespace std;

bool reweight = false;

bool normalization = true; //Turn on normalization to 1
bool log_fity = true; //Turn on log scale in Y axis
bool diff_xsection = false; //Change this to have histograms normalized to Xsection
bool flows = true; //Change this to Active/Deactivate Overflow and Underflow for all Histograms

TTree* FileReader(string inputPath, const char* treeName)
{
  // Fct to read ROOT files
  TFile* RootFile = new TFile(inputPath.c_str(),"READ");
  TTree* Tree;
  Tree = (TTree*) RootFile->Get(treeName);
  return Tree;
}

TH1D* GetHistoWeight(TTree* t, string variable, int nbins, double xmin, double xmax, string cut, string name)
{
        string sxmin, sxmax, snbins;
        stringstream ss[3];

        ss[0] << xmin;
        ss[0] >> sxmin;
        ss[1] << xmax;
        ss[1] >> sxmax;
        ss[2] << nbins;
        ss[2] >> snbins;

        string variablenew = variable + " >> h(" + snbins + "," + sxmin + "," + sxmax + ")";

        string cutnew = "1 * (" + cut + ")";
        //string cutnew = "(" + cut + ")";

//

        t->Draw(variablenew.c_str());
        TH1D *histo = (TH1D*)gDirectory->Get("h");
  if(flows)
  {
		if (histo->GetEntries()==0) return histo;

		double underflow = histo->GetBinContent(0);
		//cout << "underflow="<<underflow<<endl;
		double val = 0;
		if (underflow>0) {
			val = histo->GetBinContent(1);
      //cout<<"val= "<<val<<endl;
			histo->SetBinContent(1, val+underflow);
			 histo->SetBinContent(0, 0);
		}
		double overflow = histo->GetBinContent(nbins+1);
    //cout<<"overflow= "<<overflow<<endl;
		if (overflow>0) {
		  val = histo->GetBinContent(nbins);
		  histo->SetBinContent(nbins+1, 0);
		  histo->SetBinContent(nbins, val+overflow);
		}
  }
  //cout << "Area="<<histo->Integral()<<endl;
	//cout << "Nevents="<<histo->GetEntries()<<endl;
        histo->SetName(name.c_str());
        histo->SetTitle(name.c_str());

        return histo;
}

void Ratio_EFT_SM(TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, string variable, string EFT, int nbins, double xmin, double xmax, string selection1, string selection2, string selection3, string selection4, string selection5, string legendX, string legendY, string file_name, string lepton)
{
  TFile* file_output = new TFile(file_name.c_str(),"RECREATE");
  TTree* tree_file = new TTree("events","events");

  tree_file->Branch("");

  TH1D* Histo_SM = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection1, "Histo_SM");
  TH1D* Histo_EFT_p1 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection2, "Histo_EFT_p1");
  TH1D* Histo_EFT_p2 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection3, "Histo_EFT_p2");
  TH1D* Histo_EFT_m1 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection4, "Histo_EFT_m1");
  TH1D* Histo_EFT_m2 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection5, "Histo_EFT_m2");

  TH1D* Ratio_p2_SM = (TH1D*)Histo_EFT_p2->Clone("Ratio_p2_SM");
  TH1D* Ratio_p1_SM = (TH1D*)Histo_EFT_p1->Clone("Ratio_p1_SM");
  TH1D* Ratio_m1_SM = (TH1D*)Histo_EFT_m1->Clone("Ratio_m1_SM");
  TH1D* Ratio_m2_SM = (TH1D*)Histo_EFT_m2->Clone("Ratio_m2_SM");


  //X_sections values obtained with MG5
  double Xs_n5, Xs_n2, Xs_n1 ,Xs_SM = 35.2863;
  double Xs_p1, Xs_p2, Xs_p5;

  if(EFT == "cbwi")
  {
    Xs_n5 = 48.5315;
    Xs_n2 = 37.3699;
    Xs_n1 = 35.8251;
    Xs_p1 = 35.8091;
    Xs_p2 = 37.3955;
    Xs_p5 = 48.5184;
  }
  if(EFT == "ctwi")
  {
    Xs_n5 = 48.5683;
    Xs_n2 = 37.3838;
    Xs_n1 = 35.8029;
    Xs_p1 = 35.8126;
    Xs_p2 = 37.3893;
    Xs_p5 = 48.5945;
  }

  int N = 1000000; //sum des pts initiaux
  Ratio_p2_SM->Scale(Xs_p2/N);
  Ratio_p1_SM->Scale(Xs_p1/N);
  Ratio_m2_SM->Scale(Xs_n2/N);
  Ratio_m1_SM->Scale(Xs_n1/N);
  Histo_SM->Scale(Xs_SM/N);

  //weights_sum[5] = {SM, p2, p1, m2, m1}

  //double weights_sum[5] = {77672, 82091.3, 79047.4, 82123.4, 79057.2}; //ctwi from ctwi2p5
  //double weights_sum[5] = {77672, 81820.2, 77911.8, 81826.8, 77912.2}; //cbwi from ctwi2p5
  //nb_weights = 230054;

  //double weights_sum[5] = {83125.7, 87952.1, 84560.8, 87933.9, 84570.7}; //ctwi from cbwi_p3
  //double weights_sum[5] = {83125.7, 84795.9, 83366.5, 84795.8, 83368.8}; //cbwi from cbwi_p3
  //nb_weights = 241499;

  //double weights_sum[5] = {73926.8, 77781.1, 74899.4, 77732.2, 74872.5}; //ctwi from cbwi2p5
  //double weights_sum[5] = {73926.8, 77618.7, 74852, 77619.4, 74852.4}; //cbwi from cbwi2p5
  //int nb_weights = 226646;

  //Histo_SM->Scale(weights_sum[0]/nb_weights);
  //Ratio_p2_SM->Scale(weights_sum[1]/nb_weights);
  //Ratio_p1_SM->Scale(weights_sum[2]/nb_weights);
  //Ratio_m2_SM->Scale(weights_sum[3]/nb_weights);
  //Ratio_m1_SM->Scale(weights_sum[4]/nb_weights);


  Ratio_p2_SM->Divide(Histo_SM);
  Ratio_p1_SM->Divide(Histo_SM);
  Ratio_m2_SM->Divide(Histo_SM);
  Ratio_m1_SM->Divide(Histo_SM);
  Histo_SM->Divide(Histo_SM);



  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  //double max = (Ratio_p2_SM->GetMaximum()>Ratio_p1_SM->GetMaximum()) ? Ratio_p2_SM->GetMaximum() : Ratio_p1_SM->GetMaximum();
  double max = Ratio_p2_SM->GetMaximum();
  //Ratio_p2_SM->SetAxisRange((2-max)*1.1, max*1.25, "Y");
  Ratio_p2_SM->SetAxisRange(2-max, max*1.1, "Y");
  Ratio_p2_SM->SetXTitle(legendX.c_str());
  Ratio_p2_SM->SetYTitle(legendY.c_str());
  Ratio_p2_SM->SetLineColor(kRed);
  Ratio_p2_SM->SetLineWidth(2);
  Ratio_p2_SM->Draw("");

  Ratio_p1_SM->SetLineColor(kBlue);
  Ratio_p1_SM->SetLineWidth(2);
  Ratio_p1_SM->Draw("SAME");

  Ratio_m2_SM->SetLineColor(kOrange);
  Ratio_m2_SM->SetLineWidth(2);
  Ratio_m2_SM->Draw("SAME");

  Ratio_m1_SM->SetLineColor(kGreen);
  Ratio_m1_SM->SetLineWidth(2);
  Ratio_m1_SM->Draw("SAME");

  Histo_SM->SetLineColor(kBlack);
  Histo_SM->SetLineWidth(2);
  Histo_SM->Draw("SAME");

  double lx0 = 0.6;
  double ly0 = 0.6;
  double lx1 = 0.99;
  double ly1 = 0.99;

  string legendtitle = "Value of the EFT";

  string eft_p2_legend = EFT + "/#\\Lambda^{2} = 2 (TeV^{-2})";
  string eft_p1_legend = EFT + "/#\\Lambda^{2} = 1 (TeV^{-2})";
  string eft_m1_legend = EFT + "/#\\Lambda^{2} = -1 (TeV^{-2})";
  string eft_m2_legend = EFT + "/#\\Lambda^{2} = -2 (TeV^{-2})";
  string SM_legend = "SM";

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.05);
  legend->AddEntry(Ratio_p2_SM->GetName(), eft_p2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p1_SM->GetName(), eft_p1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m1_SM->GetName(), eft_m1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m2_SM->GetName(), eft_m2_legend.c_str(), "l");
  legend->AddEntry(Histo_SM->GetName(), SM_legend.c_str(), "l");

  legend->Draw("SAME");
  //Canvas->Print(file_name.c_str());


  //string TFile_name = "/results/ratio_madgraph/signal_proc_"+variable;
  //TFile* ratio_file = new TFile(("signal_proc_" + variable + "_" + EFT + "_" + lepton + "_TF1.root").c_str(), "RECREATE");
  TFile* ratio_file = new TFile(("results/EFT_vs_SM/Rwgt_cbwi2p5_XSect/signal_proc_" + variable + "_" + EFT + "_" + lepton + "_"+ to_string(nbins) + "bins" + "_Rwgt_cbwi2p5.root").c_str(), "RECREATE");
  ratio_file->cd();

  TGraph** ratio_Histo = new TGraph*[nbins];

  for (int i = 1 ; i<=nbins ; i++)
  {
    string number_plot = to_string(i);
    string file_name_eft = file_name;
    file_name_eft.insert(file_name.size()-4,"_"+number_plot);
    Canvas->Clear();
    //Canvas->SetLogy();
    ratio_Histo[i] = new TGraph(5);
    TF1* ratio_formula = new TF1(("bin_content_"+number_plot).c_str(),"[0]+[1]*x+[2]*x*x", -25 , 25);



    ratio_Histo[i]->SetPoint(0,-2,Ratio_m2_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(1,-1,Ratio_m1_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(2,0,Histo_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(3,1,Ratio_p1_SM->GetBinContent(i));
    ratio_Histo[i]->SetPoint(4,2,Ratio_p2_SM->GetBinContent(i));

    ratio_Histo[i]->SetMarkerStyle(kStar);

    tree_file->Fill();
    ratio_Histo[i]->Fit(ratio_formula);
    ratio_Histo[i]->Draw();
    //Canvas->Print(("test.pdf"+number_plot).c_str());


    TLegend* legend2 = new TLegend(0.6, 0.7, 0.89, 0.89, "");
    legend2->SetTextSize(0.035);
    if(variable == "PhiStar" && nbins==5)
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "0 < #phi^{*} < 1.25" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "1.25 < #phi^{*} < 2.5" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "2.5 < #phi^{*} < 4" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "4 < #phi^{*} < 5.5" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "5.5 < #phi^{*} < 6.28" ,"r");
    }

    if(variable == "cosThetaStar" && nbins==5)
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "-1 < cos(#theta^{*}) < -0.6" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.6 < cos(#theta^{*}) < -0.2" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.2 < cos(#theta^{*}) < 0.2" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.2 < cos(#theta^{*}) < 0.6" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.6 < cos(#theta^{*}) < 1" ,"r");

    }

    legend2->Draw("SAME");


    ratio_Histo[i]->GetYaxis()->SetRangeUser(0,Ratio_p2_SM->GetBinContent(i)*2);
    ratio_Histo[i]->GetYaxis()->SetTitle("EFT/SM");
    if(EFT == "ctwi") ratio_Histo[i]->GetXaxis()->SetTitle("C_{tW}^{I}");
    if(EFT == "cbwi") ratio_Histo[i]->GetXaxis()->SetTitle("C_{bW}^{I}");

    ratio_Histo[i]->SetTitle("");



    Canvas->Print(file_name_eft.c_str());

    //ratio_Histo[i]->Write(("bin_content_part_"+number_plot).c_str());

    //Canvas->Print(file_name_eft.c_str());

    ratio_formula->Write();
  }


  file_output->Write();
  ratio_file->Close();
}

void Compare_1Histos(TTree* t1, string variable, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "Histo_1");
  Histo_1->SetStats(kFALSE);

  if(normalization)
  {
    double a = Histo_1->Integral();
    Histo_1->Scale(1/a);
    Name += "_normalized";
    legendY += " normalized";
  }

  double max = Histo_1->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  double lx0, ly0, lx1, ly1;
   if (legendPlace=="legendUpLeft"){
 	 lx0 = 0.2;
 	 ly0 = 0.75;
 	 lx1 = 0.5;
 	 ly1 = 0.95;
   }
    if (legendPlace=="legendUpRight"){
 	 lx0 = 0.75;
 	 ly0 = 0.75;
 	 lx1 = 0.99;
 	 ly1 = 0.99;
   }

   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");

   legend->Draw("SAME");

   Canvas->Print(Name.c_str());
}

void Compare_2Histos(TTree* t1, TTree* t2, string variable, int nbins, double xmin, double xmax, string selection1, string selection2, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection1, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection2, "Histo_2");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle(variable.c_str());
  if(log_fity)
  {
    Canvas->SetLogy();
    Name += "_log";
    legendY += " log";
  }
  if(normalization)
  {
    double a = Histo_1->Integral();
    double b = Histo_2->Integral();
    Histo_1->Scale(1/a);
    Histo_2->Scale(1/b);
    Name += "_normalized";
    legendY += " normalized";
  }
  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  Histo_1->SetMaximum(max*1.3);
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
  legend->Draw("SAME");

  legend->SetTextSize(0.035);


  Name += ".pdf";
  Canvas->Print(Name.c_str());

  // cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  // cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;

}

void Compare_3Histos(TTree* t1, TTree* t2, TTree* t3, string variable, int nbins, double xmin, double xmax, string selection1, string selection2, string selection3, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string Name){


  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection1, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection2, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection3, "Histo_3");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle(variable.c_str());

  if(log_fity)
  {
    Canvas->SetLogy();
    Name += "_log";
    legendY += " log";
  }
  if(normalization)
  {
    Double_t integral1 = 1/Histo_1->Integral();
    Histo_1->Scale(integral1);

    Double_t integral2 = 1/Histo_2->Integral();
    Histo_2->Scale(integral2);

    Double_t integral3 = 1/Histo_3->Integral();
    Histo_3->Scale(integral3);

    Name += "_normalized";
    legendY += " normalized";
  }

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  Histo_1->SetMaximum(max*1.3 );

  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();


  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");


  Histo_3->SetLineColor(kOrange);
  Histo_3->SetLineWidth(2);
  Histo_3->Draw("SAME");


 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.7;
	 lx1 = 0.6;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.8;
	 ly0 = 0.75; //Lower Y
	 lx1 = 0.99;
	 ly1 = 0.99; //Upper Y
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.04);

  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(),"l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(),"l");
  legend->AddEntry(Histo_3->GetName(), legendEntry3.c_str(),"l");
  //legend->SetLegendSize(0.5);

  legend->Draw("SAME");

  Name += ".pdf";
  Canvas->Print(Name.c_str());
}

void Compare_4Histos(TTree* t1, TTree* t2, TTree* t3, TTree* t4, string variable, int nbins, double xmin, double xmax, string selection1, string selection2, string selection3, string selection4, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string legendEntry4, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection1, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection2, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection3, "Histo_3");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection4, "Histo_4");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);
  if(normalization)
  {
    double a = Histo_1->Integral();
    double b = Histo_2->Integral();
    double c = Histo_3->Integral();
    double d = Histo_4->Integral();

    Histo_1->Scale(1/a);
    Histo_2->Scale(1/b);
    Histo_3->Scale(1/c);
    Histo_4->Scale(1/d);
    Name += "_normalized";
    legendY += " normalized";

  }
  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  //Canvas->SetLogx();
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  Histo_3->SetLineColor(kOrange);
  Histo_3->SetLineWidth(2);
  Histo_3->Draw("SAME");

  Histo_4->SetLineColor(kGreen);
  Histo_4->SetLineWidth(2);
  Histo_4->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
  legend->AddEntry(Histo_3->GetName(), legendEntry3.c_str(), "l");
  legend->AddEntry(Histo_4->GetName(), legendEntry4.c_str(), "l");
  legend->Draw("SAME");

  Name += ".pdf";
  Canvas->Print(Name.c_str());
}

void Compare_5Histos(TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, string variable, int nbins, double xmin, double xmax, string selection1, string selection2, string selection3,string selection4, string selection5, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string legendEntry4, string legendEntry5, string Name)
{

  TH1D* Histo_1 = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection1, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection2, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection3, "Histo_3");
  TH1D* Histo_4 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection4, "Histo_4");
  TH1D* Histo_5 = GetHistoWeight(t4, variable, nbins, xmin, xmax, selection5, "Histo_5");


  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);
  Histo_4->SetStats(kFALSE);
  Histo_5->SetStats(kFALSE);

  if(normalization)
  {
    double a = Histo_1->Integral();
    double b = Histo_2->Integral();
    double c = Histo_3->Integral();
    double d = Histo_4->Integral();
    double e = Histo_5->Integral();

    Histo_1->Scale(1/a);
    Histo_2->Scale(1/b);
    Histo_3->Scale(1/c);
    Histo_4->Scale(1/d);
    Histo_5->Scale(1/e);

    Name += "_normalized";
    legendY += " normalized";
  }
  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  //Canvas->SetLogx();
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*1.1,"Y");
  if (legendX == "Top Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.0001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  if (legendX == "W Mass (GeV)")
    {
      Histo_1->SetAxisRange(0.001,max*1.1,"Y");
      Canvas->SetLogy();
    }
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  Histo_3->SetLineColor(kOrange);
  Histo_3->SetLineWidth(2);
  Histo_3->Draw("SAME");

  Histo_4->SetLineColor(kGreen);
  Histo_4->SetLineWidth(2);
  Histo_4->Draw("SAME");

  Histo_5->SetLineColor(kViolet);
  Histo_5->SetLineWidth(2);
  Histo_5->Draw("SAME");

  double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
  legend->AddEntry(Histo_3->GetName(), legendEntry3.c_str(), "l");
  legend->AddEntry(Histo_4->GetName(), legendEntry4.c_str(), "l");
  legend->AddEntry(Histo_5->GetName(), legendEntry5.c_str(), "l");
  legend->Draw("SAME");

  Name += ".pdf";
  Canvas->Print(Name.c_str());


}

void Ratio_EFT_SM_7pts(TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, TTree* t6, TTree* t7, string variable, string EFT, int nbins, double xmin, double xmax, string selection, string legendX, string legendY, string file_name, string lepton)
{
  //Fcts to Plot EFT/SM for cbwi/ctwi = [-5;5]

  TFile* file_output = new TFile(file_name.c_str(),"RECREATE");
  TTree* tree_file = new TTree("events","events");

  tree_file->Branch("");

  TH1D* Histo_SM = GetHistoWeight(t1, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p1 = GetHistoWeight(t2, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p2 = GetHistoWeight(t3, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m1= GetHistoWeight(t4, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m2 = GetHistoWeight(t5, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_p5 = GetHistoWeight(t6, variable, nbins, xmin, xmax, selection, "");
  TH1D* Histo_EFT_m5 = GetHistoWeight(t7, variable, nbins, xmin, xmax, selection, "");

  TH1D* Ratio_p2_SM = (TH1D*)Histo_EFT_p2->Clone("Ratio_p2_SM");
  TH1D* Ratio_p1_SM = (TH1D*)Histo_EFT_p1->Clone("Ratio_p1_SM");
  TH1D* Ratio_m1_SM = (TH1D*)Histo_EFT_m1->Clone("Ratio_m1_SM");
  TH1D* Ratio_m2_SM = (TH1D*)Histo_EFT_m2->Clone("Ratio_m2_SM");
  TH1D* Ratio_p5_SM = (TH1D*)Histo_EFT_p5->Clone("Ratio_p5_SM");
  TH1D* Ratio_m5_SM = (TH1D*)Histo_EFT_m5->Clone("Ratio_m5_SM");


  //Ratio_p2_SM->Divide(Histo_SM);
  //Ratio_p1_SM->Divide(Histo_SM);
  //Ratio_m2_SM->Divide(Histo_SM);
  //Ratio_m1_SM->Divide(Histo_SM);
  //Ratio_p5_SM->Divide(Histo_SM);
  //Ratio_m5_SM->Divide(Histo_SM);

  bool div_Xsection = true;
  if(div_Xsection)
  {
    double Xs_n5, Xs_n2, Xs_n1 ,Xs_SM = 35.2863;
    double Xs_p1, Xs_p2, Xs_p5;

    if(EFT == "cbwi")
    {
      Xs_n5 = 48.5315;
      Xs_n2 = 37.3699;
      Xs_n1 = 35.8251;
      Xs_p1 = 35.8091;
      Xs_p2 = 37.3955;
      Xs_p5 = 48.5184;
    }
    if(EFT == "ctwi")
    {
      Xs_n5 = 48.5683;
      Xs_n2 = 37.3838;
      Xs_n1 = 35.8029;
      Xs_p1 = 35.8126;
      Xs_p2 = 37.3893;
      Xs_p5 = 48.5945;
    }
    int N = 1000000;
    Ratio_p2_SM->Scale(Xs_p2/N);
    Ratio_p1_SM->Scale(Xs_p1/N);
    Ratio_m2_SM->Scale(Xs_n2/N);
    Ratio_m1_SM->Scale(Xs_n1/N);
    Ratio_p5_SM->Scale(Xs_p5/N);
    Ratio_m5_SM->Scale(Xs_n5/N);
    Histo_SM->Scale(Xs_SM/N);
  }

  bool ratio = true;
  if(ratio)
  {
    Ratio_p2_SM->Divide(Histo_SM);
    Ratio_p1_SM->Divide(Histo_SM);
    Ratio_m2_SM->Divide(Histo_SM);
    Ratio_m1_SM->Divide(Histo_SM);
    Ratio_p5_SM->Divide(Histo_SM);
    Ratio_m5_SM->Divide(Histo_SM);
    Histo_SM->Divide(Histo_SM);
  }
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");

  Ratio_p2_SM->SetLineColor(kRed);
  Ratio_p2_SM->SetLineWidth(2);

  Ratio_p1_SM->SetLineColor(kBlue);
  Ratio_p1_SM->SetLineWidth(2);

  Ratio_m2_SM->SetLineColor(kOrange);
  Ratio_m2_SM->SetLineWidth(2);

  Ratio_m1_SM->SetLineColor(kBlack);
  Ratio_m1_SM->SetLineWidth(2);

  Ratio_p5_SM->SetLineColor(kGreen);
  Ratio_p5_SM->SetLineWidth(2);

  Ratio_m5_SM->SetLineColor(kViolet);
  Ratio_m5_SM->SetLineWidth(2);


  Ratio_p2_SM->SetXTitle(legendX.c_str());
  Ratio_p2_SM->SetYTitle(legendY.c_str());

  double max = (Ratio_p5_SM->GetMaximum()>Ratio_m5_SM->GetMaximum()) ? Ratio_p5_SM->GetMaximum() : Ratio_m5_SM->GetMaximum();

  Ratio_p2_SM->SetAxisRange((2-max)*1.1, max*1.25, "Y");
  Ratio_p2_SM->Draw("");
  Ratio_p1_SM->Draw("SAME");
  Ratio_m1_SM->Draw("SAME");
  Ratio_m2_SM->Draw("SAME");
  Ratio_p5_SM->Draw("SAME");
  Ratio_m5_SM->Draw("SAME");

  double lx0 = 0.6;
  double ly0 = 0.6;
  double lx1 = 0.99;
  double ly1 = 0.99;
  string legendtitle = "Value of the EFT";

  string eft_p2_legend = EFT + "/#\\Lambda^{2} = 2 (TeV^{-2})";
  string eft_p1_legend = EFT + "/#\\Lambda^{2} = 1 (TeV^{-2})";
  string eft_m1_legend = EFT + "/#\\Lambda^{2} = -1 (TeV^{-2})";
  string eft_m2_legend = EFT + "/#\\Lambda^{2} = -2 (TeV^{-2})";
  string eft_p5_legend = EFT + "/#\\Lambda^{2} = 5 (TeV^{-2})";
  string eft_m5_legend = EFT + "/#\\Lambda^{2} = -5 (TeV^{-2})";



  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.05);
  legend->AddEntry(Ratio_p2_SM->GetName(), eft_p2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p1_SM->GetName(), eft_p1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m1_SM->GetName(), eft_m1_legend.c_str(), "l");
  legend->AddEntry(Ratio_m2_SM->GetName(), eft_m2_legend.c_str(), "l");
  legend->AddEntry(Ratio_p5_SM->GetName(), eft_p5_legend.c_str(), "l");
  legend->AddEntry(Ratio_m5_SM->GetName(), eft_m5_legend.c_str(), "l");

  legend->Draw("SAME");

  TGraph** ratio_Histo = new TGraph*[nbins];

    //string TFile_name = "/results/ratio_madgraph/signal_proc_"+variable;
  TFile* ratio_file = new TFile(("signal_proc_"+variable+"_"+EFT+"_"+lepton+"_TF1.root").c_str(),"RECREATE");
  ratio_file->cd();
  Canvas->Print(file_name.c_str());


  for (int i = 1 ; i<=nbins ; i++)
  {
    string number_plot = to_string(i); //Conversion of an int into a stream
    string file_name_eft = file_name;
    file_name_eft.insert(file_name.size()-4,"_"+number_plot);
    Canvas->Clear();
    //Canvas->SetLogy();
    ratio_Histo[i] = new TGraph(5);
    TF1* ratio_formula = new TF1(("bin_content_par1_"+number_plot).c_str(),"[0]+[1]*x+[2]*x*x", -3 , 3);

    bool Get_bins = false;
    if(Get_bins)
    {
      ratio_Histo[i]->SetPoint(0,-2,Histo_EFT_m2->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(1,-1,Histo_EFT_m1->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(2,0,Histo_SM->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(3,1,Histo_EFT_p1->GetBinContent(i)/Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(4,2,Histo_EFT_p2->GetBinContent(i)/Histo_SM->GetBinContent(i));
    }
    else
    {
      ratio_Histo[i]->SetPoint(0,-5,Ratio_m5_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(1,-2,Ratio_m2_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(2,-1,Ratio_m1_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(3,0,Histo_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(4,1,Ratio_p1_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(5,2,Ratio_p2_SM->GetBinContent(i));
      ratio_Histo[i]->SetPoint(6,5,Ratio_p5_SM->GetBinContent(i));
    }
    ratio_Histo[i]->SetMarkerStyle(kStar);

    tree_file->Fill();
    ratio_Histo[i]->Fit(ratio_formula);
    ratio_Histo[i]->Draw();
    //Canvas->Print(("test.pdf"+number_plot).c_str());

      TLegend* legend2 = new TLegend(0.49, 0.73, 0.75, 0.86, "");
      legend2->SetTextSize(0.035);
    if(variable == "PhiStar")
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "0 < #phi^{*} < 1.25" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "1.25 < #phi^{*} < 2.5" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "2.5 < #phi^{*} < 4" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "4 < #phi^{*} < 5.5" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "5.5 < #phi^{*} < 6.28" ,"r");
    }

    if(variable == "cosThetaStar")
    {
      if(i == 1) legend2->AddEntry(ratio_Histo[i]->GetName(), "-1 < cos(#theta^{*}) < -0.6" ,"r");
      if(i == 2) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.6 < cos(#theta^{*}) < -0.2" ,"r");
      if(i == 3) legend2->AddEntry(ratio_Histo[i]->GetName(), "-0.2 < cos(#theta^{*}) < 0.2" ,"r");
      if(i == 4) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.2 < cos(#theta^{*}) < 0.6" ,"r");
      if(i == 5) legend2->AddEntry(ratio_Histo[i]->GetName(), "0.6 < cos(#theta^{*}) < 1" ,"r");

    }

    legend2->Draw("SAME");

    ratio_Histo[i]->GetYaxis()->SetRangeUser(0,Ratio_p5_SM->GetBinContent(i)*2);
    ratio_Histo[i]->GetYaxis()->SetTitle("EFT/SM");
    if(EFT == "ctwi") ratio_Histo[i]->GetXaxis()->SetTitle("ctwi");
    if(EFT == "cbwi") ratio_Histo[i]->GetXaxis()->SetTitle("cbwi");
    ratio_Histo[i]->SetTitle("");
    Canvas->Print(file_name_eft.c_str());

    ratio_Histo[i]->Write(("bin_content_par1_"+number_plot).c_str());
    //ratio_formula->Write();
  }


  //file_output->Write();
  ratio_file->Close();




}

void Compare_9Histos(int nbins, double xmin, double xmax, string selection1, string selection2, string selection3,string selection4, string selection5, string selection6, string selection7, string selection8, string selection9, TTree* t1, TTree* t2, TTree* t3, TTree* t4, TTree* t5, TTree* t6, TTree* t7, TTree* t8, TTree* t9, string variable1, string variable2, string variable3, string variable4, string variable5, string variable6, string variable7, string variable8, string variable9, string legendX, string legendY, string legendPlace, string legendtitle, string Name){

  //Fcts to Plot all weights in the same Canvas

  TH1D *Histo_1 = GetHistoWeight(t1, variable1, nbins, xmin, xmax, selection1, "Histo_1");
  TH1D *Histo_2 = GetHistoWeight(t2, variable2, nbins, xmin, xmax, selection2, "Histo_2");
  TH1D *Histo_3 = GetHistoWeight(t3, variable3, nbins, xmin, xmax, selection3, "Histo_3");
  TH1D *Histo_4 = GetHistoWeight(t4, variable4, nbins, xmin, xmax, selection4, "Histo_4");
  TH1D *Histo_5 = GetHistoWeight(t5, variable5, nbins, xmin, xmax, selection5, "Histo_5");
  TH1D *Histo_6 = GetHistoWeight(t6, variable6, nbins, xmin, xmax, selection6, "Histo_6");
  TH1D *Histo_7 = GetHistoWeight(t7, variable7, nbins, xmin, xmax, selection7, "Histo_7");
  TH1D *Histo_8 = GetHistoWeight(t8, variable8, nbins, xmin, xmax, selection8, "Histo_8");
  TH1D *Histo_9 = GetHistoWeight(t9, variable9, nbins, xmin, xmax, selection9, "Histo_9");


  double max = Histo_1->GetMaximum();
  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle("");
  Histo_1->SetAxisRange(0,max*2.2,"Y");
  //Histo_SM->SetAxisRange(-0.3,0.5,"X");
  //Histo_SM->SetAxisRange(-10,10,"X");
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

  Histo_3->SetLineColor(kGreen);
  Histo_3->SetLineWidth(2);
  Histo_3->Draw("SAME");

  Histo_4->SetLineColor(kViolet);
  Histo_4->SetLineWidth(2);
  Histo_4->Draw("SAME");

  Histo_5->SetLineColor(kOrange);
  Histo_5->SetLineWidth(2);
  Histo_5->Draw("SAME");

  Histo_6->SetLineColor(kBlue);
  Histo_6->SetLineStyle(3);
  Histo_6->SetLineWidth(2);
  Histo_6->Draw("SAME");

  Histo_7->SetLineColor(kGreen);
  Histo_7->SetLineStyle(3);
  Histo_7->SetLineWidth(2);
  Histo_7->Draw("SAME");

  Histo_8->SetLineColor(kViolet);
  Histo_8->SetLineStyle(3);
  Histo_8->SetLineWidth(2);
  Histo_8->Draw("SAME");

  Histo_9->SetLineColor(kOrange);
  Histo_9->SetLineStyle(3);
  Histo_9->SetLineWidth(2);
  Histo_9->Draw("SAME");


  double lx0, ly0, lx1, ly1;
   if (legendPlace=="legendUpLeft"){
 	 lx0 = 0.2;
 	 ly0 = 0.75;
 	 lx1 = 0.5;
 	 ly1 = 0.95;
   }
    if (legendPlace=="legendUpRight"){
 	 lx0 = 0.75;
 	 ly0 = 0.60;
 	 lx1 = 0.99;
 	 ly1 = 0.99;
   }


   TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
   legend->SetFillColor(kWhite);
   legend->AddEntry(Histo_1->GetName(), "SM", "l");
   legend->AddEntry(Histo_2->GetName(), "C_{bw}^{I}=-2", "l");
   legend->AddEntry(Histo_3->GetName(), "C_{bw}^{I}=-1", "l");
   legend->AddEntry(Histo_4->GetName(), "C_{bw}^{I}=1", "l");
   legend->AddEntry(Histo_5->GetName(), "C_{bw}^{I}=2", "l");
   legend->AddEntry(Histo_6->GetName(), "C_{tw}^{I}=-2", "l");
   legend->AddEntry(Histo_7->GetName(), "C_{tw}^{I}=-1", "l");
   legend->AddEntry(Histo_8->GetName(), "C_{tw}^{I}=1", "l");
   legend->AddEntry(Histo_9->GetName(), "C_{tw}^{I}=2", "l");

   legend->Draw("SAME");


   Name = Name + ".pdf";
   Canvas->Print(Name.c_str());

}


void Compare_2Vars(TTree* t1, string variable1, string variable2, int nbins, double xmin, double xmax, string selection1, string selection2, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string Name){

  TH1D* Histo_1 = GetHistoWeight(t1, variable1, nbins, xmin, xmax, selection1, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t1, variable2, nbins, xmin, xmax, selection2, "Histo_2");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle(variable1.c_str());
  if(log_fity)
  {
    Canvas->SetLogy();
    Name += "_log";
    legendY += " log";
  }
  if(normalization)
  {
    double a = Histo_1->Integral();
    double b = Histo_2->Integral();
    Histo_1->Scale(1/a);
    Histo_2->Scale(1/b);
    Name += "_normalized";
    legendY += " normalized";
  }
  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  Histo_1->SetMaximum(max*1.3);
  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw();

  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME");

 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.75;
	 lx1 = 0.5;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.75;
	 ly0 = 0.75;
	 lx1 = 0.99;
	 ly1 = 0.99;
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(), "l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(), "l");
  legend->Draw("SAME");

  legend->SetTextSize(0.035);


  Name += ".pdf";
  Canvas->Print(Name.c_str());

  // cout << "Histo1 mean: "<<Histo_1->GetMean()<<endl;
  // cout << "Histo2 mean: "<<Histo_2->GetMean()<<endl;

}

void Compare_3Vars(string histoTitle, TTree* t1, string variable1, string variable2, string variable3, int nbins, double xmin, double xmax, string selection1, string selection2, string selection3, string legendX, string legendY, string legendPlace, string legendtitle, string legendEntry1, string legendEntry2, string legendEntry3, string Name){


  TH1D* Histo_1 = GetHistoWeight(t1, variable1, nbins, xmin, xmax, selection1, "Histo_1");
  TH1D* Histo_2 = GetHistoWeight(t1, variable2, nbins, xmin, xmax, selection2, "Histo_2");
  TH1D* Histo_3 = GetHistoWeight(t1, variable3, nbins, xmin, xmax, selection3, "Histo_3");

  Histo_1->SetStats(kFALSE);
  Histo_2->SetStats(kFALSE);
  Histo_3->SetStats(kFALSE);

  TCanvas* Canvas = new TCanvas("Canvas","Canvas");
  Histo_1->SetTitle(histoTitle.c_str());

  if(log_fity)
  {
    Canvas->SetLogy();
    Name += "_log";
    legendY += " log";
  }
  if(normalization)
  {
    Double_t integral1 = 1/Histo_1->Integral();
    Histo_1->Scale(integral1);

    Double_t integral2 = 1/Histo_2->Integral();
    Histo_2->Scale(integral2);

    Double_t integral3 = 1/Histo_3->Integral();
    Histo_3->Scale(integral3);

    Name += "_normalized";
    legendY += " normalized";
  }

  double max = (Histo_1->GetMaximum()>Histo_2->GetMaximum()) ? Histo_1->GetMaximum() : Histo_2->GetMaximum();
  Histo_1->SetMaximum(max*1.3 );

  Histo_1->SetXTitle(legendX.c_str());
  Histo_1->SetYTitle(legendY.c_str());
  Histo_1->SetLineColor(kRed);
  Histo_1->SetLineWidth(2);
  Histo_1->Draw("h");


  Histo_2->SetLineColor(kBlue);
  Histo_2->SetLineWidth(2);
  Histo_2->Draw("SAME h");


  Histo_3->SetLineColor(kGreen);
  Histo_3->SetLineWidth(2);
  Histo_3->Draw("SAME h");


 double lx0, ly0, lx1, ly1;
  if (legendPlace=="legendUpLeft"){
	 lx0 = 0.2;
	 ly0 = 0.7;
	 lx1 = 0.6;
	 ly1 = 0.95;
  }
   if (legendPlace=="legendUpRight"){
	 lx0 = 0.65;
	 ly0 = 0.75; //Lower Y
	 lx1 = 0.99;
	 ly1 = 0.99; //Upper Y
  }

  TLegend* legend = new TLegend(lx0, ly0, lx1, ly1, legendtitle.c_str());
  legend->SetFillColor(kWhite);
  legend->SetTextSize(0.04);

  legend->AddEntry(Histo_1->GetName(), legendEntry1.c_str(),"l");
  legend->AddEntry(Histo_2->GetName(), legendEntry2.c_str(),"l");
  legend->AddEntry(Histo_3->GetName(), legendEntry3.c_str(),"l");
  //legend->SetLegendSize(0.5);

  legend->Draw("SAME");

  Name += ".pdf";
  Canvas->Print(Name.c_str());
}
