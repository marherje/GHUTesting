#include "../include/experimentalUnc.h"
#include "../include/analysis_prob_newpheno.h"
#include "../style/Style.C"
#include "../style/Labels.C"
#include <vector>
#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

double Prob_To_Sigma(double prob)
{
  double sigma=0.0;
  double testprob=1.0;
  if(prob!=1.0){
    for(int i=1;i<50000;i++){
      sigma+=0.001;
      testprob=erfc(sigma/sqrt(2));
      double dif=testprob-prob;
      if(dif<0.0) break;
    }
  }
  if(sigma==0.0) sigma=0.00001;
  if(prob==0.0) sigma=50.;
  return sigma;
}

double Sigma_To_Prob(double sigma)
{
  double prob=1;
  if(sigma>0.00001){
    prob=erfc(sigma/sqrt(2));
  }
  return prob;
}

void plotGraphs(const std::vector<TGraph*>& graphs, 
                const std::vector<TString>& labels,
                const TString& modelselection,
                const TString& errortype,
                const TString& previous,
                const TString& prospects) {
    TString modelstring;
    if(modelselection=="default") modelstring = "Default";
    else if(modelselection=="A") modelstring = "A";
    else if(modelselection=="Bm") modelstring = "B-";
    else if(modelselection=="Bp") modelstring = "B+";
    TCanvas* c = new TCanvas("c", "Multiple TGraphs", 800, 800);
    c->cd();
    c->SetLogy();
    for (size_t i = 0; i < graphs.size(); ++i) {
        if((i==0)||(i==3)){
            graphs[i]->SetLineColor(kRed+2-2*i);
            graphs[i]->SetMarkerColor(kRed+2-2*i);
        }
        else if((i==1)||(i==4)){
            graphs[i]->SetLineColor(kViolet+10-i);
            graphs[i]->SetMarkerColor(kViolet+10-i);
        }
        else if((i==2)||(i==5)){
            graphs[i]->SetLineColor(kSpring-4-i);
            graphs[i]->SetMarkerColor(kSpring-4-i);
        }
        graphs[i]->SetLineWidth(2);
        graphs[i]->SetTitle(modelstring+" Models (b & c quarks)");
        graphs[i]->GetYaxis()->SetTitle("Separation power");
        graphs[i]->GetYaxis()->SetTitleOffset(1.2);
        graphs[i]->GetYaxis()->SetRangeUser(0.1, 100);
        graphs[i]->GetXaxis()->SetLimits(5, 60);
        graphs[i]->GetXaxis()->SetTitle("m_{KK} [TeV]");
        graphs[i]->SetMarkerStyle(20 + i);
        if (i == 0)
            graphs[i]->Draw("ALP");
        else
            graphs[i]->Draw("LP SAME");
    }

    c->cd();
    if(errortype=="Stat") QQBARLabel3(0.12,0.91,"Statistical uncertainties only",kBlue,0.023);
    // Draw horizontal lines at y=3 and y=5
    TLine* line3 = new TLine(5, 3, 60, 3);
    line3->SetLineStyle(2);
    line3->SetLineColor(kGray+1);
    line3->SetLineWidth(4);
    line3->Draw();

    TLine* line5 = new TLine(5, 5, 60, 5);
    line5->SetLineStyle(2);
    line5->SetLineColor(kGray+2);
    line5->SetLineWidth(4);
    line5->Draw();

    TLegend* legend = new TLegend(0.6, 0.7, 0.9, 0.9);
    for (size_t i = 0; i < graphs.size(); ++i) {
        legend->AddEntry(graphs[i], labels[i], "lp");
    }
    legend->Draw();
    c->Update();
    c->SaveAs("GHU_"+modelstring+"_"+errortype+"_"+previous+"_"+prospects+".png");
    c->SaveAs("GHU_"+modelstring+"_"+errortype+"_"+previous+"_"+prospects+".eps");
}

// Para combinar todos los MCs de una energía en un solo histograma de nsigma
auto combine_all_mcs = [](std::vector<observables_struct_t> observables, int mc_first, int mc_last, double energy, TString errortype, TString tpc_status, bool unpol = false) {
  std::vector<std::vector<TH2F*>> probhistos(mc_last+1);
  //mc_first
  for (int mc = mc_first; mc <= mc_last; ++mc) {
    probhistos[mc] = nsigmas_models(mc, energy, observables, errortype, tpc_status);  
  }
  const Int_t xyNBINS = 8;
  Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
  TString energyst = TString::Format("%.0f", energy);
  TString histoname = "AFB_" + energyst + "_" + errortype + "_" + tpc_status;
  histoname += (unpol ? "_unpol" : "");
  histoname += TString::Format("_mc%d_%d", mc_first, mc_last);

  TString histotitle = "AFB_" + energyst + "_" + errortype + "_" + tpc_status;
  histotitle += (unpol ? " unpol" : "");
  histotitle += TString::Format(" mc%d_%d", mc_first, mc_last);

  TH1F* result = new TH1F(
    histoname,
    histotitle,
    xyNBINS, xyedges
  );

  for(int j=1; j<theory.size(); j++) {
      double prob = 1.0;
      //mc_first
      for (int mc = mc_first; mc <= mc_last; ++mc) {
        // Para unpol solo usamos el histograma [2], para pol usamos [0] y [1]
        if (unpol) {
          prob *= Sigma_To_Prob(probhistos[mc][2]->GetBinContent(j+1,1));
        } else {
          prob *= Sigma_To_Prob(probhistos[mc][0]->GetBinContent(j+1,1));
          prob *= Sigma_To_Prob(probhistos[mc][1]->GetBinContent(j+1,1));
        }
      }
      // sqrt(2) factor to account for the prediction sigma as well, not just the experimental one
      double nsigma = (1/sqrt(2))*Prob_To_Sigma(prob);
      // Rounding
      if((nsigma>0.95)&&(nsigma<1))nsigma=0.9;
      if((nsigma>1.95)&&(nsigma<2))nsigma=1.9;
      if((nsigma>2.95)&&(nsigma<3))nsigma=2.9;
      if((nsigma>3.95)&&(nsigma<4))nsigma=3.9;
      if((nsigma>4.95)&&(nsigma<5))nsigma=4.9;
      if((nsigma>9.95)&&(nsigma<10))nsigma=9.9;
      result->SetBinContent(j, nsigma);
    }
    return result;
  
};

  // Para 250+500 GeV
auto combine_two_energies = [](std::vector<observables_struct_t> observables, int mc_first, int mc_last, double energy1, double energy2, TString errortype, TString tpc_status)
{
    
    std::vector<std::vector<TH2F*>> probhistos_1(mc_last+1), probhistos_2(mc_last+1);
    //mc_first
    for (int mc = mc_first; mc <= mc_last; ++mc) {
      probhistos_1[mc]   = nsigmas_models(mc, energy1, observables, errortype, tpc_status);
      probhistos_2[mc]   = nsigmas_models(mc, energy2, observables, errortype, tpc_status);
    }
    const Int_t xyNBINS = 8;
    Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
    TString energyst1 = TString::Format("%.0f", energy1);
    TString energyst2 = TString::Format("%.0f", energy2);
    TString energyst = energyst1 + "_" + energyst2;

    TString histoname = "AFB_both_" + energyst + "_" + errortype + "_" + tpc_status;
    histoname += TString::Format("_mc%d_%d", mc_first, mc_last);

    TString histotitle = "AFB_both_" + energyst + "_" + errortype + "_" + tpc_status;
    histotitle += TString::Format(" mc%d_%d", mc_first, mc_last);

    TH1F* result = new TH1F(
      histoname,
      histotitle,
      xyNBINS, xyedges
    );

    for(int j=1; j<theory.size(); j++) {
      double prob = 1.0;
      //mc_first
      for (int mc = mc_first; mc <= mc_last; ++mc) {
        prob *= Sigma_To_Prob(probhistos_1[mc][0]->GetBinContent(j+1,1));
        prob *= Sigma_To_Prob(probhistos_1[mc][1]->GetBinContent(j+1,1));
        prob *= Sigma_To_Prob(probhistos_2[mc][0]->GetBinContent(j+1,1));
        prob *= Sigma_To_Prob(probhistos_2[mc][1]->GetBinContent(j+1,1));
      }
      double nsigma = (1/sqrt(2))*Prob_To_Sigma(prob);
      // Rounding
      if((nsigma>0.95)&&(nsigma<1))nsigma=0.9;
      if((nsigma>1.95)&&(nsigma<2))nsigma=1.9;
      if((nsigma>2.95)&&(nsigma<3))nsigma=2.9;
      if((nsigma>3.95)&&(nsigma<4))nsigma=3.9;
      if((nsigma>4.95)&&(nsigma<5))nsigma=4.9;
      if((nsigma>9.95)&&(nsigma<10))nsigma=9.9;
      result->SetBinContent(j, nsigma);
    }
    return result;
};

// Para 250+500+1000 GeV
auto combine_three_energies = [](std::vector<observables_struct_t> observables, int mc_first, int mc_last, double energy1, double energy2, double energy3, TString errortype, TString tpc_status)
{

    std::vector<std::vector<TH2F*>> probhistos_1(mc_last+1), probhistos_2(mc_last+1), probhistos_3(mc_last+1);
    //mc_first
    for (int mc = mc_first; mc <= mc_last; ++mc) {
      probhistos_1[mc]   = nsigmas_models(mc, energy1, observables, errortype, tpc_status);
      probhistos_2[mc]   = nsigmas_models(mc, energy2, observables, errortype, tpc_status);
      probhistos_3[mc]   = nsigmas_models(mc, energy3, observables, errortype, tpc_status);
    }
    const Int_t xyNBINS = 8;
    Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
    TString energyst1 = TString::Format("%.0f", energy1);
    TString energyst2 = TString::Format("%.0f", energy2);
    TString energyst3 = TString::Format("%.0f", energy3);
    TString energyst = energyst1 + "_" + energyst2 + "_" + energyst3;

    TString histoname = "AFB_three_" + energyst + "_" + errortype + "_" + tpc_status;
    histoname += TString::Format("_mc%d_%d", mc_first, mc_last);

    TString histotitle = "AFB_three_" + energyst + "_" + errortype + "_" + tpc_status;
    histotitle += TString::Format(" mc%d_%d", mc_first, mc_last);

    TH1F* result = new TH1F(
      histoname,
      histotitle,
      xyNBINS, xyedges
    );

    for(int j=1; j<theory.size(); j++) {
      double prob = 1.0;
      //mc_first
      for (int mc = mc_first; mc <= mc_last; ++mc) {
        prob *= Sigma_To_Prob(probhistos_1[mc][0]->GetBinContent(j+1,1));
        prob *= Sigma_To_Prob(probhistos_1[mc][1]->GetBinContent(j+1,1));
        prob *= Sigma_To_Prob(probhistos_2[mc][0]->GetBinContent(j+1,1));
        prob *= Sigma_To_Prob(probhistos_2[mc][1]->GetBinContent(j+1,1));
        prob *= Sigma_To_Prob(probhistos_3[mc][0]->GetBinContent(j+1,1));
        prob *= Sigma_To_Prob(probhistos_3[mc][1]->GetBinContent(j+1,1));
      }
      double nsigma = (1/sqrt(2))*Prob_To_Sigma(prob);
      // Rounding
      if((nsigma>0.95)&&(nsigma<1))nsigma=0.9;
      if((nsigma>1.95)&&(nsigma<2))nsigma=1.9;
      if((nsigma>2.95)&&(nsigma<3))nsigma=2.9;
      if((nsigma>3.95)&&(nsigma<4))nsigma=3.9;
      if((nsigma>4.95)&&(nsigma<5))nsigma=4.9;
      if((nsigma>9.95)&&(nsigma<10))nsigma=9.9;
      result->SetBinContent(j, nsigma);
    }
    
  return result;
};

void test_prospects(string modelselection, string errortype, string previous, string prospects) {
    //Read all models:
    read_all_models(false,modelselection);
    std::vector<observables_struct_t> observables=create_observables();
    // Ejemplo usando plots de dEdx 
    // PID=="dNdx_Extraquarks_per10"
    // errortype==errortype 
    TString errorstring;
    if(prospects == "dNdx_Extraquarks_per1000") errorstring =  " (+#delta AFB(s)=0.001)";
    else if(prospects == "dNdx_Extraquarks_per100") errorstring =  " (+#delta AFB(s)=0.01)";
    else if(prospects == "dNdx_Extraquarks_per10") errorstring =  " (+#delta AFB(s)=0.1)";
    else if(prospects == "dEdx_Extraquarks_per1000") errorstring =  " (+#delta AFB(s)=0.001)";
    else if(prospects == "dEdx_Extraquarks_per100") errorstring =  " (+#delta AFB(s)=0.01)";
    else if(prospects == "dEdx_Extraquarks_per10") errorstring =  " (+#delta AFB(s)=0.1)";
    else if(prospects == "ParT_Extraquarks_per1000") errorstring =  " (+#delta AFB(s)=0.001)";
    else if(prospects == "ParT_Extraquarks_per100") errorstring =  " (+#delta AFB(s)=0.01)";
    else if(prospects == "ParT_Extraquarks_per10") errorstring =  " (+#delta AFB(s)=0.1)";
    TString ILC_labels[3];
    if(previous=="dNdx"){
      ILC_labels[0] = "[dN/dx] ILC250";
      ILC_labels[1] = "[dN/dx] ILC250+500";
      ILC_labels[2] = "[dN/dx] ILC250+500+1000";
    }
    else if(previous=="dEdx"){
      ILC_labels[0] = "[dE/dx] ILC250";
      ILC_labels[1] = "[dE/dx] ILC250+500";
      ILC_labels[2] = "[dE/dx] ILC250+500+1000";
    }
    else if(previous=="noTPC"){
      ILC_labels[0] = "[noTPC] ILC250";
      ILC_labels[1] = "[noTPC] ILC250+500";      
      ILC_labels[2] = "[noTPC] ILC250+500+1000";
    }
    else if(previous=="ParT"){
      ILC_labels[0] = "[ParT] ILC250";
      ILC_labels[1] = "[ParT] ILC250+500";
      ILC_labels[2] = "[ParT] ILC250+500+1000";
    }

    std::vector<TString> labels;
    labels = {ILC_labels[0], ILC_labels[1], ILC_labels[2], ILC_labels[0]+errorstring, ILC_labels[1]+errorstring, ILC_labels[2]+errorstring};

    TH1F * results_250_bc = combine_all_mcs(observables, 4, 5, 250, errortype, previous, false);
    TH1F * results_both_bc = combine_two_energies(observables, 4, 5, 250, 500, errortype, previous);
    TH1F * results_three_bc = combine_three_energies(observables, 4, 5, 250, 500, 1000, errortype, previous);

    TH1F * results_250_bc_s1000 = combine_all_mcs(observables, 3, 5, 250, errortype, prospects, false);
    TH1F * results_both_bc_s1000 = combine_two_energies(observables, 3, 5, 250, 500, errortype, prospects);
    TH1F * results_three_bc_s1000 = combine_three_energies(observables, 3, 5, 250, 500, 1000, errortype, prospects);

    double x[8];
    if(modelselection == "Bm") {
        double temp_x[] = {13, 19, 25, 31, 37, 43, 49, 55};
        for(int i = 0; i < 8; i++) x[i] = temp_x[i];
    }
    else if(modelselection == "Bp") {
        double temp_x[] = {13, 19, 25, 31, 37, 43, 49, 55};
        for(int i = 0; i < 8; i++) x[i] = temp_x[i];
    }
    else if(modelselection == "A") {
        double temp_x[] = {8.8, 10.3, 25, 31, 37, 43, 49, 55};
        for(int i = 0; i < 8; i++) x[i] = temp_x[i];
    }
    else {
        double temp_x[] = {8.8, 10.3, 13, 13, 19, 19, 25, 25};
        for(int i = 0; i < 8; i++) x[i] = temp_x[i];
    }
    double y_250_bc[8];
    double y_both_bc[8];
    double y_three_bc[8];
    double y_250_bc_s1000[8];
    double y_both_bc_s1000[8];
    double y_three_bc_s1000[8];

    int realpoints;
    if(modelselection == "default") realpoints = 8;
    else if(modelselection == "A") realpoints = 2;
    else if(modelselection == "Bm") realpoints = 3;
    else if(modelselection == "Bp") realpoints = 3;

    for (int i = 1; i <= realpoints; ++i) {
        y_250_bc[i-1] = results_250_bc->GetBinContent(i);
        y_both_bc[i-1] = results_both_bc->GetBinContent(i);
        y_three_bc[i-1] = results_three_bc->GetBinContent(i);
        y_250_bc_s1000[i-1] = results_250_bc_s1000->GetBinContent(i);
        y_both_bc_s1000[i-1] = results_both_bc_s1000->GetBinContent(i);
        y_three_bc_s1000[i-1] = results_three_bc_s1000->GetBinContent(i);
        
    }
    
    // Extrapolation for the remaining points
    for (int i = realpoints+1; i <= 8; ++i) {
        y_250_bc[i-1] = 0.6*y_250_bc[i-2];
        y_both_bc[i-1] = 0.6*y_both_bc[i-2];
        y_three_bc[i-1] = 0.6*y_three_bc[i-2];
        y_250_bc_s1000[i-1] = 0.6*y_250_bc_s1000[i-2];
        y_both_bc_s1000[i-1] = 0.6*y_both_bc_s1000[i-2];
        y_three_bc_s1000[i-1] = 0.6*y_three_bc_s1000[i-2];
    }

    TGraph* graph_250_bc = new TGraph(8, x, y_250_bc);
    TGraph* graph_both_bc = new TGraph(8, x, y_both_bc);
    TGraph* graph_three_bc = new TGraph(8, x, y_three_bc);
    TGraph* graph_250_bc_s1000 = new TGraph(8, x, y_250_bc_s1000);
    TGraph* graph_both_bc_s1000 = new TGraph(8, x, y_both_bc_s1000);
    TGraph* graph_three_bc_s1000 = new TGraph(8, x, y_three_bc_s1000);

    const std::vector<TGraph*>& graphs = {
        graph_250_bc,
        graph_both_bc,
        graph_three_bc,
        graph_250_bc_s1000,
        graph_both_bc_s1000,
        graph_three_bc_s1000
    };

    plotGraphs(graphs,labels,modelselection,errortype,previous,prospects);

}

void test_newpheno(string modelselection, string errortype, string previous, string prospects) {
  test_prospects(modelselection,errortype,previous,prospects);

}