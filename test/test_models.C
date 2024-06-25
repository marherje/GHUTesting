// # Copyright 2023  Adrián Irles (IFIC)
#include "../include/experimentalUnc.h"
#include "../include/analysis_prob.h"
#include "../style/Style.C"
#include "../style/Labels.C"

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

void Labels(double energy_under_test=250, TString errortype="Stat", TString tpc_status="default", int iter=8, bool pol=true, TString style="proceedings")
{
  TString energyst="250 GeV";
  if(energy_under_test==91.2) energyst="Z-Pole";
  else if(energy_under_test==380) energyst="380 GeV";
  else if(energy_under_test==500) energyst="500 GeV";
  else if(energy_under_test==1000) energyst="1 TeV";
  else if(energy_under_test==9999) energyst="250&500 GeV";  
  else if(energy_under_test==10000) energyst="250&500&1000 GeV";
  
  TString energypaper="ILC250";
  if(energy_under_test==91.2) energypaper="ILC Giga-Z";
  else if(energy_under_test==380) energypaper="380 GeV";
  else if(energy_under_test==500) energypaper="ILC500";
  else if(energy_under_test==1000) energypaper="ILC1000";
  else if(energy_under_test==9999) energypaper="ILC250+500";
  else if(energy_under_test==10000) energypaper="ILC250+500+1000*";
  
  TString lumipaper="(2000 fb^{-1})";
  if(energy_under_test==91.2) lumipaper="ILC Giga-Z";
  else if(energy_under_test==380) lumipaper="380 GeV";
  else if(energy_under_test==500) lumipaper="4000 fb^{-1}";
  else if(energy_under_test==1000) lumipaper="8000 fb^{-1}";
  else if(energy_under_test==9999) lumipaper="(2000 fb^{-1}+ 4000 fb^{-1})";
  else if(energy_under_test==10000) lumipaper="(2000 fb^{-1}+ 4000 fb^{-1}+ 8000 fb^{-1})";

  float ILCpaper=0.6;
  if(energy_under_test==91.2) ILCpaper=0.45;
  else if(energy_under_test==250) ILCpaper=0.6;
  else if(energy_under_test==380) ILCpaper=0.6;
  else if(energy_under_test==500) ILCpaper=0.6;
  else if(energy_under_test==1000) ILCpaper=0.57;
  else if(energy_under_test==9999) ILCpaper=0.45;
  else if(energy_under_test==10000) ILCpaper=0.35;
  
  TString conditions[]={TString::Format("H20, P:(-0.8,+0.3), b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b-quark, #sqrt{s}=%s",energyst.Data()),
			TString::Format("H20, P:(-0.8,+0.3), c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., c-quark, #sqrt{s}=%s",energyst.Data()),
			TString::Format("H20, P:(-0.8,+0.3), b&c quarks, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), b&c quarks, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b&c quarks, #sqrt{s}=%s",energyst.Data())};
  TString conditions_clic[]={TString::Format("CLIC-like, P:-0.8, b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, P:+0.8, b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, both  pol., b-quark, #sqrt{s}=%s",energyst.Data()),
			     TString::Format("CLIC-like, P:-0.8, c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, P:+0.8, c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, both pol., c-quark, #sqrt{s}=%s",energyst.Data()),
			     TString::Format("CLIC-like, P:-0.8, b&c quarks, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, P:+0.8, b&c quarks, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, both pol., b&c quarks. #sqrt{s}=%s",energyst.Data())};
  TString conditions_unpol[]={TString::Format("unpol., b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., b&c-quark, #sqrt{s}=%s",energyst.Data())};


  //QQBARLabel3(0.675,0.60,"ILD",kBlack,0.15);
  QQBARLabel(0.73,0.68,"",-1,0.1);
  
  if(style=="proceedings"){
    if(energy_under_test!=380){
      if(tpc_status=="noTPC")QQBARLabel2(0.3,0.225, "No PID",kRed,0.055);
      if(tpc_status=="default")QQBARLabel2(0.3,0.225, "dEdx for CM",kRed,0.055);
      if(tpc_status=="dEdx")QQBARLabel2(0.3,0.225, "dEdx for CM & FT",kRed,0.055);
      if(tpc_status=="dNdx")QQBARLabel2(0.3,0.225, "dNdx for CM & FT",kRed,0.055);
    }
    if(pol==true){
      if(energy_under_test==380)QQBARLabel2(0.20,0.15, conditions_clic[iter],kBlack,0.050);
      else QQBARLabel2(0.20,0.15, conditions[iter],kBlack,0.055);
    }
    else QQBARLabel2(0.20,0.15, conditions_unpol[iter],kBlack,0.055);
  }
  else if(style=="paper"){
    if(tpc_status=="noTPC")QQBARLabel2(0.3,0.3, "No PID",kRed,0.055);
    if(tpc_status=="dEdx")QQBARLabel2(0.3,0.3, "dEdx for CM & FT",kRed,0.055);
    if(pol==false){
      if(energy_under_test==250)energypaper="ILC250^{#diamond}";
      QQBARLabel2(ILCpaper-0.1,0.18,energypaper+" (no pol.)",kBlack,0.055);
      QQBARLabel2(ILCpaper-0.1,0.13,lumipaper,kGray+2,0.04);
    }
    else{
    QQBARLabel2(ILCpaper,0.18,energypaper,kBlack,0.055);
    QQBARLabel2(ILCpaper,0.13,lumipaper,kGray+2,0.04);
    }
  }

}


void DrawLeg(double energy_under_test=250)
{
  TString colors_st[] = {"<1 #sigma","1-2 #sigma","< 3 #sigma","3-4 #sigma","4-5 #sigma","> 5 #sigma"}; // #colors >= #levels - 1  

  Int_t colors[]={kGreen+3,kGreen+3,kGreen+3,kTeal+4,kTeal+5,kTeal+6};
  /*
  for(int i=0;i<7;i++){
    if(energy_under_test==91.2) colors[i]=colors_Z[i];
    else if(energy_under_test==250) colors[i]=colors_250[i];
    else if(energy_under_test==380) colors[i]=colors_380[i];
    else if(energy_under_test==500) colors[i]=colors_500[i];
    else if(energy_under_test==1000) colors[i]=colors_1000[i];
    else if(energy_under_test==9999) colors[i]=colors_both[i];
    else if(energy_under_test==10000) colors[i]=colors_three[i];
  }
  */
  
  gStyle->SetPalette(7,colors);
  
  TLegend *leg2 = new TLegend(0.68,0.3,1.05,0.5,"","blNDC");
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.04);
  TH2F* h_leg02[10];
  for(int ic=4; ic<6;ic++ ){
    h_leg02[ic]=new TH2F;
    h_leg02[ic]->SetFillStyle(1000);
    h_leg02[ic]->SetFillColor(colors[ic]);
    h_leg02[ic]->SetLineColor(colors[ic]);
    leg2->AddEntry(h_leg02[ic],colors_st[ic],"f");
  }
  leg2->SetFillColor(0);
  leg2->SetFillStyle(0);
  leg2->SetBorderSize(0);
  leg2->Draw();

TLegend *leg = new TLegend(0.5,0.3,0.87,0.5,"","blNDC");
leg->SetTextFont(42);
leg->SetTextSize(0.04);
TH2F* h_leg0[10];
for(int ic=2; ic<4;ic++ ){
  h_leg0[ic]=new TH2F;
  h_leg0[ic]->SetFillStyle(1000);
  h_leg0[ic]->SetFillColor(colors[ic]);
  h_leg0[ic]->SetLineColor(colors[ic]);
  leg->AddEntry(h_leg0[ic],colors_st[ic],"f");
 }
 leg->SetFillColor(0);
 leg->SetFillStyle(0);
 leg->SetBorderSize(1);
 leg->Draw();

}

void DrawTopLine()
{
  /*
  TLine *l=new TLine(-0.5,8.5,8.5,8.5);
  l->SetLineColor(kBlack);
  l->Draw();
  */
  TPad *p = new TPad("p","p",0.,0.,1.,1.);
  p->SetFillStyle(0);
  p->Draw("same");
  p->cd();
  TBox *b = new TBox(0.1,0.9,0.9,0.9);
  b->SetFillColor(kBlack);
  b->SetLineColor(kBlack);
  b->SetLineWidth(0);
  b->Draw("same");
}

/*
void DrawOfficialTitle()
{

  TPad *p = new TPad("p","p",0.,0.,1.,1.); 
  p->SetFillStyle(0); 
  p->Draw("same"); 
  p->cd();
  TBox *b = new TBox(0.1,0.91,0.9,0.99);
  b->SetFillColor(kWhite);
  b->SetLineColor(kWhite);
  b->SetLineWidth(0.0);
  b->Draw("same");

  QQBARLabel3(0.175,0.91,"Model discrimination in #sigma-level",kBlack,0.06);

}
*/

void DrawTextBins(std::vector<double> xtodraw,std::vector<double> ytodraw,std::vector<double> binvalue,Int_t color)
{
  for(int idraw=0; idraw<size(xtodraw); idraw++){ 
    double px=xtodraw.at(idraw);                                                           
    double py=ytodraw.at(idraw);      
    double sigmavalue=binvalue.at(idraw);
    //Carefully rounding so 4.96 don't turn into 5.0 and things like that
    //if((sigmavalue>0.95)and(sigmavalue<0.99))sigmavalue=0.9;
    //if((sigmavalue>1.95)and(sigmavalue<1.99))sigmavalue=1.9;
    //if((sigmavalue>2.95)and(sigmavalue<2.99))sigmavalue=2.9;
    //if((sigmavalue>3.95)and(sigmavalue<3.99))sigmavalue=3.9;
    //if((sigmavalue>4.95)and(sigmavalue<4.99))sigmavalue=4.9;
    //if((sigmavalue>9.95)and(sigmavalue<9.99))sigmavalue=9.9;
    //Moved the roundingt into the loops so "text" match
    
    if(sigmavalue>9.94)QQBARLabel3(px-0.005,py,">10",color,0.045);
    else QQBARLabel3(px,py,TString::Format("%.1f",sigmavalue),color,0.045); 
  }    
}


void test_one_energy(double energy_under_test=250, TString errortype="Stat", TString tpc_status="default", TString style="intern")
{

  TString energyst="250 GeV";
  if(energy_under_test==91.2) energyst="Z-Pole";
  else if(energy_under_test==350) energyst="350 GeV";
  else if(energy_under_test==380) energyst="380 GeV";
  else if(energy_under_test==500) energyst="500 GeV";
  else if(energy_under_test==1000) energyst="1 TeV";
  
  //cout<<"check 0"<<endl;
  read_all_models(false);
  //cout<<"check 1"<<endl;
  // stores the information for all  models (named in std::vector<TString> allmodels={"SMA","A1"};)
    // with 100% polarization options only
    // it stores the same info as the txt files
    // in a vector std::vector<model_struct_t> theory;
    // to access the info, you do: theory.at(i).model_values[ie].AFB_L[ifl]
    // indexes:
    //   i: as many as model strings
    //   ie: index for the energy: one index for every energy found in the file. If you want to know the value, you can acces to i: theory.at(i).model_values[ie].energy
    //       we can store up to 5 energy points
    //   ifl: for the flavour , 0 is empty, 1,2,3,4,5 are d,u,s,c,b

    //Content of histos: AFB_L, AFB_R, AFB_unpol, R_L, R_R, R_unpol


    // here I create the theoretical observables for the desired polarization scheme
    // I save a vector for all energies, all pol (2 in this case) and and the b and c flavours
    std::vector<observables_struct_t> observables=create_observables();
    //cout<<"check 2"<<endl;
    // I calculate the prob matrix of tested models vs hypothesis models, using th2f
    std::vector<TH2F*> probhistos_c=nsigmas_models(4, energy_under_test, observables, errortype, tpc_status) ;
    std::vector<TH2F*> probhistos_b=nsigmas_models(5, energy_under_test, observables, errortype, tpc_status) ;
    //cout<<"check 3"<<endl;
    //create final histograms
    const Int_t xyNBINS = 8;
    Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
    TH2F * AFB_results[9];                                      
    AFB_results[0] = new TH2F("AFB_b_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), b-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[1] = new TH2F("AFB_b_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[2] = new TH2F("AFB_b_comb",TString::Format("AFB_ILD_H20, b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[3] = new TH2F("AFB_c_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[4] = new TH2F("AFB_c_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[5] = new TH2F("AFB_c_comb",TString::Format("AFB_ILD_H20, c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[6] = new TH2F("AFB_b_&_AFB_c_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[7] = new TH2F("AFB_b_&_AFB_c_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[8] = new TH2F("AFB_b_&_AFB_c_comb",TString::Format("AFB_ILD_H20, b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
 
    TH2F * R_results[9];
    R_results[0] = new TH2F("R_b_L",TString::Format("R_ILD_H20 (-0.8,+0.3), b-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[1] = new TH2F("R_b_R",TString::Format("R_ILD_H20 (+0.8,-0.3), b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[2] = new TH2F("R_b_comb",TString::Format("R_ILD_H20, b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[3] = new TH2F("R_c_L",TString::Format("R_ILD_H20 (-0.8,+0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[4] = new TH2F("R_c_R",TString::Format("R_ILD_H20 (+0.8,-0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[5] = new TH2F("R_c_comb",TString::Format("R_ILD_H20, c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[6] = new TH2F("R_b_&_R_c_L",TString::Format("R_ILD_H20 (-0.8,+0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[7] = new TH2F("R_b_&_R_c_R",TString::Format("R_ILD_H20 (+0.8,-0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[8] = new TH2F("R_b_&_R_c_comb",TString::Format("R_ILD_H20, b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);

    TH2F * AFB_unpol[3];
    AFB_unpol[0] = new TH2F("AFB_b_unpol",TString::Format("AFB_ILD, unpol., b-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_unpol[1] = new TH2F("AFB_c_unpol",TString::Format("AFB_ILD, unpol., c-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_unpol[2] = new TH2F("AFB_b_c_unpol",TString::Format("AFB_ILD, unpol., b&c-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);

    TH2F * R_unpol[3];
    R_unpol[0] = new TH2F("R_b_unpol",TString::Format("R_ILD, unpol., b-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_unpol[1] = new TH2F("R_c_unpol",TString::Format("R_ILD, unpol., c-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_unpol[2] = new TH2F("R_b_c_unpol",TString::Format("R_ILD, unpol., b&c-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);

    TString conditions[]={TString::Format("H20, P:(-0.8,+0.3), b-quark #sqrt{s}=%s",energyst.Data()),TString::Format("H20, b-quark P:(+0.8,-0.3),  #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b-quark #sqrt{s}=%s",energyst.Data()),
        TString::Format("H20, P:(-0.8,+0.3), c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., c-quark #sqrt{s}=%s",energyst.Data()),
        TString::Format("H20, P:(-0.8,+0.3), b&c quarks #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), b&c quarks #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b&c quarks #sqrt{s}=%s",energyst.Data())};
    TString conditions_clic[]={TString::Format("CLIC-like, P:-0.8, b-quark #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, b-quark P:+0.8,  #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, both  pol., b-quark #sqrt{s}=%s",energyst.Data()),
			  TString::Format("CLIC-like, P:-0.8, c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, P:+0.8, c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, both pol., c-quark #sqrt{s}=%s",energyst.Data()),
			  TString::Format("CLIC-like, P:-0.8, b&c quarks #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, P:+0.8, b&c quarks #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, both pol., b&c quarks #sqrt{s}=%s",energyst.Data())};
    TString conditions_unpol[]={TString::Format("unpol., b-quark #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., b&c-quark #sqrt{s}=%s",energyst.Data())};
    TString titles_AFB[]={"AFB_{b} (e^{-}_{L}e^{+}_{R})","AFB_{b} (e^{-}_{R}e^{+}_{L})","AFB_{b} (Both pol.)","AFB_{c} (e^{-}_{L}e^{+}_{R})","AFB_{c} (e^{-}_{R}e^{+}_{L})","AFB_{c} (Both pol.)","AFB_{b} and AFB_{c} (e^{-}_{L}e^{+}_{R})","AFB_{b} and AFB_{c} (e^{-}_{R}e^{+}_{L})","AFB_{b} & AFB_{c} (Both pol.)"};
    TString titles_AFB_unpol[]={"AFB_{b} (unpol.)","AFB_{c} (unpol.)","AFB_{b} & AFB_{c} (unpol.)"};
    TString titles_R[]={"R_{b} (e^{-}_{L}e^{+}_{R})","R_{b} (e^{-}_{R}e^{+}_{L})","R_{b} (Both pol.)","R_{c} (e^{-}_{L}e^{+}_{R})","R_{c} (e^{-}_{R}e^{+}_{L})","R_{c} (Both pol.)","R_{b} and R_{c} (e^{-}_{L}e^{+}_{R})","R_{b} and R_{c} (e^{-}_{R}e^{+}_{L})","R_{b} and R_{c} (Both pol.)"};
    TString titles_R_unpol[]={"R_{b} (unpol.)","R_{c} (unpol.)","R_{b} & R_{c} (unpol.)"};

    TString title_models_x[]={"SMB","A_{1}","A_{2}","B_{1}^{-}","B_{1}^{+}","B_{2}^{-}","B_{2}^{+}","B_{3}^{-}","B_{3}^{+}"};
    TString title_models_y[]={"SMB","A_{1}","A_{2}","B_{1}^{-}","B_{1}^{+}","B_{2}^{-}","B_{2}^{+}","B_{3}^{-}","B_{3}^{+}"};

    for(int i=0; i<9; i++) {
     for(int j=1; j<theory.size(); j++) {
       AFB_results[i]->GetXaxis()->SetBinLabel(j,title_models_x[j]);
       AFB_results[i]->GetYaxis()->SetBinLabel(j,title_models_y[j]);
       R_results[i]->GetXaxis()->SetBinLabel(j,title_models_x[j]);
       R_results[i]->GetYaxis()->SetBinLabel(j,title_models_y[j]);
     }
     AFB_results[i]->GetXaxis()->SetTickLength(0.);
     AFB_results[i]->GetYaxis()->SetTickLength(0.);
     AFB_results[i]->GetYaxis()->SetLabelSize(0.075);
     AFB_results[i]->GetXaxis()->SetLabelSize(0.075);
     R_results[i]->GetXaxis()->SetTickLength(0.);
     R_results[i]->GetYaxis()->SetTickLength(0.);
     R_results[i]->GetYaxis()->SetLabelSize(0.075);
     R_results[i]->GetXaxis()->SetLabelSize(0.075);
    }  

    for(int i=0; i<3; i++) {
      for(int j=1; j<theory.size(); j++) {
	AFB_unpol[i]->GetXaxis()->SetBinLabel(j,title_models_x[j]);
	AFB_unpol[i]->GetYaxis()->SetBinLabel(j,title_models_y[j]);
	R_unpol[i]->GetXaxis()->SetBinLabel(j,title_models_x[j]);
	R_unpol[i]->GetYaxis()->SetBinLabel(j,title_models_y[j]);
      }
      AFB_unpol[i]->GetXaxis()->SetTickLength(0.);
      AFB_unpol[i]->GetYaxis()->SetTickLength(0.);
      AFB_unpol[i]->GetYaxis()->SetLabelSize(0.07);
      AFB_unpol[i]->GetXaxis()->SetLabelSize(0.07);
      R_unpol[i]->GetXaxis()->SetTickLength(0.);
      R_unpol[i]->GetYaxis()->SetTickLength(0.);
      R_unpol[i]->GetYaxis()->SetLabelSize(0.07);
      R_unpol[i]->GetXaxis()->SetLabelSize(0.07);
    }

    //Status string:
    TString errorstring;
    if(errortype=="Stat")errorstring="Stat";
    else if(errortype=="StatSyst")errorstring="StatSyst";
    else if(errortype=="StatTheoCurrent")errorstring="StatTheoCurrent";
    else if(errortype=="StatTheoGigaZ")errorstring="StatTheoGigaZ";
    else if(errortype=="StatTheo250")errorstring="StatTheo250";
    
    TString status_string;
    if(energy_under_test==91.2){
      status_string="_Zpole_"+tpc_status+"_"+errorstring;
    }
    else status_string=TString::Format("_%3.0f_"+tpc_status+"_"+errorstring,energy_under_test);

    Int_t colors[6];
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    Double_t levels[7]={0,0.99,1.99,2.99,3.99,4.99,1000.0};
    gStyle->SetPaintTextFormat(".1f");

    TCanvas *c_AFB_full = new TCanvas ("c_AFB_full"+status_string,"c_AFB_full"+status_string,800,800);
    c_AFB_full->Divide(3,3);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerSize(2.);
    std::vector<double> xWhite[9];
    std::vector<double> yWhite[9];
    std::vector<double> binWhite[9];
    std::vector<double> xBlack[9];
    std::vector<double> yBlack[9];
    std::vector<double> binBlack[9];

    //  cout<<"check 4"<<endl;
    for(int i=0; i<9; i++) {
        c_AFB_full->cd(i+1);
	int binx;
	int biny;
	for(int j=1; j<theory.size(); j++) {
            for(int k=1; k<theory.size(); k++){
	      if(j>=k) continue;
	      biny=k;
	      binx=j;
	      double nsigma=0;
		if(i<2){
		  nsigma=probhistos_b.at(i)->GetBinContent(binx+1,biny+1);
		  cout<<"Sigmas: "<<nsigma<<endl;
		}
		else if(i==2){
                  double prob1=Sigma_To_Prob(probhistos_b.at(i-2)->GetBinContent(binx+1,biny+1));
                  double prob2=Sigma_To_Prob(probhistos_b.at(i-1)->GetBinContent(binx+1,biny+1));
                  double prob=prob1*prob2;
                  nsigma=Prob_To_Sigma(prob);
                }
		else if(i<5){
		  nsigma=probhistos_c.at(i-3)->GetBinContent(binx+1,biny+1);
		}
		else if(i==5){
                  double prob1=Sigma_To_Prob(probhistos_c.at(i-5)->GetBinContent(binx+1,biny+1));
                  double prob2=Sigma_To_Prob(probhistos_c.at(i-4)->GetBinContent(binx+1,biny+1));
                  double prob=prob1*prob2;
                  nsigma=Prob_To_Sigma(prob);
                }
		else if(i<8){
		  double prob1=Sigma_To_Prob(probhistos_b.at(i-6)->GetBinContent(binx+1,biny+1));
		  double prob2=Sigma_To_Prob(probhistos_c.at(i-6)->GetBinContent(binx+1,biny+1));
		  double prob=prob1*prob2;
		  nsigma=Prob_To_Sigma(prob);
		}
		else{
		  double prob1=Sigma_To_Prob(probhistos_b.at(0)->GetBinContent(binx+1,biny+1));
                  double prob2=Sigma_To_Prob(probhistos_b.at(1)->GetBinContent(binx+1,biny+1));
		  double prob3=Sigma_To_Prob(probhistos_c.at(0)->GetBinContent(binx+1,biny+1));
                  double prob4=Sigma_To_Prob(probhistos_c.at(1)->GetBinContent(binx+1,biny+1));
                  double prob=prob1*prob2*prob3*prob4;
                  nsigma=Prob_To_Sigma(prob);
		  cout<<"Probs: "<<prob1<<", "<<prob2<<", "<<prob3<<", "<<prob4<<endl;
		  cout<<"Sigmas: "<<nsigma<<". Prob: "<<prob<<endl;
		}
		// Rounding for plots:
		if((nsigma>0.95)and(nsigma<1))nsigma=0.9;
		if((nsigma>1.95)and(nsigma<2))nsigma=1.9;
		if((nsigma>2.95)and(nsigma<3))nsigma=2.9;
		if((nsigma>3.95)and(nsigma<4))nsigma=3.9;
		if((nsigma>4.95)and(nsigma<5))nsigma=4.9;
		if((nsigma>9.95)and(nsigma<10))nsigma=9.9;


                AFB_results[i]->SetMarkerSize(2);
		AFB_results[i]->SetBinContent(j,k,nsigma);

		//For text:
		if(nsigma<3.99){
		  xWhite[i].push_back(0.2225+0.10*(j-2));
		  yWhite[i].push_back(0.235+0.10*(k-2));
		  binWhite[i].push_back(nsigma);
		}
		else{
		  xBlack[i].push_back(0.2225+0.10*(j-2));
		  yBlack[i].push_back(0.235+0.10*(k-2));
		  binBlack[i].push_back(nsigma);
		}
	    }
        }
	
        AFB_results[i]->SetTitle(titles_AFB[i]);
	AFB_results[i]->SetTitleOffset(0.8);
        AFB_results[i]->GetZaxis()->SetRangeUser(0.0,7.0);
	AFB_results[i]->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
	AFB_results[i]->Draw("col");
	//AFB_results[i]->Draw("text same");
	DrawTextBins(xBlack[i],yBlack[i],binBlack[i],kBlack);
	DrawTextBins(xWhite[i],yWhite[i],binWhite[i],kWhite);
	Labels(energy_under_test,errortype,tpc_status,i,true,style);
	DrawTopLine();
	DrawLeg(energy_under_test);
	
    }
    
    TCanvas *c_AFB_short = new TCanvas ("c_AFB_short"+status_string,"c_AFB_short"+status_string,800,800);
    c_AFB_short->Divide(2,2);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerSize(2.);
    for(int i=0;i<4;i++){
      c_AFB_short->cd(i+1);
      if(i<2){
	AFB_results[i]->Draw("col");
	//AFB_results[i]->Draw("text same");
	DrawTextBins(xBlack[i],yBlack[i],binBlack[i],kBlack);
        DrawTextBins(xWhite[i],yWhite[i],binWhite[i],kWhite);
	Labels(energy_under_test,errortype,tpc_status,i,true,style);
      }
      else{
	AFB_results[i+1]->Draw("col");
	//AFB_results[i+1]->Draw("text same");
	DrawTextBins(xBlack[i+1],yBlack[i+1],binBlack[i+1],kBlack);
        DrawTextBins(xWhite[i+1],yWhite[i+1],binWhite[i+1],kWhite);
	Labels(energy_under_test,errortype,tpc_status,i+1,true,style);
      }
      //QQBARLabel2(0.60,0.225, "PRELIMINARY",kOrange+7,0.055);
      DrawTopLine();
      DrawLeg(energy_under_test);
    }
        

    TCanvas *c_AFB_comb = new TCanvas ("c_AFB_comb"+status_string,"c_AFB_comb"+status_string,800,800);
    c_AFB_comb->cd();
    if(energy_under_test==380){
      AFB_results[8]->SetFillColor(kOrange);
      AFB_results[8]->Draw("box");
      AFB_results[8]->Draw("text same");
    }
    else{
    AFB_results[8]->Draw("col");
    //AFB_results[8]->Draw("text same");
    DrawTextBins(xBlack[8],yBlack[8],binBlack[8],kBlack);
    DrawTextBins(xWhite[8],yWhite[8],binWhite[8],kWhite);
    }
    if(style=="paper"){
      //DrawOfficialTitle();
      AFB_results[8]->SetTitle("");
      QQBARLabel3(0.11,0.91,"Between-model discrimination power (#sigma-level)",kBlack,0.0425);
    }
    Labels(energy_under_test,errortype,tpc_status,8,true,style);
    DrawTopLine();
    if(energy_under_test!=380)DrawLeg(energy_under_test);
    
    
    TCanvas *c_AFB_unpol = new TCanvas ("c_AFB_unpol"+status_string,"c_AFB_unpol"+status_string,1200,400);
    c_AFB_unpol->Divide(3,1);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerSize(2.);
    std::vector<double> xunpolWhite[9];
    std::vector<double> yunpolWhite[9];
    std::vector<double> binunpolWhite[9];
    std::vector<double> xunpolBlack[9];
    std::vector<double> yunpolBlack[9];
    std::vector<double> binunpolBlack[9];

    for(int i=0;i<3;i++){
      int binx;
      int biny;
      c_AFB_unpol->cd(i+1);
      for(int j=1; j<theory.size(); j++) {
	for(int k=1; k<theory.size(); k++){
	  if(j>=k) continue;
	  biny=k;
	  binx=j;
	  double nsigma=0;
	  if(i==0){
	    //double prob=Sigma_To_Prob(probhistos_b.at(2)->GetBinContent(binx+1,biny+1));
	    //nsigma=Prob_To_Sigma(prob);
	    nsigma=probhistos_b.at(2)->GetBinContent(binx+1,biny+1);
	  }
	  else if(i==1){
	    double nsigmas=probhistos_c.at(2)->GetBinContent(binx+1,biny+1);
            //double prob=Sigma_To_Prob(nsigmas);
            //nsigma=Prob_To_Sigma(prob);
	    nsigma=probhistos_c.at(2)->GetBinContent(binx+1,biny+1);
	  }
	  else{
	    double nsigma1;
	    double nsigma2;
	    nsigma1=probhistos_b.at(2)->GetBinContent(binx+1,biny+1);
	    nsigma2=probhistos_c.at(2)->GetBinContent(binx+1,biny+1);
	    double prob1=Sigma_To_Prob(nsigma1);
	    double prob2=Sigma_To_Prob(nsigma2);
	    double prob=prob1*prob2;
	    nsigma=Prob_To_Sigma(prob);
	  }
	  // Rounding for plots:                                                                                                                                                                              
	  if((nsigma>0.95)and(nsigma<1))nsigma=0.9;
	  if((nsigma>1.95)and(nsigma<2))nsigma=1.9;
	  if((nsigma>2.95)and(nsigma<3))nsigma=2.9;
	  if((nsigma>3.95)and(nsigma<4))nsigma=3.9;
	  if((nsigma>4.95)and(nsigma<5))nsigma=4.9;
	  if((nsigma>9.95)and(nsigma<10))nsigma=9.9;

	  AFB_unpol[i]->SetMarkerSize(2);
	  AFB_unpol[i]->SetBinContent(j,k,nsigma);
	  //For text:                        
	  if(nsigma<3.99){
	    xunpolWhite[i].push_back(0.2225+0.10*(j-2));
	    yunpolWhite[i].push_back(0.235+0.10*(k-2));
	    binunpolWhite[i].push_back(nsigma);
	  }
	  else{
	    xunpolBlack[i].push_back(0.2225+0.10*(j-2));
	    yunpolBlack[i].push_back(0.235+0.10*(k-2));
	    binunpolBlack[i].push_back(nsigma);
	  }
	}
      }
      AFB_unpol[i]->SetTitle(titles_AFB_unpol[i]);
      AFB_unpol[i]->GetZaxis()->SetRangeUser(0.0,7.0);
      AFB_unpol[i]->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      AFB_unpol[i]->Draw("col");
      
      DrawTextBins(xunpolBlack[i],yunpolBlack[i],binunpolBlack[i],kBlack);
      DrawTextBins(xunpolWhite[i],yunpolWhite[i],binunpolWhite[i],kWhite);
      Labels(energy_under_test,errortype,tpc_status,i,false,style);
      DrawTopLine();
      DrawLeg(energy_under_test);
    }

    TCanvas *c_AFB_unpol_comb = new TCanvas ("c_AFB_unpol_comb"+status_string,"c_AFB_unpol_comb"+status_string,800,800);
    
    AFB_unpol[2]->Draw("col");
    //AFB_unpol[2]->Draw("text same");
    //QQBARLabel2(0.60,0.225, "PRELIMINARY",kOrange+7,0.055);
    if(style=="paper"){
      //DrawOfficialTitle();                                                   
      AFB_unpol[2]->SetTitle("");
      QQBARLabel3(0.11,0.91,"Between-model discrimination power (#sigma-level)",kBlack,0.0425);
    }
    DrawTextBins(xunpolBlack[2],yunpolBlack[2],binunpolBlack[2],kBlack);
    DrawTextBins(xunpolWhite[2],yunpolWhite[2],binunpolWhite[2],kWhite);
    Labels(energy_under_test,errortype,tpc_status,2,false,style);
    DrawTopLine();
    DrawLeg(energy_under_test);
    
    // Save everything:
    TString savingformat[2]={".png",".eps"};
    for(int isave=0; isave<2 ;isave++){
      if(style=="proceedings"){
	//c_AFB_full->SaveAs("c_AFB_full"+status_string+savingformat[isave]);
	//c_AFB_short->SaveAs("c_AFB_short"+status_string+savingformat[isave]);
	//c_AFB_comb->SaveAs("c_AFB_comb"+status_string+savingformat[isave]);
	//c_AFB_unpol->SaveAs("c_AFB_unpol"+status_string+savingformat[isave]);
	//c_AFB_unpol_comb->SaveAs("c_AFB_unpol_comb"+status_string+savingformat[isave]);
      }
      else if(style=="paper"){
	//if(energy_under_test==91.2);
	//c_AFB_short->SaveAs("c_AFB_short"+status_string+savingformat[isave]);
	c_AFB_comb->SaveAs("c_AFB_comb"+status_string+savingformat[isave]);
	//c_AFB_unpol_comb->SaveAs("c_AFB_unpol_comb"+status_string+savingformat[isave]);
	if(energy_under_test==250)c_AFB_unpol_comb->SaveAs("c_AFB_unpol_comb"+status_string+savingformat[isave]);
      }
    }
    
}

void test_two_quarks_two_energies(TString errortype="Stat", TString tpc_status="dNdx", TString style="intern")
{
  
  //TString quarkst=TString::Format("c&b-quark";
    TString energyst="250&500 GeV";

    read_all_models(false);
    // stores the information for all  models (named in std::vector<TString> allmodels={"SMA","A1"};)
    // with 100% polarization options only
    // it stores the same info as the txt files
    // in a vector std::vector<model_struct_t> theory;
    // to access the info, you do: theory.at(i).model_values[ie].AFB_L[ifl]
    // indexes:
    //   i: as many as model strings
    //   ie: index for the energy: one index for every energy found in the file. If you want to know the value, you can acces to i: theory.at(i).model_values[ie].energy
    //       we can store up to 5 energy points
    //   ifl: for the flavour , 0 is empty, 1,2,3,4,5 are d,u,s,c,b


    // here I create the thheoretical observables for the desired polarization scheme
    // I save a vector for all energies, all pol (2 in this case) and and the b and c flavours
    std::vector<observables_struct_t> observables=create_observables();

    // I calculate the prob matrix of tested models vs hypothesis models, using th2f
    std::vector<TH2F*> probhistos_c=nsigmas_models(4, 250, observables, errortype, tpc_status) ;
    std::vector<TH2F*> probhistos_b=nsigmas_models(5, 250, observables, errortype, tpc_status) ;
    std::vector<TH2F*> probhistos_c500=nsigmas_models(4, 500, observables, errortype, tpc_status) ;
    std::vector<TH2F*> probhistos_b500=nsigmas_models(5, 500, observables, errortype, tpc_status) ;
    const Int_t xyNBINS = 8;
    Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
    TH2F * AFB_results[9];
    AFB_results[0] = new TH2F("AFB_b_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), b-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[1] = new TH2F("AFB_b_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[2] = new TH2F("AFB_b_comb",TString::Format("AFB_ILD_H20, b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[3] = new TH2F("AFB_c_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[4] = new TH2F("AFB_c_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[5] = new TH2F("AFB_c_comb",TString::Format("AFB_ILD_H20, c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[6] = new TH2F("AFB_b_&_AFB_c_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[7] = new TH2F("AFB_b_&_AFB_c_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    AFB_results[8] = new TH2F("AFB_b_&_AFB_c_comb",TString::Format("AFB_ILD_H20, b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);

    TH2F * R_results[9];
    R_results[0] = new TH2F("R_b_L",TString::Format("R_ILD_H20 (-0.8,+0.3), b-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[1] = new TH2F("R_b_R",TString::Format("R_ILD_H20 (+0.8,-0.3), b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[2] = new TH2F("R_b_comb",TString::Format("R_ILD_H20, b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[3] = new TH2F("R_c_L",TString::Format("R_ILD_H20 (-0.8,+0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[4] = new TH2F("R_c_R",TString::Format("R_ILD_H20 (+0.8,-0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[5] = new TH2F("R_c_comb",TString::Format("R_ILD_H20, c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[6] = new TH2F("R_b_&_R_c_L",TString::Format("R_ILD_H20 (-0.8,+0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[7] = new TH2F("R_b_&_R_c_R",TString::Format("R_ILD_H20 (+0.8,-0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
    R_results[8] = new TH2F("R_b_&_R_c_comb",TString::Format("R_ILD_H20, b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);

    TString conditions[]={TString::Format("H20, P:(-0.8,+0.3), b-quark #sqrt{s}=%s",energyst.Data()),TString::Format("H20, b-quark P:(+0.8,-0.3),  #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b-quark #sqrt{s}=%s",energyst.Data()),
			  TString::Format("H20, P:(-0.8,+0.3), c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., c-quark #sqrt{s}=%s",energyst.Data()),
			  TString::Format("H20, P:(-0.8,+0.3), b&c quarks #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), b&c quarks #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b&c quarks #sqrt{s}=%s",energyst.Data())};
    TString conditions_unpol[]={TString::Format("unpol., b-quark #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., b&c-quark #sqrt{s}=%s",energyst.Data())};
    TString titles_AFB[]={"AFB_{b} (e^{-}_{L}e^{+}_{R})","AFB_{b} (e^{-}_{R}e^{+}_{L})","AFB_{b} (Both pol.)","AFB_{c} (e^{-}_{L}e^{+}_{R})","AFB_{c} (e^{-}_{R}e^{+}_{L})","AFB_{c} (Both pol.)","AFB_{b} and AFB_{c} (e^{-}_{L}e^{+}_{R})","AFB_{b} and AFB_{c} (e^{-}_{R}e^{+}_{L})","AFB_{b} & AFB_{c} (Both pol.)"};
    TString titles_AFB_unpol[]={"AFB_{b} (unpol.)","AFB_{c} (unpol.)","AFB_{b} & AFB_{c} (unpol.)"};

    TString titles_R[]={"R_{b} (e^{-}_{L}e^{+}_{R})","R_{b} (e^{-}_{R}e^{+}_{L})","R_{b} (Both pol.)","R_{c} (e^{-}_{L}e^{+}_{R})","R_{c} (e^{-}_{R}e^{+}_{L})","R_{c} (Both pol.)","R_{b} and R_{c} (e^{-}_{L}e^{+}_{R})","R_{b} and R_{c} (e^{-}_{R}e^{+}_{L})","R_{b} and R_{c} (Both pol.)"};
    TString titles_R_unpol[]={"R_{b} (unpol.)","R_{c} (unpol.)","R_{b} & R_{c} (unpol.)"};

    TString title_models_x[]={"SMB","A_{1}","A_{2}","B_{1}^{-}","B_{1}^{+}","B_{2}^{-}","B_{2}^{+}","B_{3}^{-}","B_{3}^{+}"};
    TString title_models_y[]={"SMB","A_{1}","A_{2}","B_{1}^{-}","B_{1}^{+}","B_{2}^{-}","B_{2}^{+}","B_{3}^{-}","B_{3}^{+}"};

    for(int i=0; i<9; i++) {
      for(int j=1; j<theory.size(); j++) {
	AFB_results[i]->GetXaxis()->SetBinLabel(j,title_models_x[j]);
	AFB_results[i]->GetYaxis()->SetBinLabel(j,title_models_y[j]);
	R_results[i]->GetXaxis()->SetBinLabel(j,title_models_x[j]);
	R_results[i]->GetYaxis()->SetBinLabel(j,title_models_y[j]);
	
	AFB_results[i]->GetXaxis()->SetTickLength(0.);
	AFB_results[i]->GetYaxis()->SetTickLength(0.);
	AFB_results[i]->GetYaxis()->SetLabelSize(0.07);
	AFB_results[i]->GetXaxis()->SetLabelSize(0.07);
	R_results[i]->GetXaxis()->SetTickLength(0.);
	R_results[i]->GetYaxis()->SetTickLength(0.);
	R_results[i]->GetYaxis()->SetLabelSize(0.07);
	R_results[i]->GetXaxis()->SetLabelSize(0.07);
      }
    }

    //Status string:
    TString errorstring;
    if(errortype=="Stat")errorstring="Stat";
    else if(errortype=="StatSyst")errorstring="StatSyst";
    else if(errortype=="StatTheoCurrent")errorstring="StatTheoCurrent";
    else if(errortype=="StatTheoGigaZ")errorstring="StatTheoGigaZ";
    else if(errortype=="StatTheo250")errorstring="StatTheo250";
    TString status_string_both="_"+tpc_status+"_"+errorstring;

    Int_t colors[6];
    gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
    Double_t levels[7]={0,0.99,1.99,2.99,3.99,4.99,1000.0};
    gStyle->SetPaintTextFormat(".1f");
    

    TCanvas *c_AFB_full_both = new TCanvas ("c_AFB_full_both"+status_string_both,"c_AFB_full_both"+status_string_both,800,800);
    c_AFB_full_both->Divide(3,3);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerSize(2.);
    std::vector<double> xWhite[9];
    std::vector<double> yWhite[9];
    std::vector<double> binWhite[9];
    std::vector<double> xBlack[9];
    std::vector<double> yBlack[9];
    std::vector<double> binBlack[9];
    for(int i=0; i<9; i++) {
      c_AFB_full_both->cd(i+1);
      int binx;
      int biny;
      for(int j=1; j<theory.size(); j++) {
	for(int k=1; k<theory.size(); k++){
	  if(j>=k) continue;
	  biny=k;
	  binx=j;
	  double nsigma=0;
	  if(i<2){
	    double prob1=Sigma_To_Prob(probhistos_b.at(i)->GetBinContent(binx+1,biny+1));
            double prob2=Sigma_To_Prob(probhistos_b500.at(i)->GetBinContent(binx+1,biny+1));
            double prob=prob1*prob2;
            nsigma=Prob_To_Sigma(prob);
	  }
	  else if(i==2){
	    double prob1=Sigma_To_Prob(probhistos_b.at(0)->GetBinContent(binx+1,biny+1));
            double prob2=Sigma_To_Prob(probhistos_b500.at(0)->GetBinContent(binx+1,biny+1));
            double prob3=Sigma_To_Prob(probhistos_b.at(1)->GetBinContent(binx+1,biny+1));
            double prob4=Sigma_To_Prob(probhistos_b500.at(1)->GetBinContent(binx+1,biny+1));
	    double prob=prob1*prob2*prob3*prob4;
            nsigma=Prob_To_Sigma(prob);
	  }
	  else if(i<5){
	    double prob1=Sigma_To_Prob(probhistos_c.at(i-3)->GetBinContent(binx+1,biny+1));
            double prob2=Sigma_To_Prob(probhistos_c500.at(i-3)->GetBinContent(binx+1,biny+1));
            double prob=prob1*prob2;
            nsigma=Prob_To_Sigma(prob);
	  }
	  else if(i==5){
	    double prob1=Sigma_To_Prob(probhistos_c.at(0)->GetBinContent(binx+1,biny+1));
            double prob2=Sigma_To_Prob(probhistos_c500.at(0)->GetBinContent(binx+1,biny+1));
            double prob3=Sigma_To_Prob(probhistos_c.at(1)->GetBinContent(binx+1,biny+1));
            double prob4=Sigma_To_Prob(probhistos_c500.at(1)->GetBinContent(binx+1,biny+1));
            double prob=prob1*prob2*prob3*prob4;
            nsigma=Prob_To_Sigma(prob);
	  }
	  else if(i<8){
	    double prob1=Sigma_To_Prob(probhistos_b.at(i-6)->GetBinContent(binx+1,biny+1));
            double prob2=Sigma_To_Prob(probhistos_b500.at(i-6)->GetBinContent(binx+1,biny+1));
	    double prob3=Sigma_To_Prob(probhistos_c.at(i-6)->GetBinContent(binx+1,biny+1));
            double prob4=Sigma_To_Prob(probhistos_c500.at(i-6)->GetBinContent(binx+1,biny+1));
            double prob=prob1*prob2*prob3*prob4;
            nsigma=Prob_To_Sigma(prob);
	  }
	  else{
	    double prob1=Sigma_To_Prob(probhistos_b.at(0)->GetBinContent(binx+1,biny+1));
            double prob2=Sigma_To_Prob(probhistos_b500.at(0)->GetBinContent(binx+1,biny+1));
            double prob3=Sigma_To_Prob(probhistos_c.at(0)->GetBinContent(binx+1,biny+1));
            double prob4=Sigma_To_Prob(probhistos_c500.at(0)->GetBinContent(binx+1,biny+1));
	    double prob5=Sigma_To_Prob(probhistos_b.at(1)->GetBinContent(binx+1,biny+1));
            double prob6=Sigma_To_Prob(probhistos_b500.at(1)->GetBinContent(binx+1,biny+1));
            double prob7=Sigma_To_Prob(probhistos_c.at(1)->GetBinContent(binx+1,biny+1));
            double prob8=Sigma_To_Prob(probhistos_c500.at(1)->GetBinContent(binx+1,biny+1));
            double prob=prob1*prob2*prob3*prob4*prob5*prob6*prob7*prob8;
	    nsigma=Prob_To_Sigma(prob);
	    cout<<"Sigmas: "<<nsigma<<". Prob: "<<prob<<endl;
	    //cout<<"Probs.: "<<prob1<<" "<<prob2<<" "<<prob3<<" "<<prob4<<" "<<prob5<<" "<<prob6<<" "<<prob7<<" "<<prob8<<endl;
	    //cout<<"Sigma level: "<<nsigma<<". Combined prob.: "<<prob<<endl;
	  }
	  // Rounding for plots:                                                                                                                                                                              
	  if((nsigma>0.95)and(nsigma<1))nsigma=0.9;
	  if((nsigma>1.95)and(nsigma<2))nsigma=1.9;
	  if((nsigma>2.95)and(nsigma<3))nsigma=2.9;
	  if((nsigma>3.95)and(nsigma<4))nsigma=3.9;
	  if((nsigma>4.95)and(nsigma<5))nsigma=4.9;
	  if((nsigma>9.95)and(nsigma<10))nsigma=9.9;

	  AFB_results[i]->SetMarkerSize(2);
	  AFB_results[i]->SetBinContent(j,k,nsigma);
	  //For text:                                                                  
	  if(nsigma<3.99){
	    xWhite[i].push_back(0.2225+0.10*(j-2));
	    yWhite[i].push_back(0.235+0.10*(k-2));
	    binWhite[i].push_back(nsigma);
	  }
	  else{
	    xBlack[i].push_back(0.2225+0.10*(j-2));
	    yBlack[i].push_back(0.235+0.10*(k-2));
	    binBlack[i].push_back(nsigma);
	  }
	}
      }
      AFB_results[i]->SetTitle(titles_AFB[i]);
      AFB_results[i]->GetZaxis()->SetRangeUser(0.0,7.0);
      AFB_results[i]->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
      AFB_results[i]->Draw("col");
      
      //AFB_results[i]->Draw("text same");
      //QQBARLabel2(0.60,0.225, "PRELIMINARY",kOrange+7,0.055);
      DrawTextBins(xBlack[i],yBlack[i],binBlack[i],kBlack);
      DrawTextBins(xWhite[i],yWhite[i],binWhite[i],kWhite);
      Labels(9999,errortype,tpc_status,i,true,style);
      DrawTopLine();
      DrawLeg(9999);
    }

    TCanvas *c_AFB_short_both = new TCanvas ("c_AFB_short_both"+status_string_both,"c_AFB_short"+status_string_both,800,800);
    c_AFB_short_both->Divide(2,2);
    gStyle->SetOptStat(0);
    gStyle->SetMarkerSize(2.);
    for(int i=0;i<4;i++){
      c_AFB_short_both->cd(i+1);
      if(i<2){
        AFB_results[i]->Draw("col");
	//AFB_results[i]->Draw("text same");
	DrawTextBins(xBlack[i],yBlack[i],binBlack[i],kBlack);
	DrawTextBins(xWhite[i],yWhite[i],binWhite[i],kWhite);
	Labels(9999,errortype,tpc_status,i,true,style);
      }
      else{
        AFB_results[i+1]->Draw("col");
	//AFB_results[i+1]->Draw("text same");
	DrawTextBins(xBlack[i+1],yBlack[i+1],binBlack[i+1],kBlack);
	DrawTextBins(xWhite[i+1],yWhite[i+1],binWhite[i+1],kWhite);
	Labels(9999,errortype,tpc_status,i+1,true,style);
      }
      //QQBARLabel2(0.60,0.225, "PRELIMINARY",kOrange+7,0.055); 
    
      DrawTopLine();
      DrawLeg(9999);
    }
							    
    TCanvas *c_AFB_comb_both = new TCanvas ("c_AFB_comb_both"+status_string_both,"c_AFB_comb_both"+status_string_both,800,800);
    c_AFB_comb_both->cd();
    AFB_results[8]->Draw("col");
    //AFB_results[8]->Draw("text same");
    //QQBARLabel2(0.60,0.225, "PRELIMINARY",kOrange+7,0.055);
    if(style=="paper"){
      //DrawOfficialTitle();                                                                                                       
      AFB_results[8]->SetTitle("");
      QQBARLabel3(0.11,0.91,"Between-model discrimination power (#sigma-level)",kBlack,0.0425);
    }

    DrawTextBins(xBlack[8],yBlack[8],binBlack[8],kBlack);
    DrawTextBins(xWhite[8],yWhite[8],binWhite[8],kWhite);
    Labels(9999,errortype,tpc_status,8,true,style);
    DrawTopLine();
    DrawLeg(9999);
    
    //Save everything: 
    TString savingformat[2]={".png",".eps"};
    for(int isave=0; isave<2 ;isave++){
      //c_AFB_full_both->SaveAs("c_AFB_full_both"+status_string_both+savingformat[isave]);
      //c_AFB_short_both->SaveAs("c_AFB_short_both"+status_string_both+savingformat[isave]);
      c_AFB_comb_both->SaveAs("c_AFB_comb_both"+status_string_both+savingformat[isave]);
    }
    
}

void test_two_quarks_three_energies(TString errortype="Stat", TString tpc_status="dNdx", TString style="intern")
{

  TString energyst="250&500&1000 GeV";

  read_all_models(false);

  std::vector<observables_struct_t> observables=create_observables();

  // I calculate the prob matrix of tested models vs hypothesis models, using th2f
  std::vector<TH2F*> probhistos_c=nsigmas_models(4, 250, observables, errortype, tpc_status) ;
  std::vector<TH2F*> probhistos_b=nsigmas_models(5, 250, observables, errortype, tpc_status) ;
  std::vector<TH2F*> probhistos_c500=nsigmas_models(4, 500, observables, errortype, tpc_status) ;
  std::vector<TH2F*> probhistos_b500=nsigmas_models(5, 500, observables, errortype, tpc_status) ;
  std::vector<TH2F*> probhistos_c1000=nsigmas_models(4, 1000, observables, errortype, tpc_status) ;
  std::vector<TH2F*> probhistos_b1000=nsigmas_models(5, 1000, observables, errortype, tpc_status) ;
  const Int_t xyNBINS = 8;
  Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};
  TH2F * AFB_results[9];
  AFB_results[0] = new TH2F("AFB_b_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), b-quark #sqrt{s}=%s; ; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
  AFB_results[1] = new TH2F("AFB_b_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
  AFB_results[2] = new TH2F("AFB_b_comb",TString::Format("AFB_ILD_H20, b-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
  AFB_results[3] = new TH2F("AFB_c_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
  AFB_results[4] = new TH2F("AFB_c_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
  AFB_results[5] = new TH2F("AFB_c_comb",TString::Format("AFB_ILD_H20, c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
  AFB_results[6] = new TH2F("AFB_b_&_AFB_c_L",TString::Format("AFB_ILD_H20 (-0.8,+0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
  AFB_results[7] = new TH2F("AFB_b_&_AFB_c_R",TString::Format("AFB_ILD_H20 (+0.8,-0.3), b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);
  AFB_results[8] = new TH2F("AFB_b_&_AFB_c_comb",TString::Format("AFB_ILD_H20, b&c-quark #sqrt{s}=%s; ; ;",energyst.Data()),xyNBINS,xyedges,xyNBINS,xyedges);

  TString conditions[]={TString::Format("H20, P:(-0.8,+0.3), b-quark #sqrt{s}=%s",energyst.Data()),TString::Format("H20, b-quark P:(+0.8,-0.3),  #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b-quark #sqrt{s}=%s",energyst.Data()),
			TString::Format("H20, P:(-0.8,+0.3), c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., c-quark #sqrt{s}=%s",energyst.Data()),
			TString::Format("H20, P:(-0.8,+0.3), b&c quarks #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), b&c quarks #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b&c quarks #sqrt{s}=%s",energyst.Data())};
  TString conditions_unpol[]={TString::Format("unpol., b-quark #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., c-quark #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., b&c-quark #sqrt{s}=%s",energyst.Data())};
  TString titles_AFB[]={"AFB_{b} (e^{-}_{L}e^{+}_{R})","AFB_{b} (e^{-}_{R}e^{+}_{L})","AFB_{b} (Both pol.)","AFB_{c} (e^{-}_{L}e^{+}_{R})","AFB_{c} (e^{-}_{R}e^{+}_{L})","AFB_{c} (Both pol.)","AFB_{b} and AFB_{c} (e^{-}_{L}e^{+}_{R})","AFB_{b} and AFB_{c} (e^{-}_{R}e^{+}_{L})","AFB_{b} & AFB_{c} (Both pol.)"};
  TString titles_AFB_unpol[]={"AFB_{b} (unpol.)","AFB_{c} (unpol.)","AFB_{b} & AFB_{c} (unpol.)"};
  
  TString title_models_x[]={"SMB","A_{1}","A_{2}","B_{1}^{-}","B_{1}^{+}","B_{2}^{-}","B_{2}^{+}","B_{3}^{-}","B_{3}^{+}"};
  TString title_models_y[]={"SMB","A_{1}","A_{2}","B_{1}^{-}","B_{1}^{+}","B_{2}^{-}","B_{2}^{+}","B_{3}^{-}","B_{3}^{+}"};

  for(int i=0; i<9; i++) {
    for(int j=1; j<theory.size(); j++) {
      AFB_results[i]->GetXaxis()->SetBinLabel(j,title_models_x[j]);
      AFB_results[i]->GetYaxis()->SetBinLabel(j,title_models_y[j]);
      
      AFB_results[i]->GetXaxis()->SetTickLength(0.);
      AFB_results[i]->GetYaxis()->SetTickLength(0.);
      AFB_results[i]->GetYaxis()->SetLabelSize(0.075);
      AFB_results[i]->GetXaxis()->SetLabelSize(0.075);
    }
  }

  //Status string:
  TString errorstring;
  if(errortype=="Stat")errorstring="Stat";
  else if(errortype=="StatSyst")errorstring="StatSyst";
  else if(errortype=="StatTheoCurrent")errorstring="StatTheoCurrent";
  else if(errortype=="StatTheoGigaZ")errorstring="StatTheoGigaZ";
  else if(errortype=="StatTheo250")errorstring="StatTheo250";
  TString status_string_three="_"+tpc_status+"_"+errorstring;

  Int_t colors[6];
  gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
  Double_t levels[7]={0,0.99,1.99,2.99,3.99,4.99,1000.0};
  gStyle->SetPaintTextFormat(".1f");


  TCanvas *c_AFB_full_three = new TCanvas ("c_AFB_full_three"+status_string_three,"c_AFB_full_three"+status_string_three,800,800);
  c_AFB_full_three->Divide(3,3);
  gStyle->SetOptStat(0);
  gStyle->SetMarkerSize(2.);
  std::vector<double> xWhite[9];
  std::vector<double> yWhite[9];
  std::vector<double> binWhite[9];
  std::vector<double> xBlack[9];
  std::vector<double> yBlack[9];
  std::vector<double> binBlack[9];
  for(int i=0; i<9; i++) {
    c_AFB_full_three->cd(i+1);
    int binx;
    int biny;
    for(int j=1; j<theory.size(); j++) {
      for(int k=1; k<theory.size(); k++){
	if(j>=k) continue;
	biny=k;
	binx=j;
	double nsigma=0;
	if(i<2){
	  double prob1=Sigma_To_Prob(probhistos_b.at(i)->GetBinContent(binx+1,biny+1));
	  double prob2=Sigma_To_Prob(probhistos_b500.at(i)->GetBinContent(binx+1,biny+1));
	  double prob3=Sigma_To_Prob(probhistos_b1000.at(i)->GetBinContent(binx+1,biny+1));
	  double prob=prob1*prob2*prob3;
	  nsigma=Prob_To_Sigma(prob);
	}
	else if(i==2){
	  double prob1=Sigma_To_Prob(probhistos_b.at(0)->GetBinContent(binx+1,biny+1));
	  double prob2=Sigma_To_Prob(probhistos_b500.at(0)->GetBinContent(binx+1,biny+1));
	  double prob3=Sigma_To_Prob(probhistos_b1000.at(0)->GetBinContent(binx+1,biny+1));
	  double prob11=Sigma_To_Prob(probhistos_b.at(1)->GetBinContent(binx+1,biny+1));
	  double prob22=Sigma_To_Prob(probhistos_b500.at(1)->GetBinContent(binx+1,biny+1));
	  double prob33=Sigma_To_Prob(probhistos_b1000.at(1)->GetBinContent(binx+1,biny+1));
	  double prob=prob1*prob2*prob3*prob11*prob22*prob33;
	  nsigma=Prob_To_Sigma(prob);
	}
	else if(i<5){
	  double prob1=Sigma_To_Prob(probhistos_c.at(i-3)->GetBinContent(binx+1,biny+1));
	  double prob2=Sigma_To_Prob(probhistos_c500.at(i-3)->GetBinContent(binx+1,biny+1));
	  double prob3=Sigma_To_Prob(probhistos_c1000.at(i-3)->GetBinContent(binx+1,biny+1));
	  double prob=prob1*prob2*prob3;
	  nsigma=Prob_To_Sigma(prob);
	}
	else if(i==5){
	  double prob1=Sigma_To_Prob(probhistos_c.at(0)->GetBinContent(binx+1,biny+1));
	  double prob2=Sigma_To_Prob(probhistos_c500.at(0)->GetBinContent(binx+1,biny+1));
	  double prob3=Sigma_To_Prob(probhistos_c1000.at(0)->GetBinContent(binx+1,biny+1));
	  double prob11=Sigma_To_Prob(probhistos_c.at(1)->GetBinContent(binx+1,biny+1));
	  double prob22=Sigma_To_Prob(probhistos_c500.at(1)->GetBinContent(binx+1,biny+1));
	  double prob33=Sigma_To_Prob(probhistos_c1000.at(1)->GetBinContent(binx+1,biny+1));
	  double prob=prob1*prob2*prob3*prob11*prob22*prob33;
	  nsigma=Prob_To_Sigma(prob);
	}
	else if(i<8){
	  double prob1=Sigma_To_Prob(probhistos_b.at(i-6)->GetBinContent(binx+1,biny+1));
	  double prob2=Sigma_To_Prob(probhistos_b500.at(i-6)->GetBinContent(binx+1,biny+1));
	  double prob3=Sigma_To_Prob(probhistos_b1000.at(i-6)->GetBinContent(binx+1,biny+1));
	  double prob11=Sigma_To_Prob(probhistos_c.at(i-6)->GetBinContent(binx+1,biny+1));
	  double prob22=Sigma_To_Prob(probhistos_c500.at(i-6)->GetBinContent(binx+1,biny+1));
	  double prob33=Sigma_To_Prob(probhistos_c1000.at(i-6)->GetBinContent(binx+1,biny+1));
	  double prob=prob1*prob2*prob3*prob11*prob22*prob33;
	  nsigma=Prob_To_Sigma(prob);
	}
	else{
	  double prob1=Sigma_To_Prob(probhistos_b.at(0)->GetBinContent(binx+1,biny+1));
	  double prob2=Sigma_To_Prob(probhistos_b500.at(0)->GetBinContent(binx+1,biny+1));
	  double prob3=Sigma_To_Prob(probhistos_b1000.at(0)->GetBinContent(binx+1,biny+1));
	  double prob11=Sigma_To_Prob(probhistos_c.at(0)->GetBinContent(binx+1,biny+1));
	  double prob22=Sigma_To_Prob(probhistos_c500.at(0)->GetBinContent(binx+1,biny+1));
	  double prob33=Sigma_To_Prob(probhistos_c1000.at(0)->GetBinContent(binx+1,biny+1));
	  double prob111=Sigma_To_Prob(probhistos_b.at(1)->GetBinContent(binx+1,biny+1));
	  double prob222=Sigma_To_Prob(probhistos_b500.at(1)->GetBinContent(binx+1,biny+1));
	  double prob333=Sigma_To_Prob(probhistos_b1000.at(1)->GetBinContent(binx+1,biny+1));
	  double prob1111=Sigma_To_Prob(probhistos_c.at(1)->GetBinContent(binx+1,biny+1));
	  double prob2222=Sigma_To_Prob(probhistos_c500.at(1)->GetBinContent(binx+1,biny+1));
	  double prob3333=Sigma_To_Prob(probhistos_c1000.at(1)->GetBinContent(binx+1,biny+1));
	  double prob=prob1*prob2*prob3*prob11*prob22*prob33*prob111*prob222*prob333*prob1111*prob2222*prob3333;
	  nsigma=Prob_To_Sigma(prob);
	  cout<<"Sigmas: "<<nsigma<<". Prob: "<<prob<<endl;
	  cout<<"Probs: "<<prob1<<", "<<prob2<<", "<<prob3<<", "<<prob11<<", "<<prob22<<", "<<prob33<<", "<<prob111<<", "<<prob222<<", "<<prob333<<", "<<prob1111<<", "<<prob2222<<", "<<prob3333<<endl;
	  
	}
	// Rounding for plots:                                                                                                                                                                              
	if((nsigma>0.95)and(nsigma<1))nsigma=0.9;
	if((nsigma>1.95)and(nsigma<2))nsigma=1.9;
	if((nsigma>2.95)and(nsigma<3))nsigma=2.9;
	if((nsigma>3.95)and(nsigma<4))nsigma=3.9;
	if((nsigma>4.95)and(nsigma<5))nsigma=4.9;
	if((nsigma>9.95)and(nsigma<10))nsigma=9.9;

	AFB_results[i]->SetMarkerSize(2);
	AFB_results[i]->SetBinContent(j,k,nsigma);

	//For text:
	if(nsigma<3.99){
	  xWhite[i].push_back(0.2225+0.10*(j-2));
	  yWhite[i].push_back(0.235+0.10*(k-2));
	  binWhite[i].push_back(nsigma);
	}
	else{
	  xBlack[i].push_back(0.2225+0.10*(j-2));
	  yBlack[i].push_back(0.235+0.10*(k-2));
	  binBlack[i].push_back(nsigma);
	}
      }
    }
    AFB_results[i]->SetTitle(titles_AFB[i]);
    AFB_results[i]->GetZaxis()->SetRangeUser(0.0,7.0);
    AFB_results[i]->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
    AFB_results[i]->Draw("col");

    DrawTextBins(xBlack[i],yBlack[i],binBlack[i],kBlack);
    DrawTextBins(xWhite[i],yWhite[i],binWhite[i],kWhite);
    Labels(10000,errortype,tpc_status,i,true,style);
    DrawTopLine();
    DrawLeg(10000);
  }

  TCanvas *c_AFB_short_three = new TCanvas ("c_AFB_short_three"+status_string_three,"c_AFB_short"+status_string_three,800,800);
  c_AFB_short_three->Divide(2,2);
  gStyle->SetOptStat(0);
  gStyle->SetMarkerSize(2.);
  for(int i=0;i<4;i++){
    c_AFB_short_three->cd(i+1);
    if(i<2){
      AFB_results[i]->Draw("col");

      DrawTextBins(xBlack[i],yBlack[i],binBlack[i],kBlack);
      DrawTextBins(xWhite[i],yWhite[i],binWhite[i],kWhite);
      Labels(10000,errortype,tpc_status,i,true,style);
    }
    else{
      AFB_results[i+1]->Draw("col");

      DrawTextBins(xBlack[i+1],yBlack[i+1],binBlack[i+1],kBlack);
      DrawTextBins(xWhite[i+1],yWhite[i+1],binWhite[i+1],kWhite);
      Labels(10000,errortype,tpc_status,i+1,true,style);
    }

    DrawTopLine();
    DrawLeg(10000);
  }

  TCanvas *c_AFB_comb_three = new TCanvas ("c_AFB_comb_three"+status_string_three,"c_AFB_comb_three"+status_string_three,800,800);
  c_AFB_comb_three->cd();
  AFB_results[8]->Draw("col");
  //AFB_results[8]->Draw("text same");
  //QQBARLabel2(0.60,0.225, "PRELIMINARY",kOrange+7,0.055);
  if(style=="paper"){
    //DrawOfficialTitle();                                                                                                                                                 
    AFB_results[8]->SetTitle("");
    QQBARLabel3(0.11,0.91,"Between-model discrimination power (#sigma-level)",kBlack,0.0425);
  }
  DrawTextBins(xBlack[8],yBlack[8],binBlack[8],kBlack);
  DrawTextBins(xWhite[8],yWhite[8],binWhite[8],kWhite);
  Labels(10000,errortype,tpc_status,8,true,style);
  DrawTopLine();
  DrawLeg(10000);

  //Save everything:
  TString savingformat[2]={".png",".eps"};
  for(int isave=0; isave<2 ;isave++){
    //c_AFB_full_three->SaveAs("c_AFB_full_three"+status_string_three+savingformat[isave]);
    //c_AFB_short_three->SaveAs("c_AFB_short_three"+status_string_three+savingformat[isave]);
    c_AFB_comb_three->SaveAs("c_AFB_comb_three"+status_string_three+savingformat[isave]);
  }

}


void test_models(TString energy, TString errortype, string tpc_status, TString style){
  if(energy=="250") test_one_energy(250, errortype, tpc_status, style);
  else if(energy=="91.2") test_one_energy(91.2, errortype, tpc_status, style);
  else if(energy=="350") test_one_energy(350, errortype, tpc_status, style);
  else if(energy=="380") test_one_energy(380, errortype, tpc_status, style);
  else if(energy=="500") test_one_energy(500, errortype, tpc_status, style);
  else if(energy=="1000") test_one_energy(1000, errortype, tpc_status, style);
  else if(energy=="both") test_two_quarks_two_energies(errortype, tpc_status, style);
  else if(energy=="three") test_two_quarks_three_energies(errortype, tpc_status, style);
}

