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

void Labels(double energy_under_test=250, TString errortype="Stat", TString tpc_status="default", TString SM_status="nominal", int iter=8, bool pol=true, TString style="proceedings")
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
  else if(energy_under_test==10000) energypaper="ILC250+500+1000";
  
  float ILCpaper=0.6;
  if(energy_under_test==91.2) ILCpaper=0.55;
  else if(energy_under_test==250) ILCpaper=0.65;
  else if(energy_under_test==380) ILCpaper=0.6;
  else if(energy_under_test==500) ILCpaper=0.65;
  else if(energy_under_test==1000) ILCpaper=0.6;
  else if(energy_under_test==9999) ILCpaper=0.5;
  else if(energy_under_test==10000) ILCpaper=0.35;
  
  TString conditions[]={TString::Format("H20, P:(-0.8,+0.3), b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b-quark, #sqrt{s}=%s",energyst.Data()),
			TString::Format("H20, P:(-0.8,+0.3), c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., c-quark, #sqrt{s}=%s",energyst.Data()),
			TString::Format("H20, P:(-0.8,+0.3), b&c quarks, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, P:(+0.8,-0.3), b&c quarks, #sqrt{s}=%s",energyst.Data()),TString::Format("H20, both pol., b&c quarks, #sqrt{s}=%s",energyst.Data())};
  TString conditions_clic[]={TString::Format("CLIC-like, P:-0.8, b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, P:+0.8, b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, both  pol., b-quark, #sqrt{s}=%s",energyst.Data()),
			     TString::Format("CLIC-like, P:-0.8, c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, P:+0.8, c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, both pol., c-quark, #sqrt{s}=%s",energyst.Data()),
			     TString::Format("CLIC-like, P:-0.8, b&c quarks, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, P:+0.8, b&c quarks, #sqrt{s}=%s",energyst.Data()),TString::Format("CLIC-like, both pol., b&c quarks. #sqrt{s}=%s",energyst.Data())};
  TString conditions_unpol[]={TString::Format("unpol., b-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., c-quark, #sqrt{s}=%s",energyst.Data()),TString::Format("unpol., b&c-quark, #sqrt{s}=%s",energyst.Data())};


  QQBARLabel(0.75,0.70,"",-1,0.1);
  if(style=="proceedings"){
    if(energy_under_test!=380){
      if(tpc_status=="noTPC")QQBARLabel2(0.3,0.225, "No PID",kRed,0.055);
      if(tpc_status=="default")QQBARLabel2(0.3,0.225, "dEdx for CM",kRed,0.055);
      if(tpc_status=="dEdx")QQBARLabel2(0.3,0.225, "dEdx for CM & FT",kRed,0.055);
      if(tpc_status=="dNdx")QQBARLabel2(0.3,0.225, "dNdx for CM & FT",kRed,0.055);
    }
    if (SM_status=="SMA")QQBARLabel2(0.37,0.38, "sin^{2}#theta_{W}=0.23126 (SM)",kOrange+7,0.045);
    else if (SM_status=="SMB")QQBARLabel2(0.37,0.38, "sin^{2}#theta_{W}=0.23122 (SM)",kOrange+7,0.045);
    if(pol==true){
      if(energy_under_test==380)QQBARLabel2(0.20,0.15, conditions_clic[iter],kBlack,0.050);
      else QQBARLabel2(0.20,0.15, conditions[iter],kBlack,0.055);
    }
    else QQBARLabel2(0.20,0.15, conditions_unpol[iter],kBlack,0.055);
  }
  else if(style=="paper"){
    if(tpc_status=="noTPC")QQBARLabel2(0.3,0.3, "No PID",kRed,0.055);
    if(tpc_status=="dEdx")QQBARLabel2(0.3,0.3, "dEdx for CM & FT",kRed,0.055);
    if(pol==false)QQBARLabel2(0.65,0.3,"Unpolarised",kRed,0.055);
    if (SM_status=="SMA")QQBARLabel2(0.37,0.38, "sin^{2}#theta_{W}=0.23126 (SM)",kOrange+7,0.045);
    else if (SM_status=="SMB")QQBARLabel2(0.37,0.38, "sin^{2}#theta_{W}=0.23122 (SM)",kOrange+7,0.045);
    QQBARLabel2(ILCpaper,0.15,energypaper,kBlack,0.10);
  }

}


void DrawLeg()
{
  TString colors_st[] = {"<1 #sigma","1-2 #sigma","< 3 #sigma","3-4 #sigma","4-5 #sigma","> 5 #sigma"}; // #colors >= #levels - 1                              
  Int_t colors[]={kGreen+3,kGreen+3,kGreen+3,kTeal+4,kTeal+5,kTeal+6};
  
  gStyle->SetPalette(6,colors);
  
  TLegend *leg2 = new TLegend(0.,0.275,1.3,0.625,"","blNDC");
  leg2->SetTextFont(42);
  leg2->SetTextSize(0.2);
  TH2F* h_leg02[10];
  for(int ic=2; ic<6;ic++ ){
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
  TBox *b = new TBox(0.07,0.9,0.95,0.9);
  b->SetFillColor(kBlack);
  b->SetLineColor(kBlack);
  b->SetLineWidth(0);
  b->Draw("same");
}

void DrawBracket(float xpos, float ypos)
{
  TPad *p = new TPad("p","p",0.,0.,1.,1.);
  p->SetFillStyle(0);
  p->Draw("same");
  p->cd();
  TBox *b = new TBox(xpos,ypos,xpos+0.215,ypos+0.01);
  b->SetFillColor(kBlack);
  b->SetLineColor(kBlack);
  b->SetLineWidth(2);
  b->Draw("same");
  TBox *bl = new TBox(xpos,ypos,xpos+0.003,ypos+0.2);
  bl->SetFillColor(kBlack);
  bl->SetLineColor(kBlack);
  bl->SetLineWidth(2);
  bl->Draw("same");
  TBox *br = new TBox(xpos+0.215-0.003,ypos,xpos+0.215,ypos+0.2);
  br->SetFillColor(kBlack);
  br->SetLineColor(kBlack);
  br->SetLineWidth(2);
  br->Draw("same");
}

void DrawSepLine(float xpos)
{
  TPad *p = new TPad("p","p",0.,0.,1.,1.);
  p->SetFillStyle(0);
  p->Draw("same");
  p->cd();
  TBox *b = new TBox(xpos-0.001,0.05,xpos+0.001,0.9);
  b->SetFillColor(kWhite);
  b->SetLineColor(kWhite);
  b->SetLineWidth(0);
  b->Draw("same");
}

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
    
    if(sigmavalue>9.94)QQBARLabel3(px-0.006,py,">10",color,0.045);
    else QQBARLabel3(px,py,TString::Format("%.1f",sigmavalue),color,0.045);            
  }    
}

TH1F * histo_one_energy_unpol(std::vector<TH2F*> probhistos_c,std::vector<TH2F*> probhistos_b,double energy_under_test=250, TString errortype="Stat", TString tpc_status="default")
{
  const Int_t xyNBINS = 13;
  Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5};

  TString energyst;
  if(energy_under_test==91.2) energyst="Z-Pole";
  else if(energy_under_test==250) energyst="250";
  else if(energy_under_test==500) energyst="500";
  else if(energy_under_test==1000) energyst="1000";

  TH1F * results_unpol = new TH1F("AFB_"+energyst+"_"+errortype+"_"+tpc_status+"_unpol","AFB_"+energyst+"_"+errortype+"_"+tpc_status+"_unpol",xyNBINS,xyedges);

  int binx;
  int biny;
  double nsigma=0;
  double nsigma1=0;
  double nsigma2=0;
  double prob1=1;
  double prob2=1;
  double prob=0;

  for(int j=1; j<theory.size(); j++) {
    biny=0;
    binx=j;
    nsigma1=probhistos_b.at(2)->GetBinContent(binx+1,biny+1);
    nsigma2=probhistos_c.at(2)->GetBinContent(binx+1,biny+1);
    prob1=Sigma_To_Prob(nsigma1);
    prob2=Sigma_To_Prob(nsigma2);
    prob=prob1*prob2;
    nsigma=Prob_To_Sigma(prob);

    // Rounding for plots:
    if((nsigma>0.95)and(nsigma<1))nsigma=0.9;
    if((nsigma>1.95)and(nsigma<2))nsigma=1.9;
    if((nsigma>2.95)and(nsigma<3))nsigma=2.9;
    if((nsigma>3.95)and(nsigma<4))nsigma=3.9;
    if((nsigma>4.95)and(nsigma<5))nsigma=4.9;
    if((nsigma>9.95)and(nsigma<10))nsigma=9.9;

    results_unpol->SetBinContent(j,nsigma);
  }

  //delete probhistos_b;
  //  probhistos_b.at(0);
  //probhistos_c.at(0);
  
  return results_unpol;

}



TH1F * histo_one_energy(std::vector<TH2F*> probhistos_c,std::vector<TH2F*> probhistos_b,double energy_under_test=250, TString errortype="Stat", TString tpc_status="default")

{

  TString energyst="250";
  const Int_t xyNBINS = 13;
  Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5};

  TH1F * results = new TH1F("AFB_"+energyst+"_"+errortype+"_"+tpc_status,"AFB_"+energyst+"_"+errortype+"_"+tpc_status,xyNBINS,xyedges);
    
  int binx;
  int biny;
  for(int j=1; j<theory.size(); j++) {
    biny=0;
    binx=j;
    double nsigma=0;
    
    double prob1=Sigma_To_Prob(probhistos_b.at(0)->GetBinContent(binx+1,biny+1));
    double prob2=Sigma_To_Prob(probhistos_b.at(1)->GetBinContent(binx+1,biny+1));
    double prob3=Sigma_To_Prob(probhistos_c.at(0)->GetBinContent(binx+1,biny+1));
    double prob4=Sigma_To_Prob(probhistos_c.at(1)->GetBinContent(binx+1,biny+1));
    double prob=prob1*prob2*prob3*prob4;
    nsigma=Prob_To_Sigma(prob);
    
    //Rounding
    if((nsigma>0.95)and(nsigma<1))nsigma=0.9;
    if((nsigma>1.95)and(nsigma<2))nsigma=1.9;
    if((nsigma>2.95)and(nsigma<3))nsigma=2.9;
    if((nsigma>3.95)and(nsigma<4))nsigma=3.9;
    if((nsigma>4.95)and(nsigma<5))nsigma=4.9;
    if((nsigma>9.95)and(nsigma<10))nsigma=9.9;  
    results->SetBinContent(j,nsigma);
  }
  
  return results;
    
}

TH1F * histo_both_energies(std::vector<TH2F*> probhistos_c,std::vector<TH2F*> probhistos_c500,std::vector<TH2F*> probhistos_b,std::vector<TH2F*> probhistos_b500,TString errortype="Stat", TString tpc_status="default")
{
    TString energyst="250&500_GeV";

    const Int_t xyNBINS = 13;
    Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5};

    TH1F * results = new TH1F("AFB_both_"+energyst+"_"+errortype+"_"+tpc_status,"AFB_both_"+energyst+"_"+errortype+"_"+tpc_status,xyNBINS,xyedges);
    
    int binx;
    int biny;
    for(int j=1; j<theory.size(); j++) {
      biny=0;
      binx=j;
      double nsigma=0;	  
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
      
      // Rounding for plots:                                                                                                 
      if((nsigma>0.95)and(nsigma<1))nsigma=0.9;
      if((nsigma>1.95)and(nsigma<2))nsigma=1.9;
      if((nsigma>2.95)and(nsigma<3))nsigma=2.9;
      if((nsigma>3.95)and(nsigma<4))nsigma=3.9;
      if((nsigma>4.95)and(nsigma<5))nsigma=4.9;
      if((nsigma>9.95)and(nsigma<10))nsigma=9.9;
      results->SetBinContent(j,nsigma);
    }

    return results;
}

TH1F * histo_three_energies(std::vector<TH2F*> probhistos_c,std::vector<TH2F*> probhistos_c500,std::vector<TH2F*> probhistos_c1000,std::vector<TH2F*> probhistos_b,std::vector<TH2F*> probhistos_b500,std::vector<TH2F*> probhistos_b1000,TString errortype="Stat", TString tpc_status="default")
{
  TString energyst="250&500&1000 GeV";

  const Int_t xyNBINS = 13;
  Double_t xyedges[xyNBINS + 1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5, 12.5};

  TH1F * results = new TH1F("AFB_three_"+energyst+"_"+errortype+"_"+tpc_status,"AFB_three_"+energyst+"_"+errortype+"_"+tpc_status,xyNBINS,xyedges);
  int binx;
  int biny;
  for(int j=1; j<theory.size(); j++) {
    biny=0;
    binx=j;
    double nsigma=0;
    
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
    
    // Rounding for plots:                       
    if((nsigma>0.95)and(nsigma<1))nsigma=0.9;
    if((nsigma>1.95)and(nsigma<2))nsigma=1.9;
    if((nsigma>2.95)and(nsigma<3))nsigma=2.9;
    if((nsigma>3.95)and(nsigma<4))nsigma=3.9;
    if((nsigma>4.95)and(nsigma<5))nsigma=4.9;
    if((nsigma>9.95)and(nsigma<10))nsigma=9.9;
    
    results->SetBinContent(j,nsigma);
  }
  return results;
}

void test_SM_PID()
{

  //Read all models:
  read_all_models(false);
  std::vector<observables_struct_t> observables=create_observables();
  std::vector<TH2F*> probhistos_c_250_nopid=nsigmas_models(4,250, observables,"StatTheo250","noTPC");
  std::vector<TH2F*> probhistos_b_250_nopid=nsigmas_models(5,250, observables,"StatTheo250","noTPC");
  std::vector<TH2F*> probhistos_c_250_dEdx=nsigmas_models(4,250, observables,"StatTheo250","dEdx");
  std::vector<TH2F*> probhistos_b_250_dEdx=nsigmas_models(5,250, observables,"StatTheo250","dEdx");
  std::vector<TH2F*> probhistos_c_250_dNdx=nsigmas_models(4,250, observables,"StatTheo250","dNdx");
  std::vector<TH2F*> probhistos_b_250_dNdx=nsigmas_models(5,250, observables,"StatTheo250","dNdx");
  std::vector<TH2F*> probhistos_c_500_nopid=nsigmas_models(4,500, observables,"StatTheo250","noTPC");
  std::vector<TH2F*> probhistos_b_500_nopid=nsigmas_models(5,500, observables,"StatTheo250","noTPC");
  std::vector<TH2F*> probhistos_c_500_dEdx=nsigmas_models(4,500, observables,"StatTheo250","dEdx");
  std::vector<TH2F*> probhistos_b_500_dEdx=nsigmas_models(5,500, observables,"StatTheo250","dEdx");
  std::vector<TH2F*> probhistos_c_500_dNdx=nsigmas_models(4,500, observables,"StatTheo250","dNdx");
  std::vector<TH2F*> probhistos_b_500_dNdx=nsigmas_models(5,500, observables,"StatTheo250","dNdx");
  std::vector<TH2F*> probhistos_c_1000_nopid=nsigmas_models(4,1000, observables,"StatTheo250","noTPC");
  std::vector<TH2F*> probhistos_b_1000_nopid=nsigmas_models(5,1000, observables,"StatTheo250","noTPC");
  std::vector<TH2F*> probhistos_c_1000_dEdx=nsigmas_models(4,1000, observables,"StatTheo250","dEdx");
  std::vector<TH2F*> probhistos_b_1000_dEdx=nsigmas_models(5,1000, observables,"StatTheo250","dEdx");
  std::vector<TH2F*> probhistos_c_1000_dNdx=nsigmas_models(4,1000, observables,"StatTheo250","dNdx");
  std::vector<TH2F*> probhistos_b_1000_dNdx=nsigmas_models(5,1000, observables,"StatTheo250","dNdx");
  
  TH1F * results_250_unpol[3];
  results_250_unpol[0]=histo_one_energy_unpol(probhistos_c_250_nopid,probhistos_b_250_nopid,250,"StatTheo250","noTPC");
  results_250_unpol[1]=histo_one_energy_unpol(probhistos_c_250_dEdx,probhistos_b_250_dEdx,250,"StatTheo250","dEdx");
  results_250_unpol[2]=histo_one_energy_unpol(probhistos_c_250_dNdx,probhistos_b_250_dNdx,250,"StatTheo250","dNdx");
  
  TH1F * results_250[3];
  results_250[0]=histo_one_energy(probhistos_c_250_nopid,probhistos_b_250_nopid,250,"StatTheo250","noTPC");
  results_250[1]=histo_one_energy(probhistos_c_250_dEdx,probhistos_b_250_dEdx,250,"StatTheo250","dEdx");
  results_250[2]=histo_one_energy(probhistos_c_250_dNdx,probhistos_b_250_dNdx,250,"StatTheo250","dNdx");
  
  TH1F * results_both[3];
  results_both[0]=histo_both_energies(probhistos_c_250_nopid,probhistos_c_500_nopid,probhistos_b_250_nopid,probhistos_b_500_nopid,"StatTheo250","noTPC");
  results_both[1]=histo_both_energies(probhistos_c_250_dEdx,probhistos_c_500_dEdx,probhistos_b_250_dEdx,probhistos_b_500_dEdx,"StatTheo250","dEdx");
  results_both[2]=histo_both_energies(probhistos_c_250_dNdx,probhistos_c_500_dNdx,probhistos_b_250_dNdx,probhistos_b_500_dNdx,"StatTheo250","dNdx");
  
  TH1F * results_three[3];
  results_three[0]=histo_three_energies(probhistos_c_250_nopid,probhistos_c_500_nopid,probhistos_c_1000_nopid,probhistos_b_250_nopid,probhistos_b_500_nopid,probhistos_b_1000_nopid,"StatTheo250","noTPC");
  results_three[1]=histo_three_energies(probhistos_c_250_dEdx,probhistos_c_500_dEdx,probhistos_c_1000_dEdx,probhistos_b_250_dEdx,probhistos_b_500_dEdx,probhistos_b_1000_dEdx,"StatTheo250","dEdx");
  results_three[2]=histo_three_energies(probhistos_c_250_dNdx,probhistos_c_500_dNdx,probhistos_c_1000_dNdx,probhistos_b_250_dNdx,probhistos_b_500_dNdx,probhistos_b_1000_dNdx,"StatTheo250","dNdx");
  
  cout<<"sigma levels 250 unpol: "<<endl;
  for(int i=1;i<14;i++){
    cout<<results_250_unpol[0]->GetBinContent(i)<<" , "<<results_250_unpol[1]->GetBinContent(i)<<" , "<<results_250_unpol[2]->GetBinContent(i)<<endl;
  }
  cout<<"sigma levels 250: "<<endl;
  for(int i=1;i<14;i++){
    cout<<results_250[0]->GetBinContent(i)<<" , "<<results_250[1]->GetBinContent(i)<<" , "<<results_250[2]->GetBinContent(i)<<endl;
  }
  cout<<"sigma levels both: "<<endl;
  for(int i=1;i<14;i++){
    cout<<results_both[0]->GetBinContent(i)<<" , "<<results_both[1]->GetBinContent(i)<<" , "<<results_both[2]->GetBinContent(i)<<endl;
  }
  cout<<"sigma levels three: "<<endl;
  for(int i=1;i<14;i++){
    cout<<results_three[0]->GetBinContent(i)<<" , "<<results_three[1]->GetBinContent(i)<<" , "<<results_three[2]->GetBinContent(i)<<endl;
  }

  //Preparing histos

  TString title_models[]={"SMB","A_{1}","A_{2}","B_{1}^{-}","B_{1}^{+}","B_{2}^{-}","B_{2}^{+}","B_{3}^{-}","B_{3}^{+}"};
  TString title_precision[]={"O","E","N","O","E","N","O","E","N","O","E","N"};
  
  Int_t colors[6]={kGreen+3,kGreen+3,kGreen+3,kTeal+4,kTeal+5,kTeal+6};
  gStyle->SetPalette((sizeof(colors)/sizeof(Int_t)), colors);
  Double_t levels[7]={0,0.99,1.99,2.99,3.99,4.99,1000.0};
  gStyle->SetPaintTextFormat(".1f");
  
  //Prepare histo
  const Int_t xbins=12;
  const Int_t ybins=8;
  Double_t xedges[xbins+1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5, 11.5};
  Double_t yedges[ybins+1] = {-0.5, 0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5};

  TH2F * SM_comparison=new TH2F("SM_comparison","SM_comparison",xbins,xedges,ybins,yedges);
  for(int j=0; j<theory.size()-1; j++) {
    SM_comparison->GetYaxis()->SetBinLabel(j+1,title_models[j+1]);
  }
  for(int j=0; j<12; j++) {
    SM_comparison->GetXaxis()->SetBinLabel(j+1,title_precision[j]);
  }
  
  SM_comparison->SetTitle("");
  SM_comparison->GetXaxis()->SetTickLength(0.);
  SM_comparison->GetYaxis()->SetTickLength(0.);
  SM_comparison->GetYaxis()->SetLabelSize(0.07);
  SM_comparison->GetXaxis()->SetLabelSize(0.06);

  std::vector<double> xWhite;
  std::vector<double> yWhite;
  std::vector<double> binWhite;
  std::vector<double> xBlack;
  std::vector<double> yBlack;
  std::vector<double> binBlack;
  
  for(int ix=1;ix<13;ix++){
    for(int iy=1;iy<9;iy++){
      float nsigma;
      if(ix<4){
	nsigma=results_250_unpol[ix-1]->GetBinContent(iy);
      }
      else if(ix<7){
	nsigma=results_250[ix-4]->GetBinContent(iy);
      } 
      else if(ix<10){
	nsigma=results_both[ix-7]->GetBinContent(iy);
      }
      else {
	nsigma=results_three[ix-10]->GetBinContent(iy);
      }
      cout<<"sigma: "<<nsigma<<endl;      
      if((ix%3==0)and(iy==13))cout<<"_____________________"<<endl;
      SM_comparison->SetBinContent(ix,iy,nsigma);

      if(nsigma<3.99){
        xWhite.push_back(0.078+0.0738*(ix-1));
        yWhite.push_back(0.09+0.106*(iy-1));
        binWhite.push_back(nsigma);
      }
      else{
        xBlack.push_back(0.078+0.0738*(ix-1));
        yBlack.push_back(0.09+0.106*(iy-1));
        binBlack.push_back(nsigma);
      }
    }
  }
  
  TCanvas *c_SM_comparison= new TCanvas ("SM_comparison","SM_comparison",800,800);
  c_SM_comparison->cd();
  gStyle->SetOptStat(0);
  gStyle->SetMarkerSize(2.);

  c_SM_comparison->cd();
  TPad *padB = new TPad("padB", "padB", 0., 0., 0.85, 0.3);
  padB->SetTopMargin(0.065);
  padB->SetBottomMargin(0);
  padB->Draw();
  padB->cd();
  DrawBracket(0.070, 0.7);
  DrawBracket(0.2905, 0.7);
  DrawBracket(0.511, 0.7);
  DrawBracket(0.7315, 0.7);
  //DrawBracket(0.60, 0.45);
  //DrawBracket(0.80, 0.45);
  QQBARLabel2(0.082,0.50, "ILC250^{#diamond}",kBlack,0.17);
  QQBARLabel2(0.082,0.30, "(no pol.)",kBlack,0.17);
  QQBARLabel2(0.3025,0.50, "ILC250",kBlack,0.17);
  QQBARLabel2(0.5230,0.50, "ILC250",kBlack,0.17);
  QQBARLabel2(0.5230,0.30, " +500",kBlack,0.17);
  QQBARLabel2(0.7435,0.50, "ILC250",kBlack,0.17);
  QQBARLabel2(0.7435,0.30, " +500",kBlack,0.16);
  QQBARLabel2(0.7435,0.10, " +1000*",kBlack,0.16);
  
  c_SM_comparison->cd();
  TPad *padL = new TPad("padL", "padL", 0., 0.25, 0.85, 0.95);
  padL->SetTopMargin(0.1);
  padL->SetBottomMargin(0.05);
  padL->SetLeftMargin(0.065);
  padL->SetRightMargin(0.05);
  padL->Draw();
  padL->cd();
  SM_comparison->GetZaxis()->SetRangeUser(0.0,7.0);
  SM_comparison->SetContour((sizeof(levels)/sizeof(Double_t)), levels);
  SM_comparison->Draw("Col");
  DrawTextBins(xBlack,yBlack,binBlack,kBlack);
  DrawTextBins(xWhite,yWhite,binWhite,kWhite);
  DrawSepLine(0.286);
  DrawSepLine(0.507);
  DrawSepLine(0.728);
  DrawTopLine();
  QQBARLabel3(0.095,0.925,"GHU vs SM discrimination power (#sigma-level)",kBlack,0.06);

  c_SM_comparison->cd();
  TPad *padR = new TPad("padR", "padR", 0.82, 0.15, 1., 0.85);
  padR->SetTopMargin(0.1);
  padR->SetBottomMargin(0.1);
  padR->SetLeftMargin(0.001);
  padR->SetRightMargin(0.001);
  padR->Draw();
  padR->cd();
  DrawLeg();
  QQBARLabel3(0.05,0.95, "Ch. had. PID",kBlack,0.17);
  QQBARLabel2(0.1,0.88, "#bullet O: No PID",kBlack,0.16);
  QQBARLabel2(0.1,0.78, "#bullet E: #frac{dE}{dx}",kBlack,0.16);
  QQBARLabel2(0.1,0.68, "#bullet N: #frac{dN}{dx}",kBlack,0.16);
  //QQBARLabel(-0.015,0.8775,"",-1,0.1);
  c_SM_comparison->cd();
  QQBARLabel(0.845,0.9,"",-1,0.1);

  // Save everything:    
  TString savingformat[2]={".png",".eps"};
  for(int isave=0; isave<2 ;isave++){
      c_SM_comparison->SaveAs("SM_comparison_PID"+savingformat[isave]);
  }
}
