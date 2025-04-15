// # Copyright 2023  Adrián Irles (IFIC)

#include "../include/models_struct.h"
#include "TMath.h"

void read_all_models(bool debug = true)
{

    //uses std::vector<TString> allmodels = {"SMA","SMB","A1","A2","A3","B1","B2","B3","B4"}; defined in models_struct
    //std::vector<model_struct_t> theory;

    for (int i = 0; i < allmodels.size(); i++)
    {
        model_struct_t newmodel = read_model(allmodels.at(i), i, debug);
        theory.push_back(newmodel);
    }

    if (debug == true)
    {

        for (int i = 0; i < theory.size(); i++)
        {
            cout << theory.at(i).modelname << " index:" << theory.at(i).id << endl;
            for (int ie = 0; ie < 6; ie++)
                for (int ifl = 1; ifl < 6; ifl++)
                    cout << theory.at(i).model_values[ie].energy << " " << theory.at(i).model_values[ie].cross_L[ifl] << " " << theory.at(i).model_values[ie].cross_R[ifl]
                         << " " << theory.at(i).model_values[ie].AFB_L[ifl] << " " << theory.at(i).model_values[ie].AFB_R[ifl] << endl;
        }
    }
}

 std::vector<observables_struct_t> create_observables() {

    std::vector<observables_struct_t> observables;
    for (int imodel = 0; imodel < theory.size(); imodel++)
    {
      for (int ienergy = 0; ienergy < 6; ienergy++) //Z-Pole,250,350,380,500,1000
        {
	  //if (theory.at(imodel).model_values[ienergy].energy != 250)
	  //    continue;
	  
	  for (int iflav = 4; iflav < 6; iflav++)
            {
	      if(ienergy==3){ //3 is for CLIC 380
		observables_struct_t newobs;
		newobs.id = theory.at(imodel).id;
		newobs.modelname = theory.at(imodel).modelname;
		newobs.energy = theory.at(imodel).model_values[ienergy].energy;
		newobs.flav = iflav;
		newobs.Cross_L = ObsCross(imodel, ienergy, iflav, -0.8,0);
		newobs.Cross_R = ObsCross(imodel, ienergy, iflav, 0.8,0);
		newobs.Cross_unpol = ObsCross(imodel, ienergy, iflav, 0,0);
		newobs.R_L = ObsR(imodel, ienergy, iflav, -0.8,0);
		newobs.AFB_L = ObsAFB(imodel, ienergy, iflav, -0.8,0);
		newobs.R_R = ObsR(imodel, ienergy, iflav, 0.8,0);
		newobs.AFB_R = ObsAFB(imodel, ienergy, iflav,0.8,0);
		newobs.R_unpol = ObsR(imodel, ienergy, iflav, 0,0);
		newobs.AFB_unpol = ObsAFB(imodel, ienergy, iflav, 0,0);
		observables.push_back(newobs);
	      }
	      else{
		observables_struct_t newobs;
		newobs.id = theory.at(imodel).id;
		newobs.modelname = theory.at(imodel).modelname;
		newobs.energy = theory.at(imodel).model_values[ienergy].energy;
		newobs.flav = iflav;
		newobs.Cross_L = ObsCross(imodel, ienergy, iflav, -0.8,0.3);
		newobs.Cross_R = ObsCross(imodel, ienergy, iflav, 0.8,-0.3);
		newobs.Cross_unpol = ObsCross(imodel, ienergy, iflav, 0,0);
		newobs.R_L = ObsR(imodel, ienergy, iflav, -0.8,0.3);
		newobs.AFB_L = ObsAFB(imodel, ienergy, iflav, -0.8,0.3);
		newobs.R_R = ObsR(imodel, ienergy, iflav, 0.8,-0.3);
		newobs.AFB_R = ObsAFB(imodel, ienergy, iflav,0.8,-0.3);
		newobs.R_unpol = ObsR(imodel, ienergy, iflav, 0,0);
		newobs.AFB_unpol = ObsAFB(imodel, ienergy, iflav, 0,0);
		observables.push_back(newobs);
	      }
	    }
	}
    }
    return observables;

 }

std::vector<TH2F*> nsigmas_models(int quark_under_test, float energy_under_test, std::vector<observables_struct_t> observables, TString errortype, TString tpc_status) {

  
    TH2F * h_nsigmas[8];

    for(int i=0; i<8; i++) {
      h_nsigmas[i]=new TH2F(TString::Format("h_nsigmas_%i_q%i_en%i_",i,quark_under_test,int(energy_under_test))+errortype+"_"+tpc_status,TString::Format("h_nsigmas_%i_q%i_en%i_",i,quark_under_test,int(energy_under_test))+errortype+"_"+tpc_status,31,-0.5,30.5,31,-0.5,30.5);
    }

    float cross_SM_L;
    float cross_SM_R;
    float cross_SM_unpol;

    for(int iobs=0; iobs<observables.size(); iobs++) {//loop for reference model
        //this could be done with iterators, but this is good enough...

        if(observables.at(iobs).flav!=quark_under_test) continue;
        if(observables.at(iobs).energy!=energy_under_test) continue;

        int id1= observables.at(iobs).id;

        if(id1==0) { //Remember SMB
            cross_SM_L=observables.at(iobs).Cross_L;
            cross_SM_R=observables.at(iobs).Cross_R;
            cross_SM_unpol=observables.at(iobs).Cross_unpol;
        }
   
        for(int jobs=0; jobs<observables.size(); jobs++) {

            int id2=observables.at(jobs).id;
            //if(id2>id1) continue;

            float sigmas_AFB_L=0;
            float sigmas_AFB_R=0;
            float sigmas_AFB_unpol=0;

            float sigmas_R_L=0;
            float sigmas_R_R=0;
            float sigmas_R_unpol=0;

            if(observables.at(jobs).flav!=quark_under_test) continue;
            if(observables.at(jobs).energy!=energy_under_test ) continue;

            int iquark=observables.at(iobs).flav;
            if(iquark==5) iquark=6;//index to read properly the uncertainties
            int index_energy=0;
	    if(observables.at(iobs).energy==91.2) index_energy=0;
	    if(observables.at(iobs).energy==250) index_energy=1;
            if(observables.at(iobs).energy==350) index_energy=2;
	    if(observables.at(iobs).energy==380) index_energy=3;
	    if(observables.at(iobs).energy==500) index_energy=4;
	    if(observables.at(iobs).energy==1000) index_energy=5;
	    
            //----------
            // polarized Beams
            float r_stat_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*r_stat[index_energy][(iquark-4)]/100.; //renormalize of the stat. uncertainties to the cross section of the model
            float r_stat_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*r_stat[index_energy][(iquark-4)+1]/100.;
	    float afb_stat_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*afb_stat[index_energy][(iquark-4)]/100.;
            float afb_stat_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*afb_stat[index_energy][(iquark-4)+1]/100.;
	    float afb_theo_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*afb_theo_current[index_energy][(iquark-4)]/100.;
            float afb_theo_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*afb_theo_current[index_energy][(iquark-4)+1]/100.;

	    if(tpc_status=="noTPC"){
	      afb_stat_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*afb_stat_noTPC[index_energy][(iquark-4)]/100.;
	      afb_stat_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*afb_stat_noTPC[index_energy][(iquark-4+1)]/100.;
	    }
	    else if(tpc_status=="dEdx"){
              afb_stat_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*afb_stat_dEdx[index_energy][(iquark-4)]/100.;
              afb_stat_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*afb_stat_dEdx[index_energy][(iquark-4+1)]/100.;
            }
	    else if(tpc_status=="dNdx"){
	      //Temporalely changing the scaling
              afb_stat_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*afb_stat_dNdx[index_energy][(iquark-4)]/100.;
              afb_stat_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*afb_stat_dNdx[index_energy][(iquark-4+1)]/100.;
            }
	    
	    if(errortype=="StatTheoGigaZ"){
	      afb_theo_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*afb_theo_gigaZ[index_energy][(iquark-4)]/100.;
	      afb_theo_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*afb_theo_gigaZ[index_energy][(iquark-4)+1]/100.;
	    }
	    else if(errortype=="StatTheo250"){
              afb_theo_L=sqrt(cross_SM_L/observables.at(jobs).Cross_L)*afb_theo_250[index_energy][(iquark-4)]/100.;
              afb_theo_R=sqrt(cross_SM_R/observables.at(jobs).Cross_R)*afb_theo_250[index_energy][(iquark-4)+1]/100.;
            }

	    float r_syst_L=r_syst[index_energy][(iquark-4)]/100.;
            float r_syst_R=r_syst[index_energy][(iquark-4)+1]/100.;
            float afb_syst_L=afb_syst[index_energy][(iquark-4)]/100.;
            float afb_syst_R=afb_syst[index_energy][(iquark-4)+1]/100.;

            //----------
            //we assume that the unpol syst unc is the same than lef handed
            //what about the unpol stat.. uncertainty? We should assume a larger luminosity, for the moment (13thApril) 
            //let's assume a factor 4000/900 in luminosity w.r.t the L scenarios: 900fb-1 for (-0.8,0.3) 250GeV
            
            // float lum_factor=4000./900.;
            // for 500GeV, the unpolarised make no sense since this energy is not reachable by FCC
	    // note that the second index is for 350 GeV in reality
	    int iquark_u=observables.at(iobs).flav;
            float r_stat_unpolarised=r_stat_unpol[index_energy][(iquark_u-4)]/100.;
            float afb_stat_unpolarised=afb_stat_unpol[index_energy][(iquark_u-4)]/100.;
	    float r_syst_unpolarised=r_syst_unpol[index_energy][(iquark_u-4)]/100.;
	    float afb_syst_unpolarised=afb_syst_unpol[index_energy][(iquark_u-4)]/100.;
	    float afb_theo_unpolarised=afb_theo_unpol_current[index_energy][(iquark_u-4)]/100.;

	    if(tpc_status=="noTPC"){
              afb_stat_unpolarised=afb_stat_unpol_noTPC[index_energy][(iquark_u-4)]/100.;
	    }
	    else if(tpc_status=="dEdx"){
              afb_stat_unpolarised=afb_stat_unpol_dEdx[index_energy][(iquark_u-4)]/100.;
            }
	    else if(tpc_status=="dNdx"){
              afb_stat_unpolarised=afb_stat_unpol_dNdx[index_energy][(iquark_u-4)]/100.;
            }
	    
	    if(errortype=="StatTheoGigaZ"){
	      afb_theo_unpolarised=afb_theo_unpol_gigaZ[index_energy][(iquark_u-4)]/100.;
	    } 
	    else if(errortype=="StatTheo250"){
              afb_theo_unpolarised=afb_theo_unpol_250[index_energy][(iquark_u-4)]/100.;
            }
	    
	    if(errortype=="Stat"){
	      r_syst_L=0;
	      r_syst_R=0;
	      afb_syst_L=0;
	      afb_syst_R=0;
	      r_syst_unpolarised=0;
	      afb_syst_unpolarised=0;
              afb_theo_L=0;
              afb_theo_R=0;
	      afb_theo_unpolarised=0;
	    }
	    else if(errortype=="StatSyst"){
	      afb_theo_L=0;
              afb_theo_R=0;
              afb_theo_unpolarised=0;
	    }
	    else{
	      r_syst_L=0;
              r_syst_R=0;
              afb_syst_L=0;
              afb_syst_R=0;
              r_syst_unpolarised=0;
              afb_syst_unpolarised=0;
	    }

	    // ECFA ParT prospects
	    if(iquark==4){
              r_stat_L=(sqrt(0.095)/0.5)*r_stat_L;
              r_stat_R=(sqrt(0.095)/0.5)*r_stat_R;
              afb_stat_L=(sqrt(0.095)/0.5)*afb_stat_L;
              afb_stat_R=(sqrt(0.095)/0.5)*afb_stat_R;
	      afb_stat_unpolarised=(sqrt(0.095)/0.5)*afb_stat_unpolarised;
	    }
	    // Note that now iquark==6 instead of 5 to read the b values!
            if(iquark==6){
              r_stat_L=(sqrt(0.381)/0.8)*r_stat_L;
              r_stat_R=(sqrt(0.381)/0.8)*r_stat_R;
              afb_stat_L=(sqrt(0.381)/0.8)*afb_stat_L;
              afb_stat_R=(sqrt(0.381)/0.8)*afb_stat_R;
	      afb_stat_unpolarised=(sqrt(0.381)/0.8)*afb_stat_unpolarised;
	    }
      
	    //Total errors:
	    float r_total_L=sqrt(r_stat_L*r_stat_L+r_syst_L*r_syst_L);
            float r_total_R=sqrt(r_stat_R*r_stat_R+r_syst_R*r_syst_R);
            float afb_total_L=sqrt(afb_stat_L*afb_stat_L+afb_syst_L*afb_syst_L+afb_theo_L*afb_theo_L);
            float afb_total_R=sqrt(afb_stat_R*afb_stat_R+afb_syst_R*afb_syst_R+afb_theo_R*afb_theo_R);
            float r_total_unpol=sqrt(r_stat_unpolarised*r_stat_unpolarised+r_syst_unpolarised*r_syst_unpolarised);
            float afb_total_unpol=sqrt(afb_stat_unpolarised*afb_stat_unpolarised+afb_syst_unpolarised*afb_syst_unpolarised+afb_theo_unpolarised*afb_theo_unpolarised);
            	    
	    // nsigmas formula, using uncertainty of the tested model           
	    //P(P, Q) = (1/2πSPSQ) * exp[-((P-P')^2/2SP^2) - ((Q-Q')^2/2SQ^2)];
	    sigmas_AFB_L=abs(observables.at(jobs).AFB_L-observables.at(iobs).AFB_L) / abs(observables.at(iobs).AFB_L*afb_total_L);
	    sigmas_AFB_R=abs(observables.at(jobs).AFB_R-observables.at(iobs).AFB_R) / abs(observables.at(iobs).AFB_R*afb_total_R);
	    sigmas_AFB_unpol=abs(observables.at(jobs).AFB_unpol-observables.at(iobs).AFB_unpol) / abs(observables.at(iobs).AFB_unpol*afb_total_unpol);
	    sigmas_R_L=abs(observables.at(jobs).R_L-observables.at(iobs).R_L) / (observables.at(iobs).R_L*r_total_L);
	    sigmas_R_R=abs(observables.at(jobs).R_R-observables.at(iobs).R_R) / (observables.at(iobs).R_R*r_total_R);
	    sigmas_R_unpol=abs(observables.at(jobs).R_unpol-observables.at(iobs).R_unpol) / (observables.at(iobs).R_unpol*r_total_unpol);

            h_nsigmas[0]->SetBinContent(id1+1,id2+1,sigmas_AFB_L);
            h_nsigmas[1]->SetBinContent(id1+1,id2+1,sigmas_AFB_R);
	    h_nsigmas[2]->SetBinContent(id1+1,id2+1,sigmas_AFB_unpol);

            h_nsigmas[3]->SetBinContent(id1+1,id2+1,sigmas_R_L);
            h_nsigmas[4]->SetBinContent(id1+1,id2+1,sigmas_R_R);
	    h_nsigmas[5]->SetBinContent(id1+1,id2+1,sigmas_R_unpol);

        }
    }

    for(int iobs=0; iobs<observables.size(); iobs++) {//loop for reference model
        for(int jobs=0; jobs<observables.size(); jobs++) {
            for(int i=0; i<8; i++) if(h_nsigmas[i]->GetBinContent(jobs+1,iobs+1) > h_nsigmas[i]->GetBinContent(iobs+1,jobs+1))  h_nsigmas[i]->SetBinContent(iobs+1,jobs+1, h_nsigmas[i]->GetBinContent(jobs+1,iobs+1));
        }
    }

    std::vector<TH2F*> result;
    for(int i=0; i<6; i++) {
      // I'm including the unpol case, will use it for comparisons!
        result.push_back(h_nsigmas[i]);
    }

    return result;

 }
