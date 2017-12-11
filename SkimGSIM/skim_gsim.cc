#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include "TFile.h"
#include "TNtuple.h"
#include "TTree.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TMath.h"
#include "DecayBOS.h"
#include "TLorentzVector.h"
#include "TVector3.h"
#include "eloss.h"
#include "DecayTarget.h"
#include "DecayTrackCorrections.h"
#include "DecayFiducialCuts.h"


using namespace std;


// Calculates the point on each track that is the intersection between the DOCA line and each track, as well as the doca
double Calc_dtfInterDOCA(const TVector3 &locUnitDir1, const TVector3 &locUnitDir2, const TVector3 &locVertex1, const TVector3 &locVertex2, TVector3 &locInterDOCA1, TVector3 &locInterDOCA2);


void Display_Help();
/*----------------------------------- Program Execution -----------------------------------*/

int main(int argc, char *argv[]){	//argc is the number of arguments; argv is the argument
	cout << "Author: Nicholas Zachariou (nicholas@jlab.org)" << endl;
	cout << "Run with '-H' or '-h' switch for help." << endl;
	int loc_i, locMaxNumEventsPerFile = -1, locNumEvents, loc_numbbos, numbbos=0, loc_times=0; 
	vector<string> locStringInputs; //vectors (can change their size and allocate by push_back command)
	string locTempString, locTempName, locName="tree_", loc_runnum, loc_filenum, loc_tripfile="clas_";
	istringstream locIStream;
	for(loc_i = 1; loc_i < argc; loc_i++){
		if(argv[loc_i][0] != '-'){
			locStringInputs.push_back(argv[loc_i]);
			if (loc_times == 0){
				locTempName=argv[loc_i];
				locTempName=locTempName.substr(locTempName.length()-14,10);//get the run number and file extension
				locName.append(locTempName);
				locName.append(".root");
				loc_tripfile.append(locTempName);
				loc_tripfile.append(".trip");
				loc_runnum=locTempName.substr(1, 6);
				loc_filenum=locTempName.substr(8, 2);
				loc_times++;
			}
			numbbos++;
		}
		else{	//it's a switch:						
			switch(argv[loc_i][1]){	
				case 'h'://display help
					Display_Help();
					return 0; 
				case 'H'://display help
					Display_Help();
					return 0; 
				case 'M'://set MaxNumEventsPerFile
					locTempString = argv[loc_i];
					locTempString = locTempString.substr(2, locTempString.length() - 2);	//strip "-M"
					locIStream.str(locTempString);	//stores locTempString to locIStreat
					if(!(locIStream >> locMaxNumEventsPerFile))	//converts locIStream into integer and assigns it to LocMaxNumEventsPerFile
						cout << "ERROR: COMMAND LINE INPUT NOT RECOGNIZED.  MaxNumEventsPerFile NOT SET." << endl;
					break;
				default:
					break;
			}//kills switch
		}
	}//kills loop over arguments											
    if(locStringInputs.size() == 0){	//Checks if there is an input bos file
		cout << "ERROR: NEED INPUT BOS FILE." << endl;
		return 2;
	}
	

	   
    
    
    TFile *locFile = new TFile(locName.c_str(), "RECREATE");	//create a new file, change recreate to update it each time
    TTree *locTree = new TTree("h22", "Hyperon Analysis");	//name, title

    
    int loc_numofpart, loc_num_pos, loc_num_neg, loc_num_neu, loc_head_eventnum, loc_head_runnum, loc_head_time, loc_head_filenum, loc_trip_flag;
    int loc_evnt_id[3], loc_dcpb_sect[3], loc_stpb_id[3], loc_scpb_id[3], loc_evnt_charge[3];
    int loc_tagr_eid[50], loc_tagr_tid[50], loc_tagr_stat[50], loc_num_good_ph, loc_num_ph_sc, loc_num_ph_st, loc_num_ph_sc2, loc_num_ph_st2;
    int loc_index_sc[50], loc_index_st[50], loc_index_sc2[50], loc_index_st2[50];
    
    float loc_evnt_pm[3], loc_evnt_cx[3], loc_evnt_cy[3], loc_evnt_cz[3], loc_evnt_x[3], loc_evnt_y[3], loc_evnt_z[3], loc_evnt_mass[3], loc_calc_mass[3];
    float loc_evnt_be[3], loc_calc_be[3], loc_delt_be[3],  loc_scpb_t[3], loc_stpb_t[3], loc_scpb_d[3], loc_stpb_d[3], loc_tagr_tpho[50], loc_tagr_ttag[50];
    float loc_phot_e[50], loc_delt_t_sc[50], loc_delt_t_st[50], loc_miss_mass[50], loc_miss_pm[50], loc_miss_pm_corr[50], loc_miss_mass_corr[50];
    float loc_delt_t_sc2[50], loc_delt_t_st2[50], loc_miss_pm_corr2[50], loc_miss_mass_corr2[50];
    float loc_inv_mass_01, loc_inv_mass_02, loc_inv_mass_12, loc_mvrt_chi2, loc_mvrt_x, loc_mvrt_y, loc_mvrt_z, loc_beam_x, loc_beam_y, loc_beam_z;
    float loc_inv_mass_01_corr, loc_inv_mass_02_corr, loc_inv_mass_12_corr, loc_inv_mass_01_corr2, loc_inv_mass_02_corr2, loc_inv_mass_12_corr2;
    float loc_coh_edge, loc_coh_plan, loc_coh_radi;
    float loc_doca_hyperon, loc_doca_beam, loc_poca_hyperon_x, loc_poca_beam_x,loc_poca_hyperon_y, loc_poca_beam_y, loc_poca_hyperon_z, loc_poca_beam_z;
    float loc_eloss_pm[3], loc_corr_pm[3], loc_corr_phot_e[50],loc_eloss_pm2[3], loc_corr_pm2[3];
    
    bool loc_fiducial[3];

   int loc_label;
    TLorentzVector *loc_beam=new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *loc_target=new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *loc_Qp1=new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *loc_Qp2=new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *loc_Starget=new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *loc_Sbeam=new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *loc_Rp1=new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *loc_Rp2=new TLorentzVector(0.,0.,0.,0.);
    TLorentzVector *loc_Rp3=new TLorentzVector(0.,0.,0.,0.);

    
    locTree->Branch("beam","TLorentzVector",&loc_beam);
    locTree->Branch("target","TLorentzVector",&loc_target);
    locTree->Branch("qp1","TLorentzVector",&loc_Qp1);
    locTree->Branch("qp2","TLorentzVector",&loc_Qp2);
    locTree->Branch("sbeam","TLorentzVector",&loc_Sbeam);
    locTree->Branch("starget","TLorentzVector",&loc_Starget);
    locTree->Branch("rp1","TLorentzVector",&loc_Rp1);
    locTree->Branch("rp2","TLorentzVector",&loc_Rp2);
    locTree->Branch("rp3","TLorentzVector",&loc_Rp3);
    locTree->Branch("label", &loc_label, "label/I");

    
    locTree->Branch("numofpart", &loc_numofpart, "numofpart/I");
    locTree->Branch("num_pos", &loc_num_pos, "num_pos/I");
    locTree->Branch("num_neg", &loc_num_neg, "num_neg/I");
    locTree->Branch("num_neu", &loc_num_neu, "num_neu/I");
    locTree->Branch("num_good_ph", &loc_num_good_ph, "num_good_ph/I");
    locTree->Branch("num_ph_sc", &loc_num_ph_sc, "num_ph_sc/I");
    locTree->Branch("num_ph_st", &loc_num_ph_st, "num_ph_st/I");

    locTree->Branch("num_ph_sc2", &loc_num_ph_sc2, "num_ph_sc2/I");
    locTree->Branch("num_ph_st2", &loc_num_ph_st2, "num_ph_st2/I");
    
    locTree->Branch("head_eventnum", &loc_head_eventnum, "head_eventnum/I");
    locTree->Branch("head_runnum", &loc_head_runnum, "head_runnum/I");
    locTree->Branch("head_time", &loc_head_time, "head_time/I");
    locTree->Branch("evnt_id", loc_evnt_id, "evnt_id[3]/I");
    locTree->Branch("dcpb_sect", loc_dcpb_sect, "dcpb_sect[3]/I");
    locTree->Branch("stpb_id", loc_stpb_id, "stpb_id[3]/I");
    locTree->Branch("scpb_id", loc_scpb_id, "scpb_id[3]/I");
    locTree->Branch("evnt_charge", loc_evnt_charge, "evnt_charge[3]/I");
    locTree->Branch("tagr_tid", loc_tagr_tid, "tagr_tid[num_good_ph]/I");
    locTree->Branch("tagr_eid", loc_tagr_eid, "tagr_eid[num_good_ph]/I");
    locTree->Branch("tagr_stat", loc_tagr_stat, "tagr_stat[num_good_ph]/I");
    
    locTree->Branch("evnt_pm", loc_evnt_pm, "evnt_pm[3]/F");
    locTree->Branch("evnt_cx", loc_evnt_cx, "evnt_cx[3]/F");
    locTree->Branch("evnt_cy", loc_evnt_cy, "evnt_cy[3]/F");
    locTree->Branch("evnt_cz", loc_evnt_cz, "evnt_cz[3]/F");
    locTree->Branch("evnt_x", loc_evnt_x, "evnt_x[3]/F");
    locTree->Branch("evnt_y", loc_evnt_y, "evnt_y[3]/F");
    locTree->Branch("evnt_z", loc_evnt_z, "evnt_z[3]/F");
    locTree->Branch("evnt_mass", loc_evnt_mass, "evnt_mass[3]/F");
    locTree->Branch("calc_mass", loc_calc_mass, "calc_mass[3]/F");
    locTree->Branch("evnt_be", loc_evnt_be, "evnt_be[3]/F");
    locTree->Branch("calc_be", loc_calc_be, "calc_be[3]/F");
    locTree->Branch("delt_be", loc_delt_be, "delt_be[3]/F");
    locTree->Branch("scpb_t", loc_scpb_t, "scpb_t[3]/F");
    locTree->Branch("stpb_t", loc_stpb_t, "stpb_t[3]/F");
    locTree->Branch("scpb_d", loc_scpb_d, "scpb_d[3]/F");
    locTree->Branch("stpb_d", loc_stpb_d, "stpb_d[3]/F");
    locTree->Branch("tagr_tpho", loc_tagr_tpho, "tagr_tpho[num_good_ph]/F");
    locTree->Branch("tagr_ttag", loc_tagr_ttag, "tagr_ttag[num_good_ph]/F");
    locTree->Branch("phot_e", loc_phot_e, "phot_e[num_good_ph]/F");
    locTree->Branch("delt_t_sc", loc_delt_t_sc, "delt_t_sc[num_good_ph]/F");
    locTree->Branch("delt_t_st", loc_delt_t_st, "delt_t_st[num_good_ph]/F");
    locTree->Branch("delt_t_sc2", loc_delt_t_sc2, "delt_t_sc2[num_good_ph]/F");
    locTree->Branch("delt_t_st2", loc_delt_t_st2, "delt_t_st2[num_good_ph]/F");
    
    locTree->Branch("miss_mass", loc_miss_mass, "miss_mass[num_good_ph]/F");
    locTree->Branch("miss_pm", loc_miss_pm, "miss_pm[num_good_ph]/F");
    locTree->Branch("miss_pm_corr", loc_miss_pm_corr, "miss_pm_corr[num_good_ph]/F");
    locTree->Branch("miss_mass_corr", loc_miss_mass_corr, "miss_mass_corr[num_good_ph]/F");

    locTree->Branch("miss_pm_corr2", loc_miss_pm_corr2, "miss_pm_corr2[num_good_ph]/F");
    locTree->Branch("miss_mass_corr2", loc_miss_mass_corr2, "miss_mass_corr2[num_good_ph]/F");
    
    locTree->Branch("index_sc", loc_index_sc, "index_sc[num_good_ph]/I");
    locTree->Branch("index_st", loc_index_st, "index_st[num_good_ph]/I");    

    locTree->Branch("index_sc2", loc_index_sc2, "index_sc2[num_good_ph]/I");
    locTree->Branch("index_st2", loc_index_st2, "index_st2[num_good_ph]/I");
    
    locTree->Branch("inv_mass_01", &loc_inv_mass_01, "inv_mass_01/F");
    locTree->Branch("inv_mass_02", &loc_inv_mass_02, "inv_mass_02/F");
    locTree->Branch("inv_mass_12", &loc_inv_mass_12, "inv_mass_12/F");
    locTree->Branch("inv_mass_01_corr", &loc_inv_mass_01_corr, "inv_mass_01_corr/F");
    locTree->Branch("inv_mass_02_corr", &loc_inv_mass_02_corr, "inv_mass_02_corr/F");
    locTree->Branch("inv_mass_12_corr", &loc_inv_mass_12_corr, "inv_mass_12_corr/F");

    locTree->Branch("inv_mass_01_corr2", &loc_inv_mass_01_corr2, "inv_mass_01_corr2/F");
    locTree->Branch("inv_mass_02_corr2", &loc_inv_mass_02_corr2, "inv_mass_02_corr2/F");
    locTree->Branch("inv_mass_12_corr2", &loc_inv_mass_12_corr2, "inv_mass_12_corr2/F");
    
    locTree->Branch("mvrt_chi2", &loc_mvrt_chi2, "mvrt_chi2/F");
    locTree->Branch("mvrt_x", &loc_mvrt_x, "mvrt_x/F");
    locTree->Branch("mvrt_y", &loc_mvrt_y, "mvrt_y/F");
    locTree->Branch("mvrt_z", &loc_mvrt_z, "mvrt_z/F");
    locTree->Branch("beam_x", &loc_beam_x, "beam_x/F");
    locTree->Branch("beam_y", &loc_beam_y, "beam_y/F");
    locTree->Branch("beam_z", &loc_beam_z, "beam_z/F");
    locTree->Branch("doca_hyperon", &loc_doca_hyperon, "doca_hyperon/F");
    locTree->Branch("doca_beam", &loc_doca_beam, "doca_beam/F");
    locTree->Branch("poca_hyperon_x", &loc_poca_hyperon_x, "poca_hyperon_x/F");
    locTree->Branch("poca_hyperon_y", &loc_poca_hyperon_y, "poca_hyperon_y/F");
    locTree->Branch("poca_hyperon_z", &loc_poca_hyperon_z, "poca_hyperon_z/F");
    locTree->Branch("poca_beam_x", &loc_poca_beam_x, "poca_beam_x/F");
    locTree->Branch("poca_beam_y", &loc_poca_beam_y, "poca_beam_y/F");
    locTree->Branch("poca_beam_z", &loc_poca_beam_z, "poca_beam_z/F");
    locTree->Branch("eloss_pm", loc_eloss_pm, "eloss_pm[3]/F");
    locTree->Branch("eloss_pm2", loc_eloss_pm2, "eloss_pm2[3]/F");
    locTree->Branch("corr_pm", loc_corr_pm, "loc_corr_pm[3]/F");
    locTree->Branch("corr_pm2", loc_corr_pm2, "loc_corr_pm2[3]/F");

    locTree->Branch("corr_phot_e", loc_corr_phot_e, "corr_phot_e[num_good_ph]/F");
    locTree->Branch("fiducial", loc_fiducial, "fiducial[3]/O");

    string loc_EPIC_char, loc_EPIC;
    const double  c_speed=29.9792458;

    int num_evnt_row, num_tagr_row, st_stat, dc_stat, sc_stat, tagrstat, temp_fid_sec;
    int proton_pos=-1, kaon_pos=-1, negative_pos=-1, positive1_pos=-1, positive2_pos=-1;
    float temp_mass[2], phot_prop, phot_prop2, temp_fid_charge, temp_fid_vert,  temp_phot_e;
    double temp_pm_corr;
    
    const double mass_protPDG= 0.93827201, mass_pionPDG=0.13957018, mass_kaonPDG=0.493677, mass_deutPDG=1.8761238;
    const int prot_numPDG=2212, pion_numPDG=-211, kaon_numPDG=321;
    TLorentzVector Proton_4P, Kaon_4P, Pion_4P, temp_4P, Lambda_4P, Missing_4P, Photon_4P, Deuteron_4P;
    TLorentzVector Proton_4P_corr, Kaon_4P_corr, Pion_4P_corr, Photon_4P_corr;
    TLorentzVector Proton_4P_eloss, Kaon_4P_eloss, Pion_4P_eloss;
    TLorentzVector Proton_4P_eloss2, Kaon_4P_eloss2, Pion_4P_eloss2;
    TLorentzVector Proton_4P_corr2, Kaon_4P_corr2, Pion_4P_corr2;

    TVector3 old_vert1, old_vert2, new_vert1, new_vert2, old_dir1, old_dir2;

    vector<float> abs_dt_sc, abs_dt_st, abs_dt_sc2, abs_dt_st2;
    float mctkpm, mctkcx, mctkcy, mctkcz, mctkmass;
    
	DecayTarget::Get_DecayTarget().Set_RunNumber(54070);

	locNumEvents = 0;
    
    Init_dbBOS();//initialize BOS banks (calls initbos() from /packages/c_bos_io/readbos.c)
	
    
    locNumEvents=0;
    for (loc_numbbos=0; loc_numbbos<numbbos; loc_numbbos++){
        Open_dbBOSInputFile(locStringInputs[loc_numbbos].c_str());	// open first BOS input file (0)
        while(getBOS(&bcs_,1,"E")!=0){
            //while(Get_dbBOSEvent()){// get next event
            if((locMaxNumEventsPerFile <= locNumEvents) && (locMaxNumEventsPerFile > 0)){//past max
                Clean_dbBOS();	//drop banks and clean them prior to getting the next event
                break;
            }
            locNumEvents++;
            if((Get_dbEVNTbank() != NULL)&&(Get_dbSTPBbank() != NULL)&&(Get_dbTAGRbank() != NULL)&&(Get_dbSCPBbank() != NULL) && (Get_dbDCPBbank() != NULL)){// if bank exists
                cout<<locNumEvents<<endl;

                num_evnt_row = Get_dbEVNTbank()->bank.nrow;//get number of tracks in EVNT bank for that event
                num_tagr_row = Get_dbTAGRbank()->bank.nrow;//get number of tracks in EVNT bank for that event

                loc_numofpart=0;
                loc_num_good_ph=0;
                loc_num_pos=0;
                loc_num_neg=0;
                loc_num_neu=0;
                temp_mass[0]=0.0;
                temp_mass[1]=0.0;
                
                loc_beam_x=0.0;
                loc_beam_y=0.0;
                loc_beam_z=-20.0;
                
                abs_dt_sc.resize(0);
                abs_dt_st.resize(0);
                abs_dt_sc2.resize(0);
                abs_dt_st2.resize(0);
                
                loc_num_ph_sc=0;
                loc_num_ph_st=0;

                loc_num_ph_sc2=0;
                loc_num_ph_st2=0;
                
                
                for(int loc_ievt=0; loc_ievt < num_evnt_row; loc_ievt++){//loop over all tracks
                    if (Get_dbEVNTbank()->evnt[loc_ievt].status > 0){
                        if ((Get_dbEVNTbank()->evnt[loc_ievt].scstat >0) && (Get_dbEVNTbank()->evnt[loc_ievt].dcstat >0) && (Get_dbEVNTbank()->evnt[loc_ievt].ststat >0) && (Get_dbEVNTbank()->evnt[loc_ievt].charge!=0)){
                            if (Get_dbEVNTbank()->evnt[loc_ievt].charge==1){
                                loc_num_pos++;
                                if(loc_num_pos==1){
                                    temp_mass[0]=Get_dbEVNTbank()->evnt[loc_ievt].pmom*Get_dbEVNTbank()->evnt[loc_ievt].pmom*(1.0-Get_dbEVNTbank()->evnt[loc_ievt].betta*Get_dbEVNTbank()->evnt[loc_ievt].betta)/(Get_dbEVNTbank()->evnt[loc_ievt].betta*Get_dbEVNTbank()->evnt[loc_ievt].betta);
                                    positive1_pos=loc_ievt;
                                }
                                else if(loc_num_pos==2){
                                    temp_mass[1]=Get_dbEVNTbank()->evnt[loc_ievt].pmom*Get_dbEVNTbank()->evnt[loc_ievt].pmom*(1.0-Get_dbEVNTbank()->evnt[loc_ievt].betta*Get_dbEVNTbank()->evnt[loc_ievt].betta)/(Get_dbEVNTbank()->evnt[loc_ievt].betta*Get_dbEVNTbank()->evnt[loc_ievt].betta);
                                    positive2_pos=loc_ievt;
                                    if(temp_mass[0]<temp_mass[1]){
                                        proton_pos=positive2_pos;
                                        kaon_pos=positive1_pos;
                                    }
                                    else{
                                        proton_pos=positive1_pos;
                                        kaon_pos=positive2_pos;
                                    }
                                }
                            }
                            else if (Get_dbEVNTbank()->evnt[loc_ievt].charge==-1){
                                loc_num_neg++;
                                negative_pos=loc_ievt;
                            }
                        }
                        else if (((Get_dbEVNTbank()->evnt[loc_ievt].ecstat >0)||(Get_dbEVNTbank()->evnt[loc_ievt].lcstat >0)) && (Get_dbEVNTbank()->evnt[loc_ievt].charge==0))
                            loc_num_neu++;
                    }
                }
                loc_numofpart=loc_num_pos+loc_num_neg+loc_num_neu;
                if (loc_numofpart >= 3 && loc_num_pos==2 && loc_num_neg==1){
                    ///Proton///
                    loc_evnt_id[0] = Get_dbEVNTbank()->evnt[proton_pos].id;
                    loc_evnt_charge[0] = Get_dbEVNTbank()->evnt[proton_pos].charge;
                    loc_evnt_pm[0] = Get_dbEVNTbank()->evnt[proton_pos].pmom;
                    loc_evnt_be[0] = Get_dbEVNTbank()->evnt[proton_pos].betta;
                    loc_evnt_mass[0] = Get_dbEVNTbank()->evnt[proton_pos].mass;
                    loc_calc_mass[0]=loc_evnt_pm[0]*loc_evnt_pm[0]*(1.0-loc_evnt_be[0]*loc_evnt_be[0])/(loc_evnt_be[0]*loc_evnt_be[0]);
                    loc_evnt_x[0] = Get_dbEVNTbank()->evnt[proton_pos].vert.x;
                    loc_evnt_y[0] = Get_dbEVNTbank()->evnt[proton_pos].vert.y;
                    loc_evnt_z[0] = Get_dbEVNTbank()->evnt[proton_pos].vert.z;
                    loc_evnt_cx[0] = Get_dbEVNTbank()->evnt[proton_pos].dir_cos.x;
                    loc_evnt_cy[0] = Get_dbEVNTbank()->evnt[proton_pos].dir_cos.y;
                    loc_evnt_cz[0] = Get_dbEVNTbank()->evnt[proton_pos].dir_cos.z;

                    st_stat = Get_dbEVNTbank()->evnt[proton_pos].ststat;
                    sc_stat = Get_dbEVNTbank()->evnt[proton_pos].scstat;
                    dc_stat = Get_dbEVNTbank()->evnt[proton_pos].dcstat;

                    loc_scpb_t[0] = Get_dbSCPBbank()->scpb[sc_stat-1].time;
                    loc_scpb_d[0] = Get_dbSCPBbank()->scpb[sc_stat-1].path;
                    loc_stpb_t[0] = Get_dbSTPBbank()->stpb[st_stat-1].time;
                    loc_stpb_d[0] = Get_dbSTPBbank()->stpb[st_stat-1].path;
                    loc_scpb_id[0] = static_cast<int>(((Get_dbSCPBbank()->scpb[sc_stat-1].scpdht)%10000)/100);
                    loc_stpb_id[0] = static_cast<int>(((Get_dbSTPBbank()->stpb[st_stat-1].sthid)%100)/10)+(static_cast<int>((Get_dbSTPBbank()->stpb[st_stat-1].sthid)/100)-1)*4;
                    loc_dcpb_sect[0] = (Get_dbDCPBbank()->dcpb[dc_stat-1].sctr)/100;
                    loc_calc_be[0]=loc_evnt_pm[0]/sqrt(pow(loc_evnt_pm[0],2)+mass_protPDG*mass_protPDG);
                    loc_delt_be[0]=loc_evnt_be[0]-loc_calc_be[0];
                    Proton_4P.SetXYZM(loc_evnt_pm[0]*loc_evnt_cx[0], loc_evnt_pm[0]*loc_evnt_cy[0], loc_evnt_pm[0]*loc_evnt_cz[0], mass_protPDG);
                    ///Kaon///
                    loc_evnt_id[1] = Get_dbEVNTbank()->evnt[kaon_pos].id;
                    loc_evnt_charge[1] = Get_dbEVNTbank()->evnt[kaon_pos].charge;
                    loc_evnt_pm[1] = Get_dbEVNTbank()->evnt[kaon_pos].pmom;
                    loc_evnt_be[1] = Get_dbEVNTbank()->evnt[kaon_pos].betta;
                    loc_evnt_mass[1] = Get_dbEVNTbank()->evnt[kaon_pos].mass;
                    loc_calc_mass[1]=loc_evnt_pm[1]*loc_evnt_pm[1]*(1.0-loc_evnt_be[1]*loc_evnt_be[1])/(loc_evnt_be[1]*loc_evnt_be[1]);
                    loc_evnt_x[1] = Get_dbEVNTbank()->evnt[kaon_pos].vert.x;
                    loc_evnt_y[1] = Get_dbEVNTbank()->evnt[kaon_pos].vert.y;
                    loc_evnt_z[1] = Get_dbEVNTbank()->evnt[kaon_pos].vert.z;
                    loc_evnt_cx[1] = Get_dbEVNTbank()->evnt[kaon_pos].dir_cos.x;
                    loc_evnt_cy[1] = Get_dbEVNTbank()->evnt[kaon_pos].dir_cos.y;
                    loc_evnt_cz[1] = Get_dbEVNTbank()->evnt[kaon_pos].dir_cos.z;
                    
                    st_stat = Get_dbEVNTbank()->evnt[kaon_pos].ststat;
                    sc_stat = Get_dbEVNTbank()->evnt[kaon_pos].scstat;
                    dc_stat = Get_dbEVNTbank()->evnt[kaon_pos].dcstat;
                    
                    loc_scpb_t[1] = Get_dbSCPBbank()->scpb[sc_stat-1].time;
                    loc_scpb_d[1] = Get_dbSCPBbank()->scpb[sc_stat-1].path;
                    loc_stpb_t[1] = Get_dbSTPBbank()->stpb[st_stat-1].time;
                    loc_stpb_d[1] = Get_dbSTPBbank()->stpb[st_stat-1].path;
                    loc_scpb_id[1] = static_cast<int>(((Get_dbSCPBbank()->scpb[sc_stat-1].scpdht)%10000)/100);
                    loc_stpb_id[1] = static_cast<int>(((Get_dbSTPBbank()->stpb[st_stat-1].sthid)%100)/10)+(static_cast<int>((Get_dbSTPBbank()->stpb[st_stat-1].sthid)/100)-1)*4;
                    loc_dcpb_sect[1] = (Get_dbDCPBbank()->dcpb[dc_stat-1].sctr)/100;
                    loc_calc_be[1]=loc_evnt_pm[1]/sqrt(pow(loc_evnt_pm[1],2)+mass_kaonPDG*mass_kaonPDG);
                    loc_delt_be[1]=loc_evnt_be[1]-loc_calc_be[1];
                    Kaon_4P.SetXYZM(loc_evnt_pm[1]*loc_evnt_cx[1], loc_evnt_pm[1]*loc_evnt_cy[1], loc_evnt_pm[1]*loc_evnt_cz[1], mass_kaonPDG);

                    ///Pion///
                    loc_evnt_id[2] = Get_dbEVNTbank()->evnt[negative_pos].id;
                    loc_evnt_charge[2] = Get_dbEVNTbank()->evnt[negative_pos].charge;
                    loc_evnt_pm[2] = Get_dbEVNTbank()->evnt[negative_pos].pmom;
                    loc_evnt_be[2] = Get_dbEVNTbank()->evnt[negative_pos].betta;
                    loc_evnt_mass[2] = Get_dbEVNTbank()->evnt[negative_pos].mass;
                    loc_calc_mass[2]=loc_evnt_pm[2]*loc_evnt_pm[2]*(1.0-loc_evnt_be[2]*loc_evnt_be[2])/(loc_evnt_be[2]*loc_evnt_be[2]);
                    loc_evnt_x[2] = Get_dbEVNTbank()->evnt[negative_pos].vert.x;
                    loc_evnt_y[2] = Get_dbEVNTbank()->evnt[negative_pos].vert.y;
                    loc_evnt_z[2] = Get_dbEVNTbank()->evnt[negative_pos].vert.z;
                    loc_evnt_cx[2] = Get_dbEVNTbank()->evnt[negative_pos].dir_cos.x;
                    loc_evnt_cy[2] = Get_dbEVNTbank()->evnt[negative_pos].dir_cos.y;
                    loc_evnt_cz[2] = Get_dbEVNTbank()->evnt[negative_pos].dir_cos.z;
                    
                    st_stat = Get_dbEVNTbank()->evnt[negative_pos].ststat;
                    sc_stat = Get_dbEVNTbank()->evnt[negative_pos].scstat;
                    dc_stat = Get_dbEVNTbank()->evnt[negative_pos].dcstat;
                    
                    loc_scpb_t[2] = Get_dbSCPBbank()->scpb[sc_stat-1].time;
                    loc_scpb_d[2] = Get_dbSCPBbank()->scpb[sc_stat-1].path;
                    loc_stpb_t[2] = Get_dbSTPBbank()->stpb[st_stat-1].time;
                    loc_stpb_d[2] = Get_dbSTPBbank()->stpb[st_stat-1].path;
                    loc_scpb_id[2] = static_cast<int>(((Get_dbSCPBbank()->scpb[sc_stat-1].scpdht)%10000)/100);
                    loc_stpb_id[2] = static_cast<int>(((Get_dbSTPBbank()->stpb[st_stat-1].sthid)%100)/10)+(static_cast<int>((Get_dbSTPBbank()->stpb[st_stat-1].sthid)/100)-1)*4;
                    loc_dcpb_sect[2] = (Get_dbDCPBbank()->dcpb[dc_stat-1].sctr)/100;
                    loc_calc_be[2]=loc_evnt_pm[2]/sqrt(pow(loc_evnt_pm[2],2)+mass_pionPDG*mass_pionPDG);
                    loc_delt_be[2]=loc_evnt_be[2]-loc_calc_be[2];
                    Pion_4P.SetXYZM(loc_evnt_pm[2]*loc_evnt_cx[2], loc_evnt_pm[2]*loc_evnt_cy[2], loc_evnt_pm[2]*loc_evnt_cz[2], mass_pionPDG);

                    ////////////
                    
                    
                    loc_mvrt_x = 0;
                    loc_mvrt_y = 0;
                    loc_mvrt_z = 0;
                    loc_mvrt_chi2 = 0;
                    loc_head_eventnum = Get_dbHEADbank()->head[0].nevent;
                    loc_head_runnum = Get_dbHEADbank()->head[0].nrun;
                    loc_head_time = Get_dbHEADbank()->head[0].time;
                    
                    temp_4P=Proton_4P;
                    temp_fid_charge=loc_evnt_charge[0]*1.0;
                    temp_fid_sec=loc_dcpb_sect[0];
                    temp_fid_vert=loc_evnt_z[0];
                    loc_fiducial[0]=Cut_dfcTrackDirection(loc_head_runnum, temp_fid_charge, temp_fid_sec,temp_fid_vert, temp_4P);
                    
                    temp_4P=Kaon_4P;
                    temp_fid_charge=loc_evnt_charge[1]*1.0;
                    temp_fid_sec=loc_dcpb_sect[1];
                    temp_fid_vert=loc_evnt_z[1];
                    loc_fiducial[1]=Cut_dfcTrackDirection(loc_head_runnum, temp_fid_charge, temp_fid_sec,temp_fid_vert, temp_4P);
                    
                    temp_4P=Pion_4P;
                    temp_fid_charge=loc_evnt_charge[2]*1.0;
                    temp_fid_sec=loc_dcpb_sect[2];
                    temp_fid_vert=loc_evnt_z[2];
                    loc_fiducial[2]=Cut_dfcTrackDirection(loc_head_runnum, temp_fid_charge, temp_fid_sec,temp_fid_vert, temp_4P);
                    
                    temp_4P=Proton_4P+Kaon_4P;
                    loc_inv_mass_01=temp_4P.M();
                    temp_4P=Pion_4P+Kaon_4P;
                    loc_inv_mass_12=temp_4P.M();
                    Lambda_4P=Proton_4P+Pion_4P;
                    loc_inv_mass_02=Lambda_4P.M();
                    
                    Deuteron_4P.SetXYZM(0.0,0.0,0.0,mass_deutPDG);

                    loc_beam_x=0.0;
                    loc_beam_y=0.0;

                    /// Get doca vertices/////////
                    old_vert1.SetXYZ(loc_evnt_x[0], loc_evnt_y[0], loc_evnt_z[0]);
                    old_vert2.SetXYZ(loc_evnt_x[2], loc_evnt_y[2], loc_evnt_z[2]);
                    old_dir1.SetXYZ(loc_evnt_cx[0], loc_evnt_cy[0], loc_evnt_cz[0]);
                    old_dir2.SetXYZ(loc_evnt_cx[2], loc_evnt_cy[2], loc_evnt_cz[2]);
                    
                    loc_doca_hyperon=Calc_dtfInterDOCA(old_dir1, old_dir2, old_vert1, old_vert2, new_vert1, new_vert2);
                    loc_poca_hyperon_x=new_vert1.X()-(new_vert1.X()-new_vert2.X())/2.0;
                    loc_poca_hyperon_y=new_vert1.Y()-(new_vert1.Y()-new_vert2.Y())/2.0;
                    loc_poca_hyperon_z=new_vert1.Z()-(new_vert1.Z()-new_vert2.Z())/2.0;
                    
                    old_vert1.SetXYZ(loc_evnt_x[1], loc_evnt_y[1], loc_evnt_z[1]);
                    old_vert2.SetXYZ(loc_beam_x, loc_beam_y, loc_beam_z);
                    old_dir1.SetXYZ(loc_evnt_cx[1], loc_evnt_cy[1], loc_evnt_cz[1]);
                    old_dir2.SetXYZ(0, 0, 1);
                    
                    loc_doca_beam=Calc_dtfInterDOCA(old_dir1, old_dir2, old_vert1, old_vert2, new_vert1, new_vert2);
                    loc_poca_beam_x=new_vert1.X()-(new_vert1.X()-new_vert2.X())/2.0;
                    loc_poca_beam_y=new_vert1.Y()-(new_vert1.Y()-new_vert2.Y())/2.0;
                    loc_poca_beam_z=new_vert1.Z()-(new_vert1.Z()-new_vert2.Z())/2.0;
                    
                    phot_prop=(20.0+loc_poca_beam_z)/c_speed;
                    phot_prop2=(20.0+loc_evnt_z[1])/c_speed;

                    //////Correct Momenta////////////
                    Proton_4P_eloss=Proton_4P;
                    old_vert1.SetXYZ(loc_poca_hyperon_x, loc_poca_hyperon_y, loc_poca_hyperon_z);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Proton_4P_eloss, old_vert1);
                    loc_eloss_pm[0]=Proton_4P_eloss.Rho();
                    
                    Proton_4P_corr=Proton_4P;
                    temp_pm_corr=Correct_dtcTrackMomentum(54065, prot_numPDG, loc_dcpb_sect[0], Proton_4P_corr);
                    Proton_4P_corr.SetXYZM(temp_pm_corr*loc_evnt_cx[0], temp_pm_corr*loc_evnt_cy[0], temp_pm_corr*loc_evnt_cz[0], mass_protPDG);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Proton_4P_corr, old_vert1);
                    loc_corr_pm[0]=Proton_4P_corr.Rho();
                    
                    Proton_4P_eloss2=Proton_4P;
                    old_vert1.SetXYZ(loc_evnt_x[0], loc_evnt_y[0], loc_evnt_z[0]);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Proton_4P_eloss2, old_vert1);
                    loc_eloss_pm2[0]=Proton_4P_eloss2.Rho();
                    
                    Proton_4P_corr2=Proton_4P;
                    temp_pm_corr=Correct_dtcTrackMomentum(54065, prot_numPDG, loc_dcpb_sect[0], Proton_4P_corr2);
                    Proton_4P_corr2.SetXYZM(temp_pm_corr*loc_evnt_cx[0], temp_pm_corr*loc_evnt_cy[0], temp_pm_corr*loc_evnt_cz[0], mass_protPDG);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Proton_4P_corr2, old_vert1);
                    loc_corr_pm2[0]=Proton_4P_corr2.Rho();
                    
                    
                    
      
                    Pion_4P_eloss=Pion_4P;
                    old_vert1.SetXYZ(loc_poca_hyperon_x, loc_poca_hyperon_y, loc_poca_hyperon_z);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Pion_4P_eloss, old_vert1);
                    loc_eloss_pm[2]=Pion_4P_eloss.Rho();
                    
                    Pion_4P_corr=Pion_4P;
                    temp_pm_corr=Correct_dtcTrackMomentum(54065, pion_numPDG, loc_dcpb_sect[2], Pion_4P_corr);
                    Pion_4P_corr.SetXYZM(temp_pm_corr*loc_evnt_cx[2], temp_pm_corr*loc_evnt_cy[2], temp_pm_corr*loc_evnt_cz[2], mass_pionPDG);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Pion_4P_corr, old_vert1);
                    loc_corr_pm[2]=Pion_4P_corr.Rho();
                    
                    
                    Pion_4P_eloss2=Pion_4P;
                    old_vert1.SetXYZ(loc_evnt_x[2], loc_evnt_y[2], loc_evnt_z[2]);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Pion_4P_eloss2, old_vert1);
                    loc_eloss_pm2[2]=Pion_4P_eloss2.Rho();
                    
                    Pion_4P_corr2=Pion_4P;
                    temp_pm_corr=Correct_dtcTrackMomentum(54065, pion_numPDG, loc_dcpb_sect[2], Pion_4P_corr2);
                    Pion_4P_corr2.SetXYZM(temp_pm_corr*loc_evnt_cx[2], temp_pm_corr*loc_evnt_cy[2], temp_pm_corr*loc_evnt_cz[2], mass_pionPDG);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Pion_4P_corr2, old_vert1);
                    loc_corr_pm2[2]=Pion_4P_corr2.Rho();
                    
                              
                    
                    Kaon_4P_eloss=Kaon_4P;
                    old_vert1.SetXYZ(loc_poca_beam_x, loc_poca_beam_y, loc_poca_beam_z);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Kaon_4P_eloss, old_vert1);
                    loc_eloss_pm[1]=Kaon_4P_eloss.Rho();
                    
                    Kaon_4P_corr=Kaon_4P;
                    temp_pm_corr=Correct_dtcTrackMomentum(54065, kaon_numPDG, loc_dcpb_sect[1], Kaon_4P_corr);
                    Kaon_4P_corr.SetXYZM(temp_pm_corr*loc_evnt_cx[1], temp_pm_corr*loc_evnt_cy[1], temp_pm_corr*loc_evnt_cz[1], mass_kaonPDG);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Kaon_4P_corr, old_vert1);
                    loc_corr_pm[1]=Kaon_4P_corr.Rho();
                    
                    Kaon_4P_eloss2=Kaon_4P;
                    old_vert1.SetXYZ(loc_evnt_x[1], loc_evnt_y[1], loc_evnt_z[1]);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Kaon_4P_eloss2, old_vert1);
                    loc_eloss_pm2[1]=Kaon_4P_eloss2.Rho();
                    
                    Kaon_4P_corr2=Kaon_4P;
                    temp_pm_corr=Correct_dtcTrackMomentum(54065, kaon_numPDG, loc_dcpb_sect[1], Kaon_4P_corr2);
                    Kaon_4P_corr2.SetXYZM(temp_pm_corr*loc_evnt_cx[1], temp_pm_corr*loc_evnt_cy[1], temp_pm_corr*loc_evnt_cz[1], mass_kaonPDG);
                    DecayTarget::Get_DecayTarget().Correct_Momcor(Kaon_4P_corr2, old_vert1);
                    loc_corr_pm2[1]=Kaon_4P_corr2.Rho();
                    
               
                    temp_4P=Proton_4P_corr+Kaon_4P_corr;
                    loc_inv_mass_01_corr=temp_4P.M();
                    temp_4P=Pion_4P_corr+Kaon_4P_corr;
                    loc_inv_mass_12_corr=temp_4P.M();
                    temp_4P=Proton_4P_corr+Pion_4P_corr;
                    loc_inv_mass_02_corr=temp_4P.M();

                    temp_4P=Proton_4P_corr2+Kaon_4P_corr2;
                    loc_inv_mass_01_corr2=temp_4P.M();
                    temp_4P=Pion_4P_corr2+Kaon_4P_corr2;
                    loc_inv_mass_12_corr2=temp_4P.M();
                    temp_4P=Proton_4P_corr2+Pion_4P_corr2;
                    loc_inv_mass_02_corr2=temp_4P.M();
                    
                    
                    
                    /////////////////////////////////////////////////////////////
                    for (int loc_xxm=0; loc_xxm<num_tagr_row; loc_xxm++){
                        tagrstat = Get_dbTAGRbank()->tagr[loc_xxm].stat;
                        if (tagrstat == 7 || tagrstat == 15){
                            loc_tagr_tid[loc_num_good_ph] = Get_dbTAGRbank()->tagr[loc_xxm].t_id;
                            loc_tagr_eid[loc_num_good_ph] = Get_dbTAGRbank()->tagr[loc_xxm].e_id;
                            loc_tagr_ttag[loc_num_good_ph] = Get_dbTAGRbank()->tagr[loc_xxm].ttag;
                            loc_tagr_tpho[loc_num_good_ph] =  Get_dbTAGRbank()->tagr[loc_xxm].tpho;
                            loc_phot_e[loc_num_good_ph] = Get_dbTAGRbank()->tagr[loc_xxm].erg;
                            loc_tagr_stat[loc_num_good_ph] = Get_dbTAGRbank()->tagr[loc_xxm].stat;
                            
                            loc_delt_t_sc[loc_num_good_ph]=loc_scpb_t[1]-loc_scpb_d[1]/(loc_calc_be[1]*c_speed)-(loc_tagr_tpho[loc_num_good_ph]+phot_prop);
                            
                            loc_delt_t_st[loc_num_good_ph]=loc_stpb_t[1]-loc_stpb_d[1]/(loc_calc_be[1]*c_speed)-(loc_tagr_tpho[loc_num_good_ph]+phot_prop);
                            
                            loc_delt_t_sc2[loc_num_good_ph]=loc_scpb_t[1]-loc_scpb_d[1]/(loc_calc_be[1]*c_speed)-(loc_tagr_tpho[loc_num_good_ph]+phot_prop2);
                            
                            loc_delt_t_st2[loc_num_good_ph]=loc_stpb_t[1]-loc_stpb_d[1]/(loc_calc_be[1]*c_speed)-(loc_tagr_tpho[loc_num_good_ph]+phot_prop2);
                            
                            abs_dt_sc.push_back(fabs(loc_delt_t_sc[loc_num_good_ph]));
                            abs_dt_st.push_back(fabs(loc_delt_t_st[loc_num_good_ph]));
                            
                            abs_dt_sc2.push_back(fabs(loc_delt_t_sc2[loc_num_good_ph]));
                            abs_dt_st2.push_back(fabs(loc_delt_t_st2[loc_num_good_ph]));
                            
                            if (abs_dt_sc[loc_num_good_ph]<1.0)
                                loc_num_ph_sc++;
                            if (abs_dt_st[loc_num_good_ph]<1.0)
                                loc_num_ph_st++;
                            if (abs_dt_sc2[loc_num_good_ph]<1.0)
                                loc_num_ph_sc2++;
                            if (abs_dt_st2[loc_num_good_ph]<1.0)
                                loc_num_ph_st2++;
                            
                            Photon_4P.SetXYZM(0.0,0.0,loc_phot_e[loc_num_good_ph],0.0);
                            Missing_4P=Photon_4P+Deuteron_4P-Proton_4P-Kaon_4P-Pion_4P;
                            loc_miss_mass[loc_num_good_ph]=Missing_4P.M();
                            loc_miss_pm[loc_num_good_ph]=Missing_4P.Rho();
                            
                            temp_phot_e=loc_phot_e[loc_num_good_ph];
                            loc_corr_phot_e[loc_num_good_ph]=Correct_dtcPhotonEnergy(loc_head_runnum, temp_phot_e, 70);

                            Photon_4P.SetXYZM(0.0,0.0,loc_corr_phot_e[loc_num_good_ph],0.0);
                            Missing_4P=Photon_4P+Deuteron_4P-Proton_4P_corr-Kaon_4P_corr-Pion_4P_corr;
                            loc_miss_mass_corr[loc_num_good_ph]=Missing_4P.M();
                            loc_miss_pm_corr[loc_num_good_ph]=Missing_4P.Rho();

                            Missing_4P=Photon_4P+Deuteron_4P-Proton_4P_corr2-Kaon_4P_corr2-Pion_4P_corr2;
                            loc_miss_mass_corr2[loc_num_good_ph]=Missing_4P.M();
                            loc_miss_pm_corr2[loc_num_good_ph]=Missing_4P.Rho();
                            
                            
                            
                            loc_num_good_ph++;
                        }
					}

                    int temp_index_sc[loc_num_good_ph], temp_index_st[loc_num_good_ph];
                    int temp_index_sc2[loc_num_good_ph], temp_index_st2[loc_num_good_ph];

                    float temp_dt_sc[loc_num_good_ph];
                    float temp_dt_st[loc_num_good_ph];
                    
                    float temp_dt_sc2[loc_num_good_ph];
                    float temp_dt_st2[loc_num_good_ph];
                    
                    for (int loc_xxm=0; loc_xxm<loc_num_good_ph; loc_xxm++){
                        temp_dt_sc[loc_xxm]=abs_dt_sc[loc_xxm];
                        temp_dt_st[loc_xxm]=abs_dt_st[loc_xxm];
                        temp_dt_sc2[loc_xxm]=abs_dt_sc2[loc_xxm];
                        temp_dt_st2[loc_xxm]=abs_dt_st2[loc_xxm];
                    }

                    
                    TMath::Sort(loc_num_good_ph,temp_dt_sc,temp_index_sc, kFALSE);  //kFALSE for Increasing order
                    TMath::Sort(loc_num_good_ph,temp_dt_st,temp_index_st, kFALSE);  //kFALSE for Increasing order
                    TMath::Sort(loc_num_good_ph,temp_dt_sc2,temp_index_sc2, kFALSE);  //kFALSE for Increasing order
                    TMath::Sort(loc_num_good_ph,temp_dt_st2,temp_index_st2, kFALSE);  //kFALSE for Increasing order
                    
                    for (int loc_xxm=0; loc_xxm<loc_num_good_ph; loc_xxm++){
                        loc_index_sc[loc_xxm]=temp_index_sc[loc_xxm];
                        loc_index_st[loc_xxm]=temp_index_st[loc_xxm];
                        loc_index_sc2[loc_xxm]=temp_index_sc2[loc_xxm];
                        loc_index_st2[loc_xxm]=temp_index_st2[loc_xxm];
                    }
                    
                    loc_label=-1;
                    
                    if(Get_dbMCTKbank() != NULL){// if bank exists
                        mctkpm=Get_dbMCTKbank()->mctk[0].pmom;
                        mctkcx=Get_dbMCTKbank()->mctk[0].cx;
                        mctkcy=Get_dbMCTKbank()->mctk[0].cy;
                        mctkcz=Get_dbMCTKbank()->mctk[0].cz;
                        mctkmass=Get_dbMCTKbank()->mctk[0].mass;
                        loc_beam->SetXYZM(mctkpm*mctkcx,mctkpm*mctkcy,mctkpm*mctkcz,mctkpm*mctkmass);

                        mctkpm=Get_dbMCTKbank()->mctk[1].pmom;
                        mctkcx=Get_dbMCTKbank()->mctk[1].cx;
                        mctkcy=Get_dbMCTKbank()->mctk[1].cy;
                        mctkcz=Get_dbMCTKbank()->mctk[1].cz;
                        mctkmass=Get_dbMCTKbank()->mctk[1].mass;
                        loc_target->SetXYZM(mctkpm*mctkcx,mctkpm*mctkcy,mctkpm*mctkcz,mctkpm*mctkmass);
                    
                        mctkpm=Get_dbMCTKbank()->mctk[2].pmom;
                        mctkcx=Get_dbMCTKbank()->mctk[2].cx;
                        mctkcy=Get_dbMCTKbank()->mctk[2].cy;
                        mctkcz=Get_dbMCTKbank()->mctk[2].cz;
                        mctkmass=Get_dbMCTKbank()->mctk[2].mass;
                        loc_Qp1->SetXYZM(mctkpm*mctkcx,mctkpm*mctkcy,mctkpm*mctkcz,mctkpm*mctkmass);
                        
                        mctkpm=Get_dbMCTKbank()->mctk[3].pmom;
                        mctkcx=Get_dbMCTKbank()->mctk[3].cx;
                        mctkcy=Get_dbMCTKbank()->mctk[3].cy;
                        mctkcz=Get_dbMCTKbank()->mctk[3].cz;
                        mctkmass=Get_dbMCTKbank()->mctk[3].mass;
                        loc_Qp2->SetXYZM(mctkpm*mctkcx,mctkpm*mctkcy,mctkpm*mctkcz,mctkpm*mctkmass);

                        mctkpm=Get_dbMCTKbank()->mctk[4].pmom;
                        mctkcx=Get_dbMCTKbank()->mctk[4].cx;
                        mctkcy=Get_dbMCTKbank()->mctk[4].cy;
                        mctkcz=Get_dbMCTKbank()->mctk[4].cz;
                        mctkmass=Get_dbMCTKbank()->mctk[4].mass;
                        loc_Sbeam->SetXYZM(mctkpm*mctkcx,mctkpm*mctkcy,mctkpm*mctkcz,mctkpm*mctkmass);
                    
                        mctkpm=Get_dbMCTKbank()->mctk[5].pmom;
                        mctkcx=Get_dbMCTKbank()->mctk[5].cx;
                        mctkcy=Get_dbMCTKbank()->mctk[5].cy;
                        mctkcz=Get_dbMCTKbank()->mctk[5].cz;
                        mctkmass=Get_dbMCTKbank()->mctk[5].mass;
                        loc_Starget->SetXYZM(mctkpm*mctkcx,mctkpm*mctkcy,mctkpm*mctkcz,mctkpm*mctkmass);

                        mctkpm=Get_dbMCTKbank()->mctk[6].pmom;
                        mctkcx=Get_dbMCTKbank()->mctk[6].cx;
                        mctkcy=Get_dbMCTKbank()->mctk[6].cy;
                        mctkcz=Get_dbMCTKbank()->mctk[6].cz;
                        mctkmass=Get_dbMCTKbank()->mctk[6].mass;
                        loc_Rp1->SetXYZM(mctkpm*mctkcx,mctkpm*mctkcy,mctkpm*mctkcz,mctkpm*mctkmass);
                        
                        mctkpm=Get_dbMCTKbank()->mctk[7].pmom;
                        mctkcx=Get_dbMCTKbank()->mctk[7].cx;
                        mctkcy=Get_dbMCTKbank()->mctk[7].cy;
                        mctkcz=Get_dbMCTKbank()->mctk[7].cz;
                        mctkmass=Get_dbMCTKbank()->mctk[7].mass;
                        loc_Rp2->SetXYZM(mctkpm*mctkcx,mctkpm*mctkcy,mctkpm*mctkcz,mctkpm*mctkmass);
                        
                        mctkpm=Get_dbMCTKbank()->mctk[8].pmom;
                        mctkcx=Get_dbMCTKbank()->mctk[8].cx;
                        mctkcy=Get_dbMCTKbank()->mctk[8].cy;
                        mctkcz=Get_dbMCTKbank()->mctk[8].cz;
                        mctkmass=Get_dbMCTKbank()->mctk[8].mass;
                        loc_Rp3->SetXYZM(mctkpm*mctkcx,mctkpm*mctkcy,mctkpm*mctkcz,mctkpm*mctkmass);
                    
                        loc_label=Get_dbMCTKbank()->mctk[4].flag;
                    }
                 
                    if ((loc_delt_be[0]>-0.1&&loc_delt_be[0]<0.1)&&(loc_delt_be[1]>-0.1&&loc_delt_be[1]<0.1)&&(loc_delt_be[2]>-0.2&&loc_delt_be[2]<0.2)&&(loc_inv_mass_02<1.2||loc_inv_mass_12<1.2)&&(loc_label>0)){
                        locTree->Fill();

                    }
                    
                                       
                    
                    
                    
                }
            }
    
    
            Clean_dbBOS();// drop banks and clean them prior to getting the next event

        }
        Close_dbBOSInputFile();	// close the bos output file
    }

    
    
    
	locFile->Write();																			//save everything in file
	locFile->Close();																			//close file (I think this will delete all pointes to hist and ntuples)
	
	delete locFile;																				//delete pointer
    
	return 1;

}


void Display_Help(){
	cout << "DISPLAY HELP:" << endl;
	cout << "Command Line: ExecutableFile BOSInputFile" << endl;
	cout << "The optional '-M' flag is used for specifying the maximum number of events to be evaluated from each file." << endl;
}


double Calc_dtfInterDOCA(const TVector3 &locUnitDir1, const TVector3 &locUnitDir2, const TVector3 &locVertex1, const TVector3 &locVertex2, TVector3 &locInterDOCA1, TVector3 &locInterDOCA2){
    //originated from code by Jörn Langheinrich
    //you can use this function to find the DOCA to a fixed point by calling this function with locUnitDir1 and 2 parallel, and the fixed vertex as locVertex2
    double locUnitDot = locUnitDir1*locUnitDir2;
    double locDenominator = locUnitDot*locUnitDot - 1.0; /// scalar product of directions
    double locDistVertToInterDOCA1 = 0.0, locDistVertToInterDOCA2 = 0.0; //distance from vertex to DOCA point
    
    if(fabs(locDenominator) < 1.0e-15) //parallel
        locDistVertToInterDOCA1 = (locVertex2 - locVertex1)*locUnitDir2/locUnitDot; //the opposite
    else{
        double locA = (locVertex1 - locVertex2)*locUnitDir1;
        double locB = (locVertex1 - locVertex2)*locUnitDir2;
        locDistVertToInterDOCA1 = (locA - locUnitDot*locB)/locDenominator;
        locDistVertToInterDOCA2 = (locUnitDot*locA - locB)/locDenominator;
    }
    
    locInterDOCA1 = locVertex1 + locDistVertToInterDOCA1*locUnitDir1; //intersection point of DOCA line and track 1
    locInterDOCA2 = locVertex2 + locDistVertToInterDOCA2*locUnitDir2; //intersection point of DOCA line and track 2
    double locDOCA = (locInterDOCA1 - locInterDOCA2).Mag();
    return ((locVertex2.Z() > locVertex1.Z()) ? locDOCA : -1.0*locDOCA);
}


                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    
                    




















