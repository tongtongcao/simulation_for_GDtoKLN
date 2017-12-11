source /u/group/clas/builds/PRODUCTION/packages/tcl/recsis_proc.tcl
turnoff ALL;
global_section off;
turnon seb tof egn trk user pid start;
set lst_do -1;
set ltime_do -1;
setc ddl_file "/work/clas/production/claseg3/PARMS/clasbanks.ddl";

set photon_trig_type 1;

#reject all tagger photons that are outside of the st_tagger_match +/- window (ns)
set st_tagger_match 30.;
#e-t coincidence window (i assume e-counter, t-counter, in ns)
set tagger_etwin 10.;

inputfile infile;
outputfile outfile PROC 2047;
#setc chist_filename gsim.hbook;
#setc log_file_name user_ana.log;
setc outbanknames(1) "HEADHEVTEVNTDCPBECPBSCPBCCPBMCHDMCTKMCVXTAGRPARTTBTRTBERSTRESTPBTDPLTBIDSCRCECHBTGPB";
set lseb_nt_do -1;
set lall_nt_do -1;
set lmctk_nt_do -1;

#level of analysis 0: raw  2: hbt 4: tbt
set trk_level 4;
#left-right ambiguity fit chisq cut (from pawel 50, default 20, g13 20, g11 50)
set trk_lrambfit_chi2 50.;
#chisq cut in final track fit (from pawel 50, default 30, g13 40, g11 70)
set trk_tbtfit_chi2   70.;
#chisq cut in pattern recognition fit (from pawel 20, default 10, g13 40, g11 70)
set trk_prfit_chi2    70.;
#minimum number of segments with resolved left-right ambiguity (from pawel 4, default 5, g13 5)
set trk_minlramb 5;
#print out track statistics at the end of the run (g11 3)
set trk_statistics 3;

set torus_current -1500;
set mini_torus_current 0;
set poltarget_current 0;
set TargetPos(3) -20.;

setc bfield_file_name "bgrid_T67to33.fpk";
#setc prlink_file_name "prlink_eg3-1930_208M.bos";
setc prlink_file_name "prlink_tg-20pm30n.bos";
set trk_prlink_param 50;

#drift time-distance function (0 = linear (gsim))
set dc_xvst_choice 0;
#flag used to enable/disable fixed attenuation lengths (ec)
set def_atten -1;
#flag used to enable/disable default adc values (ec)
set def_adc -1;
#flag used to enable/disable default tdc values (ec)
set def_tdc -1;
#flag used to enable/disable ideal geometry (ec)
set def_geom -1;

# tell FPACK not to stop if it thinks you are running out of time
fpack "timestop -9999999999"
# do not send events to event display
set lscat $false;
set ldisplay_all $false;
# tell recsis to pause or go
setc rec_prompt "[exec whoami]_recsis> ";
go;
exit_pend;  
