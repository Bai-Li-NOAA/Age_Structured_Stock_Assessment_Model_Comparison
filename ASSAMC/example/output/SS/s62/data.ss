#V3.30
#C data file created using the SS_writedat function in the R package r4ss
#C should work with SS version: 
#C file write time: 2020-07-06 04:44:01
#
1 #_styr
30 #_endyr
1 #_nseas
12 #_months_per_seas
2 #_Nsubseasons
1.00001 #_spawn_month
-1 #_Nsexes
12 #_Nages
1 #_Nareas
2 #_Nfleets
#_fleetinfo
#_type	surveytiming	area	units	need_catch_mult	fleetname
1	-1	1	1	0	FISHERY1	#_1
3	 1	1	2	0	SURVEY1 	#_2
#_Catch data
#_year	season	fleet	catch	catch_se
 -999	1	1	  13.4845	0.00999975	#_1         
    1	1	1	 161.8137	0.00999975	#_2         
    2	1	1	 455.3107	0.00999975	#_3         
    3	1	1	 734.3635	0.00999975	#_4         
    4	1	1	 945.0887	0.00999975	#_5         
    5	1	1	 713.4732	0.00999975	#_6         
    6	1	1	1219.1395	0.00999975	#_7         
    7	1	1	1152.5528	0.00999975	#_8         
    8	1	1	2211.9228	0.00999975	#_9         
    9	1	1	1184.4242	0.00999975	#_10        
   10	1	1	1368.8869	0.00999975	#_11        
   11	1	1	1451.1586	0.00999975	#_12        
   12	1	1	1448.3339	0.00999975	#_13        
   13	1	1	1006.7889	0.00999975	#_14        
   14	1	1	1378.3061	0.00999975	#_15        
   15	1	1	1374.6971	0.00999975	#_16        
   16	1	1	1165.1279	0.00999975	#_17        
   17	1	1	2051.1590	0.00999975	#_18        
   18	1	1	1468.7139	0.00999975	#_19        
   19	1	1	1305.7008	0.00999975	#_20        
   20	1	1	1238.6223	0.00999975	#_21        
   21	1	1	1642.3032	0.00999975	#_22        
   22	1	1	1189.3158	0.00999975	#_23        
   23	1	1	1900.6640	0.00999975	#_24        
   24	1	1	1416.3633	0.00999975	#_25        
   25	1	1	1292.9600	0.00999975	#_26        
   26	1	1	1111.6343	0.00999975	#_27        
   27	1	1	1047.9727	0.00999975	#_28        
   28	1	1	1365.9105	0.00999975	#_29        
   29	1	1	 963.4537	0.00999975	#_30        
   30	1	1	1415.8173	0.00999975	#_31        
-9999	0	0	   0.0000	0.00000000	#_terminator
#_CPUE_and_surveyabundance_observations
#_Units:  0=numbers; 1=biomass; 2=F; >=30 for special types
#_Errtype:  -1=normal; 0=lognormal; >0=T
#_SD_Report: 0=no sdreport; 1=enable sdreport
#_Fleet	Units	Errtype	SD_Report
1	1	0	1	#_FISHERY1
2	0	0	1	#_SURVEY1 
#
#_CPUE_data
#_year	seas	index	obs	se_log
    1	1	2	1.533600	0.198042	#_1         
    2	1	2	1.497468	0.198042	#_2         
    3	1	2	1.581523	0.198042	#_3         
    4	1	2	1.314229	0.198042	#_4         
    5	1	2	1.243375	0.198042	#_5         
    6	1	2	1.397078	0.198042	#_6         
    7	1	2	1.402665	0.198042	#_7         
    8	1	2	1.204772	0.198042	#_8         
    9	1	2	1.284566	0.198042	#_9         
   10	1	2	1.093091	0.198042	#_10        
   11	1	2	1.179107	0.198042	#_11        
   12	1	2	0.967635	0.198042	#_12        
   13	1	2	1.002185	0.198042	#_13        
   14	1	2	1.002597	0.198042	#_14        
   15	1	2	0.884867	0.198042	#_15        
   16	1	2	1.042649	0.198042	#_16        
   17	1	2	0.927591	0.198042	#_17        
   18	1	2	0.674842	0.198042	#_18        
   19	1	2	0.878411	0.198042	#_19        
   20	1	2	0.883119	0.198042	#_20        
   21	1	2	1.013313	0.198042	#_21        
   22	1	2	0.937926	0.198042	#_22        
   23	1	2	0.971210	0.198042	#_23        
   24	1	2	0.688562	0.198042	#_24        
   25	1	2	0.556869	0.198042	#_25        
   26	1	2	0.547224	0.198042	#_26        
   27	1	2	0.599862	0.198042	#_27        
   28	1	2	0.610073	0.198042	#_28        
   29	1	2	0.579984	0.198042	#_29        
   30	1	2	0.504580	0.198042	#_30        
-9999	0	0	0.000000	0.000000	#_terminator
0 #_N_discard_fleets
#_discard_units (1=same_as_catchunits(bio/num); 2=fraction; 3=numbers)
#_discard_errtype:  >0 for DF of T-dist(read CV below); 0 for normal with CV; -1 for normal with se; -2 for lognormal
#
#_discard_fleet_info
#
#_discard_data
#
#_meanbodywt
0 #_use_meanbodywt
 #_DF_for_meanbodywt_T-distribution_like
#
#_population_length_bins
2 # length bin method: 1=use databins; 2=generate from binwidth,min,max below; 3=read vector
1 # binwidth for population size comp
10 # minimum size in the population (lower edge of first bin and size at age 0.00)
89 # maximum size in the population (lower edge of last bin)
0 #_use_lencomp
12 #_N_agebins
#
#_agebin_vector
1 2 3 4 5 6 7 8 9 10 11 12 #_agebin_vector
#
#_ageing_error
1 #_N_ageerror_definitions
#_age0	age1	age2	age3	age4	age5	age6	age7	age8	age9	age10	age11	age12
-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	#_1
 0	 0	 0	 0	 0	 0	 0	 0	 0	 0	 0	 0	 0	#_2
#
#_age_info
#_mintailcomp	addtocomp	combine_M_F	CompressBins	CompError	ParmSelect	minsamplesize
0	1e-04	1	0	0	0	1	#_FISHERY1
0	1e-04	1	0	0	0	1	#_SURVEY1 
1 #_Lbin_method: 1=poplenbins; 2=datalenbins; 3=lengths
 #_combine males into females at or below this bin number
#_Yr	Seas	FltSvy	Gender	Part	Ageerr	Lbin_lo	Lbin_hi	Nsamp	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA	NA
    1	1	1	0	0	1	-1	-1	200	0.110	0.110	0.130	0.135	0.090	0.100	0.065	0.080	0.040	0.010	0.035	0.095	#_1         
    2	1	1	0	0	1	-1	-1	200	0.070	0.080	0.105	0.125	0.110	0.105	0.065	0.070	0.055	0.040	0.030	0.145	#_2         
    3	1	1	0	0	1	-1	-1	200	0.040	0.095	0.125	0.115	0.130	0.085	0.070	0.065	0.045	0.055	0.040	0.135	#_3         
    4	1	1	0	0	1	-1	-1	200	0.095	0.080	0.110	0.130	0.120	0.095	0.060	0.045	0.075	0.030	0.030	0.130	#_4         
    5	1	1	0	0	1	-1	-1	200	0.100	0.110	0.095	0.125	0.080	0.075	0.100	0.040	0.080	0.035	0.035	0.125	#_5         
    6	1	1	0	0	1	-1	-1	200	0.090	0.125	0.130	0.125	0.075	0.080	0.080	0.025	0.065	0.060	0.015	0.130	#_6         
    7	1	1	0	0	1	-1	-1	200	0.085	0.120	0.130	0.105	0.070	0.080	0.075	0.075	0.055	0.025	0.015	0.165	#_7         
    8	1	1	0	0	1	-1	-1	200	0.110	0.110	0.170	0.115	0.095	0.050	0.065	0.055	0.025	0.030	0.050	0.125	#_8         
    9	1	1	0	0	1	-1	-1	200	0.080	0.135	0.145	0.155	0.080	0.060	0.055	0.035	0.055	0.030	0.035	0.135	#_9         
   10	1	1	0	0	1	-1	-1	200	0.105	0.165	0.115	0.135	0.130	0.070	0.065	0.020	0.030	0.025	0.020	0.120	#_10        
   11	1	1	0	0	1	-1	-1	200	0.075	0.160	0.145	0.130	0.130	0.130	0.050	0.040	0.020	0.030	0.025	0.065	#_11        
   12	1	1	0	0	1	-1	-1	200	0.080	0.135	0.220	0.145	0.095	0.080	0.065	0.045	0.035	0.025	0.010	0.065	#_12        
   13	1	1	0	0	1	-1	-1	200	0.100	0.140	0.125	0.190	0.095	0.060	0.030	0.070	0.025	0.035	0.010	0.120	#_13        
   14	1	1	0	0	1	-1	-1	200	0.075	0.125	0.150	0.105	0.140	0.110	0.090	0.055	0.055	0.025	0.010	0.060	#_14        
   15	1	1	0	0	1	-1	-1	200	0.125	0.120	0.220	0.095	0.160	0.070	0.070	0.020	0.020	0.030	0.005	0.065	#_15        
   16	1	1	0	0	1	-1	-1	200	0.095	0.170	0.130	0.170	0.100	0.075	0.055	0.055	0.040	0.015	0.025	0.070	#_16        
   17	1	1	0	0	1	-1	-1	200	0.070	0.150	0.195	0.175	0.115	0.075	0.060	0.045	0.035	0.020	0.010	0.050	#_17        
   18	1	1	0	0	1	-1	-1	200	0.160	0.105	0.170	0.180	0.140	0.080	0.020	0.045	0.015	0.030	0.015	0.040	#_18        
   19	1	1	0	0	1	-1	-1	200	0.165	0.190	0.170	0.155	0.130	0.040	0.025	0.020	0.015	0.020	0.015	0.055	#_19        
   20	1	1	0	0	1	-1	-1	200	0.265	0.200	0.135	0.055	0.125	0.070	0.040	0.020	0.010	0.020	0.005	0.055	#_20        
   21	1	1	0	0	1	-1	-1	200	0.095	0.325	0.200	0.125	0.075	0.050	0.025	0.035	0.005	0.015	0.020	0.030	#_21        
   22	1	1	0	0	1	-1	-1	200	0.175	0.135	0.270	0.150	0.100	0.035	0.050	0.055	0.005	0.010	0.005	0.010	#_22        
   23	1	1	0	0	1	-1	-1	200	0.130	0.215	0.160	0.255	0.095	0.060	0.015	0.025	0.020	0.000	0.015	0.010	#_23        
   24	1	1	0	0	1	-1	-1	200	0.110	0.185	0.315	0.115	0.155	0.050	0.030	0.005	0.010	0.005	0.010	0.010	#_24        
   25	1	1	0	0	1	-1	-1	200	0.180	0.140	0.120	0.245	0.135	0.105	0.025	0.010	0.015	0.015	0.010	0.000	#_25        
   26	1	1	0	0	1	-1	-1	200	0.125	0.200	0.210	0.155	0.125	0.060	0.070	0.040	0.005	0.010	0.000	0.000	#_26        
   27	1	1	0	0	1	-1	-1	200	0.130	0.215	0.165	0.145	0.110	0.120	0.030	0.050	0.000	0.030	0.000	0.005	#_27        
   28	1	1	0	0	1	-1	-1	200	0.130	0.215	0.170	0.170	0.130	0.045	0.075	0.010	0.020	0.005	0.015	0.015	#_28        
   29	1	1	0	0	1	-1	-1	200	0.185	0.175	0.160	0.165	0.135	0.060	0.040	0.055	0.005	0.010	0.000	0.010	#_29        
   30	1	1	0	0	1	-1	-1	200	0.195	0.275	0.165	0.140	0.080	0.065	0.020	0.015	0.020	0.005	0.010	0.010	#_30        
    1	1	2	0	0	1	-1	-1	200	0.055	0.150	0.145	0.150	0.105	0.090	0.075	0.045	0.020	0.020	0.030	0.115	#_31        
    2	1	2	0	0	1	-1	-1	200	0.065	0.165	0.100	0.110	0.095	0.095	0.060	0.055	0.040	0.040	0.040	0.135	#_32        
    3	1	2	0	0	1	-1	-1	200	0.055	0.125	0.165	0.130	0.100	0.075	0.035	0.055	0.045	0.020	0.035	0.160	#_33        
    4	1	2	0	0	1	-1	-1	200	0.055	0.105	0.130	0.150	0.080	0.100	0.050	0.065	0.040	0.035	0.030	0.160	#_34        
    5	1	2	0	0	1	-1	-1	200	0.045	0.135	0.140	0.110	0.125	0.065	0.060	0.095	0.035	0.020	0.045	0.125	#_35        
    6	1	2	0	0	1	-1	-1	200	0.060	0.135	0.140	0.075	0.110	0.090	0.110	0.075	0.060	0.020	0.020	0.105	#_36        
    7	1	2	0	0	1	-1	-1	200	0.055	0.185	0.140	0.130	0.085	0.085	0.050	0.050	0.045	0.020	0.030	0.125	#_37        
    8	1	2	0	0	1	-1	-1	200	0.050	0.155	0.225	0.135	0.095	0.075	0.035	0.040	0.045	0.015	0.015	0.115	#_38        
    9	1	2	0	0	1	-1	-1	200	0.095	0.180	0.180	0.120	0.140	0.050	0.040	0.025	0.025	0.025	0.020	0.100	#_39        
   10	1	2	0	0	1	-1	-1	200	0.090	0.210	0.130	0.150	0.110	0.050	0.060	0.035	0.030	0.035	0.020	0.080	#_40        
   11	1	2	0	0	1	-1	-1	200	0.080	0.185	0.200	0.125	0.095	0.080	0.040	0.030	0.010	0.030	0.015	0.110	#_41        
   12	1	2	0	0	1	-1	-1	200	0.065	0.150	0.210	0.200	0.100	0.080	0.065	0.045	0.000	0.015	0.000	0.070	#_42        
   13	1	2	0	0	1	-1	-1	200	0.100	0.205	0.130	0.140	0.070	0.080	0.060	0.065	0.015	0.030	0.020	0.085	#_43        
   14	1	2	0	0	1	-1	-1	200	0.065	0.185	0.200	0.165	0.125	0.080	0.070	0.045	0.030	0.010	0.010	0.015	#_44        
   15	1	2	0	0	1	-1	-1	200	0.105	0.165	0.170	0.125	0.130	0.090	0.060	0.060	0.020	0.030	0.015	0.030	#_45        
   16	1	2	0	0	1	-1	-1	200	0.130	0.230	0.130	0.125	0.085	0.045	0.070	0.045	0.035	0.010	0.015	0.080	#_46        
   17	1	2	0	0	1	-1	-1	200	0.110	0.210	0.190	0.165	0.080	0.060	0.045	0.070	0.020	0.005	0.015	0.030	#_47        
   18	1	2	0	0	1	-1	-1	200	0.080	0.160	0.245	0.175	0.070	0.090	0.045	0.040	0.020	0.010	0.025	0.040	#_48        
   19	1	2	0	0	1	-1	-1	200	0.110	0.270	0.110	0.145	0.135	0.080	0.025	0.035	0.035	0.015	0.010	0.030	#_49        
   20	1	2	0	0	1	-1	-1	200	0.155	0.300	0.210	0.095	0.075	0.060	0.035	0.025	0.020	0.000	0.005	0.020	#_50        
   21	1	2	0	0	1	-1	-1	200	0.125	0.375	0.195	0.105	0.020	0.055	0.065	0.030	0.010	0.010	0.000	0.010	#_51        
   22	1	2	0	0	1	-1	-1	200	0.095	0.165	0.375	0.165	0.055	0.005	0.025	0.050	0.020	0.020	0.000	0.025	#_52        
   23	1	2	0	0	1	-1	-1	200	0.100	0.285	0.115	0.250	0.110	0.050	0.015	0.020	0.030	0.010	0.000	0.015	#_53        
   24	1	2	0	0	1	-1	-1	200	0.060	0.210	0.325	0.120	0.155	0.040	0.030	0.025	0.015	0.010	0.005	0.005	#_54        
   25	1	2	0	0	1	-1	-1	200	0.130	0.160	0.220	0.215	0.090	0.085	0.055	0.025	0.010	0.000	0.005	0.005	#_55        
   26	1	2	0	0	1	-1	-1	200	0.105	0.305	0.210	0.065	0.145	0.045	0.065	0.015	0.025	0.020	0.000	0.000	#_56        
   27	1	2	0	0	1	-1	-1	200	0.170	0.215	0.230	0.145	0.055	0.090	0.015	0.040	0.010	0.015	0.005	0.010	#_57        
   28	1	2	0	0	1	-1	-1	200	0.090	0.305	0.195	0.135	0.085	0.065	0.075	0.010	0.035	0.000	0.000	0.005	#_58        
   29	1	2	0	0	1	-1	-1	200	0.100	0.245	0.250	0.135	0.140	0.035	0.035	0.020	0.005	0.015	0.005	0.015	#_59        
   30	1	2	0	0	1	-1	-1	200	0.190	0.270	0.220	0.140	0.040	0.065	0.040	0.000	0.015	0.000	0.010	0.010	#_60        
-9999	0	0	0	0	0	 0	 0	  0	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	0.000	#_terminator
#
#_MeanSize_at_Age_obs
0 #_use_MeanSize_at_Age_obs
0 #_N_environ_variables
0 #_N_sizefreq_methods
0 #_do_tags
0 #_morphcomp_data
0 #_use_selectivity_priors
#
999