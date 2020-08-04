#V3.30
#C data file created using the SS_writedat function in the R package r4ss
#C should work with SS version: 
#C file write time: 2020-08-03 11:38:18
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
 -999	1	1	  13.3882	0.00999975	#_1         
    1	1	1	 160.6584	0.00999975	#_2         
    2	1	1	 461.2435	0.00999975	#_3         
    3	1	1	 730.9171	0.00999975	#_4         
    4	1	1	 943.2490	0.00999975	#_5         
    5	1	1	 711.6620	0.00999975	#_6         
    6	1	1	1183.4065	0.00999975	#_7         
    7	1	1	1100.1904	0.00999975	#_8         
    8	1	1	2108.4049	0.00999975	#_9         
    9	1	1	1172.5962	0.00999975	#_10        
   10	1	1	1421.3725	0.00999975	#_11        
   11	1	1	1623.8609	0.00999975	#_12        
   12	1	1	1722.3066	0.00999975	#_13        
   13	1	1	1208.9530	0.00999975	#_14        
   14	1	1	1644.0723	0.00999975	#_15        
   15	1	1	1594.8659	0.00999975	#_16        
   16	1	1	1291.9587	0.00999975	#_17        
   17	1	1	2191.1194	0.00999975	#_18        
   18	1	1	1539.9805	0.00999975	#_19        
   19	1	1	1464.1935	0.00999975	#_20        
   20	1	1	1466.7284	0.00999975	#_21        
   21	1	1	2040.5802	0.00999975	#_22        
   22	1	1	1474.8939	0.00999975	#_23        
   23	1	1	2288.2160	0.00999975	#_24        
   24	1	1	1653.1598	0.00999975	#_25        
   25	1	1	1515.0324	0.00999975	#_26        
   26	1	1	1285.4833	0.00999975	#_27        
   27	1	1	1206.4290	0.00999975	#_28        
   28	1	1	1495.9414	0.00999975	#_29        
   29	1	1	 984.9598	0.00999975	#_30        
   30	1	1	1305.3877	0.00999975	#_31        
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
    1	1	2	1.263278	0.198042	#_1         
    2	1	2	1.235778	0.198042	#_2         
    3	1	2	1.277244	0.198042	#_3         
    4	1	2	1.339255	0.198042	#_4         
    5	1	2	1.321848	0.198042	#_5         
    6	1	2	1.160333	0.198042	#_6         
    7	1	2	1.405184	0.198042	#_7         
    8	1	2	1.068466	0.198042	#_8         
    9	1	2	1.351757	0.198042	#_9         
   10	1	2	1.304212	0.198042	#_10        
   11	1	2	1.576196	0.198042	#_11        
   12	1	2	1.335251	0.198042	#_12        
   13	1	2	1.174220	0.198042	#_13        
   14	1	2	1.072660	0.198042	#_14        
   15	1	2	0.935982	0.198042	#_15        
   16	1	2	0.875536	0.198042	#_16        
   17	1	2	0.781321	0.198042	#_17        
   18	1	2	0.939022	0.198042	#_18        
   19	1	2	0.952310	0.198042	#_19        
   20	1	2	0.952634	0.198042	#_20        
   21	1	2	0.940908	0.198042	#_21        
   22	1	2	0.925919	0.198042	#_22        
   23	1	2	0.930977	0.198042	#_23        
   24	1	2	0.756671	0.198042	#_24        
   25	1	2	0.724313	0.198042	#_25        
   26	1	2	0.562644	0.198042	#_26        
   27	1	2	0.730824	0.198042	#_27        
   28	1	2	0.517171	0.198042	#_28        
   29	1	2	0.426661	0.198042	#_29        
   30	1	2	0.478141	0.198042	#_30        
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
    1	1	1	0	0	1	-1	-1	200	0.080	0.100	0.165	0.115	0.095	0.070	0.055	0.045	0.060	0.020	0.025	0.170	#_1         
    2	1	1	0	0	1	-1	-1	200	0.070	0.140	0.115	0.125	0.135	0.085	0.085	0.045	0.025	0.040	0.030	0.105	#_2         
    3	1	1	0	0	1	-1	-1	200	0.035	0.095	0.105	0.115	0.120	0.090	0.095	0.060	0.065	0.045	0.030	0.145	#_3         
    4	1	1	0	0	1	-1	-1	200	0.055	0.085	0.105	0.140	0.105	0.060	0.090	0.060	0.055	0.030	0.050	0.165	#_4         
    5	1	1	0	0	1	-1	-1	200	0.045	0.065	0.095	0.135	0.135	0.080	0.095	0.070	0.065	0.020	0.025	0.170	#_5         
    6	1	1	0	0	1	-1	-1	200	0.145	0.050	0.105	0.115	0.135	0.070	0.055	0.060	0.050	0.015	0.035	0.165	#_6         
    7	1	1	0	0	1	-1	-1	200	0.140	0.165	0.060	0.065	0.070	0.085	0.065	0.090	0.055	0.040	0.025	0.140	#_7         
    8	1	1	0	0	1	-1	-1	200	0.055	0.260	0.165	0.040	0.080	0.060	0.070	0.060	0.040	0.050	0.025	0.095	#_8         
    9	1	1	0	0	1	-1	-1	200	0.210	0.130	0.170	0.120	0.045	0.065	0.050	0.025	0.035	0.020	0.025	0.105	#_9         
   10	1	1	0	0	1	-1	-1	200	0.140	0.235	0.095	0.110	0.145	0.040	0.020	0.020	0.025	0.020	0.025	0.125	#_10        
   11	1	1	0	0	1	-1	-1	200	0.090	0.210	0.225	0.090	0.090	0.075	0.020	0.035	0.015	0.030	0.025	0.095	#_11        
   12	1	1	0	0	1	-1	-1	200	0.040	0.065	0.265	0.285	0.045	0.115	0.050	0.010	0.030	0.005	0.010	0.080	#_12        
   13	1	1	0	0	1	-1	-1	200	0.075	0.055	0.085	0.270	0.220	0.040	0.110	0.055	0.015	0.005	0.015	0.055	#_13        
   14	1	1	0	0	1	-1	-1	200	0.125	0.105	0.080	0.065	0.180	0.240	0.030	0.070	0.025	0.000	0.015	0.065	#_14        
   15	1	1	0	0	1	-1	-1	200	0.045	0.165	0.100	0.100	0.115	0.165	0.120	0.050	0.055	0.020	0.010	0.055	#_15        
   16	1	1	0	0	1	-1	-1	200	0.090	0.110	0.165	0.145	0.055	0.050	0.125	0.125	0.040	0.025	0.025	0.045	#_16        
   17	1	1	0	0	1	-1	-1	200	0.105	0.130	0.135	0.180	0.085	0.050	0.035	0.070	0.090	0.015	0.045	0.060	#_17        
   18	1	1	0	0	1	-1	-1	200	0.220	0.160	0.150	0.090	0.110	0.045	0.050	0.050	0.065	0.025	0.015	0.020	#_18        
   19	1	1	0	0	1	-1	-1	200	0.140	0.300	0.155	0.125	0.050	0.100	0.035	0.015	0.025	0.010	0.015	0.030	#_19        
   20	1	1	0	0	1	-1	-1	200	0.225	0.210	0.220	0.135	0.065	0.015	0.040	0.030	0.010	0.000	0.020	0.030	#_20        
   21	1	1	0	0	1	-1	-1	200	0.135	0.230	0.210	0.160	0.110	0.040	0.040	0.035	0.010	0.000	0.005	0.025	#_21        
   22	1	1	0	0	1	-1	-1	200	0.040	0.235	0.250	0.175	0.125	0.045	0.045	0.015	0.035	0.010	0.005	0.020	#_22        
   23	1	1	0	0	1	-1	-1	200	0.150	0.085	0.265	0.245	0.105	0.085	0.015	0.015	0.005	0.010	0.005	0.015	#_23        
   24	1	1	0	0	1	-1	-1	200	0.080	0.305	0.065	0.150	0.135	0.120	0.090	0.020	0.010	0.005	0.015	0.005	#_24        
   25	1	1	0	0	1	-1	-1	200	0.135	0.190	0.275	0.015	0.155	0.095	0.040	0.015	0.025	0.005	0.015	0.035	#_25        
   26	1	1	0	0	1	-1	-1	200	0.145	0.135	0.160	0.265	0.080	0.050	0.065	0.020	0.035	0.015	0.005	0.025	#_26        
   27	1	1	0	0	1	-1	-1	200	0.095	0.265	0.190	0.140	0.125	0.030	0.065	0.030	0.020	0.015	0.015	0.010	#_27        
   28	1	1	0	0	1	-1	-1	200	0.075	0.120	0.315	0.180	0.080	0.105	0.045	0.040	0.010	0.015	0.010	0.005	#_28        
   29	1	1	0	0	1	-1	-1	200	0.085	0.195	0.125	0.240	0.160	0.075	0.045	0.005	0.020	0.030	0.005	0.015	#_29        
   30	1	1	0	0	1	-1	-1	200	0.405	0.065	0.180	0.050	0.145	0.045	0.015	0.050	0.005	0.005	0.010	0.025	#_30        
    1	1	2	0	0	1	-1	-1	200	0.080	0.105	0.170	0.100	0.100	0.060	0.045	0.075	0.065	0.060	0.020	0.120	#_31        
    2	1	2	0	0	1	-1	-1	200	0.040	0.145	0.135	0.100	0.145	0.105	0.070	0.060	0.025	0.040	0.035	0.100	#_32        
    3	1	2	0	0	1	-1	-1	200	0.055	0.145	0.135	0.105	0.120	0.065	0.065	0.060	0.045	0.020	0.025	0.160	#_33        
    4	1	2	0	0	1	-1	-1	200	0.025	0.090	0.100	0.175	0.130	0.060	0.065	0.055	0.055	0.030	0.025	0.190	#_34        
    5	1	2	0	0	1	-1	-1	200	0.025	0.085	0.100	0.100	0.160	0.090	0.085	0.070	0.050	0.040	0.040	0.155	#_35        
    6	1	2	0	0	1	-1	-1	200	0.095	0.080	0.110	0.110	0.100	0.085	0.040	0.045	0.060	0.055	0.040	0.180	#_36        
    7	1	2	0	0	1	-1	-1	200	0.125	0.285	0.020	0.070	0.080	0.095	0.070	0.060	0.040	0.035	0.015	0.105	#_37        
    8	1	2	0	0	1	-1	-1	200	0.055	0.275	0.185	0.060	0.060	0.045	0.080	0.065	0.035	0.035	0.015	0.090	#_38        
    9	1	2	0	0	1	-1	-1	200	0.140	0.125	0.270	0.185	0.030	0.040	0.020	0.045	0.045	0.015	0.010	0.075	#_39        
   10	1	2	0	0	1	-1	-1	200	0.130	0.345	0.075	0.155	0.105	0.020	0.030	0.005	0.035	0.010	0.005	0.085	#_40        
   11	1	2	0	0	1	-1	-1	200	0.060	0.215	0.295	0.075	0.145	0.065	0.020	0.030	0.010	0.025	0.005	0.055	#_41        
   12	1	2	0	0	1	-1	-1	200	0.040	0.100	0.270	0.215	0.060	0.115	0.065	0.020	0.015	0.015	0.015	0.070	#_42        
   13	1	2	0	0	1	-1	-1	200	0.055	0.080	0.090	0.280	0.255	0.030	0.070	0.065	0.005	0.010	0.005	0.055	#_43        
   14	1	2	0	0	1	-1	-1	200	0.115	0.105	0.105	0.060	0.275	0.165	0.000	0.035	0.050	0.015	0.015	0.060	#_44        
   15	1	2	0	0	1	-1	-1	200	0.065	0.215	0.130	0.070	0.075	0.170	0.090	0.040	0.045	0.035	0.005	0.060	#_45        
   16	1	2	0	0	1	-1	-1	200	0.055	0.140	0.265	0.100	0.070	0.090	0.070	0.095	0.010	0.025	0.035	0.045	#_46        
   17	1	2	0	0	1	-1	-1	200	0.140	0.175	0.145	0.185	0.060	0.045	0.040	0.065	0.075	0.015	0.015	0.040	#_47        
   18	1	2	0	0	1	-1	-1	200	0.155	0.250	0.145	0.090	0.115	0.090	0.010	0.020	0.040	0.040	0.005	0.040	#_48        
   19	1	2	0	0	1	-1	-1	200	0.100	0.345	0.210	0.090	0.045	0.060	0.020	0.015	0.015	0.025	0.025	0.050	#_49        
   20	1	2	0	0	1	-1	-1	200	0.140	0.265	0.265	0.110	0.075	0.035	0.040	0.010	0.005	0.005	0.030	0.020	#_50        
   21	1	2	0	0	1	-1	-1	200	0.125	0.265	0.210	0.160	0.090	0.040	0.015	0.035	0.015	0.000	0.000	0.045	#_51        
   22	1	2	0	0	1	-1	-1	200	0.040	0.235	0.280	0.135	0.145	0.065	0.055	0.010	0.010	0.000	0.005	0.020	#_52        
   23	1	2	0	0	1	-1	-1	200	0.165	0.065	0.255	0.255	0.090	0.080	0.035	0.015	0.005	0.010	0.005	0.020	#_53        
   24	1	2	0	0	1	-1	-1	200	0.080	0.325	0.105	0.130	0.200	0.055	0.070	0.015	0.005	0.000	0.000	0.015	#_54        
   25	1	2	0	0	1	-1	-1	200	0.120	0.175	0.390	0.040	0.085	0.075	0.040	0.025	0.020	0.015	0.000	0.015	#_55        
   26	1	2	0	0	1	-1	-1	200	0.085	0.245	0.205	0.225	0.035	0.065	0.050	0.035	0.030	0.015	0.000	0.010	#_56        
   27	1	2	0	0	1	-1	-1	200	0.075	0.265	0.200	0.145	0.140	0.030	0.055	0.025	0.035	0.010	0.005	0.015	#_57        
   28	1	2	0	0	1	-1	-1	200	0.115	0.135	0.310	0.175	0.080	0.085	0.010	0.025	0.030	0.015	0.005	0.015	#_58        
   29	1	2	0	0	1	-1	-1	200	0.100	0.230	0.110	0.285	0.095	0.050	0.075	0.010	0.015	0.005	0.000	0.025	#_59        
   30	1	2	0	0	1	-1	-1	200	0.275	0.125	0.215	0.095	0.135	0.040	0.035	0.040	0.015	0.005	0.000	0.020	#_60        
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