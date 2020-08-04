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
 -999	1	1	  13.396	0.00999975	#_1         
    1	1	1	 160.752	0.00999975	#_2         
    2	1	1	 461.985	0.00999975	#_3         
    3	1	1	 739.680	0.00999975	#_4         
    4	1	1	 965.994	0.00999975	#_5         
    5	1	1	 741.728	0.00999975	#_6         
    6	1	1	1282.092	0.00999975	#_7         
    7	1	1	1259.762	0.00999975	#_8         
    8	1	1	2501.417	0.00999975	#_9         
    9	1	1	1373.312	0.00999975	#_10        
   10	1	1	1651.241	0.00999975	#_11        
   11	1	1	1871.605	0.00999975	#_12        
   12	1	1	1958.734	0.00999975	#_13        
   13	1	1	1380.217	0.00999975	#_14        
   14	1	1	1938.411	0.00999975	#_15        
   15	1	1	1941.969	0.00999975	#_16        
   16	1	1	1700.161	0.00999975	#_17        
   17	1	1	3005.168	0.00999975	#_18        
   18	1	1	2094.897	0.00999975	#_19        
   19	1	1	1834.388	0.00999975	#_20        
   20	1	1	1612.335	0.00999975	#_21        
   21	1	1	1983.413	0.00999975	#_22        
   22	1	1	1271.172	0.00999975	#_23        
   23	1	1	1734.158	0.00999975	#_24        
   24	1	1	1175.152	0.00999975	#_25        
   25	1	1	1108.087	0.00999975	#_26        
   26	1	1	1140.609	0.00999975	#_27        
   27	1	1	1287.785	0.00999975	#_28        
   28	1	1	1859.528	0.00999975	#_29        
   29	1	1	1348.085	0.00999975	#_30        
   30	1	1	1905.028	0.00999975	#_31        
-9999	0	0	   0.000	0.00000000	#_terminator
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
    1	1	2	1.343482	0.198042	#_1         
    2	1	2	1.275962	0.198042	#_2         
    3	1	2	1.213569	0.198042	#_3         
    4	1	2	1.124475	0.198042	#_4         
    5	1	2	1.231031	0.198042	#_5         
    6	1	2	1.534875	0.198042	#_6         
    7	1	2	1.364297	0.198042	#_7         
    8	1	2	1.456722	0.198042	#_8         
    9	1	2	1.285481	0.198042	#_9         
   10	1	2	1.349032	0.198042	#_10        
   11	1	2	1.465800	0.198042	#_11        
   12	1	2	1.250601	0.198042	#_12        
   13	1	2	1.046086	0.198042	#_13        
   14	1	2	1.308991	0.198042	#_14        
   15	1	2	1.200348	0.198042	#_15        
   16	1	2	1.227878	0.198042	#_16        
   17	1	2	1.053231	0.198042	#_17        
   18	1	2	0.849542	0.198042	#_18        
   19	1	2	0.929783	0.198042	#_19        
   20	1	2	0.780017	0.198042	#_20        
   21	1	2	0.692622	0.198042	#_21        
   22	1	2	0.546329	0.198042	#_22        
   23	1	2	0.477687	0.198042	#_23        
   24	1	2	0.468415	0.198042	#_24        
   25	1	2	0.516000	0.198042	#_25        
   26	1	2	0.801486	0.198042	#_26        
   27	1	2	0.747796	0.198042	#_27        
   28	1	2	0.869546	0.198042	#_28        
   29	1	2	0.650470	0.198042	#_29        
   30	1	2	0.693296	0.198042	#_30        
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
    1	1	1	0	0	1	-1	-1	200	0.050	0.110	0.100	0.130	0.090	0.080	0.080	0.070	0.045	0.035	0.025	0.185	#_1         
    2	1	1	0	0	1	-1	-1	200	0.080	0.100	0.145	0.090	0.090	0.100	0.060	0.060	0.045	0.060	0.020	0.150	#_2         
    3	1	1	0	0	1	-1	-1	200	0.075	0.130	0.085	0.115	0.095	0.060	0.055	0.055	0.050	0.040	0.025	0.215	#_3         
    4	1	1	0	0	1	-1	-1	200	0.080	0.090	0.100	0.110	0.155	0.050	0.095	0.065	0.050	0.030	0.040	0.135	#_4         
    5	1	1	0	0	1	-1	-1	200	0.140	0.085	0.100	0.135	0.090	0.115	0.055	0.070	0.025	0.030	0.040	0.115	#_5         
    6	1	1	0	0	1	-1	-1	200	0.080	0.195	0.140	0.095	0.080	0.090	0.060	0.065	0.040	0.035	0.020	0.100	#_6         
    7	1	1	0	0	1	-1	-1	200	0.075	0.105	0.205	0.110	0.080	0.075	0.055	0.075	0.040	0.030	0.025	0.125	#_7         
    8	1	1	0	0	1	-1	-1	200	0.135	0.070	0.135	0.205	0.075	0.085	0.050	0.045	0.055	0.030	0.020	0.095	#_8         
    9	1	1	0	0	1	-1	-1	200	0.070	0.190	0.120	0.130	0.145	0.055	0.025	0.045	0.035	0.020	0.030	0.135	#_9         
   10	1	1	0	0	1	-1	-1	200	0.165	0.115	0.215	0.075	0.080	0.135	0.045	0.070	0.010	0.015	0.010	0.065	#_10        
   11	1	1	0	0	1	-1	-1	200	0.145	0.235	0.090	0.195	0.065	0.050	0.070	0.035	0.020	0.025	0.005	0.065	#_11        
   12	1	1	0	0	1	-1	-1	200	0.045	0.140	0.275	0.125	0.165	0.035	0.060	0.060	0.025	0.010	0.005	0.055	#_12        
   13	1	1	0	0	1	-1	-1	200	0.025	0.065	0.170	0.295	0.105	0.130	0.060	0.020	0.050	0.015	0.015	0.050	#_13        
   14	1	1	0	0	1	-1	-1	200	0.165	0.080	0.085	0.155	0.205	0.080	0.125	0.020	0.030	0.015	0.010	0.030	#_14        
   15	1	1	0	0	1	-1	-1	200	0.155	0.240	0.065	0.065	0.125	0.130	0.045	0.065	0.030	0.035	0.010	0.035	#_15        
   16	1	1	0	0	1	-1	-1	200	0.070	0.140	0.340	0.050	0.075	0.090	0.080	0.035	0.050	0.020	0.025	0.025	#_16        
   17	1	1	0	0	1	-1	-1	200	0.080	0.120	0.190	0.275	0.040	0.040	0.045	0.045	0.030	0.060	0.000	0.075	#_17        
   18	1	1	0	0	1	-1	-1	200	0.065	0.110	0.150	0.190	0.235	0.065	0.035	0.040	0.030	0.005	0.015	0.060	#_18        
   19	1	1	0	0	1	-1	-1	200	0.190	0.070	0.060	0.140	0.145	0.210	0.020	0.040	0.030	0.035	0.015	0.045	#_19        
   20	1	1	0	0	1	-1	-1	200	0.085	0.320	0.095	0.120	0.100	0.090	0.085	0.005	0.010	0.020	0.025	0.045	#_20        
   21	1	1	0	0	1	-1	-1	200	0.060	0.075	0.375	0.105	0.080	0.060	0.085	0.065	0.005	0.035	0.010	0.045	#_21        
   22	1	1	0	0	1	-1	-1	200	0.060	0.100	0.125	0.400	0.060	0.020	0.055	0.060	0.040	0.010	0.020	0.050	#_22        
   23	1	1	0	0	1	-1	-1	200	0.175	0.125	0.105	0.145	0.230	0.050	0.050	0.015	0.020	0.040	0.005	0.040	#_23        
   24	1	1	0	0	1	-1	-1	200	0.270	0.195	0.085	0.100	0.080	0.120	0.040	0.030	0.010	0.010	0.045	0.015	#_24        
   25	1	1	0	0	1	-1	-1	200	0.280	0.295	0.155	0.030	0.050	0.050	0.070	0.010	0.020	0.005	0.005	0.030	#_25        
   26	1	1	0	0	1	-1	-1	200	0.270	0.350	0.175	0.075	0.035	0.040	0.010	0.035	0.000	0.000	0.000	0.010	#_26        
   27	1	1	0	0	1	-1	-1	200	0.105	0.225	0.325	0.185	0.055	0.020	0.025	0.020	0.015	0.005	0.000	0.020	#_27        
   28	1	1	0	0	1	-1	-1	200	0.085	0.175	0.260	0.255	0.125	0.045	0.010	0.015	0.005	0.005	0.015	0.005	#_28        
   29	1	1	0	0	1	-1	-1	200	0.195	0.090	0.165	0.200	0.225	0.070	0.040	0.000	0.000	0.005	0.010	0.000	#_29        
   30	1	1	0	0	1	-1	-1	200	0.240	0.215	0.105	0.175	0.120	0.085	0.025	0.020	0.000	0.010	0.000	0.005	#_30        
    1	1	2	0	0	1	-1	-1	200	0.045	0.125	0.135	0.175	0.105	0.080	0.060	0.040	0.060	0.030	0.020	0.125	#_31        
    2	1	2	0	0	1	-1	-1	200	0.040	0.120	0.145	0.175	0.075	0.080	0.055	0.085	0.025	0.020	0.050	0.130	#_32        
    3	1	2	0	0	1	-1	-1	200	0.045	0.180	0.170	0.070	0.115	0.065	0.070	0.065	0.025	0.040	0.035	0.120	#_33        
    4	1	2	0	0	1	-1	-1	200	0.065	0.110	0.175	0.120	0.115	0.095	0.050	0.060	0.045	0.020	0.050	0.095	#_34        
    5	1	2	0	0	1	-1	-1	200	0.140	0.135	0.115	0.080	0.090	0.115	0.065	0.030	0.045	0.040	0.020	0.125	#_35        
    6	1	2	0	0	1	-1	-1	200	0.070	0.290	0.095	0.090	0.085	0.035	0.090	0.050	0.040	0.035	0.015	0.105	#_36        
    7	1	2	0	0	1	-1	-1	200	0.065	0.115	0.240	0.135	0.040	0.080	0.070	0.040	0.050	0.015	0.025	0.125	#_37        
    8	1	2	0	0	1	-1	-1	200	0.145	0.100	0.145	0.210	0.090	0.065	0.045	0.035	0.010	0.035	0.010	0.110	#_38        
    9	1	2	0	0	1	-1	-1	200	0.045	0.275	0.120	0.115	0.170	0.055	0.055	0.035	0.020	0.005	0.015	0.090	#_39        
   10	1	2	0	0	1	-1	-1	200	0.115	0.165	0.270	0.105	0.090	0.105	0.020	0.025	0.020	0.040	0.005	0.040	#_40        
   11	1	2	0	0	1	-1	-1	200	0.085	0.270	0.120	0.175	0.070	0.055	0.090	0.030	0.015	0.030	0.005	0.055	#_41        
   12	1	2	0	0	1	-1	-1	200	0.040	0.200	0.265	0.120	0.155	0.045	0.045	0.060	0.020	0.000	0.005	0.045	#_42        
   13	1	2	0	0	1	-1	-1	200	0.035	0.105	0.240	0.260	0.065	0.135	0.030	0.025	0.045	0.005	0.025	0.030	#_43        
   14	1	2	0	0	1	-1	-1	200	0.170	0.045	0.120	0.155	0.225	0.080	0.100	0.030	0.025	0.015	0.025	0.010	#_44        
   15	1	2	0	0	1	-1	-1	200	0.100	0.315	0.060	0.090	0.125	0.105	0.045	0.075	0.005	0.010	0.025	0.045	#_45        
   16	1	2	0	0	1	-1	-1	200	0.075	0.210	0.345	0.050	0.050	0.085	0.040	0.015	0.055	0.020	0.005	0.050	#_46        
   17	1	2	0	0	1	-1	-1	200	0.040	0.160	0.235	0.260	0.045	0.035	0.060	0.070	0.020	0.030	0.010	0.035	#_47        
   18	1	2	0	0	1	-1	-1	200	0.055	0.115	0.225	0.175	0.250	0.035	0.030	0.030	0.045	0.000	0.020	0.020	#_48        
   19	1	2	0	0	1	-1	-1	200	0.175	0.150	0.110	0.125	0.155	0.145	0.015	0.015	0.030	0.025	0.005	0.050	#_49        
   20	1	2	0	0	1	-1	-1	200	0.085	0.340	0.150	0.090	0.090	0.055	0.070	0.005	0.025	0.025	0.010	0.055	#_50        
   21	1	2	0	0	1	-1	-1	200	0.050	0.140	0.365	0.095	0.055	0.055	0.060	0.090	0.005	0.005	0.030	0.050	#_51        
   22	1	2	0	0	1	-1	-1	200	0.040	0.120	0.220	0.380	0.075	0.035	0.025	0.035	0.050	0.005	0.005	0.010	#_52        
   23	1	2	0	0	1	-1	-1	200	0.120	0.165	0.150	0.160	0.220	0.030	0.040	0.025	0.025	0.020	0.005	0.040	#_53        
   24	1	2	0	0	1	-1	-1	200	0.185	0.235	0.110	0.090	0.100	0.145	0.045	0.025	0.010	0.030	0.005	0.020	#_54        
   25	1	2	0	0	1	-1	-1	200	0.190	0.385	0.150	0.020	0.035	0.040	0.100	0.015	0.010	0.015	0.010	0.030	#_55        
   26	1	2	0	0	1	-1	-1	200	0.145	0.370	0.255	0.105	0.045	0.020	0.020	0.025	0.000	0.000	0.010	0.005	#_56        
   27	1	2	0	0	1	-1	-1	200	0.085	0.265	0.365	0.170	0.040	0.025	0.015	0.010	0.010	0.000	0.005	0.010	#_57        
   28	1	2	0	0	1	-1	-1	200	0.045	0.195	0.315	0.280	0.080	0.025	0.010	0.010	0.010	0.020	0.000	0.010	#_58        
   29	1	2	0	0	1	-1	-1	200	0.095	0.195	0.215	0.185	0.175	0.095	0.010	0.005	0.000	0.000	0.015	0.010	#_59        
   30	1	2	0	0	1	-1	-1	200	0.170	0.180	0.105	0.180	0.170	0.140	0.035	0.020	0.000	0.000	0.000	0.000	#_60        
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