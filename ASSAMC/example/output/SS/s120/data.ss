#V3.30
#C data file created using the SS_writedat function in the R package r4ss
#C should work with SS version: 
#C file write time: 2020-07-06 04:45:22
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
 -999	1	1	  13.449	0.00999975	#_1         
    1	1	1	 161.388	0.00999975	#_2         
    2	1	1	 460.885	0.00999975	#_3         
    3	1	1	 742.927	0.00999975	#_4         
    4	1	1	 980.186	0.00999975	#_5         
    5	1	1	 750.690	0.00999975	#_6         
    6	1	1	1297.308	0.00999975	#_7         
    7	1	1	1255.604	0.00999975	#_8         
    8	1	1	2431.867	0.00999975	#_9         
    9	1	1	1311.076	0.00999975	#_10        
   10	1	1	1530.302	0.00999975	#_11        
   11	1	1	1660.228	0.00999975	#_12        
   12	1	1	1656.866	0.00999975	#_13        
   13	1	1	1154.514	0.00999975	#_14        
   14	1	1	1592.911	0.00999975	#_15        
   15	1	1	1559.450	0.00999975	#_16        
   16	1	1	1327.291	0.00999975	#_17        
   17	1	1	2283.802	0.00999975	#_18        
   18	1	1	1598.587	0.00999975	#_19        
   19	1	1	1453.773	0.00999975	#_20        
   20	1	1	1360.367	0.00999975	#_21        
   21	1	1	1764.607	0.00999975	#_22        
   22	1	1	1194.992	0.00999975	#_23        
   23	1	1	1797.011	0.00999975	#_24        
   24	1	1	1290.578	0.00999975	#_25        
   25	1	1	1213.533	0.00999975	#_26        
   26	1	1	1068.165	0.00999975	#_27        
   27	1	1	1057.278	0.00999975	#_28        
   28	1	1	1413.713	0.00999975	#_29        
   29	1	1	1004.185	0.00999975	#_30        
   30	1	1	1409.660	0.00999975	#_31        
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
    1	1	2	1.430643	0.198042	#_1         
    2	1	2	1.339368	0.198042	#_2         
    3	1	2	1.478027	0.198042	#_3         
    4	1	2	1.486962	0.198042	#_4         
    5	1	2	1.503408	0.198042	#_5         
    6	1	2	1.406733	0.198042	#_6         
    7	1	2	1.386303	0.198042	#_7         
    8	1	2	1.262149	0.198042	#_8         
    9	1	2	1.357947	0.198042	#_9         
   10	1	2	1.158495	0.198042	#_10        
   11	1	2	1.236882	0.198042	#_11        
   12	1	2	1.065573	0.198042	#_12        
   13	1	2	1.202458	0.198042	#_13        
   14	1	2	1.168169	0.198042	#_14        
   15	1	2	0.947807	0.198042	#_15        
   16	1	2	0.991822	0.198042	#_16        
   17	1	2	0.752666	0.198042	#_17        
   18	1	2	0.856353	0.198042	#_18        
   19	1	2	0.850782	0.198042	#_19        
   20	1	2	0.879586	0.198042	#_20        
   21	1	2	0.768490	0.198042	#_21        
   22	1	2	0.711353	0.198042	#_22        
   23	1	2	0.794922	0.198042	#_23        
   24	1	2	0.599055	0.198042	#_24        
   25	1	2	0.524814	0.198042	#_25        
   26	1	2	0.679104	0.198042	#_26        
   27	1	2	0.632764	0.198042	#_27        
   28	1	2	0.574382	0.198042	#_28        
   29	1	2	0.447716	0.198042	#_29        
   30	1	2	0.581471	0.198042	#_30        
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
    1	1	1	0	0	1	-1	-1	200	0.055	0.110	0.115	0.120	0.075	0.095	0.070	0.090	0.035	0.045	0.045	0.145	#_1         
    2	1	1	0	0	1	-1	-1	200	0.075	0.130	0.150	0.110	0.120	0.085	0.065	0.060	0.020	0.055	0.050	0.080	#_2         
    3	1	1	0	0	1	-1	-1	200	0.085	0.120	0.095	0.115	0.100	0.090	0.060	0.065	0.035	0.035	0.010	0.190	#_3         
    4	1	1	0	0	1	-1	-1	200	0.075	0.125	0.100	0.115	0.065	0.105	0.100	0.065	0.045	0.030	0.025	0.150	#_4         
    5	1	1	0	0	1	-1	-1	200	0.060	0.105	0.150	0.135	0.120	0.060	0.110	0.055	0.010	0.015	0.030	0.150	#_5         
    6	1	1	0	0	1	-1	-1	200	0.045	0.170	0.115	0.120	0.075	0.105	0.065	0.065	0.035	0.055	0.050	0.100	#_6         
    7	1	1	0	0	1	-1	-1	200	0.090	0.115	0.150	0.140	0.085	0.060	0.060	0.035	0.060	0.040	0.035	0.130	#_7         
    8	1	1	0	0	1	-1	-1	200	0.060	0.155	0.085	0.185	0.150	0.065	0.065	0.030	0.050	0.035	0.015	0.105	#_8         
    9	1	1	0	0	1	-1	-1	200	0.060	0.125	0.195	0.140	0.160	0.085	0.060	0.060	0.025	0.015	0.015	0.060	#_9         
   10	1	1	0	0	1	-1	-1	200	0.090	0.140	0.145	0.230	0.075	0.120	0.045	0.025	0.020	0.015	0.015	0.080	#_10        
   11	1	1	0	0	1	-1	-1	200	0.065	0.155	0.165	0.115	0.115	0.060	0.060	0.095	0.035	0.030	0.020	0.085	#_11        
   12	1	1	0	0	1	-1	-1	200	0.070	0.150	0.180	0.195	0.090	0.095	0.035	0.055	0.015	0.030	0.025	0.060	#_12        
   13	1	1	0	0	1	-1	-1	200	0.065	0.185	0.125	0.150	0.160	0.090	0.060	0.045	0.010	0.035	0.025	0.050	#_13        
   14	1	1	0	0	1	-1	-1	200	0.115	0.085	0.175	0.115	0.085	0.110	0.065	0.075	0.030	0.040	0.020	0.085	#_14        
   15	1	1	0	0	1	-1	-1	200	0.095	0.150	0.145	0.165	0.090	0.085	0.070	0.035	0.020	0.035	0.030	0.080	#_15        
   16	1	1	0	0	1	-1	-1	200	0.065	0.115	0.220	0.135	0.140	0.090	0.085	0.040	0.005	0.015	0.015	0.075	#_16        
   17	1	1	0	0	1	-1	-1	200	0.130	0.125	0.135	0.230	0.070	0.100	0.060	0.035	0.015	0.025	0.035	0.040	#_17        
   18	1	1	0	0	1	-1	-1	200	0.220	0.130	0.125	0.080	0.145	0.040	0.050	0.055	0.025	0.040	0.030	0.060	#_18        
   19	1	1	0	0	1	-1	-1	200	0.150	0.250	0.175	0.110	0.095	0.065	0.015	0.040	0.025	0.025	0.025	0.025	#_19        
   20	1	1	0	0	1	-1	-1	200	0.160	0.165	0.220	0.110	0.095	0.075	0.045	0.045	0.030	0.015	0.010	0.030	#_20        
   21	1	1	0	0	1	-1	-1	200	0.115	0.150	0.200	0.210	0.080	0.080	0.040	0.025	0.005	0.035	0.010	0.050	#_21        
   22	1	1	0	0	1	-1	-1	200	0.060	0.155	0.285	0.150	0.115	0.050	0.050	0.035	0.030	0.005	0.010	0.055	#_22        
   23	1	1	0	0	1	-1	-1	200	0.205	0.120	0.130	0.155	0.120	0.150	0.025	0.030	0.005	0.020	0.005	0.035	#_23        
   24	1	1	0	0	1	-1	-1	200	0.100	0.295	0.150	0.140	0.115	0.095	0.045	0.020	0.015	0.010	0.010	0.005	#_24        
   25	1	1	0	0	1	-1	-1	200	0.205	0.200	0.260	0.120	0.055	0.040	0.045	0.050	0.010	0.005	0.000	0.010	#_25        
   26	1	1	0	0	1	-1	-1	200	0.220	0.260	0.120	0.180	0.060	0.035	0.050	0.035	0.025	0.005	0.005	0.005	#_26        
   27	1	1	0	0	1	-1	-1	200	0.165	0.205	0.200	0.105	0.145	0.085	0.040	0.020	0.010	0.015	0.005	0.005	#_27        
   28	1	1	0	0	1	-1	-1	200	0.155	0.120	0.305	0.210	0.055	0.070	0.035	0.015	0.015	0.015	0.005	0.000	#_28        
   29	1	1	0	0	1	-1	-1	200	0.255	0.140	0.180	0.160	0.120	0.055	0.050	0.015	0.005	0.005	0.005	0.010	#_29        
   30	1	1	0	0	1	-1	-1	200	0.235	0.200	0.205	0.100	0.090	0.070	0.005	0.030	0.025	0.005	0.020	0.015	#_30        
    1	1	2	0	0	1	-1	-1	200	0.050	0.140	0.145	0.145	0.140	0.065	0.075	0.030	0.035	0.050	0.010	0.115	#_31        
    2	1	2	0	0	1	-1	-1	200	0.040	0.145	0.150	0.110	0.140	0.055	0.065	0.045	0.045	0.045	0.035	0.125	#_32        
    3	1	2	0	0	1	-1	-1	200	0.060	0.105	0.155	0.115	0.105	0.050	0.055	0.055	0.060	0.055	0.045	0.140	#_33        
    4	1	2	0	0	1	-1	-1	200	0.100	0.140	0.145	0.155	0.105	0.085	0.055	0.055	0.050	0.035	0.010	0.065	#_34        
    5	1	2	0	0	1	-1	-1	200	0.040	0.165	0.110	0.115	0.130	0.105	0.045	0.050	0.065	0.035	0.025	0.115	#_35        
    6	1	2	0	0	1	-1	-1	200	0.070	0.195	0.120	0.165	0.070	0.080	0.020	0.065	0.040	0.025	0.020	0.130	#_36        
    7	1	2	0	0	1	-1	-1	200	0.055	0.105	0.220	0.135	0.090	0.080	0.050	0.050	0.050	0.020	0.020	0.125	#_37        
    8	1	2	0	0	1	-1	-1	200	0.090	0.210	0.125	0.180	0.085	0.050	0.085	0.015	0.010	0.015	0.030	0.105	#_38        
    9	1	2	0	0	1	-1	-1	200	0.115	0.140	0.210	0.100	0.105	0.060	0.025	0.070	0.030	0.025	0.010	0.110	#_39        
   10	1	2	0	0	1	-1	-1	200	0.085	0.250	0.130	0.155	0.080	0.095	0.070	0.020	0.040	0.005	0.010	0.060	#_40        
   11	1	2	0	0	1	-1	-1	200	0.080	0.165	0.230	0.095	0.100	0.065	0.065	0.050	0.020	0.020	0.005	0.105	#_41        
   12	1	2	0	0	1	-1	-1	200	0.105	0.165	0.165	0.145	0.080	0.075	0.065	0.050	0.030	0.030	0.005	0.085	#_42        
   13	1	2	0	0	1	-1	-1	200	0.060	0.175	0.145	0.155	0.140	0.115	0.075	0.035	0.040	0.005	0.010	0.045	#_43        
   14	1	2	0	0	1	-1	-1	200	0.100	0.150	0.245	0.105	0.090	0.100	0.055	0.045	0.040	0.035	0.010	0.025	#_44        
   15	1	2	0	0	1	-1	-1	200	0.075	0.240	0.130	0.145	0.095	0.085	0.030	0.050	0.025	0.015	0.020	0.090	#_45        
   16	1	2	0	0	1	-1	-1	200	0.080	0.175	0.235	0.105	0.135	0.040	0.025	0.040	0.020	0.055	0.025	0.065	#_46        
   17	1	2	0	0	1	-1	-1	200	0.080	0.190	0.170	0.180	0.100	0.075	0.040	0.040	0.040	0.025	0.010	0.050	#_47        
   18	1	2	0	0	1	-1	-1	200	0.100	0.195	0.155	0.115	0.160	0.045	0.085	0.045	0.030	0.025	0.010	0.035	#_48        
   19	1	2	0	0	1	-1	-1	200	0.085	0.280	0.205	0.115	0.080	0.070	0.045	0.035	0.020	0.015	0.005	0.045	#_49        
   20	1	2	0	0	1	-1	-1	200	0.145	0.190	0.250	0.120	0.075	0.065	0.060	0.025	0.000	0.030	0.010	0.030	#_50        
   21	1	2	0	0	1	-1	-1	200	0.065	0.245	0.185	0.265	0.065	0.055	0.035	0.015	0.020	0.005	0.005	0.040	#_51        
   22	1	2	0	0	1	-1	-1	200	0.060	0.180	0.305	0.170	0.140	0.055	0.010	0.015	0.030	0.010	0.005	0.020	#_52        
   23	1	2	0	0	1	-1	-1	200	0.170	0.180	0.165	0.135	0.120	0.085	0.025	0.045	0.005	0.025	0.000	0.045	#_53        
   24	1	2	0	0	1	-1	-1	200	0.060	0.385	0.130	0.095	0.135	0.075	0.065	0.015	0.005	0.010	0.010	0.015	#_54        
   25	1	2	0	0	1	-1	-1	200	0.155	0.170	0.275	0.120	0.085	0.075	0.050	0.020	0.020	0.015	0.000	0.015	#_55        
   26	1	2	0	0	1	-1	-1	200	0.125	0.245	0.155	0.200	0.075	0.045	0.055	0.045	0.025	0.015	0.000	0.015	#_56        
   27	1	2	0	0	1	-1	-1	200	0.135	0.230	0.240	0.100	0.130	0.025	0.010	0.035	0.035	0.030	0.010	0.020	#_57        
   28	1	2	0	0	1	-1	-1	200	0.115	0.215	0.305	0.170	0.035	0.085	0.035	0.010	0.005	0.010	0.000	0.015	#_58        
   29	1	2	0	0	1	-1	-1	200	0.145	0.160	0.290	0.215	0.065	0.030	0.035	0.030	0.015	0.010	0.000	0.005	#_59        
   30	1	2	0	0	1	-1	-1	200	0.160	0.225	0.165	0.155	0.120	0.085	0.025	0.030	0.000	0.005	0.015	0.015	#_60        
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