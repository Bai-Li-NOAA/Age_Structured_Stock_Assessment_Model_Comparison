#V3.30
#C data file created using the SS_writedat function in the R package r4ss
#C should work with SS version: 
#C file write time: 2020-07-06 04:46:07
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
 -999	1	1	  13.3943	0.00999975	#_1         
    1	1	1	 160.7322	0.00999975	#_2         
    2	1	1	 458.2642	0.00999975	#_3         
    3	1	1	 743.2294	0.00999975	#_4         
    4	1	1	 979.4426	0.00999975	#_5         
    5	1	1	 755.8320	0.00999975	#_6         
    6	1	1	1298.7312	0.00999975	#_7         
    7	1	1	1255.6917	0.00999975	#_8         
    8	1	1	2371.8825	0.00999975	#_9         
    9	1	1	1267.9927	0.00999975	#_10        
   10	1	1	1433.0343	0.00999975	#_11        
   11	1	1	1522.0577	0.00999975	#_12        
   12	1	1	1520.2622	0.00999975	#_13        
   13	1	1	1059.5892	0.00999975	#_14        
   14	1	1	1504.3064	0.00999975	#_15        
   15	1	1	1501.7005	0.00999975	#_16        
   16	1	1	1288.1012	0.00999975	#_17        
   17	1	1	2205.0231	0.00999975	#_18        
   18	1	1	1569.7601	0.00999975	#_19        
   19	1	1	1445.0034	0.00999975	#_20        
   20	1	1	1343.9131	0.00999975	#_21        
   21	1	1	1753.6427	0.00999975	#_22        
   22	1	1	1183.9326	0.00999975	#_23        
   23	1	1	1789.9331	0.00999975	#_24        
   24	1	1	1342.4862	0.00999975	#_25        
   25	1	1	1256.4893	0.00999975	#_26        
   26	1	1	1123.0219	0.00999975	#_27        
   27	1	1	1093.9634	0.00999975	#_28        
   28	1	1	1486.2939	0.00999975	#_29        
   29	1	1	1096.7237	0.00999975	#_30        
   30	1	1	1558.6088	0.00999975	#_31        
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
    1	1	2	1.450842	0.198042	#_1         
    2	1	2	1.635264	0.198042	#_2         
    3	1	2	1.316766	0.198042	#_3         
    4	1	2	1.475766	0.198042	#_4         
    5	1	2	1.536987	0.198042	#_5         
    6	1	2	1.377532	0.198042	#_6         
    7	1	2	1.505112	0.198042	#_7         
    8	1	2	1.374335	0.198042	#_8         
    9	1	2	1.247270	0.198042	#_9         
   10	1	2	1.302136	0.198042	#_10        
   11	1	2	1.010440	0.198042	#_11        
   12	1	2	1.208930	0.198042	#_12        
   13	1	2	1.064622	0.198042	#_13        
   14	1	2	0.970496	0.198042	#_14        
   15	1	2	1.034129	0.198042	#_15        
   16	1	2	0.863593	0.198042	#_16        
   17	1	2	0.836619	0.198042	#_17        
   18	1	2	0.859967	0.198042	#_18        
   19	1	2	0.817933	0.198042	#_19        
   20	1	2	0.711109	0.198042	#_20        
   21	1	2	0.785499	0.198042	#_21        
   22	1	2	0.640200	0.198042	#_22        
   23	1	2	0.737225	0.198042	#_23        
   24	1	2	0.696385	0.198042	#_24        
   25	1	2	0.591424	0.198042	#_25        
   26	1	2	0.669840	0.198042	#_26        
   27	1	2	0.682252	0.198042	#_27        
   28	1	2	0.627510	0.198042	#_28        
   29	1	2	0.665892	0.198042	#_29        
   30	1	2	0.640031	0.198042	#_30        
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
    1	1	1	0	0	1	-1	-1	200	0.070	0.100	0.120	0.110	0.120	0.085	0.055	0.055	0.060	0.030	0.030	0.165	#_1         
    2	1	1	0	0	1	-1	-1	200	0.105	0.115	0.110	0.105	0.125	0.045	0.080	0.060	0.035	0.045	0.040	0.135	#_2         
    3	1	1	0	0	1	-1	-1	200	0.080	0.115	0.150	0.155	0.065	0.085	0.060	0.065	0.050	0.040	0.030	0.105	#_3         
    4	1	1	0	0	1	-1	-1	200	0.055	0.170	0.105	0.105	0.060	0.110	0.095	0.080	0.045	0.040	0.025	0.110	#_4         
    5	1	1	0	0	1	-1	-1	200	0.075	0.110	0.195	0.110	0.120	0.055	0.050	0.065	0.020	0.040	0.015	0.145	#_5         
    6	1	1	0	0	1	-1	-1	200	0.095	0.090	0.120	0.125	0.135	0.060	0.050	0.050	0.060	0.050	0.020	0.145	#_6         
    7	1	1	0	0	1	-1	-1	200	0.060	0.135	0.125	0.130	0.165	0.080	0.070	0.030	0.030	0.040	0.020	0.115	#_7         
    8	1	1	0	0	1	-1	-1	200	0.065	0.110	0.180	0.145	0.090	0.095	0.075	0.065	0.045	0.035	0.005	0.090	#_8         
    9	1	1	0	0	1	-1	-1	200	0.120	0.090	0.125	0.140	0.105	0.060	0.085	0.060	0.030	0.035	0.015	0.135	#_9         
   10	1	1	0	0	1	-1	-1	200	0.090	0.140	0.090	0.140	0.110	0.100	0.080	0.050	0.035	0.030	0.020	0.115	#_10        
   11	1	1	0	0	1	-1	-1	200	0.090	0.115	0.210	0.100	0.120	0.085	0.085	0.065	0.025	0.025	0.010	0.070	#_11        
   12	1	1	0	0	1	-1	-1	200	0.120	0.165	0.150	0.170	0.100	0.075	0.055	0.025	0.015	0.045	0.020	0.060	#_12        
   13	1	1	0	0	1	-1	-1	200	0.100	0.140	0.105	0.170	0.140	0.070	0.070	0.055	0.035	0.020	0.015	0.080	#_13        
   14	1	1	0	0	1	-1	-1	200	0.120	0.125	0.210	0.175	0.060	0.095	0.030	0.040	0.035	0.025	0.005	0.080	#_14        
   15	1	1	0	0	1	-1	-1	200	0.060	0.100	0.200	0.185	0.110	0.080	0.080	0.035	0.020	0.020	0.010	0.100	#_15        
   16	1	1	0	0	1	-1	-1	200	0.110	0.080	0.145	0.155	0.160	0.095	0.065	0.070	0.015	0.005	0.015	0.085	#_16        
   17	1	1	0	0	1	-1	-1	200	0.100	0.145	0.160	0.130	0.115	0.110	0.080	0.050	0.045	0.010	0.025	0.030	#_17        
   18	1	1	0	0	1	-1	-1	200	0.095	0.115	0.210	0.100	0.090	0.135	0.090	0.045	0.055	0.015	0.010	0.040	#_18        
   19	1	1	0	0	1	-1	-1	200	0.195	0.195	0.165	0.125	0.060	0.045	0.045	0.070	0.035	0.020	0.015	0.030	#_19        
   20	1	1	0	0	1	-1	-1	200	0.115	0.230	0.195	0.155	0.070	0.045	0.025	0.055	0.030	0.045	0.005	0.030	#_20        
   21	1	1	0	0	1	-1	-1	200	0.130	0.135	0.280	0.165	0.120	0.035	0.030	0.030	0.020	0.015	0.015	0.025	#_21        
   22	1	1	0	0	1	-1	-1	200	0.115	0.165	0.220	0.185	0.100	0.065	0.070	0.020	0.015	0.010	0.020	0.015	#_22        
   23	1	1	0	0	1	-1	-1	200	0.120	0.200	0.205	0.135	0.145	0.075	0.050	0.010	0.030	0.005	0.015	0.010	#_23        
   24	1	1	0	0	1	-1	-1	200	0.170	0.210	0.225	0.090	0.070	0.120	0.040	0.015	0.030	0.010	0.015	0.005	#_24        
   25	1	1	0	0	1	-1	-1	200	0.215	0.160	0.190	0.175	0.095	0.050	0.040	0.040	0.015	0.000	0.005	0.015	#_25        
   26	1	1	0	0	1	-1	-1	200	0.145	0.220	0.230	0.140	0.125	0.030	0.045	0.015	0.015	0.020	0.005	0.010	#_26        
   27	1	1	0	0	1	-1	-1	200	0.180	0.215	0.205	0.160	0.100	0.055	0.015	0.025	0.025	0.010	0.010	0.000	#_27        
   28	1	1	0	0	1	-1	-1	200	0.245	0.260	0.175	0.085	0.095	0.055	0.030	0.025	0.015	0.005	0.000	0.010	#_28        
   29	1	1	0	0	1	-1	-1	200	0.095	0.220	0.230	0.160	0.125	0.075	0.050	0.035	0.000	0.000	0.005	0.005	#_29        
   30	1	1	0	0	1	-1	-1	200	0.145	0.185	0.245	0.135	0.115	0.090	0.020	0.025	0.025	0.000	0.000	0.015	#_30        
    1	1	2	0	0	1	-1	-1	200	0.055	0.110	0.190	0.100	0.105	0.100	0.070	0.055	0.050	0.035	0.030	0.100	#_31        
    2	1	2	0	0	1	-1	-1	200	0.035	0.140	0.140	0.105	0.110	0.070	0.075	0.045	0.055	0.040	0.025	0.160	#_32        
    3	1	2	0	0	1	-1	-1	200	0.050	0.135	0.185	0.175	0.065	0.055	0.080	0.030	0.045	0.045	0.050	0.085	#_33        
    4	1	2	0	0	1	-1	-1	200	0.075	0.185	0.135	0.125	0.100	0.070	0.060	0.045	0.035	0.020	0.025	0.125	#_34        
    5	1	2	0	0	1	-1	-1	200	0.040	0.165	0.180	0.130	0.100	0.065	0.060	0.035	0.040	0.030	0.020	0.135	#_35        
    6	1	2	0	0	1	-1	-1	200	0.065	0.130	0.090	0.195	0.125	0.075	0.050	0.055	0.030	0.025	0.035	0.125	#_36        
    7	1	2	0	0	1	-1	-1	200	0.075	0.165	0.140	0.160	0.115	0.085	0.055	0.030	0.010	0.035	0.030	0.100	#_37        
    8	1	2	0	0	1	-1	-1	200	0.050	0.195	0.165	0.110	0.095	0.090	0.040	0.040	0.030	0.055	0.025	0.105	#_38        
    9	1	2	0	0	1	-1	-1	200	0.115	0.135	0.175	0.135	0.120	0.070	0.035	0.025	0.045	0.035	0.025	0.085	#_39        
   10	1	2	0	0	1	-1	-1	200	0.105	0.185	0.110	0.130	0.095	0.100	0.060	0.055	0.030	0.030	0.015	0.085	#_40        
   11	1	2	0	0	1	-1	-1	200	0.075	0.150	0.155	0.135	0.095	0.090	0.045	0.070	0.050	0.015	0.045	0.075	#_41        
   12	1	2	0	0	1	-1	-1	200	0.105	0.250	0.135	0.155	0.045	0.095	0.050	0.045	0.025	0.040	0.005	0.050	#_42        
   13	1	2	0	0	1	-1	-1	200	0.085	0.205	0.200	0.100	0.150	0.060	0.055	0.030	0.010	0.015	0.030	0.060	#_43        
   14	1	2	0	0	1	-1	-1	200	0.050	0.135	0.250	0.175	0.065	0.105	0.035	0.040	0.040	0.020	0.015	0.070	#_44        
   15	1	2	0	0	1	-1	-1	200	0.080	0.195	0.185	0.165	0.125	0.115	0.055	0.010	0.035	0.015	0.015	0.005	#_45        
   16	1	2	0	0	1	-1	-1	200	0.060	0.145	0.130	0.185	0.150	0.095	0.070	0.060	0.025	0.025	0.020	0.035	#_46        
   17	1	2	0	0	1	-1	-1	200	0.105	0.240	0.175	0.105	0.110	0.090	0.045	0.025	0.025	0.010	0.020	0.050	#_47        
   18	1	2	0	0	1	-1	-1	200	0.135	0.190	0.155	0.100	0.100	0.070	0.090	0.040	0.040	0.025	0.010	0.045	#_48        
   19	1	2	0	0	1	-1	-1	200	0.110	0.230	0.185	0.105	0.100	0.055	0.070	0.035	0.020	0.020	0.020	0.050	#_49        
   20	1	2	0	0	1	-1	-1	200	0.110	0.250	0.270	0.070	0.120	0.050	0.030	0.040	0.015	0.005	0.000	0.040	#_50        
   21	1	2	0	0	1	-1	-1	200	0.115	0.170	0.210	0.170	0.090	0.070	0.055	0.025	0.020	0.020	0.025	0.030	#_51        
   22	1	2	0	0	1	-1	-1	200	0.165	0.200	0.210	0.235	0.085	0.025	0.030	0.020	0.000	0.000	0.000	0.030	#_52        
   23	1	2	0	0	1	-1	-1	200	0.085	0.355	0.175	0.085	0.105	0.080	0.055	0.020	0.005	0.005	0.005	0.025	#_53        
   24	1	2	0	0	1	-1	-1	200	0.110	0.195	0.255	0.145	0.075	0.110	0.050	0.025	0.005	0.005	0.000	0.025	#_54        
   25	1	2	0	0	1	-1	-1	200	0.145	0.250	0.205	0.200	0.065	0.045	0.035	0.030	0.010	0.005	0.000	0.010	#_55        
   26	1	2	0	0	1	-1	-1	200	0.145	0.300	0.210	0.130	0.090	0.060	0.020	0.020	0.000	0.005	0.000	0.020	#_56        
   27	1	2	0	0	1	-1	-1	200	0.145	0.215	0.275	0.105	0.120	0.075	0.010	0.010	0.020	0.015	0.000	0.010	#_57        
   28	1	2	0	0	1	-1	-1	200	0.150	0.305	0.165	0.170	0.060	0.050	0.070	0.005	0.010	0.010	0.000	0.005	#_58        
   29	1	2	0	0	1	-1	-1	200	0.125	0.280	0.230	0.115	0.130	0.045	0.030	0.025	0.010	0.000	0.000	0.010	#_59        
   30	1	2	0	0	1	-1	-1	200	0.135	0.180	0.250	0.200	0.085	0.070	0.025	0.010	0.020	0.010	0.010	0.005	#_60        
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