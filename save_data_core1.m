% 
% This script makes a data structure that contains all the sample mass
% depths and nuclide concentrations for core 1. Used for making plots and
% comparing model predictions to data.
%
% The samples structure contains 21 samples and include; (1) 17-OD1-C1-Surf
% (Ne-21 only), (2-7) 6 pit samples from 17-OD1-C1-PIT2-*, and (8-21) 14
% ice core samples 17-OD1-C1-*
% 
% Basically, one pastes data from the various spreadsheets into this script
% and runs it. That will store the data in a form that can be easily used
% for plotting and model-data comparison later. 
%
% Note this script includes shielding depth calculations for till and ice
% inside the script. That should make it easier to change these here
% instead of having to go back to the spreadsheet.
%
% So, if some of the observations get changed, paste in new data here and
% rerun to save revised data structure. 
%
% Greg Balco & 
% Marie Bergelin
%
% June 2019

clear all; close all;

%% Part I: data from till samples. 

% Enter all data for till samples. These data are all pasted in from
% various spreadsheets. 

% Sample depths
till_td_cm = [0  6 14 23 37 43 56];
till_bd_cm = [1 14 17 29 43 50 62];
% Density measurements 
till_measured_rho = [1.84 1.84 1.78 1.78 1.78 1.75 1.71];
% Nuclide concentrations
% Ne-21 data in Mat/g
% col 1 concentration in Mat/g, col 2 uncertainty
% There are multiple Ne-21 measurements for each sample, so column 3
% indexes measurements to samples
d21_till = [125.44	4.47 1  
    127.25	4.60 1
    141.39	4.70 1
    128.51	4.37 1
    113.47	2.11 2
    104.67	1.63 2
    105.31	1.77 3
    98.81	1.69 3	
    109.07	2.04 4
    100.82	1.71 4	
    88.62	1.90 5
    94.68	1.51 5	
    53.19	1.40 6
    57.12	1.47 6	
    36.33	1.16 7
    32.83	1.32 7];

% Al-26 data
% Only one measurement per sample, in at/g
d26_till = [NaN     NaN
    67842642	1209049
    67322924	1406815
    63433758	1333753
    58186829	1190198
    35585672	735994
    23321611	529695];

% Be-10 data 
% Only one measurement per sample 
d10_till = [NaN     NaN
    13430936	84219
    13512033	84981
    12866461	80386
    11592535	70144
    6638848     59357
    4246233     52229];

% Done entering data. Do shielding depth calculations. 

% Plot the density/depth for the till

figure; 

subplot(1,2,1); 
for a = 1:length(till_bd_cm)
    xx = [till_measured_rho(a) till_measured_rho(a)];
    yy = [till_td_cm(a) till_bd_cm(a)];
    plot(xx,yy,'r','linewidth',1); hold on;
end

set(gca,'ydir','reverse','ylim',[0 70]); grid on;
xlabel('Measured density g/cm3)'); ylabel('Depth (cm)');
title('Density measurements');

till_shielding_mass(1) = 0;
till_shielding_mass(2) = till_bd_cm(1).*till_measured_rho(1); % assume 0 to 14 cm has density of top sample
for a = 2:length(till_bd_cm);
    % then assume unsampled parts have density of sample below
    till_shielding_mass(a+1) = till_shielding_mass(a) + (till_bd_cm(a)-till_bd_cm(a-1)).*till_measured_rho(a);
end;

subplot(1,2,2);
plot(till_shielding_mass,[0 till_bd_cm],'r','linewidth',1); hold on;
plot(till_shielding_mass,[0 till_bd_cm],'ko','markerfacecolor','r');
xlabel('Cumulative mass depth (g/cm2)'); 
title('Calculated mass depth');
set(gca,'ydir','reverse','ylim',[0 70]); grid on;

% Now compute mass depths of till samples by interpolation

till_td_gcm2 = interp1([0 till_bd_cm],till_shielding_mass,till_td_cm);
till_bd_gcm2 = interp1([0 till_bd_cm],till_shielding_mass,till_bd_cm);
plot(till_td_gcm2,till_td_cm,'bo','markersize',8);
plot(till_bd_gcm2,till_bd_cm,'bo','markersize',8);

% Assemble all this stuff into a data structure
% This is a rather complicated data structure that needs to store a lot of
% information about the samples, some of which consist of multiple numbers
% per sample and not just one number. 

%sampledata = {}; % initialize empty cell array

for a = 1:length(till_td_cm)
    this_sample.td = till_td_cm(a);
    this_sample.bd = till_bd_cm(a);
    if a == 1
        this_sample.sample_name = '17-OD1-surf';  
    else
        this_sample.sample_name = ['17-OD1-PIT2-' int2str(this_sample.td) '-' int2str(this_sample.bd)];
    end
    this_sample.tdz = till_td_gcm2(a);
    this_sample.bdz = till_bd_gcm2(a);
    % Below is redundant, but will be needed for ice samples to retain info about how much
    % quartz was mixed to produce each sample
    this_sample.Mq = [1]; % Dummy value
    this_sample.Mqtd = this_sample.td;
    this_sample.Mqbd = this_sample.bd;
    this_sample.Mqtdz = this_sample.tdz;
    this_sample.Mqbdz = this_sample.bdz;
    % Attach Ne-21 data to sample
    these_N21s = find(d21_till(:,3) == a); % Find data corresponding to this sample
    this_sample.allN21s = 1e6.*d21_till(these_N21s,1)';
    this_sample.alldN21s = 1e6.*d21_till(these_N21s,2)';
    this_sample.N21 = mean(this_sample.allN21s); % Compute average of Ne-21 measurements 
    % Assign uncert to Ne-21 measurements. This will usually be the SD of
    % the measurements but can't be smaller than the best individual
    % precision. 
    this_sample.dN21 = max([std(this_sample.allN21s) min([this_sample.alldN21s])]); 
    % Attach Al-26 data
    this_sample.N26 = d26_till(a,1);
    this_sample.dN26 = d26_till(a,2);
    % Attach Be-10 data
    this_sample.N10 = d10_till(a,1);
    this_sample.dN10 = d10_till(a,2);
    
    % Place that data structure into the array of all samples
    sampledata{a} = this_sample;
end

% Result is a data structure containing all needed info about samples. We will now
% add the data for ice samples to this. 


%% Part II. Data from ice core samples.  

% 1. Core segments and sed concentration data. This includes data from all
% segments -- all segments needed to compute shielding mass. 
% From spreadsheet. Cols are top depth (cm), bot depth (cm), total wt (g),
% and sed wt (g). Zero in wt columns indicates assume clean ice density. 
d1 = [0	5	0	0
5	17	387.2	45.2
17	24	216.9	19.5
24	36	469     56.9
36	48	527.8	110.6
48	55	273.6	55
55	70	627.2	71.7
70	85	583.6	41.2
85	100	913.2	520
100	107	0	0
107	125	971.2	501.8
125	145	836.5	192
145	165	0	0
165	185	0	0
185	197	463.2	37.2
197	215	679.4	36.6
215	235	788.2	66.4
235	255	740.9	58.3
255	275	0	0
275	295	775.1	71.7
295	310	581.6	71.5
310	330	882.2	208.4
330	350	803.3	102.6
350	370	0	0
370	375	0	0
375	390	0	0
390	403	0	0
409	422	0	0
422	441	703.6	44.3
441	456	0	0
456	478	0	0
478	500	0	0
500	520	801.8	124.4
520	541	807.3	154.9
541	546	214     49.1
546	564	749.3	206.4
564	582	757.2	177.3
582	596	664.2	157.7
597	617	1084.3	585.9
617	626	359.6	87
626	649	930.9	232.8
649	672	758.8	8.8
672	692	0	0
692	710	0	0
710	721	0	0
721	739	0	0
730	756	0	0
756	781	0	0
781	800	804.7	241.3
800	819	703     163.8
819	839	860.1	210.9
839	859	768.2	126.1
859	879	754.9	130.3
879	899	807.8	152.1
899	911	403.1	79.9
911	930	782.8	225.8
930	944	651.1	161.9];

ice_td_cm = d1(:,1);
ice_bd_cm = d1(:,2);
ice_tot_wt = d1(:,3);
ice_sed_wt = d1(:,4);

% Now enter data about ice samples that actually have nuclide
% measurements.
% This is pasted in from spreadsheet and has columns of top depth (cm), bot depth
% (cm), mass of quartz mixed to make complete sample, and identifying number
% for the sample. 
d2 = [930	944	1.45	14
911	930	8.58	14
899	911	2.1     14
879	899	6.12	14
859	879	4.88	13
839	859	5.47	13
819	839	11.33	13
800	819	6.38	12
781	800	5.53	12
626	649	5.66	11
617	626	1.22	11
597	617	3.44	11
582	596	1.75	11
564	582	2.87	10
546	564	2.49	10
541	546	0.58	10
520	541	2.39	10
500	520	1.07	10
330	350	1.89	9
310	330	4.64	9
295	310	1.01	8
275	295	1.19	8
255 275 0       8
235	255	1.41	8
215	235	1.83	7
197	215	1.02	7
185	197	0.79	7
125	145	6.82	6
107	125	12.77	5
85	100	1.53	4
70	85	0.81	4
55	70	2       3
48	55	1.33	3
36	48	3.16	2
24	36	1.21	1
17	24	0.33	1
5	17	1.17	1];

% Enter data about Ne-21 concentrations
% col 1 concentration Mat/g, col 2 uncert Mat/g
% col 3 indexes to sample number
d21i = [12.51	1.31 1
10.38	0.90 1
12.82	1.17 2
11.63	1.08 2
12.49	1.36 3
16.49	1.33 3
10.83	1.16 3	
39.83	1.53 4
39.9	1.14 4	
33.8	1.13 5
44.06	1.44 5
34.13	1.15 5	
29.53	1.23 6
33.21	1.34 6
28.4	1.21 6	
13.19	0.91 7
17.16	1.03 7
11.64	0.93 7	
5.77	0.86 8
8.07	0.94 8
5.41	0.80 8	
5.49	1.01 9
15.8	1.11 9
6.3	0.98 9	
20.64	1.02 10
21.7	1.10 10
36.88	1.19 11
24.42	1.05 11 
23.72	1.07 11
49.68	1.48 12
55.53	1.42 12
61.65	1.61 13
67.44	1.70 13
65.33	1.74 14
68.72	1.63 14];

% Enter data about Al-26 concentrations
d26i = [4433605	209356
3438123	146814
2968918	147678
4224864	215432
3331110	86201
2828200	102965
1481148	124461
1173800	100649
830580	55058
450136	39671
409481	31748
489020	33816
489709	30519
470719	37952];

% Enter data about Be-10 concentrations
d10i = [708724	20399
572614	13655
572633	13291
1458947	29541
1096434	9897
871081	16570
288469	10125
160310	8999
116129	4817
147863	4089
132349	3216
418448	8103
515978	9701
481623	9124]; 

% Done entering data. Do shielding depth calculations. 

% Compute segment densities. 

rhoI = 0.917;
rhoD = 2.68;
ice_CD = ice_sed_wt./ice_tot_wt; % debris concentration
ice_CI = 1-ice_CD;
ice_rho = 1./((ice_CI./rhoI)+(ice_CD./rhoD)); % ice core segment density (g/cm3)

% That leaves NaN in cells with no sediment due to division by zero. 

% Insert ice density in NaN cells
ice_rho(find(isnan(ice_rho))) = rhoI.*ones(size(find(isnan(ice_rho))));

% Use same approach as above to compute shielding mass

ice_shielding_mass(1) = 0;
ice_shielding_mass(2) = ice_bd_cm(1).*ice_rho(1); 
for a = 2:length(ice_bd_cm)
    ice_shielding_mass(a+1) = ice_shielding_mass(a) + (ice_bd_cm(a)-ice_bd_cm(a-1)).*ice_rho(a);
end

% Make some plots

figure; subplot(1,3,1);
for a = 1:length(ice_bd_cm)
    xx = [ice_CD(a) ice_CD(a)];
    yy = [ice_td_cm(a) ice_bd_cm(a)];
    plot(xx,yy,'b','linewidth',1); hold on;
end

set(gca,'ydir','reverse','ylim',[0 1000]); grid on;
xlabel('Debris conc.'); ylabel('Depth (cm)');
title('CD'); 

subplot(1,3,2);
for a = 1:length(ice_bd_cm)
    xx = [ice_rho(a) ice_rho(a)];
    yy = [ice_td_cm(a) ice_bd_cm(a)];
    plot(xx,yy,'b','linewidth',1); hold on;
end

set(gca,'ydir','reverse','ylim',[0 1000]); grid on;
xlabel('Density g/cm3)'); ylabel('Depth (cm)');
title('Density');

subplot(1,3,3);
plot(ice_shielding_mass,[0 ice_bd_cm']','r','linewidth',1); hold on;
plot(ice_shielding_mass,[0 ice_bd_cm']','ko','markerfacecolor','r');
xlabel('Cumulative mass depth (g/cm2)'); 
title('Mass depth');
set(gca,'ydir','reverse','ylim',[0 1000]); grid on;

% Now we have a depth/shielding mass relationship for interpolating for the
% samples that were actually collected. 

% Again, load all that into the sample structure. 

% Get till linear and mass thickness in order to record all ice data in
% depth below the till surface, not the ice surface
max_till_cm = max(till_bd_cm);
max_till_gcm2 = max(till_bd_gcm2);

for a = 1:14 % there are 14 ice samples
    these_segments = find(d2(:,4) == a);
    % Obtain info about subsegments that were mixed to make complete
    % samples
    % Note converting to depth below till surface
    this_sample.Mq = d2(these_segments,3)';
    this_sample.Mqtd = max_till_cm + d2(these_segments,1)';
    this_sample.Mqbd = max_till_cm + d2(these_segments,2)';
    % Interpolate shielding depths
    this_sample.Mqtdz = max_till_gcm2 + interp1([0 ice_bd_cm']',ice_shielding_mass,d2(these_segments,1)');
    this_sample.Mqbdz = max_till_gcm2 + interp1([0 ice_bd_cm']',ice_shielding_mass,d2(these_segments,2)');
    
    % Also compute boundary values
    this_sample.td = min(this_sample.Mqtd);
    this_sample.bd = max(this_sample.Mqbd);
    this_sample.sample_name = ['17-OD1-C1-' int2str(this_sample.td-max_till_cm) '-' int2str(this_sample.bd-max_till_cm)];
    this_sample.tdz = min(this_sample.Mqtdz);
    this_sample.bdz = max(this_sample.Mqbdz);
     
    % Attach Ne-21 data to sample
    these_N21s = find(d21i(:,3) == a); % Find data corresponding to this sample
    this_sample.allN21s = 1e6.*d21i(these_N21s,1)';
    this_sample.alldN21s = 1e6.*d21i(these_N21s,2)';
    this_sample.N21 = mean(this_sample.allN21s); % Compute average of Ne-21 measurements 
    % Assign uncert to Ne-21 measurements. This will usually be the SD of
    % the measurements but can't be smaller than the best individual
    % precision. 
    this_sample.dN21 = max([std(this_sample.allN21s) min([this_sample.alldN21s])]); 
    % Attach Al-26 data
    this_sample.N26 = d26i(a,1);
    this_sample.dN26 = d26i(a,2);
    % Attach Be-10 data
    this_sample.N10 = d10i(a,1);
    this_sample.dN10 = d10i(a,2);
    
    % Place that data structure into the array of all samples
    sampledata{a+7} = this_sample;
end


% Finally, save the whole data array

save data_core1 sampledata
