%% MC_optimizer

% This script runs a monte carlo simulation for the forward model fitting
% error to the data set through optimization. The optimizer accepts a set
% of given initial data parameters x0 and objectively finds a new set of x
% values that satisfy a minimum best fit value based on the
% objective_sublimation_model.m function and given constraints.
%
% Before running this script all variable in the 'Set/define all model
% parameters' should adjusted to fit the data set and desired model
% constraints.
% 
% To run 1 Monte Carlo iteration takes ~2-3 min. If increasing the
% MC_iterations to greater than 4, replace the 'for' loop to a 'parfor'
% loop to improve running time. 
% 
% This script can be send to a supercomputer for faster results. To do this
% one will need to have set up a parallel cluster for a supercomputer. The
% script can the be in the comman window by below batch statement, where
% the pool number is the number of parallel workers running the loop
% simultaniously
% 
%   j = batch(c,'MC_optimizer','Pool',100,'CurrentFolder','.','AutoAddClientPath',false)
% 
%
% Marie Bergelin
%
% Dec 2020

clear all; close all;

%% Set/define all model parameters and dataset

% Load data to be used for analysis. If data files haven't been created,
% run preprocessor scripts to create them. 

if ~exist('data_core1.mat','file')
    save_data_core1
end

if ~exist('Pall.mat','file')
    precalculate_P_of_z
end

clear all;

load data_core1
load Pall.mat P

% Define measured till thickness (cm and g/cm2)
ztill.d = sampledata{7}.bd; % depth of the till
ztill.dz = sampledata{7}.bdz; % mass depth of the till

% Define samples to consider in model fitting
mask = [1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 0 0 0 0 0]; % Surf + Pit 2 + low ice
%mask = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 1 1 1 0 0 0 0 0];
%mask = [1 0 0 0 0 0 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0];
%mask = [1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0]; 

% Define nuclides to consider in model fitting
% nmask = [Be-21 Ne-21 Al-26]
nmask = [1 1 1]; % Define all three

% nonlcon burial mask constraint (=1 use, = 0 don't use)
burialmask = [0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1 1 1]; % all samples showing burial
%burialmask = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 0 0 0 1 1 1 1 1]; % all samples showing burial

% Number of iterations
numits = 2;
%numits = 10000;

% initial starting values for the optimizer
% X0 = [sub(m/My),erosion(m/My),Age(My),10-inh(1e5),21-inh(1e6),26-inh(1e5)]
X0 = [22.7 0.89 1.83 1.4 11.4 8.2]; % Ted Bibby's results inh is from youngest pits average conc. and age thereafter

% Upper/Lower Boundaries
minx = [0 0 0 0 0 0];
maxx = [100 2 10 20 15 50];


%%%%%% Minimum exposure settings (Lowest ice inh / zero erosion) %%%%%%%

% % X0 = [sub(m/My),erosion(m/My),Age(My),10-inh(1e5),21-inh(1e6),26-inh(1e5)]
% X0 = [22.7 0 1.83 1.2 6.4 4.1];
% 
% % Upper/Lower Boundaries
% minx = [0 0 0 1.2 6.4 4.1];
% maxx = [100 0 10 1.2 6.4 4.1];


%% Run Monte Carlo Simulation
tic

% This runs a monte carlo simulation for the error on the data

% set optimizer settings
opts = optimset('fmincon');
opts = optimset(opts,'TolFun',1e-3,'display','iter','tolX',1e-3,'MaxFunEvals',300,'MaxIter',300);

% set up empty array to save time
MCval = zeros(numits,1);
s = zeros(numits,1);
e = zeros(numits,1);
T = zeros(numits,1);
N10inh = zeros(numits,1);
N21inh = zeros(numits,1);
N26inh = zeros(numits,1);
flag = zeros(numits,1);

Cdebris = zeros(numits,1);
z_ice_loss = zeros(numits,1);


plotFlag = 0; % request plotting

% correct for Nucleogenic Ne-21 from burialmask samples

nN21 = 7e6; % Nucleogenic Ne-21 measured in Ong Valley
dnN21 = 3e6; % uncertainty Nucleogenic Ne-21

ncorr_sampledata = sampledata;

for a = 1:length(sampledata)
        
    if burialmask(a) >= 1
        % subtract Nucleogenic Ne-21 from burialmask samples only
        ncorr_sampledata{a}.N21 = sampledata{a}.N21 - nN21;
        ncorr_sampledata{a}.dN21 = sqrt((sampledata{a}.dN21)^2 + dnN21^2);
    else
    end
end


for a = 1:numits
    disp([ num2str(a),' iterations out of ',num2str(numits)])

    % Create a new random data sample set for optimizer
    randomsample = ncorr_sampledata;

    for i = 1:length(ncorr_sampledata)

        % Generate random normal distributed nuclide concentrations
        % This requires the statistics toolbox
        % randomsample{i}.N10 = normrnd(ncorr_sampledata{i}.N10,ncorr_sampledata{i}.dN10);
        % randomsample{i}.N21 = normrnd(ncorr_sampledata{i}.N21,ncorr_sampledata{i}.dN21);  
        % randomsample{i}.N26 = normrnd(ncorr_sampledata{i}.N26,ncorr_sampledata{i}.dN26);
        
        % This should do the same thing without the statistics toolbox
        randomsample{i}.N10 = ncorr_sampledata{i}.N10 + randn(1).*ncorr_sampledata{i}.dN10;
        randomsample{i}.N21 = ncorr_sampledata{i}.N21 + randn(1).*ncorr_sampledata{i}.dN21;
        randomsample{i}.N26 = ncorr_sampledata{i}.N26 + randn(1).*ncorr_sampledata{i}.dN26;
        
        
    end

    % Run optimiser
    %[optx,fval,exitflag] = fmincon(@(X) objective_sublimation_model(X,randomsample,ztill,P,mask,nmask,chi,plotFlag),X0,[],[],[],[],minx,maxx,[],opts);
    [optx,fval,exitflag] = fmincon(@(X) objective_sublimation_model(X,randomsample,ztill,P,mask,nmask,plotFlag),X0,[],[],[],[],minx,maxx,@(X)nonlcon(X,randomsample,ztill,burialmask,P),opts);
    
    % Record optimizer output
    MCval(a) = fval;
    s(a) = optx(1);
    e(a) = optx(2);
    T(a) = optx(3);
    N10inh(a) = optx(4);
    N21inh(a) = optx(5);
    N26inh(a) = optx(6);

    flag(a) = exitflag;

    % Get data on debris conc. and elevation lowering
    opt_data.s = s(a).*1e-4;
    opt_data.ET = e(a).*100.*(ztill.dz./ztill.d)./1e6;
    opt_data.T = T(a).*1e6;
    opt_data.Ztill = ztill.dz;
    
    opt_data = sublimation_model_params(opt_data);
    
    Cdebris(a) = opt_data.CD;
    z_ice_loss(a) = opt_data.z_ice_init/100;
    
    % Get burial age
    tb = get_burial_age(optx,randomsample,ztill,burialmask,P);
    
    tb_all{a,:} = tb;
    
end 

MC_results = table(MCval,s,e,T,N10inh,N21inh,N26inh,Cdebris,z_ice_loss,tb_all,flag);

%% Finish and save data 

time = toc

% % Create new folder with current date/time as its name
% cd MC_output
% foldername = datestr(datetime,30);
% mkdir(foldername)
% 
% % save workspace as a .mat file
% save(fullfile(foldername, 'output.mat'));

%% Get/Display MC Results

MC_results.flag(find(MC_results.T < 0.02)) = 0; % throw away Age < 20kyr

good = find(MC_results.flag == 1);
okay = find(MC_results.flag == 2);
bad = find(MC_results.flag == 0); 

pass = [1 2];
ok = find(ismember(MC_results.flag,pass));

% Get mean and stdev

chi = MC_results.MCval(ok);
mean_chi = mean(MC_results.MCval(ok));
std_chi = std(MC_results.MCval(ok));
s = MC_results.s(ok);
mean_s = mean(MC_results.s(ok));
std_s = std(MC_results.s(ok));
e = MC_results.e(ok);
mean_e = mean(MC_results.e(ok));
std_e = std(MC_results.e(ok));
T = MC_results.T(ok);
mean_T = mean(MC_results.T(ok));
std_T = std(MC_results.T(ok));
N10inh = MC_results.N10inh(ok)/10;
mean_N10inh = mean(MC_results.N10inh(ok))/10;
std_N10inh = std(MC_results.N10inh(ok))/10;
N21inh = MC_results.N21inh(ok);
mean_N21inh = mean(MC_results.N21inh(ok));
std_N21inh = std(MC_results.N21inh(ok));
N26inh = MC_results.N26inh(ok)/10;
mean_N26inh = mean(MC_results.N26inh(ok))/10;
std_N26inh = std(MC_results.N26inh(ok))/10;

Cdebris = MC_results.Cdebris(ok);
mean_Cdebris = mean(MC_results.Cdebris(ok));
std_Cdebris = std(MC_results.Cdebris(ok));

z_ice_loss = MC_results.z_ice_loss(ok);
mean_z_ice_loss = mean(MC_results.z_ice_loss(ok));
std_z_ice_loss = std(MC_results.z_ice_loss(ok));

ztill_dz = ztill.dz;

% Display mean and stdev stats

disp(['',newline,'Monte-Carlo simulation results:'])
disp([sprintf('out of %0.0f MC iterations',numits)])
disp(['Average and STD statistics',newline,''])
disp([sprintf('Chi-square:          %0.1f \x00B1 %0.1f',mean_chi,std_chi)])
disp([sprintf('Sublimation rate:    %0.1f \x00B1 %0.1f m/Myr',mean_s,std_s)])
disp([sprintf('Erosion rate:        %0.2f \x00B1 %0.2f m/Myr',mean_e,std_e)])
disp([sprintf('Age of ice:          %0.1f \x00B1 %0.1f Myr',mean_T,std_T)])
disp([sprintf('Be-10 inheritance:   %0.2f \x00B1 %0.2f Matoms/g',(mean_N10inh),(std_N10inh))])
disp([sprintf('Ne-21 inheritance:   %0.1f \x00B1 %0.1f Matoms/g',mean_N21inh,std_N21inh)])
disp([sprintf('Al-26 inheritance:   %0.1f \x00B1 %0.1f Matoms/g',(mean_N26inh),(std_N26inh))])
disp([sprintf('Debris Conc.:        %0.3f \x00B1 %0.3f',mean_Cdebris,std_Cdebris)])
disp([sprintf('Ice elevation loss:  %0.1f \x00B1 %0.1f m',mean_z_ice_loss,std_z_ice_loss)])


% Get CDF Statistics and display CDF Optimization stats

% stats = [median negSigma1 posSigma1]
stats = [0.5 0.16 0.84]; % 68% confidence interval
%stats = [0.5 0.025 0.975]; % 95% confidence interval

[cdfchi,sortchi] = ecdf(chi);
[cdfs,sorts] = ecdf(s);
[cdfe,sorte] = ecdf(e);
[cdfT,sortT] = ecdf(T);
[cdfN10inh,sortN10inh] = ecdf(N10inh);
[cdfN21inh,sortN21inh] = ecdf(N21inh);
[cdfN26inh,sortN26inh] = ecdf(N26inh);

[cdfCdebris,sortCdebris] = ecdf(Cdebris);
[cdfz_ice_loss,sortz_ice_loss] = ecdf(z_ice_loss);
% [cdftbmin,sorttbmin] = ecdf(tbmin);

for a = 1:length(stats)
    
    cdfchi_stats(a) = interp1(cdfchi,sortchi,stats(a));
    cdfs_stats(a) = interp1(cdfs,sorts,stats(a));   
    cdfe_stats(a) = interp1(cdfe,sorte,stats(a));
    cdfT_stats(a) = interp1(cdfT,sortT,stats(a));
    cdfN10inh_stats(a) = interp1(cdfN10inh,sortN10inh,stats(a));
    cdfN21inh_stats(a) = interp1(cdfN21inh,sortN21inh,stats(a));
    cdfN26inh_stats(a) = interp1(cdfN26inh,sortN26inh,stats(a));
    
    cdfCdebris_stats(a) = interp1(cdfCdebris,sortCdebris,stats(a));
    cdfz_ice_loss_stats(a) = interp1(cdfz_ice_loss,sortz_ice_loss,stats(a));
%     cdftbmin_stats(a) = interp1(cdftbmin,sorttbmin,stats(a));
end

disp(['',newline,'Cumulative Distribution Statistics',newline,''])
disp([sprintf('Chi-square:          %0.1f \x00B1 %0.1f(%0.1f)',cdfchi_stats(1),cdfchi_stats(3)-cdfchi_stats(1),cdfchi_stats(1)-cdfchi_stats(2))])
disp([sprintf('Sublimation rate:    %0.2f \x00B1 %0.2f(%0.2f) m/Myr',cdfs_stats(1),cdfs_stats(3)-cdfs_stats(1),cdfs_stats(1)-cdfs_stats(2))])
disp([sprintf('Erosion rate:        %0.2f \x00B1 %0.2f(%0.2f) m/Myr',cdfe_stats(1),cdfe_stats(3)-cdfe_stats(1),cdfe_stats(1)-cdfe_stats(2))])
disp([sprintf('Age of ice:          %0.2f \x00B1 %0.2f(%0.2f) Myr',cdfT_stats(1),cdfT_stats(3)-cdfT_stats(1),cdfT_stats(1)-cdfT_stats(2))])
disp([sprintf('Be-10 inheritance:   %0.2f \x00B1 %0.2f(%0.2f) Matoms/g',cdfN10inh_stats(1),(cdfN10inh_stats(3)-cdfN10inh_stats(1)),(cdfN10inh_stats(1)-cdfN10inh_stats(2)))])
disp([sprintf('Ne-21 inheritance:   %0.1f \x00B1 %0.1f(%0.1f) Matoms/g',cdfN21inh_stats(1),(cdfN21inh_stats(3)-cdfN21inh_stats(1)),(cdfN21inh_stats(1)-cdfN21inh_stats(2)))])
disp([sprintf('Al-26 inheritance:   %0.2f \x00B1 %0.2f(%0.2f) Matoms/g',cdfN26inh_stats(1),(cdfN26inh_stats(3)-cdfN26inh_stats(1)),(cdfN26inh_stats(1)-cdfN26inh_stats(2)))])
disp([sprintf('Debris Conc.:        %0.3f \x00B1 %0.3f(%0.3f)',cdfCdebris_stats(1),cdfCdebris_stats(3)-cdfCdebris_stats(1),cdfCdebris_stats(1)-cdfCdebris_stats(2))])
disp([sprintf('Ice elevation loss:  %0.1f \x00B1 %0.1f(%0.1f) m',cdfz_ice_loss_stats(1),cdfz_ice_loss_stats(3)-cdfz_ice_loss_stats(1),cdfz_ice_loss_stats(1)-cdfz_ice_loss_stats(2))])

%% Make some noise when done

S(1) = load('chirp');
S(2) = load('splat');
sound(S(1).y,S(1).Fs)
pause(1.6)
sound(S(2).y,S(2).Fs)



