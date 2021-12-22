% This script calculates production-depth profiles for spallogenic and muon
% production and saves them for later use
%
% Specifically, this script defines a series of depths between zero and
% some large values, and computes spallation, muon and total production
% rates at all depths. The production rates are stored as a structure of 4
% vectors called P that containing the depth values and the production
% rates for each of the nuclides, 10Be, 26Al, and 21Ne, and saves it as a
% .mat file called Pall. This script has no input or arguments but does
% contain parameters such as elevation and latitude specific to Ong Valley
% drill site. If used for other locations, these should be changed
% accordingly.
%
% The total production rate is the production of cosmogenic nuclide
% produced by spallation and muon interactions. Where the production due to
% spallation using stone2000.m
%
% Greg Balco and
% Marie Bergelin
% June 2019

clear all; close all;

% Define site dependent parameters

elv_m = 1678.4; % site elevation in meters.
Lat = -83.25; % site latitude in decimal degreees 

% 0. Define range of mass depths at which to calculate production rates. 

zg = [0 logspace(0,5,500)]; % Units of g/cm2 (mass depth). 
% Calculate out to max depth of 100,000 g/cm2 (1000 m of ice). Which should
% be way more than needed. Use log mesh. 

% 1. Spallogenic. 

% Define surface production rates. 
% Figure out Be-10 production rate from reference production rate and 
% scaling method at site.

pressure = antatm(elv_m); % convert to atmospheric pressure in hPa

% % Using the 'St' scaling method
% sf_St = stone2000(Lat,pressure,1);
% 
% % Note final argument of 1 suppresses calculation of muon scaling, so you
% % are only calculating spallation scaling
% 
% P10_sp_surf = sf_St.*4; % Multiply scaling factor (nondimensional) and reference production rate (atoms/g/yr)

% Using the 'LSD'n scaling method

    % CRONUS_P10 = 25.1545; % this value includes shielding factor of 0.993
    % P10_sp_surf = 25.1545./0.993
    
P10_sp_surf = 25.3318; % CRONUS, without muon and shielding

% Now also calculate Al-26 and Ne-21 production rates. You can use a
% shortcut here and just assume a fixed ratio with Be-10 production.
P26_sp_surf = 6.75.*P10_sp_surf;
P21_sp_surf = 4.03.*P10_sp_surf;

% Apply exponential scaling to produce depth profile for spallation
Lsp = 140; % units of g/cm2

topo_sf = 0.993; % measured in the field and calculated in CRONUS.

P10_sp = P10_sp_surf.*exp(-zg./Lsp)*topo_sf;
P26_sp = P26_sp_surf.*exp(-zg./Lsp)*topo_sf;
P21_sp = P21_sp_surf.*exp(-zg./Lsp)*topo_sf;

% Make some plots

figure('pos',[330   360   860   440]);
subplot(1,3,1);
plot(P10_sp,zg,'r--'); hold on;
title('^1^0Be');
xlabel('P_1_0 (atoms g^-^1 yr^-^1)');
ylabel('Mass depth (g/cm^2)');
set(gca,'ydir','reverse','xscale','log','xlim',[0.001 100],'ylim',[0 10000]);
grid on;

subplot(1,3,2);
plot(P26_sp,zg,'b--'); hold on;
title('^2^6Al');
xlabel('P_2_6 (atoms g^-^1 yr^-^1)');
%ylabel('Mass depth (g/cm2)');
set(gca,'ydir','reverse','xscale','log','xlim',[0.001 100],'ylim',[0 10000]);
grid on;

subplot(1,3,3);
plot(P21_sp,zg,'g--'); hold on;
title('^2^1Ne');
xlabel('P_2_1 (atoms g^-^1 yr^-^1)');
%ylabel('Mass depth (g/cm2)');
set(gca,'ydir','reverse','xscale','log','xlim',[0.001 100],'ylim',[0 10000]);
grid on;

drawnow;

% 2. Production due to muons

% This uses the script 'P_mu_total_alpha1.m', which implements the method
% described in two papers by Heisinger. There is some documentation of how
% the code works here: 
%
% http://hess.ess.washington.edu/math/docs/al_be_v2/al_be_fctn_desc/node21.html
%
% Unfortunately that is for a previous version so it is slightly obsolete, 
% but it gives the general idea

% First, need to define interaction likelihoods for muon production.
% 

% For Be-10
mc10.Natoms = 2.006e22;
mc10.k_neg = 0.00191 .* 0.704 .* 0.1828; % From BCO fit
mc10.sigma0 = 0.280e-30; % From BCO fit

% For Al-26
mc26.Natoms = 1.003e22;
mc26.k_neg = 0.0133 .* 0.296 .* 0.6559; % From BCO fit
mc26.sigma0 = 3.89e-30; % From BCO fit

% For Ne-21
mc21.Natoms = 1.003e22;
mc21.k_neg = 0;
mc21.sigma0 = 4.2e-30;

% Now calculate muon production. The arguments to P_mu_total_alpha1 are the
% mass depth, the atmospheric pressure at the site, and the data structure
% with the interaction constants 

P10_mu = P_mu_total_alpha1(zg,pressure,mc10);
P26_mu = P_mu_total_alpha1(zg,pressure,mc26);
P21_mu = P_mu_total_alpha1(zg,pressure,mc21);

% More plotting

subplot(1,3,1);
plot(P10_mu,zg,'r');
subplot(1,3,2);
plot(P26_mu,zg,'b');
subplot(1,3,3);
plot(P21_mu,zg,'g');


drawnow;

% 3. Calculate total production

P10 = P10_sp + P10_mu;
P26 = P26_sp + P26_mu;
P21 = P21_sp + P21_mu;

% More plotting
subplot(1,3,1);
plot(P10,zg,'r','linewidth',1.5);
subplot(1,3,2);
plot(P26,zg,'b','linewidth',1.5);
subplot(1,3,3);
plot(P21,zg,'g','linewidth',1.5);

% Save data in mat-file
P.zg = zg;
P.P = [P10' P21' P26']';

save Pall P
