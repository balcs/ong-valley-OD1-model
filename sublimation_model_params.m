function out = sublimation_model_params(in)

% This calculates various parameters of a sublimation-erosion model
% for till-covered ice. 
%
% out = sublimation_model_params(in)
% 
% in has following fields:
% 
% in.T - age of ice (yr)
% in.s - sublimation rate, i.e. rate that ice surface is lowering (cm/yr)
% in.ET - surface erosion rate of till (g/cm2/yr)
% in.Ztill - present till thickness (g/cm2)
% 
% Calculates other stuff and adds it as additional fields to the input
% structure. Returns the input structure with additional fields. 
% out.z_ice_init - implied ice thickness consumed during model time
% out.CDrhoM - implied density of debris in ice/debris mixture
% out.CD - implied weight fraction of debris in ice
% out.rhoM - implied density of ice/debris mixture
% 
% Greg Balco
% June, 2019

% Compute total ice consumed during model time

in.z_ice_init = in.s.*in.T; % ice thickness consumed during model (cm)

% Compute debris concentration in ice from sublimation rate
% and till thickness

% First, compute the quantity CD*rhoM (density of debris in mixture)
in.CDrhoM = (in.Ztill + in.T.*in.ET)./(in.T.*in.s);

% Now, compute CD (wt/wt debris concentration) and rhoM (mixture density) separately
% in.CD = (2.68*in.CDrhoM)./(0.91*2.68 + 1.77.*in.CDrhoM);
in.CD = (2.68*in.CDrhoM)./(0.91*2.68 - (0.91-2.68).*in.CDrhoM);
in.rhoM = in.CDrhoM./in.CD;

out = in;