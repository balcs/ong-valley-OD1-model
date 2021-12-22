function out = objective_sublimation_model(X,samples,ztill,P,mask,nmask,plotFlag)

% This is a wrapper for the get_misfit function that allows it to be used
% with the MATLAB optimizer fminsearch().
%
% function out = objective_sublimation_model(X,samples,ztill,P,mask,nmask);
%
% X is the parameter vector:
% X(1) sublimation rate m/Myr
% X(2) erosion rate m/Myr
% X(3) age Ma
% X(4) inherited N10 (100000 atoms)
% X(5) inherited N21 (Matoms)
% X(6) inherited N26 (100000 atoms)
%
% Greg Balco
%
% June 2019

if nargin < 7; plotFlag = 0; end

% Display input so you can see what optimizer is doing
% disp(sprintf('%0.2f %0.3f %0.2f %0.2f %0.2f %0.2f',X));

% Package input for use by get_misfit
data.s = X(1).*1e-4;
data.ET = X(2).*100.*(ztill.dz./ztill.d)./1e6;
data.T = X(3).*1e6;
data.N10inh = X(4).*1e5;
data.N21inh = X(5).*1e6;
data.N26inh = X(6).*1e5;

% Enter measured till thickness (g/cm2)
data.Ztill = ztill.dz;

% Get implied debris densities, etc. from sublimation_model_params.m
data = sublimation_model_params(data);

% Get misfit to packaged model and return to optimizer
out = get_misfit(data,samples,P,mask,nmask,0,plotFlag);

