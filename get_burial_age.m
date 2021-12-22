%% Burial Constraint to optimization model
function [tb] = get_burial_age(X,samples,ztill,burialmask,P)

% This function accepts the optimized 'X' structure and the till thickness
% structure 'ztill', obtain the sublimation parameters and set the
% inheritance parameters to zero, runs the forward model and add the
% predicted concentrations to the sampledata, finds the difference between
% the modeled and measured nuclide concentrations and caluclate the
% apparent burial age for a simple exposure and having steady-state
% erosion using the fsolver funtion.
% 
% 
% The output is a burial age structure 'tb', which includes the following
% for each nuclide ratio: Be-10/Al-26, Be-10/Ne-21, Al-26/Ne-21
% 
% tb.te_tb = [te tb] % burial age obtained from simple expore
% tb.ee_tb = [ee tb] % burial age obtained from steady-state erosion
% 
% NaN is used for samples not used for burial calculations as indictated in
% burialmaks (= 0).
%
% 
%
% Marie Bergelin
% August 2021


plotFlag = 0;

% Create structure to obtain sublimation model parameters

data.s = X(1).*1e-4;
data.ET = X(2).*100.*(ztill.dz./ztill.d)./1e6;
data.T = X(3).*1e6;
data.N10inh = 0; % Inheritance should be zero 
data.N21inh = 0;
data.N26inh = 0;

data.Ztill = ztill.dz;

data = sublimation_model_params(data);
    
% Add Predicted sample concentartions to sample data set 
samples = plot_forward_model_sample_depth(data,samples,P,plotFlag);

% subtracting model predicted nuclides from data
% this gives us the amount of nuclides present before

N10diff = zeros(1,length(samples));
N21diff = N10diff;
N26diff = N10diff;


for j = 1:length(samples)

    if burialmask(j) >= 1;

        N10diff(j) = samples{j}.N10 - samples{j}.N10p;
        N21diff(j) = samples{j}.N21 - samples{j}.N21p;
        N26diff(j) = samples{j}.N26 - samples{j}.N26p;
    else 
        N10diff(j) = NaN;
        N21diff(j) = NaN;
        N26diff(j) = NaN;
    end
end

% Define constants
% Define decay constants
l10 = 4.99e-7;
l26 = 9.83e-7;

L = 140; % g/cm2 Attenuation length (140 g/cm2 for Antarctica)

options = optimset('Display','off');

for a = 1:length(samples)

    if burialmask(a) >= 1;
        
        % define objective function using weird MATLAB syntax for function definitions

        % returns zero at correct value of X
        
        % Exposure line
        
        % N10 = (P10/l10).*(1-exp(-l10.*te)).*exp(-tb.*l10)
        % N21 = P21.*te
        % N26 = (P26/l26).*(1-exp(-l26.*te)).*exp(-tb.*l26)
        
        % vector valued function of X = [te tb]
        btefunc_N26N10 = @(X) [N10diff(a) N26diff(a)] - ([P.P(1,1) P.P(3,1)]./[l10 l26]).*(1-exp(-[l10 l26].*X(1))).*exp(-X(2).*[l10 l26]);
        btefunc_N10N21 = @(X) [N10diff(a) - (P.P(1,1)./l10).*(1-exp(-l10.*X(1))).*exp(-X(2).*l10), N21diff(a) - P.P(2,1)*X(1)];
        btefunc_N26N21 = @(X) [N26diff(a) - (P.P(3,1)./l26).*(1-exp(-l26.*X(1))).*exp(-X(2).*l26), N21diff(a) - P.P(2,1)*X(1)];

        tb{a}.te_tb_N26N10 = fsolve(btefunc_N26N10,[1e6 1e6],options);
        tb{a}.te_tb_N10N21 = fsolve(btefunc_N10N21,[1e6 1e6],options);
        tb{a}.te_tb_N26N21 = fsolve(btefunc_N26N21,[1e6 1e6],options);

        % Steady-State Erosion

        % N10 = (P10./(l10 + ee./L)).*exp(-l10.*tb);
        % N21 = P21.*L./ee
        % N26 = (P26./(l26 + ee./L)).*exp(-l26.*tb);
        
        % vector valued function of X = [ee tb]
        
        beefunc_N26N10 = @(X) [N10diff(a) N26diff(a)] - ([P.P(1,1) P.P(3,1)]./([l10 l26] + X(1)./L)).*exp(-X(2).*[l10 l26]);
        beefunc_N10N21 = @(X) [N10diff(a) - (P.P(1,1)./(l10 + X(1)./L)).*exp(-X(2).*l10), N21diff(a) - P.P(2,1).*L./X(1)];
        beefunc_N26N21 = @(X) [N26diff(a) - (P.P(3,1)./(l26 + X(1)./L)).*exp(-X(2).*l26), N21diff(a) - P.P(2,1).*L./X(1)];
        
        tb{a}.ee_tb_N26N10 = fsolve(beefunc_N26N10,[0.0001 1e6],options);
        tb{a}.ee_tb_N10N21 = fsolve(beefunc_N10N21,[0.0001 1e6],options);
        tb{a}.ee_tb_N26N21 = fsolve(beefunc_N26N21,[0.0001 1e6],options);

        %convert erosion rate from g/cm2/yr to m/Myr
        tb{a}.ee_tb_N26N10(1) = tb{a}.ee_tb_N26N10(1)./(100.*(ztill.dz./ztill.d)./1e6);
        tb{a}.ee_tb_N10N21(1) = tb{a}.ee_tb_N10N21(1)./(100.*(ztill.dz./ztill.d)./1e6);
        tb{a}.ee_tb_N26N21(1) = tb{a}.ee_tb_N26N21(1)./(100.*(ztill.dz./ztill.d)./1e6);

    else
        tb{a}.te_tb_N26N10 = [NaN NaN];
        tb{a}.te_tb_N10N21 = [NaN NaN];
        tb{a}.te_tb_N26N21 = [NaN NaN];
        
        tb{a}.ee_tb_N26N10 = [NaN NaN];
        tb{a}.ee_tb_N10N21 = [NaN NaN];
        tb{a}.ee_tb_N26N21 = [NaN NaN];       
 
    end
    
end

end

