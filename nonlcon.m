%% Burial Constraint to optimization model

function [c, ceq,samples] = nonlcon(X,samples,ztill,burialmask,P)

% This function is the non linear constraint used for optimizing the
% forward model and follows that of nonlcon in fminsearch: a Nonlinear
% constraints, specified as a function handle or function name accepts a
% vector or array x and returns two arrays, c(x) and ceq(x).

% This function accepts the optimized 'X' structure and the till thickness
% structure 'ztill', obtain the sublimation model parameters and set the
% inheritance parameters to zero, runs the forward model and add the
% predicted concentrations to the sampledata, finds the difference between
% the modeled and measured nuclide concentrations. When adding the
% concentrtion lost by decay during exposure time, the nuclide ratio for
% any of the pair cannot exceed the simple exposure line.
% 
% These ratios are used as burial constraints for the the forward model
% optimization as the optimized age cannot be older than the minimum burial age
% for any samples across all nuclide ratios.
%
% the outputs are
% c -     % Compute nonlinear inequalities at x
% ceq -   % Compute nonlinear equalities at x
%
%
% Marie Bergelin
%


plotFlag = 0;

data.s = X(1).*1e-4;
data.ET = X(2).*100.*(ztill.dz./ztill.d)./1e6;
data.T = X(3).*1e6;
data.N10inh = 0; % Inheritance should be zero 
data.N21inh = 0;
data.N26inh = 0;


data.Ztill = ztill.dz;

data = sublimation_model_params(data);

if data.CD >= 1
    
    c = data.CD - 1;
    ceq = [];
else
    
    % Add Predicted sample concentartions to sample data set 
    samples = plot_forward_model_sample_depth(data,samples,P,plotFlag);


    % subtracting model predicted nuclides from data
    % this gives us the amount of nuclides present before

    N10diff = zeros(1,length(samples));
    N21diff = N10diff;
    N26diff = N10diff;

    for j = 1:length(burialmask)

        if burialmask(j) >= 1

            N10diff(j) = samples{j}.N10 - samples{j}.N10p;
            N21diff(j) = samples{j}.N21 - samples{j}.N21p;
            N26diff(j) = samples{j}.N26 - samples{j}.N26p;
        end
    end

    % Add back what was lost to decay during burial
    tb = data.T;

    % Define decay constants
    l10 = 4.99e-7;
    l26 = 9.83e-7;

    % this would have been the concentration of nuclides in the old drift
    % before burial occured (values are normalized to production rates)
    oldN10norm = (N10diff./exp(-l10*tb))./P.P(1,1);
    oldN21norm = N21diff./P.P(2,1);
    oldN26norm = (N26diff./exp(-l26*tb))./P.P(3,1);


% The concentrations cannot exceed the exposure line in the burial plot
% in which the ratio decreases as age increases.
% 
%     % Al-26/Be-10 constraint
%     N26N10con = oldN26norm - (1/l26*(1-exp(-l26*(-1/l10*log(1-oldN10norm*l10)))));
% 
%     % Be-10/Ne-21 constraint
%     N10N21con = oldN10norm - 1/l10*(1-exp(-l10*oldN21norm));
%     
%     % Al-26/Ne-21 constraint
%     N26N21con = oldN26norm - 1/l26*(1-exp(-l26*oldN21norm));


    % Al-26/Be-10 constraint
    N26N10con = oldN26norm - (1/l26*(1-exp(-l26*(-1/l10*log(1-oldN10norm*l10)))));

    % Be-10/Ne-21 constraint
    N10N21con = oldN10norm - 1/l10*(1-exp(-l10*oldN21norm));
    
    % Al-26/Ne-21 constraint
    N26N21con = oldN26norm - 1/l26*(1-exp(-l26*oldN21norm));



    % find the maximum value of all burial ratios such that it satisfy the
    % constraint of c <= 0
    c = max([N26N10con N10N21con N26N21con]); % c <= 0
    ceq = []; % ceq = 0
    
end

end
