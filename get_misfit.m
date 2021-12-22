function out = get_misfit(data,samples,P,mask,nmask,dataDumpFlag,plotFlag)

% This function accepts the 'data' structure with sublimation model and
% inheritance parameters, runs the forward model, calculates predicted
% concentrations for selected set of samples, and returns a chi-squared
% misfit statistic. 
%
% out = get_misfit(data,samples,P,mask,nmask)
%
% data, samples, and P are upstream structures, mask is a n-element vector of
% zeros and ones that tells which samples to consider in calculating the
% misfit. 1 = use, 0 = don't use. nmask is a similar mask for
% nuclides...default is [1 1 1] meaning use all nuclides. 
%
% Greg Balco
%
% June 2019

if nargin < 7; plotFlag = 0; end
if nargin < 6; dataDumpFlag = 0; end
if nargin < 5; nmask = [1 1 1]; end

% Get prediction

p = plot_forward_model(data,P,0);

% Predict concentrations for all samples
% Initialize variables to save some time
N10p = zeros(1,length(samples));
N21p = N10p;
N26p = N10p;
M10 = N10p;
M21 = N10p;
M26 = N10p;

% Calculate cumulative sums here for trapezoidal integration
% this is based on Reimann sums
cum10 = cumsum(p.N10.*[0 diff(p.z)]);
cum21 = cumsum(p.N21.*[0 diff(p.z)]);
cum26 = cumsum(p.N26.*[0 diff(p.z)]);

if plotFlag == 1; delete(findobj('tag','misfit')); end

for a = 1:length(samples)
    if mask(a) == 1
        % OK to calculate
        this_sample = samples{a};
        segN10 = zeros(size(this_sample.Mq));
        segN21 = segN10;
        segN26 = segN10;
        % Calculate average concentrations in each interval
        % based on trapzodidal integration
        for b = 1:length(this_sample.Mq)
            zv = [this_sample.Mqtdz(b) this_sample.Mqbdz(b)];
            segN10(b) = diff(interp1(p.z,cum10,zv))./diff(zv);
            segN21(b) = diff(interp1(p.z,cum21,zv))./diff(zv);
            segN26(b) = diff(interp1(p.z,cum26,zv))./diff(zv);
        end
        % Combine based on quartz weight contributed
        N10p(a) = sum(segN10.*this_sample.Mq)./sum(this_sample.Mq);
        N21p(a) = sum(segN21.*this_sample.Mq)./sum(this_sample.Mq);
        N26p(a) = sum(segN26.*this_sample.Mq)./sum(this_sample.Mq);
        
        % Compute misfits
        % Below are error-weighted, i.e. chi-squared-like
        %M10(a) = ((N10p(a)-this_sample.N10)./this_sample.dN10); 
        %M21(a) = ((N21p(a)-this_sample.N21)./this_sample.dN21);
        %M26(a) = ((N26p(a)-this_sample.N26)./this_sample.dN26);
        
        % Below are evenly weighted by relative miss, not absolute miss
        M10(a) = ((N10p(a)-this_sample.N10)./this_sample.N10); 
        M21(a) = ((N21p(a)-this_sample.N21)./this_sample.N21);
        M26(a) = ((N26p(a)-this_sample.N26)./this_sample.N26);
        
        if plotFlag == 1
            % Update plots if requested
            figure(1); 
            subplot(1,3,1);
            plot([N10p(a) N10p(a)],[this_sample.tdz this_sample.bdz],'k','linewidth',3,'tag','misfit');
            subplot(1,3,2);
            plot([N21p(a) N21p(a)],[this_sample.tdz this_sample.bdz],'k','linewidth',3,'tag','misfit');
            %plot(N21p(a),this_sample.bdz,'ko','markerfacecolor','m','tag','misfit');
            subplot(1,3,3);
            plot([N26p(a) N26p(a)],[this_sample.tdz this_sample.bdz],'k','linewidth',3,'tag','misfit');          
           
        end
    end
    
end

if plotFlag == 1; drawnow; end

% Process 'miss' vectors to remove NaNs

M10(find(isnan(M10))) = zeros(size(find(isnan(M10))));
M21(find(isnan(M21))) = zeros(size(find(isnan(M21))));
M26(find(isnan(M26))) = zeros(size(find(isnan(M26))));

if dataDumpFlag == 1
    % Return predicted concentrations if requested
    out.N10p = N10p; out.N21p = N21p; out.N26p = N26p;
else
    % Return single misfit value for use by optimizer 
    out = nmask(1).*sum(M10(mask > 0).^2) + nmask(2).*sum(M21(mask > 0).^2) + nmask(3).*sum(M26(mask > 0).^2);

end



