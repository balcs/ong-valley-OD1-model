function out = plot_forward_model_sample_depth(data,samples,P,plotFlag)

% This generates a model prediction for sample segments identified in the
% data structure and optionally plots it on the same plots as the
% measurements.
% Run plot_core_1_data.m to get the base plots. 
% 
% out = plot_forward_model_sample_depth(data,samples,P,plotFlag)
% 
% data is data structure for sublimation model plus inheritance values
% data.N10inh, etc. 
% samples is the sampledata structure of desired sample segment
% P. is production rate structure.
% 
% returns predicted concentrations added into the sample structure
% set plotFlag = 0 to suppress plotting
%
% Greg Balco and
% Marie Bergelin
%
% July 2020

if nargin < 4; plotFlag = 0; end

% Get prediction

p = plot_forward_model(data,P,0);

% Calculate cumulative sums here for trapezoidal integration
cum10 = cumsum(p.N10.*[0 diff(p.z)]);
cum21 = cumsum(p.N21.*[0 diff(p.z)]);
cum26 = cumsum(p.N26.*[0 diff(p.z)]);

for a = 1:length(samples)

        % OK to calculate
        this_sample = samples{a};
        segN10 = zeros(size(this_sample.Mq));
        segN21 = segN10;
        segN26 = segN10;
        % Calculate average concentrations in each interval
        for b = 1:length(this_sample.Mq)
            zv = [this_sample.Mqtdz(b) this_sample.Mqbdz(b)];
            segN10(b) = diff(interp1(p.z,cum10,zv))./diff(zv);
            segN21(b) = diff(interp1(p.z,cum21,zv))./diff(zv);
            segN26(b) = diff(interp1(p.z,cum26,zv))./diff(zv);
        end
        % Combine based on quartz weight contributed
        this_sample.N10p = sum(segN10.*this_sample.Mq)./sum(this_sample.Mq);
        this_sample.N21p = sum(segN21.*this_sample.Mq)./sum(this_sample.Mq);
        this_sample.N26p = sum(segN26.*this_sample.Mq)./sum(this_sample.Mq);
        
        if plotFlag == 1
            % Update plots if requested
            figure(1); 
            subplot(1,3,1);
            plot([this_sample.N10p this_sample.N10p],[this_sample.tdz this_sample.bdz],'color',[1 0.7 0.7],'linewidth',2);
            hold on;
            
            subplot(1,3,2);
            plot([this_sample.N21p this_sample.N21p],[this_sample.tdz this_sample.bdz],'color',[0.7 0.7 1],'linewidth',2);
            hold on;
                
            subplot(1,3,3);
            plot([this_sample.N26p this_sample.N26p],[this_sample.tdz this_sample.bdz],'color',[0.7 1 0.7],'linewidth',2);          
            hold on;        
        end

        samples{a} = this_sample;
end

if plotFlag == 1
% plot model on ratio plot
    
    for a = 1:length(samples)
        nn10(a) = samples{a}.N10p./P.P(1,1);
        nn21(a) = samples{a}.N21p./P.P(2,1);
        nn26(a) = samples{a}.N26p./P.P(3,1);
    end
    
    figure(2);
    subplot(1,3,1);
    plot(nn10,nn26./nn10,'-','LineWidth',1);

    subplot(1,3,2);
    plot(nn21,nn10./nn21,'-','LineWidth',1);

    subplot(1,3,3);
    plot(nn21,nn26./nn21,'-','LineWidth',1);
end    


% Return results
out = samples;




