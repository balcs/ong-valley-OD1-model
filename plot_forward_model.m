function out = plot_forward_model(data,P,plotFlag)

% This generates a model prediction and optionally plots it on the same 
% plots as the measurements. Run plot_core_1_data.m to get the base plots. 
% 
% out = plot_forward_model(data,P,plotFlag)
% 
% data is data structure for sublimation model plus inheritance values
% data.N10inh, etc. P. is production rate structure.
% 
% returns predicted concentrations on log depth mesh
% set plotFlag = 0 to suppress plotting
%
% Greg Balco
%
% June 2019 

 

if nargin < 3; plotFlag = 1; end;

% Define mesh of depths at which to compute predicted concentrations
calczs = [0 logspace(-1,log10(1100),200)];

% Define decay constants
l10 = 4.99e-7;
l26 = 9.83e-7;

% Do forward integration at each depth
for a = 1:length(calczs)
    
    calcN10s(a) = data.N10inh.*exp(-l10.*data.T) + integral(@(t) exp(-t.*l10).*P_of_z(z_of_t(t,calczs(a),data),P,10),0,data.T);
    calcN21s(a) = data.N21inh + integral(@(t) P_of_z(z_of_t(t,calczs(a),data),P,21),0,data.T);
    calcN26s(a) = data.N26inh.*exp(-l26.*data.T) + integral(@(t) exp(-t.*l26).*P_of_z(z_of_t(t,calczs(a),data),P,26),0,data.T);

    % Without Ninh decay
%    calcN10s(a) = data.N10inh + integral(@(t) exp(-t.*l10).*P_of_z(z_of_t(t,calczs(a),data),P,10),0,data.T);
%    calcN21s(a) = data.N21inh + integral(@(t) P_of_z(z_of_t(t,calczs(a),data),P,21),0,data.T);
%    calcN26s(a) = data.N26inh + integral(@(t) exp(-t.*l26).*P_of_z(z_of_t(t,calczs(a),data),P,26),0,data.T);
   
     end


%% Plot results

if plotFlag == 1;

    % Plot on depth-concentration figure

%     figure(1); 
%     subplot(1,3,1);
%     plot(calcN10s,calczs,'-','LineWidth',1); set(gca,'ydir','reverse');
% 
%     subplot(1,3,2);
%     plot(calcN21s,calczs,'-','LineWidth',1);
% 
%     subplot(1,3,3);
%     plot(calcN26s,calczs,'-','LineWidth',1);
%     drawnow;

    figure(1); 
    subplot(1,3,1);
    plot(calcN10s,calczs,'-','color',[0 0.45 0.74],'LineWidth',1); set(gca,'ydir','reverse');
    set(gca,'FontSize',12)

    subplot(1,3,2);
    plot(calcN21s,calczs,'-','color',[0 0.45 0.74],'LineWidth',1);
    set(gca,'FontSize',12)

    subplot(1,3,3);
    plot(calcN26s,calczs,'-','color',[0 0.45 0.74],'LineWidth',1);
    set(gca,'FontSize',12)
    drawnow;


    % Plot on ratio plots

    load Pall; 
    nn10 = calcN10s./P.P(1,1);
    nn21 = calcN21s./P.P(2,1);
    nn26 = calcN26s./P.P(3,1);

    figure(2);
    subplot(1,3,1);
    plot(nn10,nn26./nn10,'-','color',[0 0.45 0.74],'LineWidth',1);

    subplot(1,3,2);
    plot(nn21,nn10./nn21,'-','color',[0 0.45 0.74],'LineWidth',1);

    subplot(1,3,3);
    plot(nn21,nn26./nn21,'-','color',[0 0.45 0.74],'LineWidth',1);
    
end

% Return results

out.z = calczs;
out.N10 = calcN10s;
out.N21 = calcN21s;
out.N26 = calcN26s;

