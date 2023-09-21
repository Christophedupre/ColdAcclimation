close all
clear all
clc

%%% ARRANGED ACCORDING TO NARRATIVE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
% CompileDataAcclimation(); %Compile all data into 1 file (run once unless data has changed)
load('data.mat') %Load compiled data
GUIAcclimation();
% OVERALL HYPOTHESES
plotHypotheses();

% 1. Is cold noxious?
HealthDecayRecovery(Pilot, Deacclimation, CstData3D); %Show decay and recovery at 22C in naive and partially acclimated animals

% 2. Does gentle cooling help?
SlowDecrease(); %Includes plotcontrols.m
plotCstGrouped(CstData3D);
plotCstAnalysis(CstData3D, Cst10C13D30Reps, Pilot);
plotCstPilot; %Optional

% 3 and 4. Is it acclimation and is it reversible?
AcclimationPersistence(Deacclimation, JustPCR, NoPCR);

% 5. Can we explain it with a model?
[Model1Data Model2Data] = GenerateODEModelData(CstData3D); % Make ODE Model data, w/o & w variab.
modelexplanation(Model1Data, Model2Data); %Cartoons to explain model
plotCstGrouped(Model1Data.HealthTruncated);
plotODEModelPerformance(Model1Data, Model2Data, CstData3D, Cst10C13D30Reps, Ramp, Pulses1_1, Pilot, Deacclimation); % Plot performance of model
Explicitsolution(Model1Data);
plotLakeData(LakeMichigan, Model1Data);

% 6. Experimental advantages
plotCaImaging_Behavior();
plotFeedingFrequency();

%%% ARRANGED ACCORDING TO DATASETS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
% CompileDataAcclimation(); %Compile all data into 1 file (run once unless data has changed)
% load('data.mat') %Load compiled data
% % PILOT EXPERIMENTS
% AcclimationPilotData(Pilot);
% % CONTROLS 
% plotcontrols(PosCtrl, NoPCR, JustPCR, Ramp, Pulses0_2, Pulses1_1, CstData3D)
% % CONSTANT TEMPERATURE, SINGLE TRACES
% plotConstantSingleTraces(CstData3D)
% % CONSTANT TEMPERATURE, GROUPED TRACES
% plotCstGrouped(CstData3D)
% % MORE ANALYSIS OF THE CONSTANT TEMPERATURE DATA (probably supplemental)
% plotCstAnalysis(CstData3D, Cst10C13D30Reps, Pilot);
% % DEACCLIMATION
% AcclimationPersistence(Deacclimation);
%%
%%% MODEL1&2: DYNAMICAL SYSTEM MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Model1Data Model2Data] = GenerateODEModel1Data(CstData3D); % Make ODE Model data, w/o & w variab.
% modelexplanation(Model2Data, CstData3D); %Cartoons to explain model
% plotModelvsExptData_SingleTraces(Model1Data, CstData3D); % Plot and compare data as single traces
% plotODEModelPerformance(Model2Data, CstData3D, Cst10C13D30Reps, Ramp, Pulses1_1, Pilot, Deacclimation); % Plot performance of model
% Explicitsolution(Model1Data);
%%
%%% LAKE TEMPERATURES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Display and analysis of lake temperatures
% plotLakeData(LakeMichigan, Model1Data);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% COMPARISON PULSES WITH STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate 13C model data
figure('Name', 'Comparison Pulses with Steps','WindowState', 'maximized')
clf;
tiledlayout(4,4);
Model13CData1 = GenerateEmptyModelStructure();
Model13CData2 = GenerateEmptyModelStructure();
cmaplines = colormap('lines');
csol.c = Model1Data.csol.c;
tspan = Model1Data.tspan;
y0 = Model1Data.y0;
f = Model1Data.TimeResolution; %Time factor (increases temporal resolution during derivation of ODE solution
NoiseAmplitude = Model1Data.NoiseAmplitude; 
imax = 15;
for theta = 5.5
    for D = 1:length(CstData3D(1,:,1)) 
        %With variability
%         cs = csol.c.*(1+(rand(3,1)-0.5)/Model1Data.cvarcoeff);% Add variability for acclimation, injury and regeneration time constants
        cs = csol.c; % No variability in coefficients
        temperature1 = linspace(0,43,44*f);
        temperature2 = [linspace(22,22,(15-(D-1))*f) linspace(theta*2+2,theta*2+2,(D-1)*f) linspace(4,4,29*f)]./22;
        temperature = [temperature1; AddNoiseToTemp(temperature2, NoiseAmplitude)]; %Variability in temperature sensed
        sol = ode23(@(t,y) AcclimfunNoise(t,y,cs,temperature), tspan,y0);
        solpts = deval(sol, tspan);
        Model13CData2.HealthFull(:,D,1) = solpts(2,:);        
        Model13CData2.HealthTruncated(:,D,1) = solpts(2,linspace(15,43,29));
        Model13CData2.AlphaFull(:,D,1) = solpts(1,:);
        Model13CData2.AlphaTruncated(:,D,1) = solpts(1,linspace(15,43,29));
        Model13CData2.csol = csol;

        %Without variability
        cs = csol.c; % No variability
        temperature = [temperature1; temperature2];
        sol = ode23(@(t,y) Acclimfun(t,y,cs,temperature), tspan,y0);
        solpts = deval(sol, tspan);
        Model13CData1.HealthFull(:,D,1) = solpts(2,:);        
        Model13CData1.HealthTruncated(:,D,1) = solpts(2,linspace(15,43,29));
        Model13CData1.AlphaFull(:,D,1) = solpts(1,:);
        Model13CData1.AlphaTruncated(:,D,1) = solpts(1,linspace(15,43,29));
        Model13CData1.csol = csol;
    end
end


for i = [5] 
    nexttile();
    h = heatmap(squeeze(CstData3D(:,:,i))');
    h.ColorbarVisible = 'off'; h.XLabel = 'Testing Time (Days)'; h.YLabel = 'Training Duration';
    h.FontSize = 8;
    % h.Title = [Labels(i) ' T= ' num2str(i*2+2) ' ^oC'];
    h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
    h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    % survival = [survival sum(CstData3D(:,:,i)>0.1,2);];
    h.Title = '13C Step - Health expt';

    nexttile()
    h2 = heatmap(squeeze(Model13CData1.HealthTruncated(:,:,1))');
    h2.ColorbarVisible = 'off'; h2.XLabel = 'Testing Time (Days)'; h2.YLabel = 'Training Duration';
    h2.FontSize = 8;
    h2.XDisplayLabels(:) = {''}; h2.YDisplayLabels(:) = {''};
    h2.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    h2.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    h2.Title = 'Health model';

    nexttile()
    h3 = heatmap(squeeze(Model13CData2.HealthTruncated(:,:,1))');
    h3.ColorbarVisible = 'off'; h3.XLabel = 'Testing Time (Days)'; h3.YLabel = 'Training Duration';
    h3.FontSize = 8;
    h3.XDisplayLabels(:) = {''}; h3.YDisplayLabels(:) = {''};
    h3.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    h3.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    h3.Title = 'Health model with noise';

    nexttile(); hold on;
    plot(linspace(0,14,15), squeeze(mean(CstData3D(:,:,i),1)), 'LineWidth',2, 'Color', cmaplines(1,:));
    plot(linspace(0,14,15), squeeze(mean(Model13CData1.HealthTruncated(:,:,1),1)), 'LineWidth',2, 'Color', cmaplines(5,:));
    plot(linspace(0,14,15), squeeze(mean(Model13CData2.HealthTruncated(:,:,1),1)), 'LineWidth',2, 'Color', cmaplines(5,:), 'LineStyle',':');
    box off; ax = gca ; ax.LineWidth = 1; ax.FontSize = 8;
    xlim([0 15]); xticks([0:5:15]); ylim([0 1]); yticks([0:1:1]);
    text(-0.2*mean(ax.XLim), mean(ax.YLim), {'avg. health'}, 'Rotation',90, 'VerticalAlignment','middle','HorizontalAlignment','center');
    % text(mean(ax.XLim),-0.4*mean(ax.YLim),'D', 'HorizontalAlignment','center');
    ax = gca; ax.XLabel.String = 'Training Duration (Days)';
    % text(ax.XLim(1),ax.YLim(2)*1.15, 'N', 'FontSize', 14);

end

colormap(cmaplist.cmapRedBlue1);

nexttile()
h4 = heatmap(Pulses1_1');
h4.ColorbarVisible = 'off'; h4.XLabel = 'Testing Time (Days)'; h4.YLabel = 'Training Duration (Day)';
h4.FontSize = 8;
h4.Title = 'Pulses - Health expt';
h4.XDisplayLabels(:) = {''}; h4.YDisplayLabels(:) = {''};
h4.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h4.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
h4.Colormap = cmaplist.cmapRedBlue1;

nexttile(); 
% Model without noise
for i = 1:imax;
    pd = i; %Pulse (i.e. training) duration (days)
%     cvar = Model2Data.csol.c.*(1+(rand(3,1)-0.5)/4);% Add variability for acclimation, injury and regeneration time constants
    cvar = Model2Data.csol.c;% Add variability for acclimation, injury and regeneration time constants
    temperature1 = linspace(0,43,44*24); 
    temperature2 = [linspace(22,22,(15-pd)*24) repmat([linspace(22,22,1) linspace(4,4,1)],1,12*pd) linspace(4,4,29*24)]./22;
    % temperature2 = AddNoiseToTemp(temperature2,NoiseAmplitude);
    temperature = [temperature1; temperature2];
    sol = ode23(@(t,y) Acclimfun(t,y,cvar,temperature), Model2Data.tspan, Model2Data.y0);
    solptstemp = deval(Model2Data.tspan,sol);    
    Model1Pulses1_1Data.HealthFull(:,1,i) = solptstemp(2,:);
    Model1Pulses1_1Data.HealthTruncated(:,1,i) = solptstemp(2,linspace(15,43,29));
    % plot(Model2Data.tspan,solptstemp(2,:), 'Color', [cmapjet(round(i/imax*256),:) 1], 'LineWidth', 2);
end
% ax = gca; ax.XLim = [-0.1 44]; ax.XTick = [0:10:44]; ax.YLim = [0 1]; ax.YTick = [0:1:1]; ax.LineWidth = 1; ax.FontSize = 8;
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Health model', 'HorizontalAlignment','center','Rotation',90);
h3 = heatmap(squeeze(Model1Pulses1_1Data.HealthTruncated(:,1,:))');
h3.ColorbarVisible = 'off'; h3.XLabel = 'Testing Time (Days)'; h3.YLabel = 'Training Duration (Day)';
h3.FontSize = 8;
h3.Title = 'Health model';
h3.XDisplayLabels(:) = {''}; h3.YDisplayLabels(:) = {''};
h3.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h3.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
h3.Colormap = cmaplist.cmapRedBlue1;
nexttile();

% Model with noise
for i = 1:imax;
    pd = i; %Pulse (i.e. training) duration (days)
%     cvar = Model2Data.csol.c.*(1+(rand(3,1)-0.5)/4);% Add variability for acclimation, injury and regeneration time constants
    cvar = Model2Data.csol.c;% Add variability for acclimation, injury and regeneration time constants
    temperature1 = linspace(0,43,44*24); 
    temperature2 = [linspace(22,22,(15-pd)*24) repmat([linspace(22,22,1) linspace(4,4,1)],1,12*pd) linspace(4,4,29*24)]./22;
    temperature2 = AddNoiseToTemp(temperature2,NoiseAmplitude);
    temperature = [temperature1; temperature2];
    sol = ode23(@(t,y) AcclimfunNoise(t,y,cvar,temperature), Model2Data.tspan, Model2Data.y0);
    solptstemp = deval(Model2Data.tspan,sol);    
    Model2Pulses1_1Data.HealthFull(:,1,i) = solptstemp(2,:);
    Model2Pulses1_1Data.HealthTruncated(:,1,i) = solptstemp(2,linspace(15,43,29));
    % plot(Model2Data.tspan,solptstemp(2,:), 'Color', [cmapjet(round(i/imax*256),:) 1], 'LineWidth', 2);
end
% ax = gca; ax.XLim = [-0.1 44]; ax.XTick = [0:10:44]; ax.YLim = [0 1]; ax.YTick = [0:1:1]; ax.LineWidth = 1; ax.FontSize = 8;
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Health model', 'HorizontalAlignment','center','Rotation',90);
h3 = heatmap(squeeze(Model2Pulses1_1Data.HealthTruncated(:,1,:))');
h3.ColorbarVisible = 'off'; h3.XLabel = 'Testing Time (Days)'; h3.YLabel = 'Training Duration (Day)';
h3.FontSize = 8;
h3.Title = 'Health model with noise';
h3.XDisplayLabels(:) = {''}; h3.YDisplayLabels(:) = {''};
h3.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h3.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
h3.Colormap = cmaplist.cmapRedBlue1;

% Training efficiency
nexttile(); hold on;
plot(linspace(0,14,15), mean(Pulses1_1,1), 'LineWidth',2, 'Color', cmaplines(1,:));
plot(linspace(0,14,15), squeeze(mean(Model1Pulses1_1Data.HealthFull(15:end,1,:),1)), 'LineWidth',2, 'Color', cmaplines(5,:));
plot(linspace(0,14,15), squeeze(mean(Model2Pulses1_1Data.HealthFull(15:end,1,:),1)), 'LineWidth',2, 'Color', cmaplines(5,:), 'LineStyle',':');
plot(linspace(0,14,15), squeeze(mean(Model13CData1.HealthTruncated(:,:,1),1)), 'LineWidth',2, 'Color', cmaplines(6,:));
    plot(linspace(0,14,15), squeeze(mean(Model13CData2.HealthTruncated(:,:,1),1)), 'LineWidth',2, 'Color', cmaplines(6,:), 'LineStyle',':');
box off; ax = gca ; ax.LineWidth = 1; ax.FontSize = 8;
xlim([0 15]); xticks([0:5:15]); ylim([0 1]); yticks([0:1:1]);
text(-0.2*mean(ax.XLim), mean(ax.YLim), {'avg. health'}, 'Rotation',90, 'VerticalAlignment','middle','HorizontalAlignment','center');
% text(mean(ax.XLim),-0.4*mean(ax.YLim),'D', 'HorizontalAlignment','center');
ax = gca; ax.XLabel.String = 'Training Duration (Days)';
% text(ax.XLim(1),ax.YLim(2)*1.15, 'N', 'FontSize', 14);


%%% COMPARISON BETWEEN PULSES OF DIFFERENT DURATION %%%%%%%%%%%%%%%%%%%%%%%
D = 15;
error = [];
theta = 5;
f = Model1Data.TimeResolution;
% f = f*2;
pdindex = [24*15 24 12 4 3 2 1];
cvar = Model1Data.csol.c;

%Steps
temperature1 = linspace(0,43,44*f);
temperature2 = [linspace(22,22,(15-(D))*f) linspace(theta*2+2,theta*2+2,(D)*f) linspace(4,4,29*f)]./22;
temperature = [temperature1; temperature2];
sol = ode23(@(t,y) Acclimfun(t,y,cvar,temperature), Model2Data.tspan, Model2Data.y0);
solpts = deval(0:0.1:43,sol);
StepsHealthFull = solpts(2,:);
StepsAlphaFull = solpts(1,:);
nexttile()
plot(0:0.1:43, StepsHealthFull); hold on;
plot(0:0.1:43, StepsAlphaFull);
plot(temperature1, temperature2);
ax = gca; ax.YLim = [0 1];
ax.Title.String = ['13C Step'];
% text(20, 0.4, ['Mean temp. = ' num2str(mean(temperature2),4)]);

for pd = pdindex;
    %Pulses
    temperature1 = linspace(0,43,44*f);
    % temperature2b = [repmat([linspace(4,4,pd) linspace(22,22,pd)],1,floor(f/(2*pd)*D)) linspace(4,4,29*f)];
    temperature2b = [repmat([linspace(4,4,pd) linspace(22,22,pd)],1,floor(f/(2*pd)*D)) linspace(4,4,29*f)];
    temperature2a = linspace(22,22,length(temperature1)-length(temperature2b));
    temperature2 = [temperature2a temperature2b]./22;
    % temperature2c = [linspace(22,22,(15-(D))*f) linspace(theta*2+2,theta*2+2,(D)*f) linspace(4,4,29*f)]./22;
    % temperature2 = [temperature2(1:100) temperature2c(101:end)];
    temperature = [temperature1; temperature2];
    sol = ode23(@(t,y) Acclimfun(t,y,cvar,temperature), Model2Data.tspan, Model2Data.y0);
    solpts = deval(0:0.01:43,sol);
    PulsesHealthFull = solpts(2,:);   
    PulsesAlphaFull = solpts(1,:);

    nexttile()
    plot(0:0.01:43, PulsesHealthFull); hold on;
    plot(0:0.01:43, PulsesAlphaFull);
    % plot(temperature1, temperature2);
    ax = gca; ax.YLim = [0 1];
    ax.Title.String = ['Pulse duration (hr) = ' num2str(pd)];
    % text(20, 0.4, ['Mean temp. = ' num2str(mean(temperature2),4)]);

    % error = [error sum((PulsesHealthFull - StepsHealthFull).^2)];
end

plot(0:0.1:43, StepsHealthFull, 'LineStyle',':','LineWidth',1); hold on;
plot(0:0.1:43, StepsAlphaFull, 'LineStyle',':','LineWidth',1);
