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
cmaplines = colormap('lines');
clf;
%%% COMPARE HEATMAPS OF MODEL VS EXPERIMENTAL DATA %%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Comparison Model vs Experiment','WindowState','maximized','Renderer','painters')
tiledlayout(6,5)
% Labels = ['BCDEFGHIJKLMNOPQ'];
Labels = ['                '];

% Experimental Data
for i = [1 2 3 5 7] 
    nexttile();
    h = heatmap(squeeze(CstData3D(:,:,i))');
    h.ColorbarVisible = 'off'; %h.XLabel = 'T (Days)'; %h.YLabel = 'Training Duration';
    h.FontSize = 8;
    % h.Title = [Labels(i) ' T= ' num2str(i*2+2) ' ^oC'];
    h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
    if i == 1
        h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
        h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
        h.XLabel = 'T (Days)'; 
        h.YLabel = 'Training Duration';
    end
    % h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    % h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    % survival = [survival sum(CstData3D(:,:,i)>0.1,2);];
end
% Model Data (without noise)
for i = [1 2 3 5 7] 
    nexttile();
    h = heatmap(squeeze(Model1Data.HealthTruncated(:,:,i))');
    h.ColorbarVisible = 'off'; %h.XLabel = 'T (Days)'; %h.YLabel = 'Training Duration';
    h.FontSize = 8;
    % h.Title = [Labels(i) ' T= ' num2str(i*2+2) ' ^oC (model)'];
    h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
    % h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    % h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    % survival = [survival sum(CstData3D(:,:,i)>0.1,2);];
end
% Model Data (with noise)
for i = [1 2 3 5 7] 
    nexttile();
    h = heatmap(squeeze(Model2Data.HealthTruncated(:,:,i))');
    h.ColorbarVisible = 'off'; %h.XLabel = 'T (Days)'; %h.YLabel = 'Training Duration';
    h.FontSize = 8;
    % h.Title = [Labels(i) ' T= ' num2str(i*2+2) ' ^oC (model)'];
    h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
    % h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    % h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    % survival = [survival sum(CstData3D(:,:,i)>0.1,2);];
end

% figure('Name','Comparison Model vs Experiment 2','WindowState','maximized','Renderer','painters')
% tiledlayout(4,5)

% Experimental Data
for i = [11 13 14 15 16]
    nexttile();
    h = heatmap(squeeze(CstData3D(:,:,i))');
    h.ColorbarVisible = 'off'; %h.XLabel = 'T (Days)'; %h.YLabel = 'Training Duration';
    h.FontSize = 8;
    % h.Title = [Labels(i) ' T= ' num2str(i*2+2) ' ^oC'];
    h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
    % h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    % h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    % survival = [survival sum(CstData3D(:,:,i)>0.1,2);];
end
% Model Data (without noise)
for i = [11 13 14 15 16]
    nexttile();
    h = heatmap(squeeze(Model1Data.HealthTruncated(:,:,i))');
    h.ColorbarVisible = 'off'; %h.XLabel = 'T (Days)'; %h.YLabel = 'Training Duration';
    h.FontSize = 8;
    % h.Title = [Labels(i) ' T= ' num2str(i*2+2) ' ^oC (model)'];
    h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
    % h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    % h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    % survival = [survival sum(CstData3D(:,:,i)>0.1,2);];
end
% Model Data (with noise)
for i = [11 13 14 15 16]
    nexttile();
    h = heatmap(squeeze(Model2Data.HealthTruncated(:,:,i))');
    h.ColorbarVisible = 'off';
    h.FontSize = 8;
    % h.Title = [Labels(i) ' T= ' num2str(i*2+2) ' ^oC (model)'];
    h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
    % if i == 16
    %     h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    %     h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    %     h.XLabel = 'T (Days)'; 
    %     h.YLabel = 'Training Duration';
    % end
    % survival = [survival sum(CstData3D(:,:,i)>0.1,2);];
end
colormap(cmaplist.cmapRedBlue1)

% % Other descriptors of the model's performance
% % Show effect of modifying training set on model performance
% % The model fails at taking variation between animals into account. 
% % The model succeeds at explaining the effect of training on acclimation.

f = gcf;
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(f, [f.Name '.pdf'])
