% function SlowDecrease()

%%% TO REMOVE%%%%%%%
% close all; clear all; clc;
load('data.mat') %Load compiled data

figureSlowDecrease = figure('Name', 'Slow Decrease', 'WindowState','maximized');
cmaplines = colormap('lines');
cmapjet = colormap('jet');
load("cmaplist.mat");
tiledlayout(4,5)

Labeltitles = 'ABCDEFGHIJKLM';
Labelcounter = 0;
Labelcounter = Labelcounter+1;

%%% NEGATIVE CONTROLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Negative controls
% nexttile()
% plot([linspace(0,14,15) 14], [linspace(22,22,15) 4],'color',[cmaplines(1,:) 0.8], 'LineWidth',2); hold on;
% plot(linspace(14,28,15), linspace(4,4,15),'color',[cmaplines(2,:) 0.8], 'LineWidth',2);
% xlim([0 28]); xticks([14 28]); xticklabels({'TD';'TD+28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 22])
% yticks([0 22]); yticklabels({'0';'22'})
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% text(ax.XLim(1),ax.YLim(2)*1.15, Labeltitles(Labelcounter), 'FontSize', 14);
nexttile(); hold on;
StepDuration = 1:2:14;
for i = 1:length(StepDuration)
    plot([0 StepDuration(i) StepDuration(i) 42], [22 22 4 4],'Color',cmaplines(round(256/7*i),:),'LineWidth',1)    
end
ax = gca; ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'Temperature (C)';
Labelcounter = Labelcounter+1;

nexttile(); 
h = heatmap(NoPCR');
h.ColorbarVisible = 'off'; h.XLabel = 'T (Days)'; h.YLabel = 'TD';
h.FontSize = 8;
h.Title = [Labeltitles(Labelcounter) '. Negative Controls (NoPCR)']; Labelcounter = Labelcounter+1;
h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};

nexttile();
h = heatmap(JustPCR');
h.ColorbarVisible = 'off'; h.XLabel = 'T (Days)'; h.YLabel = 'TD';
h.FontSize = 8;
h.Title = [Labeltitles(Labelcounter) '. Negative Controls (JustPCR)']; Labelcounter = Labelcounter+1;
h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};

%%% POSITIVE CONTROLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile()
% plot([linspace(0,14,15) 14], [linspace(4,4,15) 4],'color',[cmaplines(1,:) 0.8], 'LineWidth',2); hold on;
% plot(linspace(14,28,15), linspace(4,4,15),'color',[cmaplines(2,:) 0.8], 'LineWidth',2);
% xlim([0 28]); xticks([14 28]); xticklabels({'D';'D+28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 22])
% yticks([0 4 22]); yticklabels({'0';'4';'22'})
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% text(ax.XLim(1),ax.YLim(2)*1.15, Labeltitles(Labelcounter), 'FontSize', 14); Labelcounter = Labelcounter+1;
plot([0 42], [4 4], 'Color',cmaplines(1,:),'LineWidth',1);
ax = gca; ax.YLim = [0 22]; ax.Box = 'off'; ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'Temperature (C)';


nexttile(); 
h = heatmap(PosCtrl');
h.ColorbarVisible = 'off'; h.XLabel = 'T (Days)'; h.YLabel = 'TD';
h.FontSize = 8;
h.Title = [Labeltitles(Labelcounter) '. Positive Controls (=acclimated)']; Labelcounter = Labelcounter + 1;
h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
h.ColorLimits = [0 0.1]; %Otherwise color is normalized and because they all have same health it looks like they all have health of 0.5;

%%% PLOT CONTROLS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Time = linspace(0,length(CstData3D(:,1,1))-1,length(CstData3D(:,1,1)));
p = [];

%%% RAMPS OF VARYING SLOPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot display properties

nexttile(); hold on;
RampDuration = [13 11 9 7 3 1];
for i = 1:length(RampDuration)
    plot([0 14-RampDuration(i) 42], [22 4 4], 'Color',cmaplines(round(256/6*i),:),'LineWidth',1);    
end
ax = gca;
ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'Temperature (C)';
% ax = gca; ax.LineWidth = 1; ax.FontSize = 10; ax.YLim = [0 22]; ax.YTick = [4 22]; 
% ax.YTickLabel = {'4';'22'};ax.XLim = [0 28]; ax.XTick = [0 14 21 28]; ax.XTickLabel = {'0'; '14'; '[...]'; '42'};
% text(ax.XLim(1),ax.YLim(2)*1.15, Labeltitles(Labelcounter), 'FontSize', 14); Labelcounter = Labelcounter+1;


%%% DATA RAMPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(); hold on;
GradientDuration = fliplr([10/0.5*16/24, 10/0.625*16/24, 10/0.75*16/24, 10/1*16/24, 10/2*16/24, 0.5/18*16/24]); % Duration(days) =  slope*deltatheta*(day/hour) = 10/0.5 * 16 * (1/24);
Pilot.dataGradientAll = [Pilot.NotHabituated(:,2) Pilot.AcclimGradient2c(:,2) Pilot.AcclimGradient1c(:,2) Pilot.AcclimGradient0_75c(:,2) Pilot.AcclimGradient0_625c(1:10,2) Pilot.AcclimGradient0_5c(1:10,2)];

HealthMeans = [mean(Pilot.NotHabituated(:,2)), mean(Pilot.AcclimGradient2c(:,2)),...
    mean(Pilot.AcclimGradient1c(:,2)),mean(Pilot.AcclimGradient0_75c(:,2)),...
    mean(Pilot.AcclimGradient0_625c(:,2)), mean(Pilot.AcclimGradient0_5c(:,2))];
HealthSem = [std(Pilot.NotHabituated(:,2))./sqrt(length(Pilot.NotHabituated(:,2))),...
    std(Pilot.AcclimGradient2c(:,2))./sqrt(length(Pilot.AcclimGradient2c(:,2))),...
    std(Pilot.AcclimGradient1c(:,2))./sqrt(length(Pilot.AcclimGradient1c(:,2))),...
    std(Pilot.AcclimGradient0_75c(:,2))./sqrt(length(Pilot.AcclimGradient0_75c(:,2))),...
    std(Pilot.AcclimGradient0_625c(:,2))./sqrt(length(Pilot.AcclimGradient0_625c(:,2))),...
    std(Pilot.AcclimGradient0_5c(:,2))./sqrt(length(Pilot.AcclimGradient0_5c(:,2)))];

GradientDurationS = repmat(GradientDuration,size(Pilot.dataGradientAll,1),1);
GradientDurationS = GradientDurationS + (rand(size(GradientDurationS,1),size(GradientDurationS,2))-0.5).*2;
Data = Pilot.dataGradientAll + (rand(size(Pilot.dataGradientAll,1),size(Pilot.dataGradientAll,2))-0.5)./10;

for i = 1:size(GradientDurationS,2)
    s1 = scatter(GradientDurationS(:,i), Data(:,i),'MarkerFaceColor',...
        cmaplines(round(256/6*i),:), 'MarkerEdgeColor','none','SizeData',5);
end

[xData, yData] = prepareCurveData(GradientDuration(1,:),HealthMeans);
ft = fittype( '1./(1+exp(-(x-of)/tau))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = [200 43]; % [of tau]
[fitresultgrad, gof] = fit( xData, yData, ft, opts );
% plot(GradientDuration(1,:), HealthMeans)
b = bar(GradientDuration(1,:), HealthMeans);
b.EdgeColor = 'none'; b.FaceColor = cmaplines(1,:); b.FaceAlpha = 0.5; b.BarWidth = 0.5;

% fill([GradientDuration(1,:), fliplr(GradientDuration(1,:))], [HealthMeans+HealthSem fliplr(HealthMeans-HealthSem)],cmaplines(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'w');

ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10; box off;
ax.YLim = [0 1.2]; ax.YTick = [0:1:1]; ax.YLabel.String = 'Health at T=14'; ax.XLim = [-1 16];
ax.XTick = 0:5:15; ax.XLabel.String = 'D (days)';
text(mean(ax.XLim)*0.5,1.1, ['n = ' num2str(length(Pilot.dataGradientAll(:,1))) ' each'], 'FontSize', 10, 'HorizontalAlignment','left');
% text(ax.XLim(1),ax.YLim(2)*1.15, Labeltitles(Labelcounter), 'FontSize', 14);

% nexttile(); %ax = gca; ax.Visible = 'off';
% plot(GradientDuration, survival,'color',[cmaplines(1,:) 0.8], 'LineWidth',2);
% ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10; box off;
% ax.YLim = [0 110]; ax.XTick = 0:5:15; ax.XLabel.String = 'D (days)';
% ax.YLabel.String = 'Survival (%)';
% text(ax.XLim(1),ax.YLim(2)*1.15, Labeltitles(Labelcounter), 'FontSize', 14); Labelcounter = Labelcounter+1;

%%% TRUNCATED RAMPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ramp
nexttile(); hold on;
RampDuration1 = 1:2:14;
RampDuration2 = 22-RampDuration1.*18/RampDuration1(end);
for i = 1:length(RampDuration1)
    plot([0 RampDuration1(i) RampDuration1(i) 14 42], [22 RampDuration2(i) 4 4 4],'LineWidth',1);
end
ax = gca; ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'Temperature (C)';
% 
% nexttile()
% plot([linspace(0,7,15)], [linspace(22,13,15)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2); hold on;
% plot([linspace(7,14,8)], [linspace(13,4,8)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2, 'LineStyle',':');
% % plot([linspace(7,14,8)], [linspace(4,4,8)],'color',[cmap(2,:) 0.8], 'LineWidth',2);
% plot([linspace(7,7,8)], [linspace(13,4,8)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2);
% plot(linspace(7,21,15), linspace(4,4,15),'color',[cmaplines(2,:) 0.8], 'LineWidth',2);
% plot(linspace(21,28,8), linspace(4,4,8),'color',[cmaplines(2,:) 0.8], 'LineWidth',2, 'LineStyle',':');
% xlim([0 28]); xticks([7 21]); xticklabels({'D';'D+28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 22])
% yticks([4 22]); yticklabels({'4';'22'})
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% text(ax.XLim(1),ax.YLim(2)*1.15, Labeltitles(Labelcounter), 'FontSize', 14); Labelcounter = Labelcounter+1;

nexttile()
h = heatmap(Ramp');
h.ColorbarVisible = 'off'; h.XLabel = 'T (Days)'; h.YLabel = 'TD';
h.FontSize = 8;
h.Title = [Labeltitles(Labelcounter) '. Ramp']; Labelcounter = Labelcounter + 1;
h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};

%%% PULSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nexttile()
% plot([linspace(0,3.5,3) 3.5 linspace(3.5,7,4) linspace(7,10.5,3) 10.5 linspace(10.5,14,4)], [linspace(22,22,3) 4 linspace(4,4,4) linspace(22,22,3) 4 linspace(4,4,4)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2); hold on;
% plot(linspace(14,28,15), linspace(4,4,15),'color',[cmaplines(2,:) 0.8], 'LineWidth',2);
% xlim([0 28]); xticks([14 28]); xticklabels({'D';'D+28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 22])
% yticks([4 22]); yticklabels({'4';'22'})
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% text(ax.XLim(1),ax.YLim(2)*1.15, Labeltitles(Labelcounter), 'FontSize', 14); Labelcounter = Labelcounter + 1;

nexttile(); hold on;
PulsesHighTimePoints = (2:2:24*14)/24;
PulsesLowTimePoints = ((2:2:24*14)-2)/24;
PulsesTimePoints = sort([PulsesHighTimePoints PulsesLowTimePoints]);
PulsesTempPoints = repmat([4 4 22 22],1,length(PulsesTimePoints)./4);
for i=1:48:length(PulsesTimePoints)-48
    plot(PulsesTimePoints(1,i:i+48), PulsesTempPoints(1,i:i+48));
end
plot([12 42], [4 4], 'Color', cmaplines(1,:));
plot([20 20 35],[20 10 10],'Color','k');
plot([20 22 22 24 24 26 26 28 28 30 30 32 32],[10 10 20 20 10 10 20 20 10 10 20 20 10],'Color','k');
text(19,8,'0','FontSize',8)
text(31,8,'0.5','FontSize',8)
ax = gca; ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'Temperature (C)';

nexttile()
h = heatmap(Pulses0_2');
h.ColorbarVisible = 'off'; h.XLabel = 'T (Days)'; h.YLabel = 'TD';
h.FontSize = 8;
h.Title = [Labeltitles(Labelcounter) '. T_{Pulse} = 0']; Labelcounter = Labelcounter + 1;
h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};

nexttile()
h = heatmap(Pulses1_1');
h.ColorbarVisible = 'off'; h.XLabel = 'T (Days)'; h.YLabel = 'TD';
h.FontSize = 8;
h.Title = [Labeltitles(Labelcounter) '. T_{Pulse} = 2']; Labelcounter = Labelcounter + 1;
h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
h.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
h.ColorbarVisible = 'on';

nexttile(); ax = gca; ax.Visible = 'off';
nexttile(); ax = gca; ax.Visible = 'off';


colormap(cmaplist.cmapRedBlue1);

f = gcf;
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(f, [f.Name '.pdf'])

% end