%%% TO REMOVE%%%%%%%
close all; clear all; clc;
load('data.mat') %Load compiled data

figureRegenVsAcclim = figure('Name', 'Cst Pilot', 'WindowState','maximized');
cmaplines = colormap('lines');
tiledlayout(4,5)

%%% STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile()
time = [linspace(0,0.5,4) linspace(0.5,11,11) linspace(11,28,14)];
temperature = [linspace(22,22,4) linspace(11,11,11) linspace(4,4,14)]./22;
plot(time, temperature, 'Color',cmaplines(1,:), 'LineWidth',2); hold on;
ax = gca; ax.Box = 'off'; ax.YLim(2) = 1.1;
ax.XTick = [0.5 11 22]; ax.XTickLabel = {'0','D', '14'};
ax.YTick = [0 0.5 1]; ax.YTickLabel = {'0','\theta', '22'}; ax.XLim = [0 28];
text(14,0.45, 'Measure health'); text(21.8,0.3,'|');
ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;

nexttile();
ConstantData = [Pilot.Constant12C(:,2) Pilot.Constant15C(:,2) Pilot.Constant18C(:,2) Pilot.Constant22C(:,2)];
h = heatmap(ConstantData');
h.ColorbarVisible = 'off'; h.XLabel = 'D'; h.YLabel = '\theta';
h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''}; h.FontSize = 8;
h.XDisplayLabels([1 6 11]) = {'0', '5', '10'};
h.YDisplayLabels([1 2 3 4]) = {'12', '15', '18', '22'};
% h.Title = 'Health at end of training';

nexttile()
TimeVector = 1:length(ConstantData(:,1));
p1 = plot(TimeVector,ConstantData(:,1), 'Color', [0 0 0], 'LineStyle', ':'); hold on;
p2 = plot(TimeVector,ConstantData(:,2), 'Color', [0 0 0], 'LineStyle', '-.');
p3 = plot(TimeVector,ConstantData(:,3), 'Color', [0 0 0], 'LineStyle', '--');
p4 = plot(TimeVector,ConstantData(:,4), 'Color', [0 0 0]);
ax = gca; ax.Box = 'off'; ax.YLabel.String = 'Health at T=14'; ax.XLabel.String = 'D';
l = legend([p1 p2 p3 p4], {'\theta = 12', '\theta = 15', '\theta = 18', '\theta = 22'});
l.Box = 'on'; l.EdgeColor = 'none'; l.Location = "northwest";

nexttile(); ax = gca; ax.Visible = 'off';
nexttile(); ax = gca; ax.Visible = 'off';

f = gcf;
exportgraphics(f, [f.Name '.pdf'])