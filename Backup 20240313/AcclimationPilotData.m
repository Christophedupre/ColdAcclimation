function AcclimationPilotData(Pilot);
% %%% TO REMOVE%%%%%%%
% close all; clear all; clc;
% load('data.mat') %Load compiled data
%%%%%%%%%%%%%%%%%%%%
figurePilot = figure('Name', 'Pilot Data', 'WindowState','maximized');
cmaplines = colormap('lines');
tiledlayout(4,5)

%%% RAMPS OF VARYING SLOPE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile()
hold on; 
time = linspace(0,672,672)./24; %28 days
temp = [ceil(linspace(22,4,256)/0.625)*(0.625) linspace(4,4,672-256)];
plot(time,temp, 'Color',cmaplines(1,:),'LineWidth',2)

% Plot display properties
ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10; ax.YLim = [0 25]; ax.XLim = [0 max(time)];
ax.YLabel.String = 'Temperature (C)'; ax.XLabel.String = 'Time (days)'; ax.Box = 'off';
ax.XTick = [0 11 21.8]; ax.XTickLabel = {'0','D','14'}; ax.YTick = [0 22]; 
text(14,10, 'Measure health'); text(21.8,6.6,'|');

%%% DATA RAMPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile()
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

[xData, yData] = prepareCurveData(GradientDuration,HealthMeans);
ft = fittype( '1./(1+exp(-(x-of)/tau))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.StartPoint = [200 43]; % [of tau]
[fitresultgrad, gof] = fit( xData, yData, ft, opts );

errorbar(GradientDuration,HealthMeans,HealthSem,...
    'Color','k','Linewidth',1,'LineStyle','none'); hold on;
bar(GradientDuration,HealthMeans,0.4, 'FaceColor', 'k', 'FaceAlpha',0.2, 'EdgeColor', 'none');
for i=1:length(GradientDuration)
    s = scatter(GradientDuration(i), Pilot.dataGradientAll(:,i),'filled', 'MarkerFaceColor', [0 0 0],...
            'AlphaData', 0.1*ones(1,length(GradientDuration(i))), 'MarkerFaceAlpha','flat', 'SizeData', 15, 'jitter','on', 'jitterAmount',0.17);
end
p = plot(linspace(0,15,16),feval(fitresultgrad,linspace(0,15,16)),'color',cmaplines(1,:),'LineWidth',2);

ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10; box off;
ax.YLim = [0 1.2]; ax.YTick = [0:1:1]; ax.YLabel.String = 'Health at T=14'; ax.XLim = [-1 16];
ax.XTick = 0:5:15; ax.XLabel.String = 'D (days)';
text(mean(ax.XLim)*0.32,1.1, ['n = ' num2str(length(Pilot.dataGradientAll(:,1))) ' for each TR'], 'FontSize', 10);

nexttile(); ax = gca; ax.Visible = 'off';
nexttile(); ax = gca; ax.Visible = 'off';
nexttile(); ax = gca; ax.Visible = 'off';

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

%%% PULSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile()
first = 30; last = 53;
TemperatureSteps = repmat([linspace(20,20,5) linspace(4,4,6)],1,6);
TemperatureSteps2 = linspace(4,4,length(TemperatureSteps(first:last)));
TemperatureSteps(first:last) = NaN;
TimeSteps = linspace(0,15,length(TemperatureSteps));
TimeSteps2 = TimeSteps(first:last);
TemperatureProbe = linspace(4,4,24);
TimeProbe = linspace(15,34,length(TemperatureProbe));
plot(TimeSteps,TemperatureSteps,'color',[cmaplines(1,:) 0.8],'LineWidth',2); hold on;
plot(TimeProbe, TemperatureProbe,'color',[cmaplines(2,:) 0.8],'LineWidth',2);
plot(TimeSteps2,TemperatureSteps2, 'k--','LineWidth',1); hold on;
xticks([15 34]); xticklabels({'0';'1'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)'])
box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
legend('Training', 'Test', 'box', 'off', 'FontSize', 8)

nexttile();
iterations = 10;
TimeColdShocks = 1:1:48*iterations;
TemperatureColdShocks24hr = repmat([linspace(20,20,24) linspace(4,4,24)],1,iterations);
TemperatureColdShocks8hr = repmat([linspace(22,22,16) linspace(4,4,8)],1,iterations*2);
% TemperatureColdShocks1hr = repmat([linspace(22,22,1) linspace(4,4,1)],1,iterations*24);
plot(TimeColdShocks, TemperatureColdShocks24hr, 'color', [cmaplines(1,:) 0.8],'LineWidth', 2); hold on;
ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10; ax.Box = 'off';
ax.XLim = [0 96]; ax.XLabel.String = 'Time (hrs)';
ax.YTick = [0 4 20]; ax.YLabel.String = ['T (' char(176) 'C)'];
ax.XTick = [24 48 72]; ax.XTickLabel = {0; 'TP'; '2TP'};

nexttile()
hold on;
plot(Pilot.AcclimLevelCRS24hr(:,1)*2,Pilot.AcclimLevelCRS24hr(:,2), 'color', [cmaplines(1,:) 0.6], 'LineWidth', 2);
plot(Pilot.AcclimLevelCRS8hr(:,1)*1,Pilot.AcclimLevelCRS8hr(:,2), 'color', [cmaplines(2,:) 0.6], 'LineWidth', 2);
plot(Pilot.AcclimLevelCRS1hr(:,1)*2/24,Pilot.AcclimLevelCRS1hr(:,2), 'color', [cmaplines(3,:) 0.6], 'LineWidth', 2);
plot(Pilot.AcclimLevelCRS0_5hr(:,1)*1/24,Pilot.AcclimLevelCRS0_5hr(:,2), 'color', [cmaplines(4,:) 0.6], 'LineWidth', 2);
plot(Pilot.AcclimLevelCRS0_1hr(:,1)*1/120,Pilot.AcclimLevelCRS0_1hr(:,2), 'color', [cmaplines(5,:) 0.6], 'LineWidth', 2);
lgd = legend('24', '8', '1', '0.5', '0.1', 'Interpreter', 'latex', 'FontSize', 7,... 
    'EdgeColor', 'none','NumColumns',1, 'box', 'on');
lgd.Title.String = '$TP$(hr)';
lgd.Location = 'east';
% lgd.Position = [0.84 0.655 0.07 0.27];
set(lgd.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));
ylabel(['Health']); yticks([0:1]); box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
xlabel('Training duration (days)'); xlim([0 25])
lgd.Position(3) = lgd.Position(3) - 0.03;
lgd.Position(1) = lgd.Position(1) + 0.03;

nexttile()
plot(Pilot.AcclimPulsesSym(:,1), Pilot.AcclimPulsesSym(:,2),'o', 'color', [cmaplines(1,:) 0.1], 'LineWidth', 2, 'Marker', '*');
hold on; box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10; ylim([0 15]); xlim([0 25]);
ylabel('Time to H=1 (D)'); xlabel('TP (hrs)')
ax.YLabel.FontSize = 10;

nexttile()
ax = gca; ax.Visible = 'off';

%%% ASYMMETRIC PULSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nexttile()
first = 30; last = 53;
TemperatureSteps = repmat([linspace(20,20,5) linspace(4,4,6)],1,6);
TemperatureSteps2 = linspace(4,4,length(TemperatureSteps(first:last)));
TemperatureSteps(first:last) = NaN;
TimeSteps = linspace(0,15,length(TemperatureSteps));
TimeSteps2 = TimeSteps(first:last);
TemperatureProbe = linspace(4,4,24);
TimeProbe = linspace(15,34,length(TemperatureProbe));
plot(TimeSteps,TemperatureSteps,'color',[cmaplines(1,:) 0.8],'LineWidth',2); hold on;
plot(TimeProbe, TemperatureProbe,'color',[cmaplines(2,:) 0.8],'LineWidth',2);
plot(TimeSteps2,TemperatureSteps2, 'k--','LineWidth',1); hold on;
xticks([15 34]); xticklabels({'0';'1'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)'])
box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
legend('Training', 'Test', 'box', 'off', 'FontSize', 10)

nexttile()
iterations = 2;
TemperatureColdShocks1hr = repmat([linspace(20,20,80) linspace(4,4,40)],1,iterations);
TimeColdShocks = 1:1:length(TemperatureColdShocks1hr);
plot(TimeColdShocks, TemperatureColdShocks1hr, 'color', [cmaplines(1,:) 0.8],'LineWidth', 2); hold on;
ylabel(['T (' char(176) 'C)']); xticks([80 120 200]); xticklabels({0; 'TP'; 2}); yticks([0 4 20]); box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
xlabel('Time (hrs)');

nexttile()
hold on;

plot(Pilot.AcclimLevelCRAS1hr1hr(:,1)*2/24,Pilot.AcclimLevelCRAS1hr1hr(:,2), 'color', [cmaplines(1,:) 0.6], 'LineWidth', 2);
plot(Pilot.AcclimLevelCRAS0_5hr1_5hr(:,1)*2/24,Pilot.AcclimLevelCRAS0_5hr1_5hr(:,2), 'color', [cmaplines(2,:) 0.6], 'LineWidth', 2);
plot(Pilot.AcclimLevelCRAS0_2hr1_8hr(:,1)*2/24,Pilot.AcclimLevelCRAS0_2hr1_8hr(:,2), 'color', [cmaplines(3,:) 0.6], 'LineWidth', 2);
plot(Pilot.AcclimLevelCRAS0_1hr1_9hr(:,1)*2/24,Pilot.AcclimLevelCRAS0_1hr1_9hr(:,2), 'color', [cmaplines(4,:) 0.6], 'LineWidth', 2);
plot(Pilot.AcclimLevelCRAS0_05hr1_95hr(:,1)*2/24,Pilot.AcclimLevelCRAS0_05hr1_95hr(:,2), 'color', [cmaplines(5,:) 0.6], 'LineWidth',2);
plot(Pilot.AcclimLevelCRAS0_025hr1_975hr(:,1)*2/24,Pilot.AcclimLevelCRAS0_025hr1_975hr(:,2), 'color', [cmaplines(6,:) 0.6], 'LineWidth',2);
plot(Pilot.AcclimLevelCRAS0hr2hr(:,1)*2/24,Pilot.AcclimLevelCRAS0hr2hr(:,2), 'color', [cmaplines(7,:) 0.6], 'LineWidth',2);

ylabel(['Health']); yticks([0:1]); box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
xlim([0 12]); ylim([0 1.2]);
% xticks([linspace(0,14,15)]);
xlabel('Training duration (days)')
lgdPulsesAsym = legend('1','0.5', '0.2', '0.1', '0.05', '0.025', '0', 'Interpreter', 'latex', 'FontSize', 6,... 
    'EdgeColor', 'none','NumColumns',1, 'box', 'on');
lgdPulsesAsym.Title.String = '$TP$(hr)';
% lgdPulsesAsym.Position(1) = lgdPulsesAsym.Position(1)*2;
% lgdPulsesAsym.Position = [0.83 0.665 0.07 0.25];
lgdPulsesAsym.Location = 'east';
set(lgdPulsesAsym.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

lgdPulsesAsym.Position(3) = lgdPulsesAsym.Position(3) - 0.03;
lgdPulsesAsym.Position(1) = lgdPulsesAsym.Position(1) + 0.03;

nexttile()
hold on;
A = Pilot.AcclimLevelCRAS0hr2hrExtTest(1:end-4,:); B = Pilot.AcclimLevelCRAS0hr2hrExtTest(end-3:end,:);
[h,p] = ttest2(A,B);
p1 = plot(0:5,mean(A),'color', [cmaplines(7,:) 0.6],'LineWidth',2);
scatter(0:5,A','filled', 'MarkerFaceColor', cmaplines(7,:),...
        'AlphaData', 0.1*ones(1,length(0:5)), 'MarkerFaceAlpha',...
        'flat', 'SizeData', 15, 'jitter','on', 'jitterAmount',0.12);
p2 = plot(0:5,mean(B),'color', [cmaplines(1,:) 0.6],'LineWidth',2);
scatter(0:5,B','filled', 'MarkerFaceColor', cmaplines(1,:),...
        'AlphaData', 0.1*ones(1,length(0:5)), 'MarkerFaceAlpha',...
        'flat', 'SizeData', 15, 'jitter','on', 'jitterAmount',0.12);
text(4.5,1,'*p<0.05','color','k','FontSize',10);
plot([1 2 3 4 5], [5.8 5.8 5 5 5],'k*','MarkerSize', 4)
% text(1,6.5,'0 is not fully acclimated')

ylabel(['Health']); yticks([0:1]); box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
xlim([0 6]); xticks(0:5); ylim([0 1.2]);
% xticks([linspace(0,14,15)]);
xlabel('Extended test duration (days)')
lgdPulsesAsym5D = legend([p1 p2],'0 (n=12)', '+Ctrl (n=4)', 'Interpreter', 'latex', 'FontSize', 8,... 
    'EdgeColor', 'none','NumColumns',1, 'box', 'on');
lgdPulsesAsym5D.Location = 'southwest';
% lgdPulsesAsym5D.Title.String = '$TP$(hr)';
set(lgdPulsesAsym5D.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));


a = gcf;
exportgraphics(a, [a.Name '.pdf'])

%%% GRAPHS FOR TALK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% figGraphsForTalk = figure('Name', 'Graphs for talk', 'renderer', 'painters', 'WindowState', 'maximized');
% % figGraphsForTalk.Position = [10 0 750 780];
% nrows = 3; ncols = 4;
% tiledlayout(nrows,ncols)
% 
% % Theta = 10C; D = 0
% nexttile(2+(1-1)*ncols)
% TemperatureSteps = [linspace(0,14,141)];
% AImax = 1; offset = 0; tau = 1;
% AI = exp(-(TemperatureSteps-offset)./tau);
% plot(TemperatureSteps, AI, 'LineWidth',1, 'Color', [0 0 0])
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% xlabel('Time (Days)')
% ylabel('Health')
% ylabel({'\bf{\theta = 10C}','\rm {Health}'})
% ylim([0 AImax])
% title('D = 0')
% 
% % Theta = 10C; D = 14
% nexttile(3+(1-1)*ncols)
% TemperatureSteps = [linspace(0,14,141)];
% AImax = 1; offset = 0; tau = 10000;
% AI = exp(-(TemperatureSteps-offset)./tau);
% plot(TemperatureSteps, AI, 'LineWidth',1, 'Color', [0 0 0])
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% xlabel('Time (Days)')
% ylabel('Health')
% ylim([0 AImax])
% title('D = 14')
% 
% % Theta = 10C, AI vs D
% nexttile(4+(1-1)*ncols)
% TemperatureSteps = [linspace(0,14,141)];
% AImax = 1; offset = 7.5; tau = 1;
% AI = AImax./(1+exp(-(TemperatureSteps-offset)./tau));
% plot(TemperatureSteps, AI, 'LineWidth',1, 'Color', [0 0 0])
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% xlabel('D')
% ylabel('\alpha')
% ylim([0 AImax])
% title('\alpha vs D')
% 
% % Theta = 22C; D = 0
% nexttile(2+(2-1)*ncols)
% TemperatureSteps = [linspace(0,14,141)];
% AImax = 1; offset = 0; tau = 1;
% AI = exp(-(TemperatureSteps-offset)./tau);
% plot(TemperatureSteps, AI, 'LineWidth',1, 'Color', [0 0 0])
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% xlabel('Time (Days)')
% ylabel({'\bf{\theta = 22C}','\rm {Health}'})
% ylim([0 AImax])
% 
% % Theta = 22C; D = 14
% nexttile(3+(2-1)*ncols)
% TemperatureSteps = [linspace(0,14,141)];
% AImax = 1; offset = 0; tau = 1;
% AI = exp(-(TemperatureSteps-offset)./tau);
% plot(TemperatureSteps, AI, 'LineWidth',1, 'Color', [0 0 0])
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% xlabel('Time (Days)')
% ylabel('Health')
% ylim([0 AImax])
% 
% % Theta = 22C, AI vs D
% nexttile(4+(2-1)*ncols)
% TemperatureSteps = [linspace(0,14,141)];
% AImax = 1; offset = 75; tau = 1;
% AI = AImax./(1+exp(-(TemperatureSteps-offset)./tau));
% plot(TemperatureSteps, AI, 'LineWidth',1, 'Color', [0 0 0])
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% xlabel('D')
% ylabel('\alpha')
% ylim([-0.015 AImax])
% 
% % Theta = 18C; D = 0
% nexttile(2+(3-1)*ncols)
% TemperatureSteps = [linspace(0,14,141)];
% AImax = 1; offset = 0; tau = 1;
% AI = exp(-(TemperatureSteps-offset)./tau);
% plot(TemperatureSteps, AI, 'LineWidth',1, 'Color', [0 0 0])
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% xlabel('Time (Days)')
% % ylabel({'$\bf{\theta = 10C}$','A.I.'},'interpreter','latex')
% ylabel({'\bf{\theta = 18C}','\rm {Health}'})
% ylim([0 AImax])
% 
% % Theta = 18C; D = 14
% nexttile(3+(3-1)*ncols)
% TemperatureSteps = [linspace(0,14,141)];
% AImax1 = 1; offset1 = 6; tau1 = 1.5;
% AImax2 = 1; offset2 = 8; tau2 = 1.5;
% % AI = exp(-(TemperatureSteps-offset1)./tau1) + AImax./(1+exp(-(TemperatureSteps-offset2)./tau2));
% AI1 = AImax1./(1+exp((TemperatureSteps-offset1)./tau1));
% AI2 = AImax2./(1+exp(-(TemperatureSteps-offset2)./tau2));
% AI = AI1 + AI2;
% % plot(TemperatureSteps, AI1, 'LineWidth',1, 'Color', [1 0 0]); hold on;
% % plot(TemperatureSteps, AI2, 'LineWidth',1, 'Color', [0 1 0]);
% plot(TemperatureSteps, AI, 'LineWidth',1, 'Color', [0 0 0])
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% xlabel('Time (Days)')
% ylabel('Health')
% ylim([0 AImax])
% 
% % Theta = 18C, AI vs D
% nexttile(4+(3-1)*ncols)
% TemperatureSteps = [linspace(0,14,141)];
% AImax = 0.8; offset = 10; tau = 2;
% AI = AImax./(1+exp(-(TemperatureSteps-offset)./tau));
% plot(TemperatureSteps, AI, 'LineWidth',1, 'Color', [0 0 0])
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% xlabel('D')
% ylabel('\alpha')
% ylim([0 1])
% 
% % Experiment diagram
% nexttile(1+(1-1)*ncols)
% TemperatureCst = [linspace(12,12,99) 4];
% TemperatureTest = linspace(4,4,100);
% TimeCst = [linspace(1,length(TemperatureCst),length(TemperatureCst)-1) length(TemperatureCst)];
% TimeTest = linspace(length(TemperatureCst),200,100);
% p1 = plot(TimeCst,TemperatureCst,'color',[cmaplines(1,:) 0.8],'LineWidth',2); hold on;
% p2 = plot(TimeTest,TemperatureTest,'color',[cmaplines(2,:) 0.8],'LineWidth',2);
% xlim([0 200])
% xticks([100 200]); xticklabels({'0';'28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 20])
% yticks([0 12 20]); yticklabels({'0';'\theta';'20'});
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
% legend([p1 p2],'Training', 'Test', 'box', 'off', 'FontSize', 10)
% 
% nexttile(1+(3-1)*ncols)
% box off;axis off
% % xlabel('$\alpha(D,\theta) = \frac{1}{1+e^{-\frac{D(\theta)+off}{tau(\theta)}}}$','interpreter','latex', 'box', 'off')
% text(-0.5,1.5,{'Acclimation? Yes','How? Temperature exposure','$\theta = $Training Temp.','$D = $Training Dur.','$\alpha(D,\theta) = \frac{1}{1+e^{-\frac{D-of(\theta)}{\tau(\theta)}}}$' ...
%     ,'$of(\theta) = 6$','$\tau(\theta)\epsilon[2,+\infty]$', '$\tau(\theta)=1.6+\frac{\theta-10}{12},\theta \epsilon [10,16]$','Process?','$\alpha(t+1)=\alpha(t) + k_pe(t)+k_i\int_0^t e(\tau)d\tau$','$k_p$ = proportionality constant','$k_i$ = integration constant','$e(t)$ = error = $\alpha_{target} - \alpha(t)$', '$\alpha_{target} = 1-\frac{\theta-4}{18}, \theta \epsilon [4,22]$'},'interpreter','latex', 'FontSize',15)
% 
% a = gcf;
% exportgraphics(a, [a.Name '.pdf'])



end