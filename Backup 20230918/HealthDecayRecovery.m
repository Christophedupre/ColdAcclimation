function HealthDecayRecovery(Pilot, Deacclimation, CstData3D)
%%
% close all; clear all; clc;
load('data.mat')
%%% LOAD DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
HRN = Pilot.DecayRegen.HRN; %Health Recovery Naive
HRNmean = mean(HRN); HRNsem = std(HRN)./sqrt(size(HRN,1));
HDN = Pilot.DecayRegen.HDN; %Health Decay Naive
HDNmean = mean(HDN); HDNsem = std(HDN)./sqrt(size(HDN,1));
HRA = Pilot.DecayRegen.HRA; %Health Recovery Acclimated
HRAmean = mean(HRA); HRAsem = std(HRA)./sqrt(size(HRA,1));
HDA = Pilot.DecayRegen.HDA; %Health Decay Acclimated
HDAmean = mean(HDA); HDAsem = std(HDA)./sqrt(size(HDA,1));
HDRN = [HDN(:,1:end-1) HRN]; % Health Decay and Recovery Naive
HDRNmean = mean(HDRN); HDRNsem = std(HDRN)./sqrt(length(HDRN));
HDRA = [HDA(:,1:end-1) HRA]; % Health Decay and Recovery Acclimated
HDRAmean = mean(HDRA); HDRAsem = std(HDRA)./sqrt(length(HDRA));

% figtotal = figure('Name', 'Cold-induced injury', 'WindowState', 'maximized','Renderer','painters');
figtotal = figure('Name', 'Cold-induced injury', 'WindowState', 'maximized');
cmaplines = colormap('lines');
tiledlayout(3,4)

%%% HEALTH INDICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Labels = [{'H'} {'I'} {'J'} {'K'} {'L'} {'M'}];
% for j=0:5
%     i = 5-j;
%     nexttile();
%     cropfactor = 4;
%     I = imread(['Health' num2str(i) '.jpg']);
%     imin = round(size(I,1)./cropfactor);
%     imax = round(size(I,1)./cropfactor*(cropfactor-1));    
%     I = I(imin:imax,:);
%     imshow(I);
%     text(1,100,num2str(i/5))
%     ax = gca;
%     text(ax.XLim(1),-ax.YLim(2)*0.1, Labels(j+1), 'FontSize', 14);
% end

%%% SURVIVAL OF NAIVE ANIMALS MOVED TO 4C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cartoon
nexttile(); hold on;
yyaxis right; ax = gca; ax.YColor = [0 0 0];
TempProtocol = [22 22 4 4];
p0 = plot([-6 0 0 28], TempProtocol,'color','k','LineWidth',2);
yticks([4 22]); yticklabels({['4'];['22']}); %xticklabels([0 15 30]);
xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([-2 24]);
box off; ax2 = gca; ax2.LineWidth = 1.5; ax2.FontSize = 10; ax2.Color = [1 1 1 0.5];

% Data
yyaxis left; ax = gca; ax.YColor = [0 0 0];;

A = Deacclimation(29:34,:);
B = CstData3D(:,:,10); B(end,:) = [];
C = [A;B';JustPCR(1:end-1,:)';NoPCR(1:end-1,:)']; %Pool data from naive animals and animals from step experiment trained at 22C;
HealthNaive = C;
HealthNaiveAverage = mean(HealthNaive(:,:));
HealthNaiveSEM = std(HealthNaive(:,:))./sqrt(size(HealthNaive,1));

TimeHealth = 0:1:27;
TimeSEM = [TimeHealth fliplr(TimeHealth)];

inBetweenSEM = [HealthNaiveAverage + HealthNaiveSEM fliplr(HealthNaiveAverage - HealthNaiveSEM)];

survival = sum(HealthNaive>0,1);
survival = survival./max(survival)*100;
fill(TimeSEM, inBetweenSEM(:,:), cmaplines(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'w');

p2 = plot(TimeHealth, HealthNaiveAverage(:,:), 'Color', cmaplines(1,:), 'LineWidth', 2, 'LineStyle','-');
plot([-6 0], [1 1], 'Color', cmaplines(1,:), 'LineWidth', 2); %Add health for the last 5 days before the experiment

HealthNaive = [ones(51,5) HealthNaive]; %Add health for the last 5 days before the 
HealthNaive = HealthNaive + (rand(size(HealthNaive))-0.5)./4;
TimeHealthNaive = -5:1:27;
TimeHealthNaive = repmat(TimeHealthNaive,size(HealthNaive,1),1);
TimeHealthNaive = TimeHealthNaive + rand(size(TimeHealthNaive))./1;
p4 = scatter(TimeHealthNaive,HealthNaive,'filled', 'MarkerFaceColor', cmaplines(1,:),...
        'SizeData', 1);
p1 = plot(TimeHealth, survival/100, 'Color', cmaplines(2,:), 'LineWidth', 2, 'LineStyle','-');
plot([-10 40], [0 0], 'LineStyle','--', 'Color','k', 'LineWidth',1);


ylim([-0.25 1.2]); yticks(0:1:1); 
xlabel('Time (days)');
ylabel('Health, Survival'); box('off')
ax = gca; ax.XLim = [-6.2 30]; %ax.XTick = 0:10:30;
ax.LineWidth = 1.5; ax.FontSize = 10;

lgd1 = legend([p0 p1 p4(1) p2],'Temperature' ,'Survival','Raw','Average', 'box', 'off', 'FontSize',8);
lgd1.NumColumns = 2; lgd1.FontSize = 6;
% lgd1.Position(3) = lgd1.Position(3)/2; lgd1.Position(1) = lgd1.Position(1) + 0.01;
text(ax.XLim(1),ax.YLim(2)*1.1, 'O', 'FontSize', 14);

% nexttile(); ax = gca; ax.Visible = 'off';
% nexttile(); ax = gca; ax.Visible = 'off';

%%% PLOT DECAY AND RECOVERY CURVES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cartoon
nexttile(); hold on;
% ax2 = axes('Position',[.3 .28 .31 .13]);
% yyaxis right; ax = gca; ax.YLabel.String = 'Survival (%)'; ax.YColor = [0 0 0];
% TempProtocol = [22 22 4 4 22 22];
% plot([-1 0 0 5 5 10], TempProtocol,'color','k','LineWidth',2);
% xticks([0 5 10]); yticks([4 22]); yticklabels({['4'];['22']}); xticklabels([0 5 10]);
% xlim([-1 10]); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([2 22]);
% box off; ax2 = gca; ax2.LineWidth = 1.5; ax2.FontSize = 10; ax2.Color = [1 1 1 0.5];
% text(ax2.XLim(1),ax2.YLim(2)*1.1, 'P', 'FontSize', 14);

% nexttile(); hold on;
ColorID = cmaplines(1,:);
p1 = plot(linspace(0,size(HDNmean,2)-1,size(HDNmean,2)),HDNmean,'color',...
    [ColorID 0.5],'LineWidth',2);

Time = 0:1:size(HDRNmean,2)-1; TimeSEM = [Time fliplr(Time)];
inBetweenSEM = [HDRNmean + HDRNsem fliplr(HDRNmean - HDRNsem)];
fill(TimeSEM, inBetweenSEM, ColorID, 'FaceAlpha', 0.5, 'EdgeColor', 'none');

HDRN = HDRN + (rand(size(HDRN))-0.5)./4;
TimeHDRN = 0:1:size(HDRN,2)-1;
TimeHDRN = repmat(TimeHDRN,size(HDRN,1),1);
TimeHDRN = TimeHDRN + rand(size(TimeHDRN))./1;

scatter(TimeHDRN,HDRN,'filled', 'MarkerFaceColor', ColorID,...
        'SizeData', 1);

%Add health for the last 5 days before the experiment
plot([-2 0], [1 1], 'Color', cmaplines(1,:), 'LineWidth', 2);
scatter([-2 -2 -2], [1 1 1],'filled', 'MarkerFaceColor', cmaplines(1,:),...
        'AlphaData', 0.1*ones(1,3), 'MarkerFaceAlpha',...
        'flat', 'SizeData', 25, 'jitter','on', 'jitterAmount',0.12);

box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10; ylabel('Health');
ax.YTick = 0:1:1; ax.XTick = 0:5:10; ax.XLim = [-2.2 10.2];
% xlabel (['Time (days)']);xlim([-0.1 6]); xticks(0:1:5)
ylim([-0.25 1.2]);
plot([-10 40], [0 0], 'LineStyle','--', 'Color','k', 'LineWidth',1);

box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10; ylabel('Health'); xlabel('Time (days)')


yyaxis right; ax = gca; ax.YColor = [0 0 0];
TempProtocol = [22 22 4 4 22 22];
plot([-1 0 0 5 5 10], TempProtocol,'color','k','LineWidth',2);
xticks([0 5 10]); yticks([4 22]); yticklabels({['4'];['22']}); xticklabels([0 5 10]);
xlim([-1 10]); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([-2 24]);
box off; ax2 = gca; ax2.LineWidth = 1.5; ax2.FontSize = 10; ax2.Color = [1 1 1 0.5];
text(ax2.XLim(1),ax2.YLim(2)*1.1, 'P', 'FontSize', 14);


%%% MEASURE HEALTH REGENERATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Reg22C = (Pilot.DecayRegen.HRN');

% nexttile()
% % plot(HRN');hold on;
% hold on;
% %Plot representative data
% Reg4C = [CstData3D(end-16:end,[6 11],4) CstData3D(end-16:end,[6 7],5) CstData3D(end-16:end,[8 9 10 11],6)];
% Reg4Cmean =  mean(Reg4C,2);
% Reg4Cstd = std(Reg4C');
% Reg22Cmean = mean(Reg22C,2);
% Reg22Cstd = std(Reg22C');
% 
% p1 = plot(Reg4Cmean, 'Color',cmaplines(1,:), 'LineWidth',2);
% % plot(Reg4C(:,4), 'Color', cmaplines(1,:), 'LineWidth', 2);
% p2 = plot(Reg22Cmean, 'Color', cmaplines(2,:), 'LineWidth', 2);
% % plot(Reg22C(:,1), 'Color', [cmaplines(2,:) 1], 'LineWidth', 2);
% 
% Xvalues = [1:length(Reg22Cmean)  fliplr(1:length(Reg22Cmean))];
% Yvalues = [Reg22Cmean' + Reg22Cstd fliplr(Reg22Cmean' - Reg22Cstd)];
% fill(Xvalues, Yvalues, 'c', 'FaceColor', cmaplines(2,:), 'FaceAlpha', 0.2, ...
%     'EdgeColor','none')
% 
% Xvalues = [1:length(Reg4Cmean)  fliplr(1:length(Reg4Cmean))];
% Yvalues = [Reg4Cmean' + Reg4Cstd fliplr(Reg4Cmean' - Reg4Cstd)];
% fill(Xvalues, Yvalues, 'c', 'FaceColor', cmaplines(1,:), 'FaceAlpha', 0.2, ...
%     'EdgeColor','none')
% 
% %Measure slopes of regeneration
% [ValMin4C IndexMin4C] = min(Reg4C);
% [ValMax4C IndexMax4C] = max(Reg4C);
% Slopes4C = (ValMax4C-ValMin4C)./(IndexMax4C-IndexMin4C);
% 
% [ValMin22C IndexMin22C] = min(Reg22C);
% [ValMax22C IndexMax22C] = max(Reg22C);
% Slopes22C = (ValMax22C-ValMin22C)./(IndexMax22C-IndexMin22C);
% Slopes22C(31) = []; % Exclude animal that died
% 
% ax = gca; ax.YLim = [0 1]; ax.LineWidth = 1.5; ax.FontSize = 10;
% ax.YLabel.String = 'Health'; ax.XLabel.String = 'Time (days)';
% text(0*ax.XLim(2),ax.YLim(2)*1.15, ' ')
% l = legend([p1(1) p2(1)], {'4C (n=8)','22C (n=32)'}); l.Box = 'off'; l.Location = 'southeast';
% text(ax.XLim(1),ax.YLim(2)*1.1, 'S', 'FontSize', 14);



%%% TEMPERATURE THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cartoon
% nexttile(); hold on;
% % ax2 = axes('Position',[.65 .28 .22 .15]); 
% TempProtocol = [linspace(20,20,50) linspace(4,4,24*5) linspace(20,20,50)];
% plot(-50:(size(TempProtocol,2)-51), TempProtocol,'color','k','LineWidth',2);
% xticks([0 120]); yticks([4 20]); yticklabels({['\theta'];['22']}); xticklabels([0 1]);
% xlim([-50 170]); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)'])
% box off; ax2 = gca; ax2.LineWidth = 1.5; ax2.FontSize = 10; ax2.Color = [1 1 1 0.5];
% text(ax2.XLim(1),ax2.YLim(2)*1.1, 'U', 'FontSize', 14);

nexttile(); hold on;
% ax2 = axes('Position',[.3 .28 .31 .13]); 
TempProtocol = [linspace(22,22,50) linspace(4,4,24*1) linspace(22,22,24*9)];
plot(-50:(size(TempProtocol,2)-51), TempProtocol,'color','k','LineWidth',2);
xticks([0 24 120 240]); yticks([4 22]); yticklabels({['\theta'];['22']}); xticklabels([0 1 5 10]);
xlim([-50 24*10]); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']);
box off; ax2 = gca; ax2.LineWidth = 1.5; ax2.FontSize = 10; ax2.Color = [1 1 1 0.5];
text(ax2.XLim(1),ax2.YLim(2)*1.1, 'P', 'FontSize', 14);

nexttile(); hold on;
A = sortrows(Pilot.TempThreshold,2); A = rmmissing(A);
D.Temp = A(:,2);
D.Health_after = A(:,4);
ivalues = unique(D.Temp);
meanhealth = []; semhealth = [];
% figtempthresh = figure('Name', 'Temperature Threshold', 'renderer', 'painters', 'Position', [452 301 350 400]);
% figure(figtempthresh)
% clf(figtempthresh); hold on;
for i = 1:size(ivalues,1)
    meanhealth = [meanhealth mean(D.Health_after(D.Temp == ivalues(i)))];
    semhealth = [semhealth std(D.Health_after(D.Temp == ivalues(i)))/sqrt(size(D.Health_after(D.Temp == ivalues(i)),1))];
end
% errorbar(ivalues,meanhealth,semhealth,semhealth,'Linewidth',2,'LineStyle','none');

Time = ivalues'; TimeSEM = [Time fliplr(Time)];
inBetweenSEM = [meanhealth + semhealth fliplr(meanhealth - semhealth)];
% p2 = fill(TimeSEM, inBetweenSEM, cmaplines(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');

HealthAfter = D.Health_after;
TempHealthAfter = D.Temp; %+ rand(size(D.Temp,1),1)./2;
HealthAfter = HealthAfter + (rand(size(HealthAfter))-0.5)./10;
p4 = scatter(TempHealthAfter,HealthAfter,'filled', 'MarkerFaceColor', cmaplines(1,:),...
    'SizeData', 4);

b = bar(ivalues,meanhealth);
b.EdgeColor = 'none';b.FaceColor = cmaplines(1,:); b.FaceAlpha = 0.4;

% %Complete with data from Cst Experiments
% HealthCst = squeeze(CstData3D(1,1:3,end-6:end));
% HealthCstR = HealthCst + (rand(size(HealthCst,1),size(HealthCst,2))-0.5)*0.2;
% TempCst = repmat([22:2:34],3,1);
% TempCst = TempCst + (rand(size(TempCst,1), size(TempCst,2))-0.5)*0.1;
% scatter(TempCst, HealthCstR,'filled', 'MarkerFaceColor', cmaplines(1,:),'SizeData', 4);
% b2 = bar(TempCst(1,:), mean(HealthCst));
% b2.EdgeColor = 'none'; b2.FaceColor = cmaplines(2,:); b2.FaceAlpha = 0.4; b2.BarWidth = 0.5;


% plot(ivalues,meanhealth,'color',cmaplines(1,:),'LineWidth',2);
box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
xlim([3.5 20.5]); ylim([0 1.2]); ylabel('Health after 1 day at \theta'); xlabel('\theta ({\circ}C)');
% xticks([5 10 15 20]);

dim = [.32 .82 .27 .07];
str = 'n = 3 for each \theta';
% text(mean(ax.XLim)./2,ax.YLim(2)*1.1,str);
ax.Title.String = str;
text(ax.XLim(1),ax.YLim(2)*1.1, 'F', 'FontSize', 14);

%%% OTHER VERSION OF TEMPERATURE THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(); hold on;
%Complete with data from Cst Experiments
HealthCst = squeeze(CstData3D(1,:,:));
HealthCstR = HealthCst + (rand(size(HealthCst,1),size(HealthCst,2))-0.5)*0.2;
TempCst = repmat([4:2:34],size(HealthCst,1),1);
TempCst = TempCst + (rand(size(TempCst,1), size(TempCst,2))-0.5)*0.1;

scatter(TempCst(:,1:8), HealthCstR(:,1:8),'filled', 'MarkerFaceColor', cmaplines(1,:),'SizeData', 4);
b4 = bar(TempCst(1,1:8), mean(HealthCst(:,1:8)));
b4.EdgeColor = 'none'; b4.FaceColor = cmaplines(1,:); b4.FaceAlpha = 0.4; b4.BarWidth = 0.5;

scatter(TempCst(:,9:16), HealthCstR(:,9:16),'filled', 'MarkerFaceColor', cmaplines(2,:),'SizeData', 4);
b3 = bar(TempCst(1,9:16), mean(HealthCst(:,9:16)));
b3.EdgeColor = 'none'; b3.FaceColor = cmaplines(2,:); b3.FaceAlpha = 0.4; b3.BarWidth = 0.5;

% plot(ivalues,meanhealth,'color',cmaplines(1,:),'LineWidth',2);
box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 10;
xlim([3.5 34.5]); ylim([0 1.2]); ylabel('Health after 1-14 days at T'); xlabel('T ({\circ}C)');
% xticks([5 10 15 20]);

dim = [.32 .82 .27 .07];
% str = 'n = 3 for each \theta';
% text(mean(ax.XLim)./2,ax.YLim(2)*1.1,str);
% ax.Title.String = str;
text(ax.XLim(1),ax.YLim(2)*1.1, 'F', 'FontSize', 14);

%%% Temperature of dish when moved between bench and fridge %%%%%%%%%%%%%%%
nexttile(); hold on;
time = Pilot.ATTM.Data(:,1)';
meantemperature = mean(Pilot.ATTM.Data(:,2:4),2)';
p1 = plot(time, meantemperature);
DataSEM = std(Pilot.ATTM.Data(:,2:4),[],2)';
DataSEM = [meantemperature+DataSEM fliplr(meantemperature - DataSEM)];
fill([time fliplr(time)], DataSEM, cmaplines(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');

meantemperature = mean(Pilot.ATTM.Data(:,5:7),2)';
p2 = plot(time, meantemperature);
DataSEM = std(Pilot.ATTM.Data(:,2:4),[],2)';
DataSEM = [meantemperature+DataSEM fliplr(meantemperature - DataSEM)];
fill([time fliplr(time)], DataSEM, cmaplines(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');

ax = gca; ax.YLabel.String = 'Temperature'; ax.XLabel.String = 'Time (min)'; ax.LineWidth = 1.5;
ax.FontSize = 10; ax.YTick = [4 22];
% l = legend([p1 p2], {'Bench to Fridge', 'Fridge to Bench'});
l = legend([p1 p2], {'22 to 4', '4 to 22'});
l.Box = 'off';
text(45,23,'n=3')
text(ax.XLim(1),ax.YLim(2)*1.1, 'Q', 'FontSize', 14);

%%% QUANTIFICATION OF HEALING SPEED AT DIFFERENT TEMPERATURES %%%%%%%%%%%%%
nexttile()
% plot(HRN');hold on;
hold on;
%Plot representative data
Reg4C = [CstData3D(end-16:end,[6 11],4) CstData3D(end-16:end,[6 7],5) CstData3D(end-16:end,[8 9 10 11],6)];
Reg4Cmean =  mean(Reg4C,2);
Reg4Cstd = std(Reg4C');
Reg22Cmean = mean(Reg22C,2);
Reg22Cstd = std(Reg22C');

p1 = plot(Reg4Cmean, 'Color',cmaplines(1,:), 'LineWidth',2);
% plot(Reg4C(:,4), 'Color', cmaplines(1,:), 'LineWidth', 2);
p2 = plot(Reg22Cmean, 'Color', cmaplines(2,:), 'LineWidth', 2);
% plot(Reg22C(:,1), 'Color', [cmaplines(2,:) 1], 'LineWidth', 2);

Xvalues = [1:length(Reg22Cmean)  fliplr(1:length(Reg22Cmean))];
Yvalues = [Reg22Cmean' + Reg22Cstd fliplr(Reg22Cmean' - Reg22Cstd)];
fill(Xvalues, Yvalues, 'c', 'FaceColor', cmaplines(2,:), 'FaceAlpha', 0.2, ...
    'EdgeColor','none')

Xvalues = [1:length(Reg4Cmean)  fliplr(1:length(Reg4Cmean))];
Yvalues = [Reg4Cmean' + Reg4Cstd fliplr(Reg4Cmean' - Reg4Cstd)];
fill(Xvalues, Yvalues, 'c', 'FaceColor', cmaplines(1,:), 'FaceAlpha', 0.2, ...
    'EdgeColor','none')

%Measure slopes of regeneration
[ValMin4C IndexMin4C] = min(Reg4C);
[ValMax4C IndexMax4C] = max(Reg4C);
Slopes4C = (ValMax4C-ValMin4C)./(IndexMax4C-IndexMin4C);

[ValMin22C IndexMin22C] = min(Reg22C);
[ValMax22C IndexMax22C] = max(Reg22C);
Slopes22C = (ValMax22C-ValMin22C)./(IndexMax22C-IndexMin22C);
Slopes22C(31) = []; % Exclude animal that died

ax = gca; ax.YLim = [0 1]; ax.LineWidth = 1.5; ax.FontSize = 10;
ax.YLabel.String = 'Health'; ax.XLabel.String = 'Time (days)';
text(0*ax.XLim(2),ax.YLim(2)*1.15, ' ')
l = legend([p1(1) p2(1)], {'4C (n=8)','22C (n=32)'}); l.Box = 'off'; l.Location = 'southeast';
text(ax.XLim(1),ax.YLim(2)*1.1, 'S', 'FontSize', 14);

nexttile(); hold on;
s1 = scatter(1, Slopes4C, 'MarkerFaceColor', cmaplines(1,:), 'MarkerEdgeColor','none', ...
    'MarkerFaceAlpha', 0.5, 'jitter','on', 'jitterAmount',0.1);
s2 = scatter(4, Slopes22C, 'MarkerFaceColor', cmaplines(2,:), 'MarkerEdgeColor','none', ...
    'MarkerFaceAlpha', 0.5, 'jitter','on', 'jitterAmount',0.1);
% l = legend([s1(1) s2(1)], {'4C','22C'}); l.Box = 'off'; l.Location = 'northwest';

ax = gca; ax.XLim = [0 5]; ax.YLim = [0 0.5]; ax.LineWidth = 1.5;
ax.XTick = [1 4]; ax.XTickLabel = {'4C' '22C'}; ax.FontSize = 10;
ax.YTick = [0:0.25:0.5];
ax.YLabel.String = '$\Delta H/\Delta T$'; ax.XLabel.String = 'Temperature';
ax.YLabel.Interpreter = 'latex';
text(ax.XLim(1),ax.YLim(2)*1.1, 'T', 'FontSize', 14);


f = gcf;
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(f, [f.Name '.pdf'])

end