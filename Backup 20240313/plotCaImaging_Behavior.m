function plotCaImaging_Behavior();
close all
clear all
clc

figure('Name','Calcium imaging and Behavior','WindowState','maximized')
tiledlayout(4,6);

load('data.mat');
cmaplines = colormap('lines');
nexttile();
for i=1:5
    plot(Behavior.LCtrain4Ctime./60, Behavior.LCtrain4C(i,:)'+i*1.2,'Color',cmaplines(1,:)); hold on;
end
ax = gca; ax.YTick = []; ax.XLabel.String = 'Time (min)'; ax.Box = 'off'; ax.YLabel.String = 'Events'; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*1.1,'A. Longitudinal contractions events');

nexttile()
for i=1:4
    plot(Behavior.LCtrain22Ctime./60, Behavior.LCtrain22C(i,:)'+i*1.2,'Color',cmaplines(2,:)); hold on;
end
ax = gca; ax.YTick = []; ax.XLabel.String = 'Time (min)'; ax.Box = 'off'; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*1.1,'B.');

nexttile(); hold on;
LC4Cfrequency = sum(Behavior.LCtrain4C,2)./Behavior.LCtrain4Ctime(end)*3600;
LC22Cfrequency = sum(Behavior.LCtrain22C,2)./Behavior.LCtrain22Ctime(end)*3600;
scatter(1,LC4Cfrequency,10,cmaplines(1,:),'filled', ...
    'MarkerFaceAlpha',1,'jitter','on','jitterAmount',0.05);

scatter(2,LC22Cfrequency,10,cmaplines(2,:),'filled', ...
    'MarkerFaceAlpha',1,'jitter','on','jitterAmount',0.05);
% scatter(xpos,HealthData,'filled', 'MarkerFaceColor', 'k',...
%         'SizeData', 25, 'jitter','on', 'jitterAmount',0.15);
b = bar(1, mean(mean(LC4Cfrequency)));
b.BarWidth = 0.2; b.EdgeColor = 'none'; b.FaceColor = cmaplines(1,:); b.FaceAlpha = 0.2;
b = bar(2, mean(mean(LC22Cfrequency)));
b.BarWidth = 0.2; b.EdgeColor = 'none'; b.FaceColor = cmaplines(2,:); b.FaceAlpha = 0.2;
ax = gca; ax.XLim = [0 3]; ax.XTick = [1 2]; ax.XTickLabel = {'4C','22C'}; ax.LineWidth = 1;
ax.YLabel.String = 'Frequency (1/hr)';
t = text(0,ax.YLim(2)*1.2,'C.');

%Display mean and sem of LC frequency
disp(['LC4C frequency = ' num2str(mean(LC4Cfrequency)) '+/-' num2str(std(LC4Cfrequency)./sqrt(length(LC4Cfrequency)))])
disp(['LC22C frequency = ' num2str(mean(LC22Cfrequency)) '+/-' num2str(std(LC22Cfrequency)./sqrt(length(LC22Cfrequency)))])

nexttile()
for i=1:4
    plot(Behavior.SomVel4CTime(i,:), Behavior.SomVel4C(i,:)'+i*1.2,'Color',cmaplines(1,:)); hold on;
end

ax = gca; ax.YTick = []; ax.XLabel.String = 'Time (sec)'; ax.Box = 'off'; ax.XLim = [0 400]; ax.YLabel.String = 'Velocity (a.u.)'; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*1.1,'D. Somersaulting duration');

nexttile()
for i=1:3
    plot(Behavior.SomVel22CTime(i,:), Behavior.SomVel22C(i,:)'+i*1.2,'Color',cmaplines(2,:)); hold on;
end
ax = gca; ax.YTick = []; ax.XLabel.String = 'Time (sec)'; ax.Box = 'off'; ax.XLim = [0 400]; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*1.1,'E.');

nexttile()
scatter(1,Behavior.SomDur4C,10,cmaplines(1,:),'filled', ...
    'MarkerFaceAlpha',1,'jitter','on','jitterAmount',0.05);
hold on;
scatter(2,Behavior.SomDur22C,10,cmaplines(2,:),'filled', ...
    'MarkerFaceAlpha',1,'jitter','on','jitterAmount',0.05);

b = bar(1, mean(mean(Behavior.SomDur4C)));
b.BarWidth = 0.2; b.EdgeColor = 'none'; b.FaceColor = cmaplines(1,:); b.FaceAlpha = 0.2;
b = bar(2, mean(mean(Behavior.SomDur22C)));
b.BarWidth = 0.2; b.EdgeColor = 'none'; b.FaceColor = cmaplines(2,:); b.FaceAlpha = 0.2;

%Display mean and sem of Somersaulting duration
disp(['mean SomDur4C = ' num2str(mean(Behavior.SomDur4C)) '+/-' num2str(std(Behavior.SomDur4C)./sqrt(length(Behavior.SomDur4C)))])
disp(['mean SomDur22C = ' num2str(mean(Behavior.SomDur22C)) '+/-' num2str(std(Behavior.SomDur22C)./sqrt(length(Behavior.SomDur22C)))])

ax = gca; ax.XLim = [0 3]; ax.XTick = [1 2]; ax.XTickLabel = {'4C','22C'}; ax.LineWidth = 1;
ax.YLabel.String = 'Duration (sec)';
t = text(0,ax.YLim(2)*1.1,'F.');

%%% CALCIUM IMAGING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(); hold on;
for i=1:2
    plot(CaImaging.CB4Ctimevector(1,:), normalize(CaImaging.CB4C(i,:),'range')+i*1.2-1,'Color',cmaplines(1,:));
end
ax = gca; ax.YTick = []; ax.YLabel.String = 'Fluor. (a.u.)';
ax.XLabel.String = 'Time (sec)'; ax.Box = 'off'; ax.XLim = [0 600]; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*1.1,'G. Contraction bursts');

%CBs
nexttile(); hold on;
for i=1:3
    plot(CaImaging.CB22Ctimevector(1,:), normalize(CaImaging.CB22C(i,:),'range')+i*1.2-1,'Color',cmaplines(2,:));
end
ax = gca; ax.YTick = [];
ax.XLabel.String = 'Time (sec)'; ax.Box = 'off'; ax.XLim = [0 30]; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*1.1,'H.');
ax1 = gca;

nexttile();hold on;
scatter(1,CaImaging.IntervalsCB4C,10,'filled','MarkerFaceColor',cmaplines(1,:),'jitter','on', 'jitteramount', 0.05)
scatter(2,CaImaging.IntervalsCB22C,10,'filled','MarkerFaceColor',cmaplines(2,:),'jitter','on', 'jitteramount', 0.05)
b = bar(1, mean(mean(CaImaging.IntervalsCB4C)));
b.BarWidth = 0.2; b.EdgeColor = 'none'; b.FaceColor = cmaplines(1,:); b.FaceAlpha = 0.2;
b = bar(2, mean(mean(CaImaging.IntervalsCB22C)));
b.BarWidth = 0.2; b.EdgeColor = 'none'; b.FaceColor = cmaplines(2,:); b.FaceAlpha = 0.2;
ax = gca; ax.XLim = [0 3]; ax.XTick = [1 2]; ax.XTickLabel = {'4C','22C'}; ax.YLabel.String = {'Interval (sec)'}; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*0.9,'I.');

%Display mean and sem of CB intervals
disp(['mean IntervalsCB4C = ' num2str(mean(CaImaging.IntervalsCB4C)) '+/-' num2str(std(CaImaging.IntervalsCB4C)./sqrt(length(CaImaging.IntervalsCB4C)))])
disp(['mean IntervalsCB22C = ' num2str(mean(CaImaging.IntervalsCB22C)) '+/-' num2str(std(CaImaging.IntervalsCB22C)./sqrt(length(CaImaging.IntervalsCB22C)))])

%RPs
nexttile(); hold on;
for i=1:2
    plot(CaImaging.RP4Ctimevector(1,:), normalize(CaImaging.RP4C(i,:),'range')+i*1.2-1,'Color',cmaplines(1,:));
end
ax = gca; ax.YTick = []; ax.YLabel.String = 'Fluor. (a.u.)'; 
ax.XLabel.String = 'Time (sec)'; ax.Box = 'off'; ax.XLim = [0 350]; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*1.1,'J. Rhythmic potentials');

nexttile(); hold on;
for i=1:3
    plot(CaImaging.RP22Ctimevector(1,:), normalize(CaImaging.RP22C(i,:),'range')+i*1.2-1,'Color',cmaplines(2,:));
end
ax = gca; ax.YTick = [];
ax.XLabel.String = 'Time (sec)'; ax.Box = 'off'; ax.XLim = [0 30]; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*1.1,'K.');
ax2 = gca;

nexttile();hold on;
scatter(1,CaImaging.IntervalsRP4C,10,'filled','MarkerFaceColor',cmaplines(1,:),'jitter','on', 'jitteramount', 0.05)
scatter(2,CaImaging.IntervalsRP22C,10,'filled','MarkerFaceColor',cmaplines(2,:),'jitter','on', 'jitteramount', 0.05)
b = bar(1, mean(mean(CaImaging.IntervalsRP4C)));
b.BarWidth = 0.2; b.EdgeColor = 'none'; b.FaceColor = cmaplines(1,:); b.FaceAlpha = 0.2;
b = bar(2, mean(mean(CaImaging.IntervalsRP22C)));
b.BarWidth = 0.2; b.EdgeColor = 'none'; b.FaceColor = cmaplines(2,:); b.FaceAlpha = 0.2;
ax = gca; ax.XLim = [0 3]; ax.XTick = [1 2]; ax.XTickLabel = {'4C','22C'}; ax.YLabel.String = {'Interval (sec)'}; ax.LineWidth = 1;
t = text(0,ax.YLim(2)*1.1,'L.');

%Display mean and sem of RP intervals
disp(['mean IntervalsRP4C = ' num2str(mean(CaImaging.IntervalsRP4C)) '+/-' num2str(std(CaImaging.IntervalsRP4C)./sqrt(length(CaImaging.IntervalsRP4C)))])
disp(['mean IntervalsRP22C = ' num2str(mean(CaImaging.IntervalsRP22C)) '+/-' num2str(std(CaImaging.IntervalsRP22C)./sqrt(length(CaImaging.IntervalsRP22C)))])

% nexttile(); 
axposition = ax1.Position;
axes('Position',[axposition(1)+axposition(3)/2 axposition(2)+axposition(4)/1.7 axposition(3)/1.25 axposition(4)/2]);
hold on;
meanF = mean(CaImaging.IndividualTransientsCB4C,1);
semF = std(CaImaging.IndividualTransientsCB4C,1);
signaly = [(meanF + semF) fliplr(meanF - semF)];
signalx = [CaImaging.IndividualTransientsCB4CTime fliplr(CaImaging.IndividualTransientsCB4CTime)];
plot(CaImaging.IndividualTransientsCB4CTime, mean(CaImaging.IndividualTransientsCB4C,1), 'Color',cmaplines(1,:))
fill(signalx, signaly, cmaplines(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');

meanF = mean(CaImaging.IndividualTransientsCB22C,1);
semF = std(CaImaging.IndividualTransientsCB22C,1);
signaly = [(meanF + semF) fliplr(meanF - semF)];
signalx = [CaImaging.IndividualTransientsCB22CTime fliplr(CaImaging.IndividualTransientsCB22CTime)];
plot(CaImaging.IndividualTransientsCB22CTime, mean(CaImaging.IndividualTransientsCB22C,1), 'Color',cmaplines(2,:))
fill(signalx, signaly, cmaplines(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
ax = gca; ax.LineWidth = 1; ax.XLim = [0 6];
ax.YTick = []; %ax.YLabel.String = 'Fluor. (a.u.)';ax.XLabel.String = 'Time (sec)';
% t = text(0,ax.YLim(2)*1.1,'M. Contraction bursts - average');

% nexttile(); 
axposition = ax2.Position;
axes('Position',[axposition(1)+axposition(3)/2.8 axposition(2)+axposition(4)/1.7 axposition(3)/1.25 axposition(4)/2]);
hold on;
meanF = mean(CaImaging.IndividualTransientsRP4C,1);
semF = std(CaImaging.IndividualTransientsRP4C,1);
signaly = [(meanF + semF) fliplr(meanF - semF)];
signalx = [CaImaging.IndividualTransientsRP4CTime fliplr(CaImaging.IndividualTransientsRP4CTime)];
plot(CaImaging.IndividualTransientsRP4CTime, mean(CaImaging.IndividualTransientsRP4C,1), 'Color',cmaplines(1,:))
fill(signalx, signaly, cmaplines(1,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');

meanF = mean(CaImaging.IndividualTransientsRP22C,1);
semF = std(CaImaging.IndividualTransientsRP22C,1);
signaly = [(meanF + semF) fliplr(meanF - semF)];
signalx = [CaImaging.IndividualTransientsRP22CTime fliplr(CaImaging.IndividualTransientsRP22CTime)];
plot(CaImaging.IndividualTransientsRP22CTime, mean(CaImaging.IndividualTransientsRP22C,1), 'Color',cmaplines(2,:))
fill(signalx, signaly, cmaplines(2,:), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
ax = gca; ax.LineWidth = 1; ax.XLim = [0 6];
ax.YTick = []; %ax.YLabel.String = 'Fluor. (a.u.)';ax.XLabel.String = 'Time (sec)';
% t = text(0,ax.YLim(2)*1.1,'N. Rhythmic potentials - average');

f = gcf;
exportgraphics(f, [f.Name '.jpg'])

end