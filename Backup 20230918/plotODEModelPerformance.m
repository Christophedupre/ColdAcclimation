function plotODEModelPerformance(Model1Data, Model2Data, CstData3D, Cst10C13D30Reps, Ramp, Pulses1_1, Pilot, Deacclimation)

figModelPerformance = figure('Name', 'Model Performance', 'renderer', 'painters');
figModelPerformance.WindowState = 'maximized';
clf(figModelPerformance)
tiledlayout(4,5)
cmaplines = colormap("lines");
cmapjet = colormap('jet');
cmapturbo = colormap('turbo');
load("cmaplist.mat");

%%% HEATMAP OF TRAINING EFFICIENCY (EXPERIMENTAL DATA) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Heatmap of training efficiency
nexttile()
hold off;
h = heatmap(squeeze(mean(CstData3D,1))');
h.ColorbarVisible = 'off'; h.XLabel = 'Training Duration (Days)'; h.YLabel = 'Training Temp. (C)';
h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''}; h.FontSize = 8;
h.XDisplayLabels([1 5 10 15]) = {'0', '5', '10', '15'};
h.YDisplayLabels([1 6 11 16]) = {'4', '14', '24', '34'};
h.Title = 'Training efficiency (expt data)';

%%% HEATMAP OF TRAINING EFFICIENCY (MODEL DATA W/O NOISE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile();
h0 = heatmap(squeeze(mean(Model1Data.HealthTruncated(8:12,:,:),1))');
h0.ColorbarVisible = 'off'; h0.XLabel = 'Training Duration (Days)'; h0.YLabel = 'Training Temp. (C)';
h0.XDisplayLabels(:) = {''}; h0.YDisplayLabels(:) = {''}; h0.FontSize = 8;
h0.XDisplayLabels([1 5 10 15]) = {'0', '5', '10', '15'};
h0.YDisplayLabels([1 6 11 16]) = {'4', '14', '24', '34'};
cmaphisto = h0.Colormap;
h0.Colormap = cmaphisto;
h0.Title = 'Training efficiency (model data w/o noise)';

%%% HEATMAP OF TRAINING EFFICIENCY (MODEL DATA WITH NOISE) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile();
h0 = heatmap(squeeze(mean(Model2Data.HealthTruncated(8:12,:,:),1))');
h0.ColorbarVisible = 'off'; h0.XLabel = 'Training Duration (Days)'; h0.YLabel = 'Training Temp. (C)';
h0.XDisplayLabels(:) = {''}; h0.YDisplayLabels(:) = {''}; h0.FontSize = 8;
h0.XDisplayLabels([1 5 10 15]) = {'0', '5', '10', '15'};
h0.YDisplayLabels([1 6 11 16]) = {'4', '14', '24', '34'};
cmaphisto = h0.Colormap;
h0.Colormap = cmaphisto;
h0.Title = 'Training efficiency (model data w noise)';

%%% Error and survival %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
metric.error = zeros(1,15,16); 
metric.SurvivalModel = 15*ones(1,16); 
metric.SurvivalData = 15*ones(1,16); 

for theta = 1:length(CstData3D(1,1,:))
    for D = 1:length(CstData3D(1,:,1))
        metric.error(1,D,theta) = [sum(sum((Model2Data.HealthTruncated(:,D,theta)-CstData3D(:,D,theta)').^2))];        
        % Measure survival rate
        if (min(Model2Data.HealthTruncated(:,D,theta)-0.01)<=0)
            metric.SurvivalModel(theta) = metric.SurvivalModel(theta) - 1;
        end
        if (min(CstData3D(:,D,theta))<=0)
            metric.SurvivalData(theta) = metric.SurvivalData(theta) - 1;
        end
    end
end

%Plot animal survival
nexttile(); hold on;
YData = zeros(1,length([metric.SurvivalData])*2);
YData(1:2:length([metric.SurvivalData])*2) = [metric.SurvivalData];
YData(2:2:length([metric.SurvivalData])*2) = [metric.SurvivalData];
XData = [4:1:35];
p1 = plot(XData, YData,'LineWidth',2,'Color',cmaplines(1,:));
YData = zeros(1,length([metric.SurvivalModel])*2);
YData(1:2:length([metric.SurvivalModel])*2) = [metric.SurvivalModel];
YData(2:2:length([metric.SurvivalModel])*2) = [metric.SurvivalModel];
XData = [4:1:35];
p2 = plot(XData, YData, 'LineWidth',2,'Color', cmaplines(5,:));
ax = gca; ax.XTick = [4:6:34]; ax.Box = 'off'; ax.LineWidth = 1; ax.FontSize = 8; ax.YTick = [0:15:15]; ax.YLim = [0 16];
ax.XLim = [0 36];
text(mean(ax.XLim), -0.5*mean(ax.YLim),'Training temperature (C)','HorizontalAlignment','center'); 
text(-0.2*mean(ax.XLim), mean(ax.YLim),'Survival', 'HorizontalAlignment','center','Rotation',90);
% text(ax.XLim(1),ax.YLim(2)*1.15, 'D', 'FontSize', 14);
l = legend([p1 p2], [{'Expt'}, {'Model'}]);
l.Box = 'off';

%%% VARIABILITY OVER 30 REPEATS AT D = 13 AND THETA = 4 %%%%%%%%%%%%%%%%%%%

nexttile(); hold on;
solpts = zeros(44,30);
f = Model1Data.TimeResolution;
NoiseAmplitude = Model1Data.NoiseAmplitude;

% Generate model data
for i = 1:30
%     cvar = Model2Data.csol.c.*(1+(rand(3,1)-0.5)/Model2Data.cvarcoeff);% Add variability for acclimation, injury and regeneration time constants    
    cvar = Model2Data.csol.c; % No var in coefficients
    temperature1 = linspace(0,43,44*f);
    temperature2 = [linspace(22,22,(15-(14-1))*f) linspace(4*2+2,4*2+2,(14-1)*f) linspace(4,4,29*f)]./22; %Duration = 13 days, training temperature = 10C
    temperature = [temperature1; AddNoiseToTemp(temperature2,NoiseAmplitude)];    
    sol = ode23(@(t,y) AcclimfunNoise(t,y,cvar,temperature), Model2Data.tspan, Model2Data.y0);
    solptstemp = deval(Model2Data.tspan,sol);
    solpts(:,i) = solptstemp(2,:);
end

%Plot experimental data
p1 = plot(0:1:28, mean(Cst10C13D30Reps,2), 'Color', [0 0 0], 'LineWidth', 2);
sd = std(Cst10C13D30Reps');%./sqrt(size(Cst10C13D30Reps,1));
fi = fill([0:1:28 fliplr(0:1:28)], [mean(Cst10C13D30Reps,2)'-sd fliplr(mean(Cst10C13D30Reps,2)'+sd)], cmaplines(1,:));
fi.FaceAlpha = 0.5; fi.EdgeAlpha = 0;

XData = repmat([1:1:29],30,1); XData = XData + (rand(30,29)-0.5)./1;
YData = Cst10C13D30Reps'; YData = YData + (rand(30,29)-0.5)./5;
s = scatter(XData,YData,'MarkerFaceColor',...
        cmaplines(1,:), 'MarkerEdgeColor','none','SizeData',2);

%Plot model data
p2 = plot(0:1:29, mean(solpts(15:44,:),2), 'Color', cmaplines(5,:), 'LineWidth', 2);
sd2 = std(solpts');%./sqrt(size(Cst10C13D30Reps,1));
sd2 = sd2(:,15:44);
fi2 = fill([0:1:29 fliplr(0:1:29)], [mean(solpts(15:44,:),2)'-sd2 fliplr(mean(solpts(15:44,:),2)'+sd2)], cmaplines(5,:));
fi2.FaceAlpha = 0.5; fi2.EdgeAlpha = 0;

XData = repmat([0:1:29],30,1); XData = XData + (rand(30,30)-0.5)./1;
YData = solpts(15:44,:)'; YData = YData + (rand(30,30)-0.5)./5;
s = scatter(XData,YData,'MarkerFaceColor',...
        cmaplines(5,:), 'MarkerEdgeColor','none','SizeData',2);

yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 8;
xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1.2]); yticks([0:1]);
text(-0.2*mean(ax.XLim), mean(ax.YLim),'Health', 'HorizontalAlignment','center','Rotation',90);
text(ax.XLim(2)*0.02,ax.YLim(2)*0.2,{'TT = 10 C, TD = 13 Days, 30 reps'}, 'BackgroundColor', [1 1 1 0.75], 'FontSize',9);
text(mean(ax.XLim), -0.5*mean(ax.YLim),'T (days)','HorizontalAlignment','center'); 
% text(ax.XLim(1),ax.YLim(2)*1.15, 'E', 'FontSize', 14);

% nexttile(); axis 'off';
% nexttile(); axis 'off';

%Plot statistics about experimental and model data
% % Mean health distribution
% nexttile()
% h1 = histogram(mean(Cst10C13D30Reps)); hold on;
% h2 = histogram(mean(solpts(15:44,:)));
% h1.Normalization = 'probability'; h2.Normalization = 'probability';
% h1.EdgeColor = 'none'; h1.FaceAlpha = 0.5; h1.BinEdges = linspace(0,1,14);
% h2.EdgeColor = 'none'; h2.FaceAlpha = 0.5; h2.BinEdges = linspace(0,1,14); h2.FaceColor = cmaplines(5,:);
% ax = gca; ax.Box = 'off'; ax.LineWidth = 1; ax.FontSize = 8; ax.XLim = [0 1]; ax.YTick = [0:0.7:0.7]; ax.YLim = [0 0.7];
% text(mean(ax.XLim), -0.5*mean(ax.YLim),'Avg health', 'HorizontalAlignment','center');
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Probability', 'HorizontalAlignment','center','Rotation',90);
% 
% % Min health distribution
% nexttile()
% h1 = histogram(min(Cst10C13D30Reps)); hold on;
% h2 = histogram(min(solpts(15:44,:)));
% h1.Normalization = 'probability'; h2.Normalization = 'probability';
% h1.EdgeColor = 'none'; h1.FaceAlpha = 0.5; h1.BinEdges = linspace(0,1,14);
% h2.EdgeColor = 'none'; h2.FaceAlpha = 0.5; h2.BinEdges = linspace(0,1,14); h2.FaceColor = cmaplines(5,:);
% ax = gca; ax.Box = 'off'; ax.LineWidth = 1; ax.FontSize = 8; ax.XLim = [0 1]; ax.YTick = [0:0.7:0.7]; ax.YLim = [0 0.7];
% text(mean(ax.XLim), -0.5*mean(ax.YLim),'Min health', 'HorizontalAlignment','center');
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Probability', 'HorizontalAlignment','center','Rotation',90);

% Final health distribution
% nexttile()
% h1 = histogram(Cst10C13D30Reps(end,:)); hold on;
% h2 = histogram(solpts(end,:));
% h1.Normalization = 'probability'; h2.Normalization = 'probability';
% h1.EdgeColor = 'none'; h1.FaceAlpha = 0.5; h1.BinEdges = linspace(0,1,14);
% h2.EdgeColor = 'none'; h2.FaceAlpha = 0.5; h2.BinEdges = linspace(0,1,14); h2.FaceColor = cmaplines(5,:);
% ax = gca; ax.Box = 'off'; ax.LineWidth = 1; ax.FontSize = 8; ax.XLim = [0 1]; ax.YTick = [0:1:1]; ax.YLim = [0 1];
% text(mean(ax.XLim), -0.5*mean(ax.YLim),'Final health', 'HorizontalAlignment','center');
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Probability', 'HorizontalAlignment','center','Rotation',90);

% Max health distribution
% nexttile()
% h1 = histogram(max(Cst10C13D30Reps)); hold on;
% h2 = histogram(max(solpts(15:44,:)));
% h1.EdgeColor = 'none'; h1.FaceAlpha = 0.5; h1.BinEdges = linspace(0,1,14);
% h2.EdgeColor = 'none'; h2.FaceAlpha = 0.5; h2.BinEdges = linspace(0,1,14); h2.FaceColor = cmaplines(5,:);
% ax = gca;ax.Box = 'off'; ax.LineWidth = 1; ax.FontSize = 8;
% text(ax.XLim(1),ax.YLim(2)*1.1,{'Max health distribution'});

% nexttile
% ax = gca; ax.Visible = 'off';


%%% TRUNCATED RAMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Truncated Ramps Protocol Diagram
% nexttile(); hold on;
% TrainingDuration = 0:1:14;
% imax = 15;
% for i=[1 2 5 10 15]
%     % Time and Temperature Vectors for Training Phase
%     X = []; Y = [];
%     X = [X linspace(0,14,15)]; %Training Time
%     Y = [Y linspace(22, 22, 15-TrainingDuration(i))]; 
%     slope = (22-4)/14;
%     Y = [Y linspace(22,22-slope*TrainingDuration(i),TrainingDuration(i))];%Training Temperature
%     idx = length(X); %index of end of training
%     % Test phase
%     X = [X linspace(14,43,29)];
%     Y = [Y linspace(4,4,29)];
% %     plot(X,Y, 'Color', cmaplines(1,:), 'LineStyle',':', 'LineWidth',2);
%     plot(X,Y, 'Color', [cmapjet(round(i/imax*256),:) 1], 'LineStyle',':', 'LineWidth',2);
%     if(i == 10)
% %         plot(X,Y, 'Color', [cmapjet(round(i/imax*256),:) 1], 'LineStyle','-', 'LineWidth',2);
%     end
%     plot(X(idx+1:end), Y(idx+1:end), 'Color', cmaplines(2,:), 'LineWidth',2);
% end
% 
% ax = gca; ax.LineWidth = 1; ax.FontSize = 8; ax.XLim = [0 28]; 
% ax.XTick = [0 14 21 28]; ax.XTickLabel = {'0' '14' '[...]' '42'}; 
% ylim([0 22]); yticks([4 22]); yticklabels({'4';'22'})
% text(ax.XLim(2)*0.05,ax.YLim(2)*0.9,{'Truncated Ramp'}, 'BackgroundColor', [1 1 1 0.75]);
% text(mean(ax.XLim), -0.5*mean(ax.YLim) ,'T (days)', 'HorizontalAlignment','center');
% text(-0.2*mean(ax.XLim), mean(ax.YLim) , ['T (' char(176) 'C)'], 'HorizontalAlignment','center','Rotation',90);

% %Other diagram for truncated ramp
% nexttile()
% plot([linspace(0,7,15)], [linspace(22,13,15)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2); hold on;
% plot([linspace(7,14,8)], [linspace(13,4,8)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2, 'LineStyle',':');
% % plot([linspace(7,14,8)], [linspace(4,4,8)],'color',[cmap(2,:) 0.8], 'LineWidth',2);
% plot([linspace(7,7,8)], [linspace(13,4,8)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2);
% plot(linspace(7,21,15), linspace(4,4,15),'color',[cmaplines(2,:) 0.8], 'LineWidth',2);
% plot(linspace(21,28,8), linspace(4,4,8),'color',[cmaplines(2,:) 0.8], 'LineWidth',2, 'LineStyle',':');
% xlim([0 28]); xticks([0 7 21]); xticklabels({'0';'D';'D+28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 22])
% yticks([4 22]); yticklabels({'4';'20'})
% box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 8;
% text(ax.XLim(1),ax.YLim(2)*1.15, 'F', 'FontSize', 14);

%Other diagram for truncated ramp
nexttile(); hold on;
RampDuration1 = 1:4:14;
RampDuration2 = 22-RampDuration1.*18/RampDuration1(end);
for i = 1:length(RampDuration1)
    plot([0 RampDuration1(i) RampDuration1(i) 14 42], [22 RampDuration2(i) 4 4 4]+i/4,'LineWidth',1);
end
ax = gca; ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'T (C)';

% Experimental Data
nexttile(); 
% RampExtended = [NaN(15,15); Ramp];
h2 = heatmap(Ramp');
h2.ColorbarVisible = 'off'; h2.XLabel = 'Testing Time (Days)'; h2.YLabel = 'Training Duration (Day)';
h2.FontSize = 8;
h2.Title = 'Health expt';
h2.XDisplayLabels(:) = {''}; h2.YDisplayLabels(:) = {''};
h2.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h2.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
h2.Colormap = cmaphisto;

RampDuration = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14];

% Model without noise
for i = 1:length(RampDuration);
%     cvar = Model2Data.csol.c.*(1+(rand(3,1)-0.5)/4);% Add variability for acclimation, injury and regeneration time constants
    cvar = Model2Data.csol.c;% no variability in c
    temperature1 = linspace(0,43,44*f);
    temperature2 = [linspace(22,22,(15-RampDuration(i))*f) linspace(22,22-((22-4)/14)*RampDuration(i),RampDuration(i)*f) linspace(4,4,29*f)]./22;
    % temperature2 = AddNoiseToTemp(temperature2,NoiseAmplitude);
    temperature = [temperature1; temperature2];
    sol = ode23(@(t,y) Acclimfun(t,y,cvar,temperature), Model2Data.tspan, Model2Data.y0);
    solpts = deval(sol,Model2Data.tspan);
    Model1DataTR.HealthFull(:,1,i) = solpts(2,:);
    Model1DataTR.HealthTruncated(:,1,i) = solpts(2,linspace(15,43,29));
    Model1DataTR.AlphaFull(:,1,i) = solpts(1,:);
end

nexttile(); 
h1 = heatmap(squeeze(Model1DataTR.HealthTruncated(:,1,:))');
h1.ColorbarVisible = 'off'; h1.XLabel = 'Testing Time (Days)'; h1.YLabel = 'Training Duration (Day)';
h1.FontSize = 8;
h1.Title = 'Health model';
h1.XDisplayLabels(:) = {''}; h1.YDisplayLabels(:) = {''};
h1.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h1.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
h1.Colormap = cmaphisto;


% Model with noise
for i = 1:length(RampDuration);
%     cvar = Model2Data.csol.c.*(1+(rand(3,1)-0.5)/4);% Add variability for acclimation, injury and regeneration time constants
    cvar = Model2Data.csol.c;% no variability in c
    temperature1 = linspace(0,43,44*f);
    temperature2 = [linspace(22,22,(15-RampDuration(i))*f) linspace(22,22-((22-4)/14)*RampDuration(i),RampDuration(i)*f) linspace(4,4,29*f)]./22;
    temperature2 = AddNoiseToTemp(temperature2,NoiseAmplitude);
    temperature = [temperature1; temperature2];
    sol = ode23(@(t,y) AcclimfunNoise(t,y,cvar,temperature), Model2Data.tspan, Model2Data.y0);
    solpts = deval(sol,Model2Data.tspan);
    Model2DataTR.HealthFull(:,1,i) = solpts(2,:);
    Model2DataTR.HealthTruncated(:,1,i) = solpts(2,linspace(15,43,29));
    Model2DataTR.AlphaFull(:,1,i) = solpts(1,:);
end

nexttile(); 
h1 = heatmap(squeeze(Model2DataTR.HealthTruncated(:,1,:))');
h1.ColorbarVisible = 'off'; h1.XLabel = 'Testing Time (Days)'; h1.YLabel = 'Training Duration (Day)';
h1.FontSize = 8;
h1.Title = 'Health model with noise';
h1.XDisplayLabels(:) = {''}; h1.YDisplayLabels(:) = {''};
h1.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h1.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
h1.Colormap = cmaphisto;

nexttile(); hold on;
p1 = plot(linspace(0,14,15), mean(Ramp,1), 'LineWidth',2, 'Color', cmaplines(1,:));
p2 = plot(linspace(0,14,15), squeeze(mean(Model1DataTR.HealthFull(15:end,1,:),1)), 'LineWidth',2, 'Color', cmaplines(5,:));
p3 = plot(linspace(0,14,15), squeeze(mean(Model2DataTR.HealthFull(15:end,1,:),1)), 'LineWidth',2, 'Color', cmaplines(5,:),'LineStyle',':');
l = legend([p1 p2 p3],[{'expt'},{'model'}, {'model with noise'}]);
l.Box = 'off'; l.Location = "northwest";
box off; ax = gca ; ax.LineWidth = 1; ax.FontSize = 8;
xlim([0 15]); xticks([0:5:15]); ylim([0 1]); yticks([0:1:1]);
text(-0.2*mean(ax.XLim), mean(ax.YLim), {'avg. health'}, 'Rotation',90, 'VerticalAlignment','middle','HorizontalAlignment','center');
% text(mean(ax.XLim),-0.4*mean(ax.YLim),'D', 'HorizontalAlignment','center');
ax = gca; ax.XLabel.String = 'Training Duration (Days)';
% text(ax.XLim(1),ax.YLim(2)*1.15, 'I', 'FontSize', 14);

% nexttile(); hold on;
% %Survival
% SurvivalRamp = ceil(min(Ramp(:,:)));
% SurvivalModel = ceil(squeeze(min(Model3Data.HealthFull(:,1,:)))-0.01);
% 
% 
% b1 = bar([0:1:14],[SurvivalRamp], 'EdgeColor', 'none','FaceColor',cmaplines(1,:), 'FaceAlpha',0.5); hold on;
% b2 = bar([0:1:14],[SurvivalModel], 'EdgeColor','none','FaceColor',cmaplines(5,:), 'FaceAlpha',0.5);
% ax = gca; ax.XTick = [0:5:15]; ax.Box = 'off'; ax.LineWidth = 1; ax.FontSize = 8; ax.YTick = [0:2:2]; ax.YLim = [0 2]; ax.XLim = [0 15];
% % text(mean(ax.XLim), -0.4*mean(ax.YLim),'D','HorizontalAlignment','center');
% ax = gca; ax.XLabel.String = 'Training Duration (Days)';
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Survival', 'HorizontalAlignment','center','Rotation',90);
% % text(ax.XLim(1),ax.YLim(2)*1.15, 'J', 'FontSize', 14);
% % l = legend([b1 b2],'Data','Model'); l.Box = 'off'; l.FontSize = 7; %l.Position(1:2) = [0.86 0.89];

%%% PULSES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cartoon for pulses diagram
nexttile(); hold on;
StepDuration = [1:4:14];
for i = 1:length(StepDuration)
    plot([0 StepDuration(i) StepDuration(i) 42], [13 13 4 4]+i/2-1, 'Color',cmaplines(round(256/length(StepDuration)*i),:),'LineWidth',1);    
    if i == 1
        r = rectangle('Position',[0 4+i/4-1 StepDuration(i) 13-4], 'FaceColor', [cmaplines(round(256/length(StepDuration)*i),:),0.5]);
    else
        r = rectangle('Position',[StepDuration(i-1) 4+i/4-1 4 13-4], 'FaceColor', [cmaplines(round(256/length(StepDuration)*i),:),0.5]);
    end
    r.EdgeColor = 'none';
end
l = line([StepDuration(i-1) 18],[13+i/2-1 20],'Color',cmaplines(4,:));
l = line([StepDuration(i-1)+4 18],[13+i/2-1 13],'Color',cmaplines(4,:));

plot([20 20 35],[20 8 8]+2,'Color','k');
plot([20 22 22 24 24 26 26 28 28 30 30 32 32 34],[10 10 20 20 10 10 20 20 10 10 20 20 10 10]+2,'Color','k');
text(19,8,'0','FontSize',8)
text(31,8,'0.5','FontSize',8)
ax = gca; ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'T (C)';


ax = gca;
% ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'Train. Temp (C)';
ax.YLim = [0 25];
% ax.XTick = []; ax.YTick = [];

% generate model data
nexttile(); 
% hold on; 
clear solpts solptstemp;
imax = 15;

% Experimental Data (Pulses)
% hold on;
% Time = linspace(15,length(Pulses1_1(:,1))-1+15,length(Pulses1_1(:,1)));
% imax = length(Pulses1_1(1,:));
% for i = 1:imax
%     plot(Time, Pulses1_1(:,i), 'Color',[cmapjet(round(i/imax*256),:) 1], 'LineWidth', 2);
% end
% ax = gca; ax.XLim = [-0.1 44]; ax.XTick = [0:10:44]; ax.YLim = [0 1]; ax.YTick = [0:1:1]; ax.LineWidth = 1; ax.FontSize = 8;
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Health expt', 'HorizontalAlignment','center','Rotation',90);

% Pulses1_1Extended = [NaN(15,15); Pulses1_1];
% h4 = heatmap(Pulses1_1Extended');
h4 = heatmap(Pulses1_1');
h4.ColorbarVisible = 'off'; h4.XLabel = 'Testing Time (Days)'; h4.YLabel = 'Training Duration (Day)';
h4.FontSize = 8;
h4.Title = 'Health expt';
h4.XDisplayLabels(:) = {''}; h4.YDisplayLabels(:) = {''};
h4.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
h4.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
h4.Colormap = cmaphisto;

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
h3.Colormap = cmaphisto;
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
h3.Colormap = cmaphisto;

% Training efficiency
nexttile(); hold on;
plot(linspace(0,14,15), mean(Pulses1_1,1), 'LineWidth',2, 'Color', cmaplines(1,:));
plot(linspace(0,14,15), squeeze(mean(Model1Pulses1_1Data.HealthFull(15:end,1,:),1)), 'LineWidth',2, 'Color', cmaplines(5,:));
plot(linspace(0,14,15), squeeze(mean(Model2Pulses1_1Data.HealthFull(15:end,1,:),1)), 'LineWidth',2, 'Color', cmaplines(5,:), 'LineStyle',':');
box off; ax = gca ; ax.LineWidth = 1; ax.FontSize = 8;
xlim([0 15]); xticks([0:5:15]); ylim([0 1]); yticks([0:1:1]);
text(-0.2*mean(ax.XLim), mean(ax.YLim), {'avg. health'}, 'Rotation',90, 'VerticalAlignment','middle','HorizontalAlignment','center');
% text(mean(ax.XLim),-0.4*mean(ax.YLim),'D', 'HorizontalAlignment','center');
ax = gca; ax.XLabel.String = 'Training Duration (Days)';
% text(ax.XLim(1),ax.YLim(2)*1.15, 'N', 'FontSize', 14);

% % Survival
% nexttile(); hold on;
% SurvivalPulses = ceil(min(Pulses1_1(:,:)));
% SurvivalModel = ceil(squeeze(min(Model2Pulses1_1Data.HealthFull(:,1,:)))-0.01);
% 
% b1 = bar([0:1:14],[SurvivalPulses], 'EdgeColor', 'none','FaceColor',cmaplines(1,:), 'FaceAlpha',0.5); hold on;
% b2 = bar([0:1:14],[SurvivalModel], 'EdgeColor','none','FaceColor',cmaplines(5,:), 'FaceAlpha',0.5);
% ax = gca; ax.XTick = [0:5:15]; ax.Box = 'off'; ax.LineWidth = 1; ax.FontSize = 8; ax.YTick = [0:2:2]; ax.YLim = [0 2]; ax.XLim = [0 15];
% % text(mean(ax.XLim),-0.4*mean(ax.YLim),'D', 'HorizontalAlignment','center');
% ax = gca; ax.XLabel.String = 'Training Duration (Days)';
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Survival', 'HorizontalAlignment','center','Rotation',90);
% % text(ax.XLim(1),ax.YLim(2)*1.15, 'O', 'FontSize', 14);
% % l = legend([b1 b2],'Data','Model'); l.Box = 'off'; l.FontSize = 7; %l.Position(1:2) = [0.86 0.89];

% %Plot animal survival
% nexttile(); hold on;
% SurvivalPulses = ceil(min(Pulses1_1(:,:)));
% SurvivalModel = ceil(squeeze(min(Model2Pulses1_1Data.HealthFull(:,1,:)))-0.01);
% YData = zeros(1,length([metric.SurvivalData])*2);
% YData(1:2:length([SurvivalPulses])*2) = [SurvivalPulses];
% YData(2:2:length([SurvivalPulses])*2) = [SurvivalPulses];
% XData = [4:1:35];
% p1 = plot(XData, YData,'LineWidth',2,'Color',cmaplines(1,:));
% YData = zeros(1,length([metric.SurvivalModel])*2);
% YData(1:2:length([SurvivalModel])*2) = [SurvivalModel];
% YData(2:2:length([SurvivalModel])*2) = [SurvivalModel];
% XData = [4:1:35];
% p2 = plot(XData, YData, 'LineWidth',2,'Color', cmaplines(5,:));
% ax = gca; ax.XTick = [4:6:34]; ax.Box = 'off'; ax.LineWidth = 1; ax.FontSize = 8; ax.YTick = [0:15:15]; ax.YLim = [0 2];
% ax.XLim = [0 36];
% text(mean(ax.XLim), -0.5*mean(ax.YLim),'Training Duration (Days)','HorizontalAlignment','center'); 
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Survival', 'HorizontalAlignment','center','Rotation',90);
% % text(ax.XLim(1),ax.YLim(2)*1.15, 'D', 'FontSize', 14);
% l = legend([p1 p2], [{'Expt'}, {'Model'}]);
% l.Box = 'off';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Complete ramps with various slopes
% RampDuration = [13 11 9 7 3 1/48];
RampDuration = [13 11 9 7 3 1];
xpos = fliplr(round([10/0.5*16/24, 10/0.625*16/24, 10/0.75*16/24, 10/1*16/24, 10/2*16/24, 1/24],3));

%Cartoon
% nexttile(); hold on;
% imax = length(RampDuration);
% for i = 1:length(RampDuration)
%     Xtraining = [linspace(0,14-i,15-RampDuration(i)) linspace(15-i,14,RampDuration(i))];
%     Xtesting = linspace(14,43,29);
%     Ytraining = [linspace(22,22,15-RampDuration(i)) linspace(22,4,RampDuration(i))];
%     Ytesting = linspace(4,4,length(Xtesting));
% %     plot(Xtraining, Ytraining, 'Color',cmaplines(1,:), 'LineWidth',2, LineStyle=':')
%     plot(Xtraining, Ytraining, 'Color',[cmapjet(round(i/imax*256),:) 1], 'LineWidth',2, LineStyle=':')    
%     plot(Xtesting, Ytesting, 'Color',cmaplines(2,:), 'LineWidth',2)
%     if i == 9
% %         plot(Xtraining, Ytraining, 'Color',[cmapjet(round(i/imax*256),:) 1], 'LineWidth',2)
%     end
% end
% ax = gca; ax.LineWidth = 1; ax.FontSize = 8; ax.YLim = [0 22]; ax.YTick = [4 22]; 
% ax.YTickLabel = {'4';'22'};ax.XLim = [0 28]; ax.XTick = [0 14 21 28]; ax.XTickLabel = {'0'; '14'; '[...]'; '42'};
% 
% text(ax.XLim(2)*0.05,ax.YLim(2)*0.9,{'Ramp with variable slope'}, 'BackgroundColor', [1 1 1 0.75]);
% text(-0.2*mean(ax.XLim), mean(ax.YLim) , ['T (' char(176) 'C)'], 'HorizontalAlignment','center','Rotation',90);
% text(0.9*mean(ax.XLim), -0.5*mean(ax.YLim) ,'T (days)', 'HorizontalAlignment','center');
% text(ax.XLim(1),ax.YLim(2)*1.15, 'P', 'FontSize', 14);

%Other cartoon for ramps with variable slope
nexttile(); hold on;
RampDuration = [13 11 9 7 3 1];
for i = 1:length(RampDuration)
    plot([0 14-RampDuration(i) 42], [22 4 4]+i/5, 'Color',cmaplines(round(256/6*i),:),'LineWidth',1);    
end
ax = gca;
ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'T (C)';

nexttile(); 
% hold on;
RampDurationInv = fliplr(RampDuration);
for R = 1:10
    for i = 1:length(RampDurationInv);
%         cvar = Model2Data.csol.c.*(1+(rand(3,1)-0.5)/4);% Add variability for acclimation, injury and regeneration time constants
        cvar = Model2Data.csol.c;% Add variability for acclimation, injury and regeneration time constants
        temperature1 = linspace(0,43,44*f);
        temperature2 = [linspace(22,22,(15-RampDurationInv(i))*f) linspace(22,4,RampDurationInv(i)*f) linspace(4,4,29*f)]./22;
        temperature2 = AddNoiseToTemp(temperature2,NoiseAmplitude);
        temperature = [temperature1; temperature2];
        sol = ode23(@(t,y) AcclimfunNoise(t,y,cvar,temperature), Model2Data.tspan, Model2Data.y0);
        solpts = deval(sol,Model2Data.tspan);
        ModelRampsData.HealthFull(:,R,i) = solpts(2,:);
        ModelRampsData.AlphaFull(:,R,i) = solpts(1,:);
    end
end

% imax = length(RampDuration);
% tmax = length(ModelRampsData.HealthFull(:,1,1));
% Time = linspace(0,tmax-1,tmax);
% 
% for i = 1:length(RampDuration)
%     plot(Time, squeeze(ModelRampsData.HealthFull(:,1,imax+1-i)),'LineWidth',2, 'Color', [cmapjet(round(i/imax*256),:) 1]);
% end
% 
% ax = gca; ax.YLim = [0 1]; ax.YTick = 0:1:1; ax.XLim = [0 44];
% box off; ax.LineWidth = 1; ax.FontSize = 8;
% text(-0.2*mean(ax.XLim), mean(ax.YLim),'Health model', 'HorizontalAlignment','center','Rotation',90);

h4 = heatmap(squeeze(ModelRampsData.HealthFull(:,1,:))');
h4.ColorbarVisible = 'off'; h4.XLabel = 'Testing Time (Days)'; h4.YLabel = 'Training Duration (Day)';
h4.FontSize = 8;
h4.Title = 'Q. Health model';
h4.XDisplayLabels(:) = {''}; h4.YDisplayLabels(:) = {''};
h4.XDisplayLabels([1 5 10 15 20 25 30 35 40]) = {'0', '5', '10', '15', '20', '25', '30', '35', '40'};
h4.YDisplayLabels([1 2 3 4 5 6]) = {'1', '3', '7', '9', '11', '13'};
% h4.Colormap = cmaphisto;

nexttile(); hold on;
RampsModel2Data = squeeze(ModelRampsData.HealthFull(end,:,:));
[StdRampsModel2Data, MeanRampsModel2Data] = std(RampsModel2Data);
% plot(xpos,MeanRampsModel2Data,'Color',cmaplines(5,:),'LineWidth',2);
% f1 = fill([xpos fliplr(xpos)], [MeanRampsModel2Data + StdRampsModel2Data fliplr(MeanRampsModel2Data - StdRampsModel2Data)], cmaplines(5,:));
% f1.FaceAlpha = 0.5; f1.EdgeAlpha = 0;
xData = repmat(xpos,10,1) + (rand(10,6)-0.5)*2;
yData = RampsModel2Data + (rand(10,6)-0.5)/5;
s1 = scatter(xData,yData,'filled', 'MarkerFaceColor', cmaplines(5,:),...
        'SizeData', 2);

[StdRampsExptData, MeanRampsExptData] = std(Pilot.dataGradientAll);
% plot(xpos,MeanRampsExptData,'Color',cmaplines(1,:),'LineWidth',2);
% f1 = fill([xpos fliplr(xpos)], [MeanRampsExptData + StdRampsExptData fliplr(MeanRampsExptData - StdRampsExptData)], cmaplines(1,:));
% f1.FaceAlpha = 0.5; f1.EdgeAlpha = 0;
ax = gca; ax.YLim = [0 1]; ax.XLim = [0 14]; ax.XTick = [0:2:14];
xData = repmat(xpos,10,1) + (rand(10,6)-0.5)*2;
yData = Pilot.dataGradientAll + (rand(10,6)-0.5)/5;
s1 = scatter(xData,yData,'filled', 'MarkerFaceColor', cmaplines(1,:),...
        'SizeData', 2);
yticks([0:1:1]); box off; ax.LineWidth = 1; ax.FontSize = 8;
text(-0.2*mean(ax.XLim), mean(ax.YLim), {'Final Health'}, 'Rotation',90, 'VerticalAlignment','middle','HorizontalAlignment','center');
% text(mean(ax.XLim),-0.4*mean(ax.YLim),'Training Duration', 'HorizontalAlignment','center');
ax.XLabel.String = 'Training Duration (Day)';

% text(ax.XLim(1),ax.YLim(2)*1.15, 'R', 'FontSize', 14);

b = bar(xpos, [MeanRampsExptData; MeanRampsModel2Data]);
b(1).EdgeColor = 'none'; b(1).FaceColor = cmaplines(1,:); b(1).FaceAlpha = 0.5;
b(2).EdgeColor = 'none'; b(2).FaceColor = cmaplines(5,:); b(2).FaceAlpha = 0.5;

nexttile(); hold on;
SurvivalRamps = sum(ceil(Pilot.dataGradientAll(:,:)),1);
SurvivalModel = sum(ceil(squeeze(ModelRampsData.HealthFull(end,:,:)-0.01)));
b1 = bar([1:1:6],[SurvivalRamps], 'EdgeColor', 'none','FaceColor',cmaplines(1,:), 'FaceAlpha',0.5); hold on;
b2 = bar([1:1:6],[SurvivalModel], 'EdgeColor','none','FaceColor',cmaplines(5,:), 'FaceAlpha',0.5);

ax = gca; ax.XTick = [0:1:6]; ax.XTickLabel = {'', '1', '3', '7', '9', '11', '13'}; ax.Box = 'off'; ax.LineWidth = 1; ax.FontSize = 8; ax.YTick = [0:10:10]; ax.YLim = [0 10];
ax.XLim(1) = 0;
% text(mean(ax.XLim), -0.5*mean(ax.YLim),'Training Duration (Days)','HorizontalAlignment','center'); 
ax = gca; ax.XLabel.String = 'Training Duration (Days)';
text(-0.2*mean(ax.XLim), mean(ax.YLim),'Survival', 'HorizontalAlignment','center','Rotation',90);
% text(ax.XLim(1),ax.YLim(2)*1.15, 'S', 'FontSize', 14);

%Reinitialize first heatmap colormap
cmaphisto = colormap(cmaplist.cmapRedBlue1);
h0.Colormap = cmaphisto;
h1.Colormap = cmaphisto;
h2.Colormap = cmaphisto;
h3.Colormap = cmaphisto;
h4.Colormap = cmaphisto;


% % Heatmap Model Data
% nexttile()
% hold off;
% h = heatmap(squeeze(mean(Model2Data.HealthTruncated(8:12,:,:),1))');
% h.ColorbarVisible = 'off'; h.XLabel = 'D'; h.YLabel = 'TRT';
% h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''}; h.FontSize = 8;
% h.XDisplayLabels([1 5 10 15]) = {'0', '5', '10', '15'};
% h.YDisplayLabels([1 6 11 16]) = {'4', '14', '24', '34'};

fig = gcf;
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(fig, [fig.Name '.pdf'])

%%% SECOND FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figModelPerformance2 = figure('Name', 'Model Performance 2', 'renderer', 'opengl');
figModelPerformance2.WindowState = 'maximized';
clf(figModelPerformance2)
tiledlayout(4,5)

% Deacclimation cartoon
nexttile(); hold on;
cmapturbo = colormap('turbo');
cmaplines = colormap('lines');
StepDuration = (2:6:28);
for i = 1:length(StepDuration)
    plot([-40 -34 -34 -20 -20 0 0 StepDuration(i) StepDuration(i)], [22 22 12 12 4 4 22 22 4]+i/3,'Color',cmapturbo(round(256/length(StepDuration)*i),:),'LineWidth',1)    
    plot([StepDuration(i) 56], [4 4]+i/3,'Color',cmapturbo(round(256/length(StepDuration)*i),:),'LineWidth',1, 'LineStyle','--');    
end
ax = gca; ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'T (C)';
ax.LineWidth = 1; ax.YLim = [0 28];
% ax.YLim = [0 22]; ax.YTick = [4 22]; 
ax.XTick = [-40:20:60]; ax.FontSize = 10; ax.XLim = [-40 56]; ax.XTickLabel = [{'0'},{'20'},{'0'},{'20'},{'40'},{'60'}];
% ax.Title.String = 'A';

% Generate Deacclimation data
cvar = Model2Data.csol.c; % No var in coefficients
tspan = 0:1:99;
temperature1 = linspace(0,99,100*f);
health1 = [];
alpha1 = [];
health2 = [];
alpha2 = [];

% Without noise
% Acclimated - deaclimated animals
for i = 1:28    
    temperature2 = [linspace(22,22,i*f)]./22; %Duration = 13 days, training temperature = 10C
    % if i == 18
    %     temperature2 = [linspace(4,4,i*f)]./22; %Duration = 13 days, training temperature = 10C
    % end
    temperature2 = [temperature2 linspace(4,4,length(temperature1)-length(temperature2))./22];
    temperature = [temperature1; temperature2]; %AddNoiseToTemp(temperature2)
    y0 = [4/22 1]; %Start with a fully acclimated animal
    sol = ode23(@(t,y) Acclimfun(t,y,cvar,temperature), tspan, y0);
    solpts = deval(tspan,sol);
    alpha1 = [alpha1; solpts(1,i:i+27)];
    health1 = [health1; solpts(2,i:i+27)];    
end

% Naive animals
for i = 1:6
    temperature2 = [linspace(22,22,0*f)]./22; %Duration = 13 days, training temperature = 10C
    temperature2 = [temperature2 linspace(4,4,length(temperature1)-length(temperature2))./22];
    temperature = [temperature1; temperature2]; %AddNoiseToTemp(temperature2)
    y0 = [1 1]; %Start with a fully deacclimated animal
    sol = ode23(@(t,y) Acclimfun(t,y,cvar,temperature), tspan, y0);
    solpts = deval(tspan,sol);
    alpha1 = [alpha1; solpts(1,1:1+27)];
    health1 = [health1; solpts(2,1:1+27)];
end

% With noise
% Acclimated - deaclimated animals
for i = 1:28    
    temperature2 = [linspace(22,22,i*f)]./22; %Duration = 13 days, training temperature = 10C
    temperature2 = [temperature2 linspace(4,4,length(temperature1)-length(temperature2))./22];
    temperature = [temperature1; AddNoiseToTemp(temperature2,NoiseAmplitude)]; %AddNoiseToTemp(temperature2)
    y0 = [4/22 1]; %Start with a fully acclimated animal
    sol = ode23(@(t,y) AcclimfunNoise(t,y,cvar,temperature), tspan, y0);
    solpts = deval(tspan,sol);
    alpha2 = [alpha2; solpts(1,i:i+27)];
    health2 = [health2; solpts(2,i:i+27)];    
end

% Naive animals
for i = 1:6
    temperature2 = [linspace(22,22,0*f)]./22; %Duration = 13 days, training temperature = 10C
    temperature2 = [temperature2 linspace(4,4,length(temperature1)-length(temperature2))./22];
    temperature = [temperature1; AddNoiseToTemp(temperature2,NoiseAmplitude)]; %AddNoiseToTemp(temperature2)
    y0 = [1 1]; %Start with a fully deacclimated animal
    sol = ode23(@(t,y) AcclimfunNoise(t,y,cvar,temperature), tspan, y0);
    solpts = deval(tspan,sol);
    alpha2 = [alpha2; solpts(1,1:1+27)];
    health2 = [health2; solpts(2,1:1+27)];
end

% Heatmap - experimental data
nexttile();
h0 = heatmap(Deacclimation);
h0.ColorbarVisible = 'off'; h0.XLabel = 'Testing Time (days)'; h0.YLabel = 'Training Duration (Days)';
h0.FontSize = 8;
h0.XDisplayLabels(:) = {''}; h0.YDisplayLabels(:) = {''};
h0.XDisplayLabels([1 10 20]) = {'1', '10', '20'};
h0.YDisplayLabels([1 10 20 29]) = {'1', '10', '20', 'N'};
cmaphisto = colormap(cmaplist.cmapRedBlue1);
h0.Colormap = cmaphisto;

% Heatmap - model without noise
nexttile();
h0 = heatmap(health1);
h0.ColorbarVisible = 'off'; h0.XLabel = 'Testing Time (days)'; h0.YLabel = 'Training Duration (Days)';
h0.FontSize = 8;
h0.XDisplayLabels(:) = {''}; h0.YDisplayLabels(:) = {''};
h0.XDisplayLabels([1 10 20]) = {'1', '10', '20'};
h0.YDisplayLabels([1 10 20 29]) = {'1', '10', '20', 'N'};
% cmaphisto = h0.Colormap;
cmaphisto = colormap(cmaplist.cmapRedBlue1);
h0.Colormap = cmaphisto;

% nexttile(); hold on;
% for i = 1:size(health1,1)
%     plot(health1(i,:));
%     pause()
% end

% Heatmap - model with noise
nexttile();
h0 = heatmap(health2);
h0.ColorbarVisible = 'off'; h0.XLabel = 'Testing Time (days)'; h0.YLabel = 'Training Duration (Days)';
h0.FontSize = 8;
h0.XDisplayLabels(:) = {''}; h0.YDisplayLabels(:) = {''};
h0.XDisplayLabels([1 10 20]) = {'1', '10', '20'};
h0.YDisplayLabels([1 10 20 29]) = {'1', '10', '20', 'N'};
% cmaphisto = h0.Colormap;
cmaphisto = colormap(cmaplist.cmapRedBlue1);
h0.Colormap = cmaphisto;


% Plot - can we infer acclimation index?
nexttile(); hold on;
plot(mean(Deacclimation(:,1:end),2),'Color',cmaplines(1,:),'LineWidth',2)
plot(mean(health1(:,1:end),2), 'Color',cmaplines(5,:),'LineWidth',2)
plot(mean(health2(:,1:end),2), 'Color',cmaplines(5,:),'LineWidth',2, 'LineStyle', ':')
ax = gca; ax.LineWidth = 1; ax.FontSize = 8; ax.YLim = [0 1]; ax.YTick = [0 1]; 
ax.XLim = [0 28];
text(-0.2*mean(ax.XLim), mean(ax.YLim) , ['Average Health'], 'HorizontalAlignment','center','Rotation',90);
text(mean(ax.XLim), -0.5*mean(ax.YLim) ,'D', 'HorizontalAlignment','center');

fig = gcf;
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(fig, [fig.Name '.pdf'])

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

%%% THIRD FIGURE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figModelPerformance3 = figure('Name', 'Model Performance 3', 'WindowState', 'maximized', 'renderer', 'opengl');
tiledlayout(4,5)

% Cartoon for step at 12C for 14 Days
nexttile(); hold on;
plot([0 14 14 42], [12 12 4 4], 'Color', [0 0 0], 'LineWidth',2)

ax = gca; ax.XLabel.String = 'Time (Days)'; ax.YLabel.String = 'T (C)';
ax.YLim = [0 30];

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
    h.Title = 'Health expt';

    nexttile()
    h2 = heatmap(squeeze(Model1Data.HealthTruncated(:,:,i))');
    h2.ColorbarVisible = 'off'; h2.XLabel = 'Testing Time (Days)'; h2.YLabel = 'Training Duration';
    h2.FontSize = 8;
    h2.XDisplayLabels(:) = {''}; h2.YDisplayLabels(:) = {''};
    h2.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    h2.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    h2.Title = 'Health model';

    nexttile()
    h3 = heatmap(squeeze(Model2Data.HealthTruncated(:,:,i))');
    h3.ColorbarVisible = 'off'; h3.XLabel = 'Testing Time (Days)'; h3.YLabel = 'Training Duration';
    h3.FontSize = 8;
    h3.XDisplayLabels(:) = {''}; h3.YDisplayLabels(:) = {''};
    h3.XDisplayLabels([1 5 10 15 20 25]) = {'0', '5', '10', '15', '20', '25'};
    h3.YDisplayLabels([1 6 11 15]) = {'0', '5', '10', '14'};
    h3.Title = 'Health model with noise';

    nexttile(); hold on;
    plot(linspace(0,14,15), squeeze(mean(CstData3D(:,:,i),1)), 'LineWidth',2, 'Color', cmaplines(1,:));
    plot(linspace(0,14,15), squeeze(mean(Model1Data.HealthTruncated(:,:,i),1)), 'LineWidth',2, 'Color', cmaplines(5,:));
    plot(linspace(0,14,15), squeeze(mean(Model2Data.HealthTruncated(:,:,i),1)), 'LineWidth',2, 'Color', cmaplines(5,:), 'LineStyle',':');
    box off; ax = gca ; ax.LineWidth = 1; ax.FontSize = 8;
    xlim([0 15]); xticks([0:5:15]); ylim([0 1]); yticks([0:1:1]);
    text(-0.2*mean(ax.XLim), mean(ax.YLim), {'avg. health'}, 'Rotation',90, 'VerticalAlignment','middle','HorizontalAlignment','center');
    % text(mean(ax.XLim),-0.4*mean(ax.YLim),'D', 'HorizontalAlignment','center');
    ax = gca; ax.XLabel.String = 'Training Duration (Days)';
    % text(ax.XLim(1),ax.YLim(2)*1.15, 'N', 'FontSize', 14);

end
colormap(cmaplist.cmapRedBlue1)

f = gcf;
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(f, [f.Name '.pdf'])

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
nexttile()
plot(0:0.1:43, StepsHealthFull); hold on;
plot(temperature1, temperature2);
ax = gca; ax.YLim = [0 1];
ax.Title.String = ['13C Step'];
% text(20, 0.4, ['Mean temp. = ' num2str(mean(temperature2),4)]);

for pd = pdindex;
    %Pulses
    temperature1 = linspace(0,43,44*f);
    temperature2b = [repmat([linspace(22,22,pd) linspace(4,4,pd)],1,floor(f/(2*pd)*D)) linspace(4,4,29*f)];
    temperature2a = linspace(22,22,length(temperature1)-length(temperature2b));
    temperature2 = [temperature2a temperature2b]./22;
    temperature = [temperature1; temperature2];
    sol = ode23(@(t,y) Acclimfun(t,y,cvar,temperature), Model2Data.tspan, Model2Data.y0);
    solpts = deval(0:0.1:43,sol);
    PulsesHealthFull = solpts(2,:);   
    
    nexttile()
    plot(0:0.1:43, PulsesHealthFull); hold on;
    % plot(temperature1, temperature2);
    ax = gca; ax.YLim = [0 1];
    ax.Title.String = ['Pulse duration (hr) = ' num2str(pd)];
    % text(20, 0.4, ['Mean temp. = ' num2str(mean(temperature2),4)]);
    
    % error = [error sum((PulsesHealthFull - StepsHealthFull).^2)];
end

end