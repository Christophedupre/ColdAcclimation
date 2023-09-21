function modelexplanation(Model1Data, Model2Data);
load('data.mat')
figCst = figure('Name', 'Model Explanation', 'renderer', 'painters'); 
figCst.WindowState = 'maximized';
tiledlayout(2,2)
cmaplines = colormap('lines');
% Current model
nexttile()
tspan = Model1Data.tspan;
y0 = [1 1];
cs = Model1Data.csol.c;
NoiseAmplitude = Model1Data.NoiseAmplitude;
f = Model1Data.TimeResolution;
theta = 10;
D = 14; %Training duration goes from 1:15, therefore D = 14 corresponds to training duration of 13 days.

%Without variability
temperature1 = linspace(0,43,44*f);
temperature2 = [linspace(22,22,(15-(D-1))*f) linspace(theta*2+2,theta*2+2,(D-1)*f) linspace(4,4,29*f)]./22;
temperature = [temperature1; temperature2]; %Variability in temperature sensed
sol = ode23(@(t,y) Acclimfun(t,y,cs,temperature), tspan,y0);
solpts = deval(sol, tspan);
p1 = plot(tspan,solpts(1,:), 'LineWidth',2,'Color',cmaplines(2,:)); hold on;
p2 = plot(tspan,solpts(2,:), 'LineWidth',2,'Color',cmaplines(5,:));
p3 = plot(temperature1,temperature2, 'LineWidth',2,'Color',[0 0 0]); 
p4 = plot(14:1:42, mean(CstData3D(:,7:15,5),2), 'LineWidth',2,'Color',cmaplines(1,:));
% p4 = plot(15:1:43, CstData3D(:,15,5), 'LineWidth',2,'Color',cmaplines(1,:));
ax = gca; ax.XLim = [tspan(1) tspan(end)]; ax.YLim = [0 1.1]; ax.Box = 'off';
ax.LineWidth = 2; ax.FontSize = 15;
text(-0.1*ax.XLim(2),mean(ax.YLim),'Health, Acclimation, Temperature','HorizontalAlignment','center','Rotation',90,'FontSize', 15);
text(mean(ax.XLim),-0.2*ax.YLim(2),'Time (Days)','HorizontalAlignment','center','Rotation',0,'FontSize',15);
str = '$\frac{d\alpha}{dt} = (\theta - \alpha)*c_1$; $\frac{dH}{dt} = (\theta - \alpha)*c_2 + (H* \theta)*c_3$';
text(0.01*ax.XLim(2),1.1*ax.YLim(2),str, 'Interpreter','latex','FontSize',15)

l = legend([p1 p2 p3 p4],{'\alpha (Acclimation)','H (Health)','\theta (Temperature)','Raw Data'},'Location','east');
l.EdgeColor = 'none'; l.Color = [1 1 1 0.1];

%With variability
nexttile()
temperature1 = linspace(0,43,44*f);
temperature2 = [linspace(22,22,(15-(D-1))*f) linspace(theta*2+2,theta*2+2,(D-1)*f) linspace(4,4,29*f)]./22;
temperature2wnoise = AddNoiseToTemp(temperature2, NoiseAmplitude);
temperature = [temperature1; temperature2wnoise]; %Variability in temperature sensed
sol = ode23(@(t,y) AcclimfunNoise(t,y,cs,temperature), tspan,y0);
solpts = deval(sol, tspan);
p1 = plot(tspan,solpts(1,:), 'LineWidth',2,'Color',cmaplines(2,:)); hold on;
p2 = plot(tspan,solpts(2,:), 'LineWidth',2,'Color',cmaplines(5,:));
temp2 = AddNoiseToTemp(temperature2, NoiseAmplitude/10);
temp2 = temp2(round(linspace(1,length(temp2),length(tspan))));
p3 = plot(tspan,temp2, 'LineWidth',2,'Color',[0 0 0]); 
% p3 = plot(temperature1,temperature2wnoise, 'LineWidth',2,'Color',[0 0 0]); 
p4 = plot(14:1:42, mean(CstData3D(:,7:15,5),2), 'LineWidth',2,'Color',cmaplines(1,:));
% p4 = plot(15:1:43, CstData3D(:,15,5), 'LineWidth',2,'Color',cmaplines(1,:));
ax = gca; ax.XLim = [tspan(1) tspan(end)]; ax.YLim = [0 1.1]; ax.Box = 'off';
ax.LineWidth = 2; ax.FontSize = 15;
text(-0.1*ax.XLim(2),mean(ax.YLim),'Health, Acclimation, Temperature','HorizontalAlignment','center','Rotation',90,'FontSize',15);
text(mean(ax.XLim),-0.2*ax.YLim(2),'Time (Days)','HorizontalAlignment','center','Rotation',0,'FontSize',15);
str = 'Idem, With Noise in Temperature Sensing';
text(0.01*ax.XLim(2),1.1*ax.YLim(2),str, 'Interpreter','latex','FontSize',15)

% nexttile()
% ax = gca; ax.Visible = 'off';

fig = gcf;
% f.Color = [1 1 1];
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(fig, [fig.Name '.pdf'])

fig = gcf;
% f.Color = [1 1 1];
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(fig, [fig.Name '.pdf'])


figModelExpl = figure('Name', 'Model Explanation 2', 'renderer', 'painters'); 
figModelExpl.WindowState = 'maximized';
tiledlayout(4,4)
%%% PLOT ALPHA VS HEALTH （PRELIM. SLOPE FIELD) %%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = Model1Data.tspan;
nexttile(); hold on;
cmapjet = colormap('jet');
cmaplines = colormap('lines');
imax = max(max(max((Model1Data.AlphaFull))));
for theta = 1:length(CstData3D(1,1,:))
    for D = 1:length(CstData3D(1,:,1))        
        if(Model1Data.HealthFull(end,D,theta) <0.1)
            if (D == 14) 
                plot(Model1Data.AlphaFull(15:end,D,theta),Model1Data.HealthFull(15:end,D,theta), 'Color', cmaplines(2,:),'LineWidth',0.5);                
            end
            p1 = scatter(Model1Data.AlphaFull(15,D,theta),Model1Data.HealthFull(15,D,theta),...
                'MarkerFaceColor',cmaplines(2,:),'MarkerFaceAlpha',0.5, 'MarkerEdgeColor','none');
        else
            if (D == 14) 
                plot(Model1Data.AlphaFull(15:end,D,theta),Model1Data.HealthFull(15:end,D,theta), 'Color', cmaplines(1,:),'LineWidth',0.5);
            end
            p2 = scatter(Model1Data.AlphaFull(15,D,theta),Model1Data.HealthFull(15,D,theta),...
                'MarkerFaceColor',cmaplines(1,:),'MarkerFaceAlpha',0.5, 'MarkerEdgeColor','none');
        end
        box off; ax = gca; ax.LineWidth = 1;
    end
end
p3 = scatter([-1],[-1],'MarkerFaceColor',[0 0 0],'MarkerFaceAlpha',0.5, 'MarkerEdgeColor','none');
p4 = plot(-1,-1,'Color', [0 0 0],'LineWidth',0.5);
ax = gca; ax.YTick = [0 1] ;ax.XLabel.String = '\alpha'; ax.YLabel.String = 'Health'; ax.XLim(1) = 0; ax.YLim(1) = 0;
text(ax.XLim(1),ax.YLim(2)*1.15, 'B', 'FontSize', 14);
l = legend([p4(1) p3(1) p1(1) p2(1)],{'Test phase','T.P. onset','Died','Survived'}); l.Box = 'off'; l.Location = 'southeast';

%%% FULL SLOPE FIELD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Show effect of dropping an animal to 4C given a specific initial contition
nexttile(); hold on;
cmaplines = colormap('lines');

VectorSurvived = [];
VectorDied = [];
tspan = Model1Data.tspan;
for AlphaInitial = 0.2:0.1:1.5;
    for HealthInitial = 0.2:0.1:1;
        y0 = [AlphaInitial,HealthInitial];
        cs = Model1Data.csol.c;
        temperature1 = linspace(0,43,44);
        temperature2 = [linspace(4,4,44)]./22;
        temperature = [temperature1; temperature2]; %Can add Variability in temperature sensed with AddNoiseToTemp(temperature2)
        sol = ode23(@(t,y) Acclimfun(t,y,cs,temperature), tspan,y0);
        solpts = deval(sol, tspan);
        DeltaAlpha = solpts(1,2) - solpts(1,1);
        DeltaHealth = solpts(2,2) - solpts(2,1);        
        if (abs(AlphaInitial - 0.6) <0.1 && abs(HealthInitial-1) <0.1)
            AlphaFullSurvived = solpts(1,:);
            HealthFullSurvived = solpts(2,:);
        end
        if (abs(AlphaInitial - 1) <0.1 && abs(HealthInitial-1) <0.1)
            AlphaFullDied = solpts(1,:);
            HealthFullDied = solpts(2,:);
        end

        if(solpts(2,:)>0.2)
            VectorSurvived = [VectorSurvived; AlphaInitial DeltaAlpha HealthInitial DeltaHealth];
        else
            VectorDied = [VectorDied; AlphaInitial DeltaAlpha HealthInitial DeltaHealth];
        end
    end
end

% Plot
quiver(VectorSurvived(:,1),VectorSurvived(:,3),VectorSurvived(:,2),VectorSurvived(:,4),'Color',cmaplines(1,:),'LineWidth',1);
% plot(Model1Data.AlphaFull(15:end,14,4),Model1Data.HealthFull(15:end,14,4), 'Color', cmaplines(1,:),'LineWidth',2);
plot(AlphaFullSurvived, HealthFullSurvived, 'Color', cmaplines(1,:),'LineWidth',2);
quiver(VectorDied(:,1),VectorDied(:,3),VectorDied(:,2),VectorDied(:,4),'Color',cmaplines(2,:),'LineWidth',1);
% plot(Model1Data.AlphaFull(15:end,14,12),Model1Data.HealthFull(15:end,14,12), 'Color', cmaplines(2,:),'LineWidth',2);
plot(AlphaFullDied, HealthFullDied, 'Color', cmaplines(2,:),'LineWidth',2);

ax = gca; ax.YLim = [0 1]; ax.YTick = [0 1]; ax.LineWidth = 1;
text(ax.XLim(1),ax.YLim(2)*1.15, 'C', 'FontSize', 14);
% ax.XLim = [0 1];
text(mean(ax.XLim), -0.5*mean(ax.YLim) ,'\alpha', 'HorizontalAlignment','center');
text(-0.2*mean(ax.XLim), mean(ax.YLim) , 'Health', 'HorizontalAlignment','center','Rotation',90);

nexttile()
% timepoints = 
plot([0 1 1 14 14 42 42 60], [22 22 10 10 4 4 22 22])


% GENERATE 30 REPEATS WITH AND WITHOUT NOISE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tspan = Model1Data.tspan;
y0 = [1 1];
cs = Model1Data.csol.c;
NoiseAmplitude = Model1Data.NoiseAmplitude;
f = Model1Data.TimeResolution;
theta = 4;
D = 14;

% Generate model data with noise
solpts = zeros(44,30);
for i = 1:30
    temperature1 = linspace(0,43,44*f);
    temperature2 = [linspace(22,22,(15-(D-1))*f) linspace(theta*2+2,theta*2+2,(D-1)*f) linspace(4,4,29*f)]./22;
    temperature = [temperature1; AddNoiseToTemp(temperature2, NoiseAmplitude)]; %Variability in temperature sensed
    sol = ode23(@(t,y) AcclimfunNoise(t,y,cs,temperature), tspan,y0);
    solptstemp = deval(sol, tspan);
    solpts(:,i) = solptstemp(2,:);
end

ModelDataWithNoise = solpts;

% Generate model data without noise
solpts = zeros(44,30);
for i = 1:30
    temperature1 = linspace(0,43,44*f);
    temperature2 = [linspace(22,22,(15-(D-1))*f) linspace(theta*2+2,theta*2+2,(D-1)*f) linspace(4,4,29*f)]./22;
    temperature = [temperature1; temperature2]; %Variability in temperature sensed
    sol = ode23(@(t,y) Acclimfun(t,y,cs,temperature), tspan,y0);
    solptstemp = deval(sol, tspan);
    solpts(:,i) = solptstemp(2,:);
end

ModelDataWithoutNoise = solpts;

nexttile(); hold on;
%Plot experimental data
p1 = plot(0:1:28, mean(Cst10C13D30Reps,2), 'Color', [0 0 0], 'LineWidth', 2);
sd = std(Cst10C13D30Reps');%./sqrt(size(Cst10C13D30Reps,1));
fi = fill([0:1:28 fliplr(0:1:28)], [mean(Cst10C13D30Reps,2)'-sd fliplr(mean(Cst10C13D30Reps,2)'+sd)], cmaplines(1,:));
fi.FaceAlpha = 0.5; fi.EdgeAlpha = 0;

XData = repmat([1:1:29],30,1); XData = XData + (rand(30,29)-0.5)./1;
YData = Cst10C13D30Reps'; YData = YData + (rand(30,29)-0.5)./5;
s = scatter(XData,YData,'MarkerFaceColor',...
        cmaplines(1,:), 'MarkerEdgeColor','none','SizeData',2);

yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 10;
xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1.2]); yticks([0:1]);
text(-0.2*mean(ax.XLim), mean(ax.YLim),'Health', 'HorizontalAlignment','center','Rotation',90);
text(ax.XLim(2)*0.05,ax.YLim(2)*0.9,{'T = 10C, Training Dur. = 13D, 30 reps'}, 'BackgroundColor', [1 1 1 0.75]);
text(mean(ax.XLim), -0.5*mean(ax.YLim),'Testing Time (days)','HorizontalAlignment','center'); 
text(ax.XLim(1),ax.YLim(2)*1.15, 'E', 'FontSize', 14);

nexttile(); hold on;
%Plot model data without noise
p2 = plot(0:1:29, mean(ModelDataWithoutNoise(15:44,:),2), 'Color', cmaplines(5,:), 'LineWidth', 2);
sd2 = std(ModelDataWithoutNoise');%./sqrt(size(Cst10C13D30Reps,1));
sd2 = sd2(:,15:44);
fi2 = fill([0:1:29 fliplr(0:1:29)], [mean(ModelDataWithoutNoise(15:44,:),2)'-sd2 fliplr(mean(ModelDataWithoutNoise(15:44,:),2)'+sd2)], cmaplines(5,:));
fi2.FaceAlpha = 0.5; fi2.EdgeAlpha = 0;

XData = repmat([0:1:29],30,1); %XData = XData + (rand(30,30)-0.5)./1;
YData = ModelDataWithoutNoise(15:44,:)'; %YData = YData + (rand(30,30)-0.5)./5;
s = scatter(XData,YData,'MarkerFaceColor',...
        cmaplines(5,:), 'MarkerEdgeColor','none','SizeData',2);

yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 10;
xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1.2]); yticks([0:1]);
text(-0.2*mean(ax.XLim), mean(ax.YLim),'Health', 'HorizontalAlignment','center','Rotation',90);
% text(ax.XLim(2)*0.05,ax.YLim(2)*0.9,{'\theta = 10, TD = 14, 30 reps'}, 'BackgroundColor', [1 1 1 0.75]);
text(mean(ax.XLim), -0.5*mean(ax.YLim),'Testing Time (days)','HorizontalAlignment','center'); 
text(ax.XLim(1),ax.YLim(2)*1.15, 'E', 'FontSize', 14);

nexttile(); hold on;
%Plot model data with noise
p2 = plot(0:1:29, mean(ModelDataWithNoise(15:44,:),2), 'Color', cmaplines(5,:), 'LineWidth', 2);
sd2 = std(ModelDataWithNoise');%./sqrt(size(Cst10C13D30Reps,1));
sd2 = sd2(:,15:44);
fi2 = fill([0:1:29 fliplr(0:1:29)], [mean(ModelDataWithNoise(15:44,:),2)'-sd2 fliplr(mean(ModelDataWithNoise(15:44,:),2)'+sd2)], cmaplines(5,:));
fi2.FaceAlpha = 0.5; fi2.EdgeAlpha = 0;

XData = repmat([0:1:29],30,1); XData = XData + (rand(30,30)-0.5)./1;
YData = ModelDataWithNoise(15:44,:)'; YData = YData + (rand(30,30)-0.5)./5;
s = scatter(XData,YData,'MarkerFaceColor',...
        cmaplines(5,:), 'MarkerEdgeColor','none','SizeData',2);

yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 10;
xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1.2]); yticks([0:1]);
text(-0.2*mean(ax.XLim), mean(ax.YLim),'Health', 'HorizontalAlignment','center','Rotation',90);
% text(ax.XLim(2)*0.05,ax.YLim(2)*0.9,{'\theta = 10, TD = 14, 30 reps'}, 'BackgroundColor', [1 1 1 0.75]);
text(mean(ax.XLim), -0.5*mean(ax.YLim),'Testing Time (days)','HorizontalAlignment','center'); 
text(ax.XLim(1),ax.YLim(2)*1.15, 'F', 'FontSize', 14);


f = gcf;
% f.Color = [1 1 1];
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(f, [f.Name '.pdf'])

end