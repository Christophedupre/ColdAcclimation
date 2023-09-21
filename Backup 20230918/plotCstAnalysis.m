function plotCstAnalysis(CstData3D, Cst10C13D30Reps, Pilot)

figCstDataAnalysis = figure('Name', 'Constant Data Analysis', 'renderer', 'painters');
figCstDataAnalysis.WindowState = 'maximized';
clf(figCstDataAnalysis)
nrows = 4; ncols = 5;
tiledlayout(nrows, ncols);

cmaplines = colormap("lines");
cmapjet = colormap("jet");


% Variability: 30 animals at T = 13, theta = 10
nexttile(); hold on;
Time = 0:1:28;
for j = 1:size(Cst10C13D30Reps,2) % Over all animals
    s = scatter(Time, Cst10C13D30Reps(:,j));
    s.MarkerFaceColor = cmapjet(1,:); s.MarkerFaceAlpha = 1; s.MarkerEdgeColor = 'none';
    s.SizeData = 10; s.XJitter = 'rand';
    yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 10;
    xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1]); yticks([0:0.5:1]);  
end
   
p1 = plot(Time, mean(Cst10C13D30Reps,2), 'Color', [0 0 0], 'LineWidth', 2);
sd = std(Cst10C13D30Reps');%./sqrt(size(Cst10C13D30Reps,1));
fi = fill([Time fliplr(Time)], [mean(Cst10C13D30Reps,2)'-sd fliplr(mean(Cst10C13D30Reps,2)'+sd)], cmaplines(1,:));
fi.FaceAlpha = 0.5; fi.EdgeAlpha = 0;
text(0*ax.XLim(2),ax.YLim(2)*1.15, 'Variability for \theta = 10, D = 13 (n = 30)');


%fit sigmoid to acclimation performance
TrainingDurations = 0:size(CstData3D,2)-1;
TrainingPerformance = squeeze(mean(CstData3D,1))';
tau = []; of = [];
for i = 1:size(TrainingPerformance,1) %to each training temperature
    %fit tau and offset
    str = '1./(1+exp(-(x-of)/tau))';    
    ft = fittype(str, 'independent', 'x', 'dependent', 'y'); %fit function (sigmoid)
    [xData, yData] = prepareCurveData(TrainingDurations,TrainingPerformance(i,:)); %Data to fit
    opts = fitoptions('Method', 'NonlinearLeastSquares');
    opts.StartPoint = [6 13]; % [of tau]
    opts.Lower = [0 0]; %Cannot have negative taus or offset
    opts.Upper = [inf inf]; 
    [fitresultconstant, gof] = fit(xData, yData, ft, opts); %compute the fit
    %extract parameters tau and offset
    tau = [tau fitresultconstant.tau]; 
    of = [of fitresultconstant.of];     
end

%%% TEST FOR DIFFERENCE IN TRAINING EFFICIENCY BETWEEN PROTOCOLS %%%%%%%%%%
nexttile(); hold on;
for i = 1:length(CstData3D(1,1,:))
    xpos = linspace(i*2+2,i*2+2,15);
    s = scatter(xpos,mean(CstData3D(:,:,i),1));
    s.MarkerFaceColor = cmapjet(1,:); s.MarkerFaceAlpha = 1; s.MarkerEdgeColor = 'none';
    s.SizeData = 10; %s.XJitter = 'rand';
end
plot(linspace(4,34,16),squeeze(mean(mean(CstData3D(:,:,:),1))), 'Color',cmaplines(1,:),'LineWidth',2);
text(0*ax.XLim(2),ax.YLim(2)*1.15, 'Difference in training efficiency?');
ax = gca; ax.XLabel.String = '\theta'; 
ax.YLabel.String = 'Avg health';

ttestresults = zeros(length(CstData3D(1,1,:)),length(CstData3D(1,1,:)));
for i = 1:length(CstData3D(1,1,:))
    for j = 1:length(CstData3D(1,1,:))
        [H,ttestresults(i,j)] = ttest(mean(CstData3D(:,:,i),1),mean(CstData3D(:,:,j),1));
    end
end

nexttile()
h = heatmap(ttestresults);
h.ColorbarVisible = 'off'; 
h.XLabel = '\theta'; h.YLabel = '\theta';
h.FontSize = 8;
h.Title = ['Ttest (P-values)'];
h.XDisplayLabels(:) = {''}; h.YDisplayLabels(:) = {''};
h.XDisplayLabels([1 6 11 16]) = {'4', '14', '24', '34'};
h.YDisplayLabels([1 6 11 16]) = {'4', '14', '24', '34'};

nexttile(); ax = gca; ax.Visible = "off";
nexttile(); ax = gca; ax.Visible = "off";

%Plot fit results
for i = [4 1];
    nexttile(); hold on;
    plot(TrainingDurations, ft(of(i), tau(i), TrainingDurations), 'Color', cmapjet(round(i/length(CstData3D(1,1,:))*256),:), 'LineWidth', 2);
    plot(TrainingDurations, TrainingPerformance(i,:), "*", 'Color', cmapjet(round(i/length(CstData3D(1,1,:))*256),:))
    box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 10;
    xlabel('D'); ylabel('Performance'); 
    ax = gca; ax.YLim = [0 1.2];
    text(ax.XLim(2)*0.05,ax.YLim(2)*0.92,{['\theta = ' num2str(i*2+2)...
        ';  \tau = ' num2str(tau(i), 1) ';  of = ' num2str(of(i),1)]}); 
    
    if i == 4
        text(0*ax.XLim(2),ax.YLim(2)*1.15, 'Acclimation process is sigmoid');
    end

end

%Plot tau and offset
nexttile()
stot = [];
for i=1:length(tau)
    s = scatter(tau(i), of(i)); hold on;
    s.MarkerFaceColor = cmapjet(round(i/length(tau)*256),:);
    s.MarkerFaceAlpha = 0.5; s.MarkerEdgeColor = 'none';
    box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 10;
    stot = [stot s];
end
ylabel('Offset'); xlabel('tau'); ylim([0 15]); xlim([0 15]);

%Plot tau and offset (log scale)
nexttile()
stot = [];
for i=1:length(tau)
    s = scatter(log(tau(i))/log(10), log(of(i))/log(10)); hold on;
    s.MarkerFaceColor = cmapjet(round(i/length(tau)*256),:);
    s.MarkerFaceAlpha = 0.5; s.MarkerEdgeColor = 'none';
    box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 10;
    stot = [stot s];
end
ylabel('log(Offset)'); xlabel('log(tau)'); %ylim([0 15]); xlim([0 15]);


% Plot legend box
str = num2str(((1:1:length(tau))*2+2)');
l = legend(stot, str, 'NumColumns', 4, 'Box','off');
l.Position = [0.8244    0.0905    0.1351    0.1490];
l.Title.String = 'Training Temperatures';

f = gcf;
exportgraphics(f, [f.Name '.pdf'])


clear i s l str stot xData yData ax tau of ft opts fitresultconstant gof nrows ncols TrainingDurations TrainingPerformance Time

end
