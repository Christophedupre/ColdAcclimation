function plotHypotheses()
    load('data.mat');
    %%%
    fig = figure('Name','4 main questions','WindowState', 'maximized');
    cmaplines = colormap('lines');
    tiledlayout(4,4)

    nexttile(); hold on;
    A = Deacclimation(29:34,:);
    B = CstData3D(:,:,10); B(end,:) = [];
    C = [A;B';JustPCR(1:end-1,:)';NoPCR(1:end-1,:)']; %Pool data from naive animals and animals from step experiment trained at 22C;
    HealthNaive = C;    
    survival = sum(HealthNaive>0,1);
    survival = survival./max(survival);
    p1 = plot([0 14 14 42], [22 22 4 4], 'Color','k', 'LineWidth', 2);
    ax = gca; ax.YLim = [2 22*1.05]; ax.Box = 'off'; ax.Title.String = ['Cold is lethal (n=' num2str(size(HealthNaive,1)) ')']; ax.FontSize = 9;
    ax.XLabel.String = 'Time (days)'; ax.YLabel.String = 'Temperature (C)';
    yyaxis right; ax = gca; ax.YLabel.String = 'Survival (%)'; ax.YColor = [0 0 0];
    p2 = plot([0 (1:1:28)+14], [survival(1) survival], 'Color',cmaplines(2,:), 'LineWidth', 2);

    TimeHealthNaive = (0:1:size(HealthNaive,2)-1)+13;
    TimeHealthNaive = repmat(TimeHealthNaive,size(HealthNaive,1),1);
    TimeHealthNaive = TimeHealthNaive + rand(size(TimeHealthNaive))./1;
    HealthNaive = HealthNaive + (rand(size(HealthNaive))-0.5)./4;
    p4 = scatter(TimeHealthNaive,HealthNaive,'filled', 'MarkerFaceColor', cmaplines(1,:),...
        'SizeData', 1);
    p5 = plot(TimeHealthNaive(1,:),mean(HealthNaive,1), 'Color',cmaplines(1,:),'LineWidth',2, 'LineStyle','-');

    ax.YLim = [0 max(survival)*1.05]; ax.LineWidth = 2;
    l = legend([p1 p2], {'Temp.', 'Survival'}); l.Box = 'off';
    text(ax.XLim(1),ax.YLim(2)*1.1, 'D', 'FontSize', 12);

    % Acclimation is possible, without jitter
    nexttile(); hold on;
    % survival = sum(CstData3D(:,:,5)>0,2);
    survival = sum(Cst10C13D30Reps'>0,2);
    survival = survival./max(survival);
    plot([0 1 1 14 14 42], [22 22 10 10 4 4], 'Color','k', 'LineWidth', 2)
    ax = gca; ax.YLim = [2 22*1.05]; ax.Box = 'off'; ax.Title.String = 'Acclimation is possible (n=30)'; ax.FontSize = 9;
    ax.XLabel.String = 'Time (days)'; ax.YLabel.String = 'Temperature (C)';

    yyaxis right; ax = gca; ax.YLabel.String = 'Survival (%)'; ax.YColor = [0 0 0];
    plot([0 (0:1:28)+14], survival, 'Color',cmaplines(2,:), 'LineWidth', 2)

    TimeCst10C13D = (0:1:size(Cst10C13D30Reps,1)-1)+13;
    TimeCst10C13D = repmat(TimeCst10C13D,size(Cst10C13D30Reps,2),1);
    TimeCst10C13D = TimeCst10C13D + rand(size(TimeCst10C13D))./1;
    Cst10C13D30Reps = Cst10C13D30Reps + (rand(size(Cst10C13D30Reps))-0.5)./4;
    p4 = scatter(TimeCst10C13D,Cst10C13D30Reps','filled', 'MarkerFaceColor', cmaplines(1,:),...
        'SizeData', 1);
    p5 = plot(TimeCst10C13D(1,:),mean(Cst10C13D30Reps,2), 'Color',cmaplines(1,:),'LineWidth',2, 'LineStyle','-');

    ax.YLim = [0 max(survival)*1.05]; ax.LineWidth = 2;
    text(ax.XLim(1),ax.YLim(2)*1.1, 'E', 'FontSize', 12);

    nexttile(); hold on;
    survival = [];
    survival = sum(Deacclimation(1:5,:)>0,1);
    survival = survival./max(survival);    
    plot([0 9 9 14 14 42], [4 4 22 22 4 4], 'Color','k', 'LineWidth', 2)
    ax = gca; ax.YLim = [2 22*1.05]; ax.Box = 'off'; ax.Title.String = 'Acclimation is persistent (n=5)'; ax.FontSize = 9;
    ax.XLabel.String = 'Time (days)'; ax.YLabel.String = 'Temperature (C)';
    yyaxis right; ax = gca; ax.YLabel.String = 'Survival (%)'; ax.YColor = [0 0 0];
    plot([0 (1:1:28)+14],[survival(1) survival], 'Color',cmaplines(2,:), 'LineWidth', 2)

    TimeDeacclimation5D = (0:1:size(Deacclimation(1:5,:),2)-1)+13;
    TimeDeacclimation5D = repmat(TimeDeacclimation5D,size(Deacclimation(1:5,:),1),1);
    TimeDeacclimation5D = TimeDeacclimation5D + rand(size(TimeDeacclimation5D))./1;
    HealthDeacclimation(1:5,:) = Deacclimation(1:5,:) + (rand(size(Deacclimation(1:5,:)))-0.5)./10;
    p4 = scatter(TimeDeacclimation5D,HealthDeacclimation(1:5,:),'filled', 'MarkerFaceColor', cmaplines(1,:),...
        'SizeData', 1);
    p5 = plot(TimeDeacclimation5D(1,:),mean(Deacclimation(1:5,:),1), 'Color',cmaplines(1,:),'LineWidth',2, 'LineStyle','-');

    ax.YLim = [0 max(survival)*1.05]; ax.LineWidth = 2;
    text(ax.XLim(1),ax.YLim(2)*1.1, 'F', 'FontSize', 12);

    nexttile(); hold on;
    survival = [];
    survival = sum(Deacclimation(19:28,:)>0,1); %Naive animals are in Deacclimation(29:34,:)
    survival = survival./max(survival);
    plot([0 2 2 14 14 42], [4 4 22 22 4 4], 'Color','k', 'LineWidth', 2)
    ax = gca; ax.YLim = [2 22*1.05]; ax.Box = 'off'; ax.Title.String = 'Acclimation is reversible (n=10)'; ax.FontSize = 9;
    ax.XLabel.String = 'Time (days)'; ax.YLabel.String = 'Temperature (C)';
    yyaxis right; ax = gca; ax.YLabel.String = 'Survival (%)'; ax.YColor = [0 0 0];
    plot([0 (1:1:28)+14],[survival(1) survival], 'Color',cmaplines(2,:), 'LineWidth', 2)
    
    TimeDeacclimation19D = (0:1:size(Deacclimation(19:28,:),2)-1)+13;
    TimeDeacclimation19D = repmat(TimeDeacclimation19D,size(Deacclimation(19:28,:),1),1);
    TimeDeacclimation19D = TimeDeacclimation19D + rand(size(TimeDeacclimation19D))./1;
    HealthDeacclimation(19:28,:) = Deacclimation(19:28,:) + (rand(size(Deacclimation(19:28,:)))-0.5)./10;
    p4 = scatter(TimeDeacclimation19D,HealthDeacclimation(19:28,:),'filled', 'MarkerFaceColor', cmaplines(1,:),...
        'SizeData', 1);
    p5 = plot(TimeDeacclimation19D(1,:),mean(Deacclimation(19:28,:),1), 'Color',cmaplines(1,:),'LineWidth',2, 'LineStyle','-');

    ax.YLim = [0 max(survival)*1.05]; ax.LineWidth = 2;
    text(ax.XLim(1),ax.YLim(2)*1.1, 'G', 'FontSize', 12);

    f = gcf;
    % exportgraphics(f, [f.Name '.jpg'])
    exportgraphics(f, [f.Name '.pdf'])

%%
end