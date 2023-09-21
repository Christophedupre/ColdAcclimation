function plotConstantSingleTraces(CstData3D)
figCstSingleTraces = figure('Name', 'Constant Temp Single Traces', 'renderer', 'painters'); 
figCstSingleTraces.WindowState = 'maximized';
clf(figCstSingleTraces)
ncols = 15; nrows = 16; %Number of columns and rows in this tiledlayout
tiledlayout(nrows,ncols)
Time = linspace(0,length(CstData3D(:,1,1))-1,length(CstData3D(:,1,1)));

for j = 1:size(CstData3D,3)
    for i=1:size(CstData3D,2)
        nexttile()
        plot(Time, CstData3D(:,i,j), 'LineWidth', 2)
        box off;
        xlim([-0.1 28]); xticks([]); xticklabels([]); ylim([0 1]); yticks([]);
        if(j == 1) text(15,1.4,[num2str(i-1)]); end %Training duration
        if(i == 1) text(-15,0.6,[num2str(j*2+2) 'C']); end %Training temp
    end    
end

t1 = annotation('textbox', 'EdgeColor', 'none', 'String',...
    'Training Temperature (C)', 'Position', [0.05 0.35 0.1 0.1],...
    'FontSize',15,'Rotation', 90);

t2 = annotation('textbox', 'EdgeColor', 'none', 'String',...
    'Training Duration (Days)', 'Position', [0.42 0.9 0.1 0.1],...
    'FontSize',15,'Rotation', 0);

f = gcf;
exportgraphics(f, [f.Name '.pdf'])
clear f i j ncols nrows t1 t2 Time

end