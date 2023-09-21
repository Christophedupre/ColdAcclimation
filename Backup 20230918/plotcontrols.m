function plotcontrols(PosCtrl, NoPCR, JustPCR, Ramp, Pulses0_2, Pulses1_1, CstData3D)
figControls = figure('Name', 'Controls', 'renderer', 'painters'); 
figControls.WindowState = 'maximized';
cmaplines = colormap("lines");
cmapjet = colormap("jet");

figure(figControls)
clf(figControls)
ncols = 5; nrows = 4;
tiledlayout(nrows, ncols)
Time = linspace(0,length(CstData3D(:,1,1))-1,length(CstData3D(:,1,1)));
p = [];

%Positive controls
nexttile()
plot([linspace(0,14,15) 14], [linspace(4,4,15) 4],'color',[cmaplines(1,:) 0.8], 'LineWidth',2); hold on;
plot(linspace(14,28,15), linspace(4,4,15),'color',[cmaplines(2,:) 0.8], 'LineWidth',2);
xlim([0 28]); xticks([14 28]); xticklabels({'0';'28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 22])
yticks([0 4]); yticklabels({'0';'4'})
box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 12;

nexttile(); hold on;
for j = 1:size(PosCtrl,2) % Over all training durations
    p1 = plot(Time, PosCtrl(:,j), 'Color',[cmapjet(round(j/size(PosCtrl,2)*256),:) 1], 'LineWidth', 2);
    yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 12;
    xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1]); yticks([0:0.5:1]);
    p = [p p1];    
end
xlabel('Time (days)'); ylabel('A.I.'); text(1, 1.1, 'Positive Controls (=acclimated)');

%Negative controls
nexttile()
plot([linspace(0,14,15) 14], [linspace(22,22,15) 4],'color',[cmaplines(1,:) 0.8], 'LineWidth',2); hold on;
plot(linspace(14,28,15), linspace(4,4,15),'color',[cmaplines(2,:) 0.8], 'LineWidth',2);
xlim([0 28]); xticks([14 28]); xticklabels({'0';'28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 22])
yticks([0 22]); yticklabels({'0';'22'})
box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 12;

nexttile(); hold on;
for j = 1:size(PosCtrl,2) % Over all training durations
    p1 = plot(Time, NoPCR(:,j), 'Color',[cmapjet(round(j/size(NoPCR,2)*256),:) 1], 'LineWidth', 2);
    yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 12;
    xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1]); yticks([0:0.5:1]);    
end
text(1, 1.1, 'Negative Controls (NoPCR)');

nexttile(); hold on;
for j = 1:size(PosCtrl,2) % Over all training durations
    p1 = plot(Time, JustPCR(:,j), 'Color',[cmapjet(round(j/size(JustPCR,2)*256),:) 1], 'LineWidth', 2);
    yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 12;
    xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1]); yticks([0:0.5:1]);    
end
text(1, 1.1, 'Negative Controls (JustPCR)');

%Ramp
nexttile()
plot([linspace(0,7,15)], [linspace(22,13,15)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2); hold on;
plot([linspace(7,14,8)], [linspace(13,4,8)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2, 'LineStyle',':');
% plot([linspace(7,14,8)], [linspace(4,4,8)],'color',[cmap(2,:) 0.8], 'LineWidth',2);
plot([linspace(7,7,8)], [linspace(13,4,8)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2);
plot(linspace(7,21,15), linspace(4,4,15),'color',[cmaplines(2,:) 0.8], 'LineWidth',2);
plot(linspace(21,28,8), linspace(4,4,8),'color',[cmaplines(2,:) 0.8], 'LineWidth',2, 'LineStyle',':');
xlim([0 28]); xticks([7 21]); xticklabels({'0';'28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 22])
yticks([4 22]); yticklabels({'4';'20'})
box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 12;

nexttile(); hold on;
for j = 1:size(PosCtrl,2) % Over all training durations
    p1 = plot(Time, Ramp(:,j), 'Color',[cmapjet(round(j/size(NoPCR,2)*256),:) 1], 'LineWidth', 2);
    yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 12;
    xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1]); yticks([0:0.5:1]);    
end
text(1, 1.1, 'Ramp');

%Pulses
nexttile(11)
plot([linspace(0,3.5,3) 3.5 linspace(3.5,7,4) linspace(7,10.5,3) 10.5 linspace(10.5,14,4)], [linspace(22,22,3) 4 linspace(4,4,4) linspace(22,22,3) 4 linspace(4,4,4)],'color',[cmaplines(1,:) 0.8], 'LineWidth',2); hold on;
plot(linspace(14,28,15), linspace(4,4,15),'color',[cmaplines(2,:) 0.8], 'LineWidth',2);
xlim([0 28]); xticks([14 28]); xticklabels({'0';'28'}); xlabel('Time (days)'); ylabel(['T (' char(176) 'C)']); ylim([0 22])
yticks([4 22]); yticklabels({'4';'20'})
box off; ax = gca; ax.LineWidth = 1.5; ax.FontSize = 12;

nexttile(12); hold on;
for j = 1:size(PosCtrl,2) % Over all training durations
    p1 = plot(Time, Pulses0_2(:,j), 'Color',[cmapjet(round(j/size(NoPCR,2)*256),:) 1], 'LineWidth', 2);
    yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 12;
    xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1]); yticks([0:0.5:1]);    
end
text(1, 1.1, 'T_{Pulse} = 0');

nexttile(13); hold on;
for j = 1:size(PosCtrl,2) % Over all training durations
    p1 = plot(Time, Pulses1_1(:,j), 'Color',[cmapjet(round(j/size(NoPCR,2)*256),:) 1], 'LineWidth', 2);
    yticks([0:5]); box off; ax = gca; ax.LineWidth = 1; ax.FontSize = 12;
    xlim([-0.1 28]); xticks([0:5:28]); ylim([0 1]); yticks([0:0.5:1]);    
end
text(1, 1.1, 'T_{Pulse} = 2');

legendstring = num2str((1:length(p))');
lgdCtrls = legend([p],legendstring,'Interpreter', 'latex', 'FontSize', 9,... 
    'EdgeColor', 'none','NumColumns',5, 'box', 'on');
lgdCtrls.Title.String = 'Training duration (Days)';
lgdCtrls.Position = [0.78 0.11 0.1229 0.1100];
set(lgdCtrls.BoxFace, 'ColorType','truecoloralpha', 'ColorData',uint8(255*[1;1;1;.8]));

nexttile() %Just to remove a bug caused by the legend display
h = heatmap(squeeze(mean(CstData3D,1))'); h.ColorbarVisible = 'off'; 
h.Visible = 'off';

f = gcf;
exportgraphics(f, [f.Name '.pdf'])
clear p ax i j p1 legendstring h f lgdCtrls ncols nrows;

end