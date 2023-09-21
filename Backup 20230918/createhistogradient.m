% Creates histogram for all the acclimation gradient data
function createhistogradient(data)
cmap = colormap("lines");
cmap2 = colormap("jet");

% xpos = [1 2 3 4 5 6];
xpos = fliplr(round([10/0.5*16/24, 10/0.625*16/24, 10/0.75*16/24, 10/1*16/24, 10/2*16/24, 1/24],3));
errorbar(xpos,mean(data(:,:)),std(data(:,:))./sqrt(size(data,1)),std(data(:,:))./sqrt(size(data,1)),...
    'Color','k','Linewidth',1,'LineStyle','none')
hold on;
bar(xpos,mean(data(:,:)),0.4, 'FaceColor', 'k', 'FaceAlpha',0.2, 'EdgeColor', 'w');
duration = fliplr(round([10/0.5*16/24, 10/0.625*16/24, 10/0.75*16/24, 10/1*16/24, 10/2*16/24, 1/24]));
% XvaluesRT = [1:6].*ones(size(data,1),size(data,2));
XvaluesRT = xpos.*ones(size(data,1),size(data,2));
for i = 1:size(data,2)
    s = scatter(XvaluesRT(:,i), data(:,i),'filled', 'MarkerFaceColor', [0 0 0],...
        'AlphaData', 0.1*ones(1,length(XvaluesRT)), 'MarkerFaceAlpha','flat', 'SizeData', 15, 'jitter','on', 'jitterAmount',0.17);
end


% xticks(xpos);
% labelsnames = {[num2str(10/0.5*16/24,3)]; [num2str(10/0.625*16/24,3)]; [num2str(10/0.75*16/24,2)];...
%     [num2str(10/1*16/24,2)]; [num2str(10/2*16/24,2)]; [num2str(round(0.5/18*16/24,2))]};
% xticklabels(fliplr(labelsnames'));

% xtickangle(30);
ylim([0 1.2])
yticks([0:0.2:1])
ylabel('A.I.')
xlim([-1 16])
xticklabels([])
% xlabel('Ramp duration (days)')
text(2.8,1.3, ['n = ' num2str(length(data(:,1))) ' per ramp duration'], 'FontSize', 14);
ax = gca; ax.LineWidth = 1.5; ax.FontSize = 14; box off;
