% Creates histogram for all the acclimation constant data
function createhistoconstant(data)
cmap = colormap("jet");
xpos = 0:size(data,2)-1;

ydata = mean(data(:,:),1);
hold on;
bar(xpos,ydata,0.4, 'FaceColor', 'k', 'FaceAlpha',0.2, 'EdgeColor', 'w');
% Xvalues = [1:size(data,2)].*ones(size(data,1),size(data,2));
Xvalues = xpos.*ones(size(data,1),size(data,2));
for i = 1:size(data,2)
    s = scatter(Xvalues(:,i), data(:,i), 'filled', 'MarkerFaceColor',...
        [0 0 0], 'AlphaData', 0*ones(1,length(data(:,i))), 'MarkerFaceAlpha','flat', 'SizeData', 15, 'jitter','on', 'jitterAmount',0.01);
end

% xticks(xpos);
% labelsnames = {(xpos)};
% xticklabels(fliplr(labelsnames'));

ylim([0 1.2])
yticks([0 1])
ylabel('A.I.')
xlim([-1 16])
% xticks([0:5:15])
% xlabel('Step duration (days)')
% text(4.2,5.6, ['n = ' num2str(length(data(:,1))) ' per step duration'], 'FontSize', 14);
ax = gca; ax.LineWidth = 1.5; ax.FontSize = 14; box off;
