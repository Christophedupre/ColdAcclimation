cmaplines = colormap('lines');
figure('Name','Feeding Intervals','WindowState','maximized'); 
tiledlayout(2,5);
nexttile()
hold on;

FI4C = Pilot.FeedingInterval4C;
FI22C = Pilot.FeedingInterval22C;
b = bar(1,mean(FI4C),'FaceColor',cmaplines(1,:),'BarWidth',0.5);
b.EdgeColor = 'none';
b = bar(2,mean(FI22C),'FaceColor',cmaplines(2,:),'BarWidth',0.5);
b.EdgeColor = 'none';
scatter(1,FI4C,15,[0 0 0],'filled','MarkerFaceAlpha',1,'jitter','on','jitteramount',0.05);
scatter(2,FI22C,15,[0 0 0],'filled','MarkerFaceAlpha',1,'jitter','on','jitteramount',0.05);

ax = gca;
ax.YLabel.String = 'Feeding intervals (days)'; ax.XTick = [1 2]; ax.LineWidth = 1;
ax.FontSize = 12; ax.XTickLabel = {'4C','22C'};

disp(['mean feeding interval at 4C = ' num2str(mean(FI4C),2) '+/-' num2str(std(FI4C)./sqrt(length(FI4C)),2)]);
disp(['mean feeding interval at 22C = ' num2str(mean(FI22C),2) '+/-' num2str(std(FI22C)./sqrt(length(FI22C)),2)]);

l = line([1 2],[190 190]);
l.LineWidth = 2; l.Color = [0 0 0];
t = text(1.5,192,'***'); t.HorizontalAlignment = "center"; t.FontSize = round(t.FontSize*1.5);

f = gcf;
exportgraphics(f, [f.Name '.jpg'])