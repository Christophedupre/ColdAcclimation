function plotLakeData(LakeMichigan, ModelData)
figlake = figure('Name','Lake Temperature','WindowState','maximized','Renderer','painters');
cmaplines = colormap('lines');
cmapjet = colormap('jet');
cmapturbo = colormap('turbo');
tiledlayout(3,4)
load("cmaplist.mat");

%%% DISPLAY LAKE TEMPERATURE AND FOURRIER TRANSFORM OVER MULTIPLE YEARS %%%
nexttile(); hold on;
cmaplines = colormap('lines');
signal = LakeMichigan.temperatureSurfaceAllYears;
Fs = 1/(60*60);     % Sampling frequency. It is 1/hour, therefore 1/60*60sec   
T = 1/Fs;           % Sampling period (1 hour)       
L = length(signal); % Length of signal (#samples)
t = (0:L-1)/24;     % Time vector (days)

p1 = plot(t,signal, 'Color',cmaplines(1,:), 'LineWidth',1); 
ax = gca; ax.YLim = [0 26]; ax.LineWidth = 1; ax.FontSize = 10; ax.XLim = [0 max(t)];
ax.XLabel.String = 'time (days)'; ax.YLabel.String = 'Temperature (C)';
ax.XTick = floor(linspace(0,max(t),15));
ax.XTickLabel = num2str((0:1:14)'); %Indicate time in month of the year instead of day
ax.XLabel.String = 'time (years)'; ax.Box = 'off';
% ax.Title.String = 'Temp at depth = 9m';
text(ax.XLim(1),ax.YLim(2)*1.15, 'A.', 'FontSize', 14);
% l = legend('Depth (m) = 9'); l.Box = 'off';
% l.Title.String = ('Depth (m)');

text(ax.XLim(2)*0.7,ax.YLim(2)*0.9,{'Depth = 9m'}, 'BackgroundColor', [1 1 1 0.75], 'FontSize',9);

nexttile()
Y = fft(signal);
P2 = abs(Y); % 2-sided spectrum
P1 = P2(1:floor(L/2+1)); % single-sided spectrum P1
P1(2:end-1) = 2*P1(2:end-1);
% P1 = P1./max(P1); % Normalize
f2 = Fs*(0:L/2)/L*60*60*24; %Frequency, in 1/days
plot(f2,P1, 'Color',cmaplines(1,:), 'LineWidth',1); 
ax = gca; ax.XLabel.String = 'frequency (1/day)'; ax.YLabel.String = 'Amplitude (a.u.)';
ax.YLim = [0 4000]; ax.XLim = [0 2];
ax.YTick = [];
ax.LineWidth = 1; ax.FontSize = 10; ax.Box = 'off';
text(ax.XLim(1),ax.YLim(2)*1.15, 'B', 'FontSize', 14);

% ax2 = axes('Position',[0.65 0.65 0.28 0.28]);
nexttile()
f2 = Fs*(0:L/2)/L*60*60*24*365; %Frequency, in 1/year
plot(f2,P1, 'Color',cmaplines(1,:), 'LineWidth',1); 
ax = gca; ax.XLabel.String = 'frequency (1/year)'; ax.YLabel.String = 'Amplitude (a.u.)';
ax.XLim = [0 2]; ax.YTick = [];
ax.LineWidth = 1; ax.FontSize = 10; ax.Box = 'off';
text(ax.XLim(1),ax.YLim(2)*1.24, 'C', 'FontSize', 14);
% text(ax.XLim(1),(ax.YLim(2)-ax.YLim(1))*1.15+ax.YLim(1), 'C', 'FontSize', 14)

%%% DISPLAY TEMPERATURE AT DIFFERENT DEPTHS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(); hold on;
cmapjet = colormap('jet');
p = [];
t = 0:1:length(LakeMichigan.temperature(1,:))-1;
for i = 1:3:length(LakeMichigan.z)
    p1 = plot(t,LakeMichigan.temperature(i,:), 'Color',[cmapjet(floor(i*256/length(LakeMichigan.z)),:) 0.5], 'LineWidth',1); 
%     p1 = plot(t,LakeMichigan.temperature(i,:), 'Color',[cmaplines(i,:) 0.5], 'LineWidth',1); 
    p = [p p1];
end
ax = gca; ax.YLim = [0 24]; ax.LineWidth = 1; ax.FontSize = 10;
ax.XLabel.String = 'time (days)'; ax.YLabel.String = 'Temperature (C)';
ax.XTick = floor(linspace(0,max(t),12));
ax.XTickLabel = {'S','O','N','D','J','F','M','A','M','J','J','A'}; %Indicate time in month of the year instead of day
ax.XLabel.String = 'time (months)';
ax.XLim = [0 length(t)];
text(ax.XLim(1),ax.YLim(2)*1.15, 'D', 'FontSize', 14);

l = legend(p,num2str(LakeMichigan.z(1:3:end)));
l.Box = 'off';
l.Position(1) = ax.Position(1)+ax.Position(3)./5;
l.Position(2) = ax.Position(2)+ax.Position(4)./3;
l.Position(4) = l.Position(4)./2; l.Position(3) = l.Position(3)./3;
l.Location = 'north';
l.Title.String = ('Depth (m)');

%%% DISPLAY SOME STATISTICS AT DIFFERENT DEPTHS %%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile()
hold on;
% ax = gca; ax.Visible = 'off';
% Compute min, max and average temperatures
mintemperature = zeros(1,length(LakeMichigan.z));
maxtemperature = mintemperature; meantemperature = mintemperature;
stdevtemperature = mintemperature;
for i = 1:length(LakeMichigan.z)
    mintemperature(i) = min(LakeMichigan.temperature(i,:));
    maxtemperature(i) = max(LakeMichigan.temperature(i,:));
    meantemperature(i) = mean(LakeMichigan.temperature(i,:));
    stdevtemperature(i) = std(LakeMichigan.temperature(i,:));
end
% p1 = plot(LakeMichigan.z, mintemperature, 'Color',cmaplines(1,:), 'LineWidth',1);
% p2 = plot(LakeMichigan.z, maxtemperature, 'Color',cmaplines(2,:), 'LineWidth',1);
p3 = plot(LakeMichigan.z, meantemperature, 'Color',cmaplines(1,:), 'LineWidth',1);
% p4 = plot(LakeMichigan.z, stdevtemperature, 'Color',cmaplines(4,:), 'LineWidth',1);
x = [LakeMichigan.z' fliplr(LakeMichigan.z')];
y = [meantemperature+stdevtemperature fliplr(meantemperature-stdevtemperature)];
f = fill(x,y,cmaplines(1,:),'FaceColor',cmaplines(1,:),'FaceAlpha',0.2,'EdgeColor','none');

l = legend([p3 f], {'Mean','Stdev'});
l.Box = 'off';
ax = gca; ax.XLabel.String = 'Depth (meters)'; ax.YLabel.String = 'Temperature (C)';
ax.LineWidth = 1; ax.FontSize = 10; ax.Box = 'off';
text(ax.XLim(1),ax.YLim(2)*1.15, 'E', 'FontSize', 14);

% %%% DISPLAY LAKE TEMPERATURE OVER 1 YEAR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nexttile()
% signal = LakeMichigan.temperature(1,:);
% Fs = 1/(60*60);     % Sampling frequency. It is 1/hour, therefore 1/60*60sec   
% T = 1/Fs;           % Sampling period (1 hour)       
% L = length(signal); % Length of signal (#samples)
% t = (0:L-1)/24;     % Time vector (days)
% % s = [9:12 1:8];      % Time vector (seasons)
% 
% p1 = plot(t,signal, 'Color',cmaplines(1,:), 'LineWidth',1); 
% ax = gca; ax.YLim = [0 24]; ax.LineWidth = 1; ax.FontSize = 10; ax.XLim = [0 max(t)];
% ax.XLabel.String = 'time (days)'; ax.YLabel.String = 'Temperature (C)';
% ax.XTick = floor(linspace(0,max(t),12));
% ax.XTickLabel = {'S','O','N','D','J','F','M','A','M','J','J','A'}; %Indicate time in month of the year instead of day
% ax.XLabel.String = 'time (months)'; ax.Box = 'off';
% text(ax.XLim(1),ax.YLim(2)*1.15, 'F', 'FontSize', 14);
% % rectangle('Position',[180 3.5 10 1], 'LineWidth', 2) %Indicate where zoom in was made

%%% ZOOM IN A DATA PERIOD ILLUSTRATING DAILY FLUCTUATIONS %%%%%%%%%%%%%%%%%
% nexttile()
% signal = LakeMichigan.temperature(1,:);
% plot(t,signal, 'Color',cmaplines(1,:), 'LineWidth',1); 
% ax = gca; ax.YLim = [3.5 4.5]; ax.XLim = [180 200];
% ax.LineWidth = 1; ax.FontSize = 10; ax.Box = 'off';
% ax.XLabel.String = 'time (days)'; ax.YLabel.String = 'Temperature (C)';
% signaltruncated = signal(200*24:240*24);
% 
% % Plot FFT of signal
% Fs = 1/(60*60); % Sampling frequency. It is 1/hour, therefore 1/60*60sec   
% T = 1/Fs; % Sampling period (1 hour)       
% L = length(signaltruncated); % Length of signal (#samples)
% t = (0:L-1)/24;        % Time vector (days)
% 
% nexttile()
% Y = fft(signaltruncated);
% P2 = abs(Y); % 2-sided spectrum
% P1 = P2(1:floor(L/2+1)); % single-sided spectrum P1
% P1(2:end-1) = 2*P1(2:end-1);
% f2 = Fs*(0:L/2)/L*60*60*24; %Frequency, in 1/days
% plot(f2,P1, 'Color',cmaplines(1,:), 'LineWidth',1); 
% ax = gca; ax.XLabel.String = 'frequency (1/day)'; ax.YLabel.String = 'Amplitude';
% ax.LineWidth = 1; ax.FontSize = 10; ax.Box = 'off';
% ax.YLim = [0 400]; ax.XLim = [0 2];

%%% MEASUREMENT OF TEMPERATURE GRADIENTS OVER VARIOUS DURATIONS %%%%%%%%%%%
signal = LakeMichigan.temperature(1,:);
nexttile(); hold on;
DeltaTs = [1 12 24 24*10 24*30 24*90 24*180];
DeltaTLabels = {'1h', '12h', '1D', '10D', '1M', '3M', '6M'};
RampGradients = -[0.5 0.625 0.75 1 2]./10;
maxdiff = 0; mindiff = 0;
% Plot Delta\theta
for i = 1:length(DeltaTs)
    %shift signal
    deltaT = DeltaTs(i);
    sample = signal;
    sample1 = [ones(1,deltaT)*sample(1) sample];
    sample2 = [sample ones(1,deltaT)*sample(end)];
    diff = (sample2 - sample1);
    n = floor(length(diff)*1/10);
    ri = ceil(rand(1,n)*length(diff)); %Random indices
    s = scatter(linspace(i,i,n),diff(1,ri),'Marker', '.','XJitter','rand',...
        'XJitterWidth',0.5,'MarkerEdgeColor',cmaplines(i,:));
    
    % Add "violin" plot
    [a,edges] = histcounts(diff);
    pdf = a./max(a)/4;
    x = [pdf fliplr(-pdf)]+i; 
    y = [linspace(min(diff),max(diff),length(a)) fliplr(linspace(min(diff),max(diff),length(a)))];
    fill(x,y,cmaplines(1,:), 'FaceColor',cmaplines(i,:),'FaceAlpha',0.2,'EdgeColor',[0 0 0]);

    % Look at max and min difference
    if (max(diff)>maxdiff) maxdiff = max(diff); end
    if (min(diff)<mindiff) mindiff = min(diff); end
end


ax = gca; ax.XLim = [0 length(DeltaTs)+1];
ax.XTick = 1:length(DeltaTs); ax.XTickLabel = DeltaTLabels;
ax.YLabel.String = 'Temperature Difference (C)'; ax.XLabel.String = 'Time interval';
ax.LineWidth = 1; ax.FontSize = 10; ax.Box = 'off';
text(ax.XLim(1),(ax.YLim(2)-ax.YLim(1))*1.15+ax.YLim(1), 'G', 'FontSize', 14)

% % Plot Normalized Delta\theta (in C/hr)
% nexttile(); hold on;
% maxdiff = 0; mindiff = 0;
% for i = 1:length(DeltaTs)
%     %shift signal
%     deltaT = DeltaTs(i);
%     sample = signal;
%     sample1 = [ones(1,deltaT)*sample(1) sample];
%     sample2 = [sample ones(1,deltaT)*sample(end)];
%     diff = (sample2 - sample1)/deltaT;
%     n = floor(length(diff)*1/1);
%     ri = ceil(rand(1,n)*length(diff)); %Random indices
%     s = scatter(linspace(i,i,n),diff(1,ri),'Marker', '.','XJitter','rand',...
%         'XJitterWidth',0.5,'MarkerEdgeColor',cmaplines(i,:));
% 
%     % Add "violin" plot
%     [a,edges] = histcounts(diff);
%     pdf = a./max(a)/4;
%     x = [pdf fliplr(-pdf)]+i; 
%     y = [linspace(min(diff),max(diff),length(a)) fliplr(linspace(min(diff),max(diff),length(a)))];
%     fill(x,y,cmaplines(1,:), 'FaceColor',cmaplines(i,:),'FaceAlpha',0.2,'EdgeColor', [0 0 0])
% 
%      % Look at max and min difference
%     if (max(diff)>maxdiff) maxdiff = max(diff); end
%     if (min(diff)<mindiff) mindiff = min(diff); end
% end
% disp(['max slope = ' num2str(maxdiff) 'C/hr'])
% disp(['min slope = ' num2str(mindiff) 'C/hr'])
% disp(['max temp = ' num2str(max(signal)) 'C'])
% disp(['min temp = ' num2str(min(signal)) 'C'])
% % maxdiff
% % mindiff
% 
% for i =1:length(RampGradients)
%     % Add lines showing the gradients used in the ramp experiments
%     p = plot([1 7],ones(1,2).*RampGradients(i),'Color',[0 0 0],'LineStyle','--');
% end
% l = legend([p],'Ramp experiments');
% l.Box = 'off';
% 
% ax = gca; ax.XLim = [0 length(DeltaTs)+1];
% ax.XTick = 1:length(DeltaTs); ax.XTickLabel = DeltaTLabels;
% ax.YLabel.String = 'd\theta/dt (C/Hr)'; ax.XLabel.String = 'Time interval';
% ax.YLim = [-0.25 0.15];
% ax.LineWidth = 1; ax.FontSize = 10; ax.Box = 'off';
% text(ax.XLim(1),(ax.YLim(2)-ax.YLim(1))*1.15+ax.YLim(1), 'H', 'FontSize', 14)

%%% TEST WHETHER MODEL PREDICTS SURVIVAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile(); hold on;
cmaplines = colormap('lines');
solpts = zeros(2,length(LakeMichigan.temperature));

% Generate model data
cvar = ModelData.csol.c; % No var in coefficients
temperature1 = linspace(0,365,length(LakeMichigan.temperature(1,:)));
temperature2 = [LakeMichigan.temperature(1,:)]./22; %Duration = 13 days, training temperature = 10C
temperature = [temperature1; temperature2];
sol = ode23(@(t,y) Acclimfun(t,y,cvar,temperature), temperature1, ModelData.y0);
solptstemp = deval(temperature1,sol);
solpts(1,:) = solptstemp(1,:);
solpts(2,:) = solptstemp(2,:);

% % Generate model data - Lake Superior
% cvar = ModelData.csol.c; % No var in coefficients
% load('AcclimationLakeSuperiorData.mat');
% solpts = zeros(2,length(LakeSuperior.temperature));
% temperature1 = linspace(0,365,length(LakeSuperior.temperature(1,:)));
% temperature2 = [LakeSuperior.temperature(1,:)]./22; %Duration = 13 days, training temperature = 10C
% temperature = [temperature1; temperature2];
% sol = ode23(@(t,y) Acclimfun(t,y,cvar,temperature), temperature1, ModelData.y0);
% solptstemp = deval(temperature1,sol);
% solpts(1,:) = solptstemp(1,:);
% solpts(2,:) = solptstemp(2,:);

%Plot model data
p1 = plot(temperature1, solpts(2,:), 'Color', cmaplines(2,:), 'LineWidth', 2);
p2 = plot(temperature1, solpts(1,:), 'Color', cmaplines(1,:), 'LineWidth', 2);
p3 = plot(temperature1, temperature2, 'Color', [0 0 0 0.5], 'LineWidth', 1);
% the data used for zoom in starts at time point 5285 and ends at 5500;

% p4 = plot()

ax = gca; ax.LineWidth = 1; ax.FontSize = 10; ax.XLim = [0 max(temperature1)];
ax.XLabel.String = 'time (days)'; ax.YLabel.String = 'Normalized value';
ax.XTick = floor(linspace(0,max(temperature1),12));
ax.XTickLabel = {'S','O','N','D','J','F','M','A','M','J','J','A'}; %Indicate time in month of the year instead of day
ax.XLabel.String = 'time (months)'; ax.Box = 'off';
l = legend([p1 p2 p3], {'Health','Acclimation', 'Temperature'});
l.Box = 'off'; l.Location = 'north';
text(ax.XLim(1),ax.YLim(2)*1.15, 'I', 'FontSize', 14);

%Plot close-up
nexttile(); hold on;
p1 = plot(temperature1, solpts(2,:), 'Color', cmaplines(2,:), 'LineWidth', 2);
p2 = plot(temperature1, solpts(1,:), 'Color', cmaplines(1,:), 'LineWidth', 2);
p3 = plot(temperature1, temperature2, 'Color', [0 0 0 0.5], 'LineWidth', 1);

ax = gca; ax.LineWidth = 1; ax.FontSize = 10; ax.XLim = [0 max(temperature1)];
ax.XLabel.String = 'time (days)'; ax.YLabel.String = 'Normalized value';
ax.XTick = floor(linspace(0,max(temperature1),12));
ax.XTickLabel = {'S','O','N','D','J','F','M','A','M','J','J','A'}; %Indicate time in month of the year instead of day
ax.XLabel.String = 'time (months)'; ax.Box = 'off';
l = legend([p1 p2 p3], {'Health','Acclimation', 'Temperature'});
l.Box = 'off'; l.Location = 'north';
ax.XLim = [214 245];
ax.YLim = [0.15 0.25];
text(ax.XLim(1),(ax.YLim(2)-ax.YLim(1))*1.15+ax.YLim(1), 'J', 'FontSize', 14);

%%% HEALTH 28 DAYS AFTER END OF RAMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ramps that have high slope and amplitude give less acclimation and 
% are therefore more lethal if they are followed by 50 days of cold
nexttile()
amax = 20;
% dmax = linspace(24,24*20,24*20-23);
dmax = linspace(1,61,61);
% InitialTemp = amax+1;
InitialTemp = 22;
health = zeros(amax,length(dmax));
label = [];
Delay = 10;
for a=1:amax
    for di=1:length(dmax)
        d = dmax(di);
        if (a == 1) 
            label = [label {num2str(d)}];
        end
        % tspan = 0:0.1:d/24;
        tspan = linspace(0,d+Delay,d+Delay+1);
        amplitude = a/22;
        y0 = [InitialTemp/22 1];
        c = ModelData.csol.c;
        temperature1 = tspan;
        temperature2 = [linspace(InitialTemp/22,InitialTemp/22-amplitude,d) linspace(InitialTemp/22-amplitude,InitialTemp/22-amplitude,Delay+1)];
        temperature = [temperature1;temperature2];
        sol = ode23(@(t,y) Acclimfun(t,y,c,temperature) , tspan, y0);
        solpts = deval(sol,tspan);
        health(a,di) = [(solpts(2,end))];
    end
end
% plot(health);
h = heatmap(health);
h.YLabel = 'Ramp Amplitude (C)'; h.XLabel = 'Ramp Duration (days)';
h.XDisplayLabels(:) = {''};
h.XDisplayLabels(5:5:end) = label(5:5:end);
temp = h.YDisplayLabels(:);
h.YDisplayLabels(:) = {''};
h.YDisplayLabels(5:5:end) = temp(5:5:end);
h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';
h.Title = ['K. Health ' num2str(Delay) ' days after end of ramp'];

h.Colormap = cmaplist.cmapRedBlue1;

%%% RAMP AMPLITUDE VS RAMP DURATION IN ENVIRONMENTAL DATA %%%%%%%%%%%%%%%%%

nexttile()
signal = LakeMichigan.temperature(1,:);
Fs = 1/(60*60);     % Sampling frequency. It is 1/hour, therefore 1/60*60sec   
T = 1/Fs;           % Sampling period (1 hour)       
L = length(signal); % Length of signal (#samples)
t = (0:L-1)/24;     % Time vector (days)
% s = [9:12 1:8];      % Time vector (seasons)

p1 = plot(t,signal, 'Color',cmaplines(1,:), 'LineWidth',1);
deltamax = 24*60; 

HData = [];
% plot(delta(1,:));
edges = 0:1:20;
deltaT = 1:24*1:deltamax;
for i = deltaT
    H = histcounts(signal(1:end-i)-signal(i+1:end),edges,'Normalization','probability');
    HData = [HData; H];
end

YLabels = [];
for i = 1:length(edges) YLabels = [YLabels {num2str(i)}]; end

XLabels = [];
for i = deltaT XLabels = [XLabels {num2str(ceil(i./24))}]; end


h = heatmap(HData');
h.YDisplayLabels(:) = {''};
h.XDisplayLabels(:) = {''};
h.YDisplayLabels(5:5:end) = YLabels(5:5:end-1);
h.XDisplayLabels(5:5:end) = XLabels(5:5:end);
h.YLabel = 'Ramp Amplitude (C)'; h.XLabel = 'Ramp Duration (days)';
% h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';
h.Title = 'L. Ramp ampl. and dur. probab.';

% ax = gca; ax.YLim = [0 24]; ax.LineWidth = 1; ax.FontSize = 10; ax.XLim = [0 max(t)];
% ax.XLabel.String = 'time (days)'; ax.YLabel.String = 'Temperature (C)';
% ax.XTick = floor(linspace(0,max(t),12));
% ax.XTickLabel = {'S','O','N','D','J','F','M','A','M','J','J','A'}; %Indicate time in month of the year instead of day
% ax.XLabel.String = 'time (months)'; ax.Box = 'off';

% colormap(cmaplist.cmapRedBlue1);
% colormap(cmaplist.cmapDarkSky);
h.Colormap = cmapturbo;

nexttile()
ax = gca; p = plot([0], [1], 'Color',cmapturbo(1,:));
ax.Colormap = cmapturbo;
c = colorbar('west');
% c.Location = "east";
c.Ticks = [0 1];
c.Label.String = 'Probability';
ax.Visible = 'off';



f = gcf;
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(f, [f.Name '.pdf'])

%%% OVERNIGHT VS SEASONAL DROPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name', 'Simulations', 'WindowState', 'maximized')
tiledlayout(3,4)
%%% OVERNIGHT VS SEASONAL RAMPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From overnight (i.e. ~10-20 hours) to seasonal (i.e. >50 days) ramps, 
% there is no ramp which by itself is lethal. A constant temperature 
% decrease alone is not harmful. However, once the temperature has
% decreased, staying at a low temperature can be harmful.

nexttile()
amax = 20;
% dmax = linspace(24,24*20,24*20-23);
dmax = linspace(1,61,61);
% InitialTemp = amax+1;
InitialTemp = 22;
health = zeros(amax,length(dmax));
acclim = zeros(amax,length(dmax));
label = [];
for a=1:amax
    for di=1:length(dmax)
        d = dmax(di);
        if (a == 1) 
            label = [label {num2str(d)}];
        end
        % tspan = 0:0.1:d/24;
        tspan = linspace(0,d,100);
        amplitude = a/22;
        y0 = [InitialTemp/22 1];
        c = ModelData.csol.c;
        temperature1 = tspan;
        temperature2 = [linspace(InitialTemp/22,InitialTemp/22-amplitude,length(tspan))];
        temperature = [temperature1;temperature2];
        sol = ode23(@(t,y) Acclimfun(t,y,c,temperature) , tspan, y0);
        solpts = deval(sol,tspan);
        health(a,di) = [(solpts(2,end))];
        acclim(a,di) = [(solpts(1,end))];
    end
end
% plot(health);
h = heatmap(health);
h.YLabel = 'Ramp Amplitude (C)'; h.XLabel = 'Ramp Duration (days)';
h.XDisplayLabels(:) = {''};
h.XDisplayLabels(5:5:end) = label(5:5:end);
temp = h.YDisplayLabels(:);
h.YDisplayLabels(:) = {''};
h.YDisplayLabels(5:5:end) = temp(5:5:end);
h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';
h.Title = 'Health at end of ramp';

nexttile()
h = heatmap(acclim);
h.YLabel = 'Ramp Amplitude (C)'; h.XLabel = 'Ramp Duration (days)';
h.XDisplayLabels(:) = {''};
h.XDisplayLabels(5:5:end) = label(5:5:end);
temp = h.YDisplayLabels(:);
h.YDisplayLabels(:) = {''};
h.YDisplayLabels(5:5:end) = temp(5:5:end);
h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';
h.Title = '\alpha at end of ramp';

%%% HEALTH 28 DAYS AFTER END OF RAMP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ramps that have high slope and amplitude give less acclimation and 
% are therefore more lethal if they are followed by 50 days of cold
nexttile()
amax = 20;
% dmax = linspace(24,24*20,24*20-23);
dmax = linspace(1,61,61);
% InitialTemp = amax+1;
InitialTemp = 22;
health = zeros(amax,length(dmax));
label = [];
for a=1:amax
    for di=1:length(dmax)
        d = dmax(di);
        if (a == 1) 
            label = [label {num2str(d)}];
        end
        % tspan = 0:0.1:d/24;
        tspan = linspace(0,d+Delay,d+Delay+1);
        amplitude = a/22;
        y0 = [InitialTemp/22 1];
        c = ModelData.csol.c;
        temperature1 = tspan;
        temperature2 = [linspace(InitialTemp/22,InitialTemp/22-amplitude,d) linspace(InitialTemp/22-amplitude,InitialTemp/22-amplitude,Delay+1)];
        temperature = [temperature1;temperature2];
        sol = ode23(@(t,y) Acclimfun(t,y,c,temperature) , tspan, y0);
        solpts = deval(sol,tspan);
        health(a,di) = [(solpts(2,end))];
    end
end
% plot(health);
h = heatmap(health);
h.YLabel = 'Ramp Amplitude (C)'; h.XLabel = 'Ramp Duration (days)';
h.XDisplayLabels(:) = {''};
h.XDisplayLabels(5:5:end) = label(5:5:end);
temp = h.YDisplayLabels(:);
h.YDisplayLabels(:) = {''};
h.YDisplayLabels(5:5:end) = temp(5:5:end);
h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';
h.Title = 'Health 28 days after end of ramp';

%%% RAMP AMPLITUDE VS RAMP DURATION IN ENVIRONMENTAL DATA %%%%%%%%%%%%%%%%%

nexttile()
signal = LakeMichigan.temperature(1,:);
Fs = 1/(60*60);     % Sampling frequency. It is 1/hour, therefore 1/60*60sec   
T = 1/Fs;           % Sampling period (1 hour)       
L = length(signal); % Length of signal (#samples)
t = (0:L-1)/24;     % Time vector (days)
% s = [9:12 1:8];      % Time vector (seasons)

p1 = plot(t,signal, 'Color',cmaplines(1,:), 'LineWidth',1);
deltamax = 24*60; 

HData = [];
% plot(delta(1,:));
edges = 0:1:20;
deltaT = 1:24*1:deltamax;
for i = deltaT
    H = histcounts(signal(1:end-i)-signal(i+1:end),edges,'Normalization','probability');
    HData = [HData; H];
end

YLabels = [];
for i = 1:length(edges) YLabels = [YLabels {num2str(i)}]; end

XLabels = [];
for i = deltaT XLabels = [XLabels {num2str(ceil(i./24))}]; end


h = heatmap(HData');
h.YDisplayLabels(:) = {''};
h.XDisplayLabels(:) = {''};
h.YDisplayLabels(5:5:end) = YLabels(5:5:end-1);
h.XDisplayLabels(5:5:end) = XLabels(5:5:end);
h.YLabel = 'Ramp Amplitude (C)'; h.XLabel = 'Ramp Duration (days)';
% h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';
h.Title = 'Ramp ampl. and dur. probab.';

% ax = gca; ax.YLim = [0 24]; ax.LineWidth = 1; ax.FontSize = 10; ax.XLim = [0 max(t)];
% ax.XLabel.String = 'time (days)'; ax.YLabel.String = 'Temperature (C)';
% ax.XTick = floor(linspace(0,max(t),12));
% ax.XTickLabel = {'S','O','N','D','J','F','M','A','M','J','J','A'}; %Indicate time in month of the year instead of day
% ax.XLabel.String = 'time (months)'; ax.Box = 'off';

%%% VAR RAMPS FOLLOWED BY 4C %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile();
tspan = 0:1:60;
y0 = [22/22 1];
c = ModelData.csol.c;
numslopes = 40;
health = zeros(numslopes,length(tspan));
acclim = zeros(numslopes,length(tspan));

for i = 1:numslopes
    temperature1 = tspan;
    temperature2 = [linspace(22/22,4/22,i) linspace(4/22,4/22,length(tspan)-i)];
    temperature = [temperature1;temperature2];
    sol = ode23(@(t,y) Acclimfun(t,y,c,temperature) , tspan, y0);
    solpts = deval(sol,tspan);
    health(i,:) = solpts(2,:);
    acclim(i,:) = solpts(1,:);
end

Xlabels = [];
for i=1:length(tspan) Xlabels = [Xlabels {(num2str(tspan(i)+1))}]; end
Ylabels = [];
for i=1:numslopes Ylabels = [Ylabels {(num2str(i))}]; end

h = heatmap(health);
h.YLabel = '18C Ramp Duration (Days)'; h.XLabel = 'Time (Days)';
h.XDisplayLabels(:) = {''};
h.XDisplayLabels(5:5:end) = Xlabels(5:5:end);
temp = h.YDisplayLabels(:);
h.YDisplayLabels(:) = {''};
h.YDisplayLabels(5:5:end) = temp(5:5:end);
h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';
h.Title = 'Ramp duration affects acclimation (H)';

nexttile()
h = heatmap(acclim);
h.YLabel = '18C Ramp Duration (Days)'; h.XLabel = 'Time (Days)';
h.XDisplayLabels(:) = {''};
h.XDisplayLabels(5:5:end) = Xlabels(5:5:end);
temp = h.YDisplayLabels(:);
h.YDisplayLabels(:) = {''};
h.YDisplayLabels(5:5:end) = temp(5:5:end);
h.ColorLimits = [0 1];
h.ColorbarVisible = 'off';
h.Title = 'Ramp duration affects acclimation (\alpha)';

load("cmaplist.mat");
colormap(cmaplist.cmapRedBlue1);


% %%% VAR RAMPS FOLLOWED BY VAR 4C DURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nexttile();
% y0 = [22/22 1];
% c = ModelData.csol.c;
% numslopes = 40;
% MaxDuration4C = 40;
% health = zeros(numslopes,MaxDuration4C);
% 
% for i = 1:numslopes
%     for j = 1:MaxDuration4C
%         tspan = 0:1:(i+j-1);
%         temperature1 = tspan;
%         temperature2 = [linspace(22/22,4/22,i) linspace(4/22,4/22,j)];
%         temperature = [temperature1;temperature2];
%         sol = ode23(@(t,y) Acclimfun(t,y,c,temperature) , tspan, y0);
%         solpts = deval(sol,tspan);
%         health(i,j) = solpts(2,end);
%     end
% end
% 
% Xlabels = [];
% for i=1:MaxDuration4C Xlabels = [Xlabels {(num2str(i))}]; end
% Ylabels = [];
% for i=1:numslopes Ylabels = [Ylabels {(num2str(i))}]; end
% 
% h = heatmap(health);
% h.YLabel = 'Ramp Duration (Days)'; h.XLabel = 'Time after ramp (Days)';
% h.XDisplayLabels(:) = {''};
% h.XDisplayLabels(5:5:end) = Xlabels(5:5:end);
% temp = h.YDisplayLabels(:);
% h.YDisplayLabels(:) = {''};
% h.YDisplayLabels(5:5:end) = temp(5:5:end);
% h.ColorLimits = [0 1];
% h.ColorbarVisible = 'off';
% h.Title = 'Effect of 4C duration after ramps';
% %%% OVERNIGHT VS SEASONAL DROPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % The duration of the period after which the temperature has dropped is
% % important
% nexttile()
% amax = 22;
% % dmax = linspace(24,24*20,24*20-23);
% dmax = linspace(1,61,61);
% % InitialTemp = amax+1;
% InitialTemp = amax+1;
% InitialTemp = 22;
% health = zeros(amax,length(dmax));
% label = [];
% for a=1:amax
%     for di=1:length(dmax)
%         d = dmax(di);
%         if (a == 1) 
%             label = [label {num2str(d)}];
%         end
%         % tspan = 0:0.1:d/24;
%         tspan = linspace(0,d,100);
%         amplitude = a/22;
%         y0 = [InitialTemp/22 1];
%         c = ModelData.csol.c;
%         temperature1 = tspan;
%         temperature2 = [linspace(InitialTemp/22-amplitude,InitialTemp/22-amplitude,length(tspan))];
%         temperature = [temperature1;temperature2];
%         sol = ode23(@(t,y) Acclimfun(t,y,c,temperature) , tspan, y0);
%         solpts = deval(sol,tspan);
%         health(a,di) = [(solpts(2,end))];
%     end
% end
% % plot(health);
% h = heatmap(health);
% h.YLabel = 'Drop Amplitude (C)'; h.XLabel = 'Time after drop (days)';
% h.XDisplayLabels(:) = {''};
% h.XDisplayLabels(5:5:end) = label(5:5:end);
% temp = h.YDisplayLabels(:);
% h.YDisplayLabels(:) = {''};
% h.YDisplayLabels(5:5:end) = temp(5:5:end);
% h.ColorLimits = [0 1];
% h.ColorbarVisible = 'off';
% h.Title = 'Drop ampl. is important';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f = gcf;
% exportgraphics(f, [f.Name '.jpg'])
exportgraphics(f, [f.Name '.pdf'])

end
