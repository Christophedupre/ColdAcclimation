function Explicitsolution(Model1Data)

figModelPredictions = figure('Name','Model Predictions','WindowState','maximized');
figure(figModelPredictions)
clf;
tiledlayout(4,4)
cmaplines = colormap('lines');
%%% ACCLIMATION - SHOW THAT THE EXPLICIT SOLUTION IS CORRECT %%%%%%%%%%%%%%
nexttile(); hold on;

c1 = Model1Data.csol.c(1);
theta = 4/22;%Test temperature
alpha0 = 1;%Initial acclimation level 
tspan = linspace(0,14,15);

%ODE, rewritten independently from acclimfun
% tspan = linspace(0,14,15);
% sol = ode45(@(t,alpha) (4/22-alpha)*c1,tspan,1);
% solpts = deval(sol,tspan);

%ODE using acclimfun
cs = Model1Data.csol.c;
y0 = Model1Data.y0;
temperature1 = tspan;
temperature2 = linspace(theta,theta,15);
temperature = [temperature1; temperature2];
sol = ode45(@(t,y) Acclimfun(t,y,cs,temperature), tspan,y0);
solptsODEacclimfun = deval(sol, tspan);

%Explicit solution
A = (alpha0 - theta);
alpha = theta + A*exp(-c1*tspan);

%Prediction
alphat = theta*1.1; %Target Acclimation level
t = (-1/c1)*log((alphat-theta)/A);

p0 = plot(tspan,solptsODEacclimfun(1,:), 'Color',cmaplines(1,:)); % ODE
% p1 = plot(tspan,solpts, 'Color',cmaplines(1,:)); % ODE
p2 = plot(tspan, alpha, 'Marker', '+', 'MarkerEdgeColor', cmaplines(2,:),'Color','none'); % Solution
p3 = plot(t,alphat,'Marker','*','MarkerEdgeColor', cmaplines(1,:),'Color','none'); %Prediction

legend([p0 p2 p3],{'ODE','Explicit Solution',['Prediction for \alpha = ' num2str(alphat)]},'box','off')
ax = gca; ax.Box = 'off'; ax.LineWidth = 1;
ax.YLabel.String = '\alpha';
ax.XLabel.String = 'Time (Days)';

%%% PREDICTIONS FOR VARIOUS STARTING ALPHA0 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha0 = linspace(alphat,1,20);
A = (alpha0 - theta);
t = (-1/c1)*log((alphat-theta)./A);
nexttile()
plot(alpha0,t);
ax = gca; ax.Box = 'off'; ax.LineWidth = 1;
ax.YLabel.String = ['Time to \alpha = ' num2str(alphat,2) ' (days)'];
ax.XLabel.String = '\alpha0';

%%% HEALTH - SHOW THAT THE EXPLICIT SOLUTION IS CORRECT %%%%%%%%%%%%%%%%%%%
c1 = Model1Data.csol.c(1);
c2 = Model1Data.csol.c(2);
c3 = Model1Data.csol.c(3);
alpha0 = 1;%Initial acclimation level
h0 = 1;%Initial health level
A = (alpha0 - theta);
B = h0 - (A*c2)/(c1+theta*c3);

%ODE, rewritten independently from acclimfun
% tspan = linspace(0,14,15);
% sol = ode45(@(t,h) (-A*exp(-c1*t)*c2 + c3*theta*h),tspan,1);
% solptsODE1 = deval(sol,tspan);

%Explicit solution, piecewise-defined
health = ones(1,length(tspan));
h0 = 1; i0 = 0;
for i = 1:length(tspan)    
    health(i) = ((A*c2)/(c1+c3*theta))*exp(-c1*tspan(i)) + B*exp(theta*c3*tspan(i));
    if(health(i)<0.2) %Death
        if(h0 == 1) h0 = health(i); i0 = i; end
        health(i) = h0*exp(-tspan(i-i0+1));
    end
    if(health(i)>1) %Max health reached        
        health(i) = 1;
    end
end

%Prediction
nexttile(); hold on;
% p1 = plot(tspan, solptsODE1, 'Color',cmaplines(1,:)); % ODE
p2 = plot(tspan, solptsODEacclimfun(2,:), 'Color',cmaplines(1,:)); % ODE
p3 = plot(tspan, health, '+', 'MarkerEdgeColor', cmaplines(2,:),'Color','none'); % Solution

legend([p2 p3], {'ODE','Explicit Solution'},'Box','off')
ax = gca; ax.Box = 'off'; ax.LineWidth = 1; ax.YLim(1) = 0; ax.YLim(2) = 1;
ax.YLabel.String = 'H';
ax.XLabel.String = 'Time (Days)';