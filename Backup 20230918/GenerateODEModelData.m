%This code optimizes the c1, c2 and c3 parameters of the ODE model
%Optimization method taken from https://www.mathworks.com/help/optim/ug/fit-ode-problem-based-least-squares.html

function [Model1Data Model2Data] = GenerateODEModelData(CstData3D)

%%% CREATE MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ExperimentalData = CstData3D; %Can remove and clean up

%Make optimization problem
tspan = [linspace(0,14,15) linspace(15,43,29)]; %Time vectors for training and testing
y0 = [1 1]; %Initial conditions for acclimation and health
c = optimvar('c',3,'LowerBound',[0.1 0 0],'UpperBound',[0.4 10 10]); %Optimization variable c1,c2,c3
myfcn = fcn2optimexpr(@CtoODE,c,tspan,y0); %Optimization expression
myfcn = myfcn(linspace(15,43,29),:,:); %Only take the time points that correspond to the experimental data
obj = sum(sum(sum((ExperimentalData(:,:,:) - myfcn(:,:,:)).^2))); %objective function. Here we can select a subset of the data if needed
prob = optimproblem('Objective',obj);

% %Solve optimization problem
% c0.c = [0.32 0.38 0.13]; %Initial conditions for c1, c2, c3
% tic
% [csol sumsq] = solve(prob,c0);
% toc

% Skip optimization and use previous results
disp('Optimization skipped, using best results from previous optimizations')
% csol.c = [1.00; 1.12; 0.11]; %Fast c1
% csol.c = [0.37; 0.39; 0.13]; %Intermediate c1
% csol.c = [0.12; 0.18; 0.13]; %Slow c1
load('csFromGUI.mat'); csol.c = cs'; %Get the cs that were found using the GUI

clear c0 sumsq prob obj myfcn c

fprintf('Results: c1 = %1.2f, c2 = %1.2f, c3 =  %1.2f \n',csol.c); %Display results

%%% USE MODEL TO GENERATE DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Model1Data = GenerateEmptyModelStructure();
Model2Data = GenerateEmptyModelStructure();
f = Model1Data.TimeResolution; %Time factor (increases temporal resolution during derivation of ODE solution
NoiseAmplitude = Model1Data.NoiseAmplitude; 
for theta = 1:length(CstData3D(1,1,:))
    for D = 1:length(CstData3D(1,:,1)) 
        %With variability
%         cs = csol.c.*(1+(rand(3,1)-0.5)/Model1Data.cvarcoeff);% Add variability for acclimation, injury and regeneration time constants
        cs = csol.c; % No variability in coefficients
        temperature1 = linspace(0,43,44*f);
        temperature2 = [linspace(22,22,(15-(D-1))*f) linspace(theta*2+2,theta*2+2,(D-1)*f) linspace(4,4,29*f)]./22;
        temperature = [temperature1; AddNoiseToTemp(temperature2, NoiseAmplitude)]; %Variability in temperature sensed
        sol = ode23(@(t,y) AcclimfunNoise(t,y,cs,temperature), tspan,y0);
        solpts = deval(sol, tspan);
        Model2Data.HealthFull(:,D,theta) = solpts(2,:);        
        Model2Data.HealthTruncated(:,D,theta) = solpts(2,linspace(15,43,29));
        Model2Data.AlphaFull(:,D,theta) = solpts(1,:);
        Model2Data.AlphaTruncated(:,D,theta) = solpts(1,linspace(15,43,29));
        Model2Data.csol = csol;

        %Without variability
        cs = csol.c; % No variability
        temperature = [temperature1; temperature2];
        sol = ode23(@(t,y) Acclimfun(t,y,cs,temperature), tspan,y0);
        solpts = deval(sol, tspan);
        Model1Data.HealthFull(:,D,theta) = solpts(2,:);        
        Model1Data.HealthTruncated(:,D,theta) = solpts(2,linspace(15,43,29));
        Model1Data.AlphaFull(:,D,theta) = solpts(1,:);
        Model1Data.AlphaTruncated(:,D,theta) = solpts(1,linspace(15,43,29));
        Model1Data.csol = csol;
    end
end

Model1Data.name = 'ODE without variation';
Model2Data.name = 'ODE with variation';
end