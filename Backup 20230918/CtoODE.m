function solpts = CtoODE(c,tspan,y0)
solpts = zeros(44,15,16);
for j = 1:15;%Over all training durations
    for i=1:16 %Over all training temperatures
        temperature1 = linspace(0,43,44);
        temperature2 = [linspace(22,22,15-(j-1)) linspace(i*2+2,i*2+2,j-1) linspace(4,4,29)]./22;
        temperature = [temperature1; temperature2];
        sol = ode45(@(t,y) Acclimfun(t,y,c,temperature), tspan, y0);
        solptsinterm = deval(sol,tspan);
        solptsinterm(1,:) = []; %only health is used for optimization
        solpts(:,j,i) = solptsinterm;
    end
end
end