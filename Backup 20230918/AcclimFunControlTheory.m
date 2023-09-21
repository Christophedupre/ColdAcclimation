function [Healthlevel Acclimlevel Atarget] = AcclimFunControlTheory(Temperature);

ki = 0.01; %Integration constant
HealthDecaySpeed = 5; %Speed at which health decays if temperature is too low
Acclimlevel = zeros(1,length(Temperature))+1;
Atarget = zeros(1,length(Temperature))+1;
DeltaA = zeros(1,length(Temperature));
Healthlevel = zeros(1,length(Temperature))+1;
Time = linspace(0,length(Temperature),length(Temperature));

for i = 1:length(Temperature)
%     Atarget(i) = (22-(Temperature(i)-4))./22;
    Atarget(i) = Temperature(i)./22;
    if(i>5) %integrator can only start after 5th timepoint
        DeltaA(i) = ki*sum(Atarget(i-5:i)-Acclimlevel(i-5:i));
    else
        DeltaA(i) = ki*sum(Atarget(1)-Acclimlevel(1))*5;
    end
    
    if i == 1
        Acclimlevel(i) = 1 + DeltaA(i);
        Healthlevel(i) = 1 + DeltaA(i)*HealthDecaySpeed;
    else
        Acclimlevel(i) = Acclimlevel(i-1) + DeltaA(i);
        Healthlevel(i) = Healthlevel(i-1) + DeltaA(i)*HealthDecaySpeed;
    end
    % Regeneration
    if (Healthlevel(i) < 1)
        Healthlevel(i) = Healthlevel(i) + (1-Healthlevel(i))*0.2;
    end
    if(Healthlevel(i) > 1)
        Healthlevel(i) = 1; % Healthlevel cannot exceed 1;
    end
end
% Atarget(i+1) = Atarget(i); %for the last datapoint

end