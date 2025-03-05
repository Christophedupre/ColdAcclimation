function dydt = AcclimfunNoise(t,y,c,temperature)
temp = interp1(temperature(1,:),temperature(2,:),t); %selects the right temperature

dydt = [(temp-y(1))*c(1);((temp-y(1)))*c(2) + y(2)*temp*c(3)]; %regeneration depends on health and temp.

%If health goes below threshold, the animal dies
if(y(2)<=0.2) 
    dydt(2) = -y(2);
end

%If health goes above one, delta health can only be negative
if(y(2)+dydt(2)>=1.3 && dydt(2)>0) 
        dydt(2) = 1-y(2); %y(2) will reach 1 but won't exceed it
end

end


