function dydt = Acclimfun2(t,y,c,temperature)
temp = interp1(temperature(1,:),temperature(2,:),t); %selects the right temperature

if(y(1)<temp)
    c1f = c(1)+0.5;
    % c1f = c(1);
else
    c1f = c(1);
end

dydt = [(temp-y(1))*c1f;((temp-y(1)))*c(2) + y(2)*temp*c(3)]; %regeneration depends on health and temp.

if(y(2)<=0.2) %If health goes below threshold, the animal dies
    dydt(2) = -y(2);
%     dydt(2) = 0;
end

if(y(2)+dydt(2)>=1 && dydt(2)>0) %If health goes above one, delta health can only be negative
        dydt(2) = 1-y(2); %y(2) will reach 1 but won't exceed it
end

end


