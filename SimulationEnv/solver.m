function [ynew,ymod] = solver(t,y,dt,parms)

%Runge-Kutta method 3rd order

if isfield(parms, 'OL')
    % Open loop simulation
    [s1,ymod] = equ_ol(t,y,dt,parms);
    s2 = equ_ol(t+dt,ymod+s1*dt,dt,parms);
    s3 = equ_ol(t+(1/2)*dt,ymod+(s1+s2)/4*dt,dt,parms);
    s = (s1 + 4*s3 + s2)/6;
else
    % closed loop sim.
    [s1,ymod] = equ(t,y,dt,parms);
    s2 = equ(t+dt,ymod+s1*dt,dt,parms);
    s3 = equ(t+(1/2)*dt,ymod+(s1+s2)/4*dt,dt,parms);
    s = (s1 + 4*s3 + s2)/6;
end

ynew = ymod + s * dt;
