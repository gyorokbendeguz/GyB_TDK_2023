function [yprime,ymod] = equ_ol(t,y,~,parms)
% Brush tyre with sliding at zero speed
% y = [Mz; qi]
%
%

x = parms.x;
dx = parms.dx;
l = parms.l;
k = parms.k;
% J = parms.J;
n = parms.n;
epsilon = parms.epsilon;
% Psi = parms.Psi;
Time = parms.Time;
Fz = parms.Fz;
Omega = parms.Omega;
mu = parms.mu;
mu_slip = parms.mu_slip;
a = parms.a;

% psi = interp1(Time, Psi, t);
omega = interp1(Time, Omega, t);
fz = interp1(Time, Fz, t);
q = y(2:end);

qcrp = 3*fz*mu/(4*a^3*k)*(a^2-x.^2);
q_slip = 3*fz*mu_slip/(4*a^3*k)*(a^2-x.^2);

qd = zeros(1,n+1);

for i=1:n+1
    if abs(q(i))< qcrp(i)
        qd(i) = (l-x(i))*omega;
    else
        qd(i) = 0;
        q(i) = (q_slip(i)-epsilon)*sign(q(i));
    end
end

defmom=0;
defmomd=0;
for i=2:n+1
    defmom=defmom+(l-(x(i)+x(i-1))/2)*(q(i)+q(i-1))/2*dx;
    defmomd=defmomd+(l-(x(i)+x(i-1))/2)*(qd(i)+qd(i-1))/2*dx;
end

% if -k*defmom > 6.3552 % from meas.
%     defmom = sign(defmom)*6.3552/k;
%     defmomd = 0;
% end
% if -k*defmom < -5.8531
%     defmom = sign(defmom)*5.8531/k;
%     defmomd = 0;
% end

ymod = [-k*defmom, q];
yprime = [-k*defmomd, qd];

end
