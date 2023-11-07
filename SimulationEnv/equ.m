function [yprime,ymod] = equ(t,y,~,parms)
% Brush tyre with sliding at zero speed
% y = [psi; psid; qi]
%
%

x = parms.x;
dx = parms.dx;
l = parms.l;
k = parms.k;
J = parms.J;
n = parms.n;
qcrp = parms.qcrp;
q_slip = parms.q_slip;
epsilon = parms.epsilon;
kp = parms.kp;
kd = parms.kd;
ki = parms.ki;
psidesprofile = parms.psidesprofile;

psi = y(1);
ome = y(2);
q = y(3:n+3);
z = y(n+4);

qd = zeros(1,n+1);

for i=1:n+1
    if abs(q(i))< qcrp(i)
        qd(i) = (l-x(i))*ome;
    else
        qd(i) = 0;
        q(i) = (q_slip(i)-epsilon)*sign(q(i));
    end
end

defmom=0;
%defmomd=0;
for i=2:n+1
    defmom=defmom+(l-(x(i)+x(i-1))/2)*(q(i)+q(i-1))/2*dx;
    %defmomd=defmomd+(l-(x(i)+x(i-1))/2)*(qd(i)+qd(i-1))/2*dx;
end

psides = [];
%psides = interp1(psidesprofile(:,1),psidesprofile(:,2),t,'linear','extrap');
for i = 1:length(psidesprofile(:,1))
    if t <= psidesprofile(i,1)
        psides = psidesprofile(i,2);
        break;
    end
end

if isempty(psides)
    psides = psidesprofile(end,2);
end

M = -kp*(psi-psides)-kd*ome-ki*z;

zd = psi-psides;

ymod = [psi, ome, q, z];
yprime = [ome, -k/J*defmom + M/J, qd, zd];


end

