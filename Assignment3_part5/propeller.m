function [T, Qm, dn] = propeller(n, Qd, Qm, D, rho, KT, KQ, h)

% Constants
Km = 0.6; Tm = 10; Im = 100000;
Q = KQ * (rho*D^5*abs(n)*n);

% Gains
Y = Qd / Km;

dQm = (1 / Tm) * ( - Qm + Y * Km ); 
Qm = Qm + dQm*h;


% Calculate updated n
dn = 1/Im * (Qm - Q);

n = n + h*dn;


% Q = KQ * (rho*D^5*abs(n)*n);
T = KT * (rho*D^4*abs(n)*n);


end

