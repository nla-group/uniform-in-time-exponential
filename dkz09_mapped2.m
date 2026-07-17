function [pts,r,rho,pls] = dkz09_mapped2(n, c, d)
% Mapped Zolotarev roots/poles from DKZ09
% n returned pts lie on [0, inf)
% nodal rational function r with these 
% roots and poles on negative axis
% is maximal on [c,d] and minimal on [0,inf)
% Here, c < d < 0 !!!

% first map to symmetric condenser

bb = c+ sqrt(c*(c-d));
dd = 2*c - bb;
tr = @(z) (z - bb)./(z - dd); % this maps [c,d] U [0,+inf) onto [-1,-a],[a,1]

lmin = -tr(d);
lmax = -tr(c);

% now get DKZ points on [-1,-a],[a,1]

del = lmin/lmax;
mu = (1-sqrt(del))/(1+sqrt(del)); mu = mu^2;

%K1 = ellipke(mu^2); K2 = ellipke(1-mu^2);
L = -log(mu)/pi;
[K1,K2] = util_ellipkkp(L);
rho = exp(-pi/4*K2/K1);

%Kp = ellipke(1-del^2);
L = -log(del)/pi;
[~,Kp] = util_ellipkkp(L);

t = (1:2:2*n-1)/(2*n)*Kp;

%[sn,cn,dn] = ellipj(t, 1-del^2);
[sn,cn,dn] = util_ellipjc(1i*t,L);
cn = 1./cn; dn = dn.*cn; sn = -1i*sn.*cn;

pts = dn*lmax;

% now map back to [c,d] U [0,+inf)
trinv = @(x) (bb-dd*x)./(1-x); 

pls = trinv(-pts); % poles on [c,d]
pts = trinv(pts);  % roots on [0,+inf)

%pls = -pts;

r = @(x) 1;
for j = 1:n
    r = @(x) r(x) .* (x - pts(j))./(x - pls(j));
end

end