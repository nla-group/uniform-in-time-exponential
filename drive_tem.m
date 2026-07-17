% This codes creates the plots in Figure 7.
% It requires the Rational Krylov Toolbox 
% and the Advanpix Multiprecision Toolbox 
% to be in the MATLAB path.

% Since this codes require solving linear systems of size 181302 X 181302 with 
% (non-shared) complex shifts computed by CF approximants, it can take
% several hours to run. 

clear all, close all

myaddpath;

mp.Digits(128) % set working precision to 128 digits

load('femdata181302.mat','C','M','q');  % load the matrices C, M and vector q
C = C + tril(C,-1)';
M = M + tril(M,-1)';
% Me'(t) + Ce(t) = 0, for all t, with Me(0) = q.
e0 = M\q;  % initial vector e0 = e(0)
Mnorm = sqrt((e0.')*M*e0);  % normalized const

t = logspace(-6,-3,31);     % time 
N = 1e3; Z = logspace(-12,12,N); Z = Z(:); % discretization 
f = exp(-t.*Z);             % functions to be fitted

nCF = 13; % degree of CF approximation
nZ = 21; % degree of Zolotarev approximation
nR = 18; % degree of RKFIT approximation

%% Compute the CF near best approximants as the reference 
%load('cf_for_exp.mat', 'cfapprx')

[po, re] = modified_cf(nCF);  % poles (po) and residues (re) for approximating exp(-z)

for j = 1:length(t)  % poles and residues scaled for exp(-tz)
    pole = po/t(j);
    pole = pole(:);
    pole = pole.';
    resi = re/t(j);
    for k = 1:length(pole)
        vCF(:,k) = (C-pole(k)*M)\q; 
    end
    cfapprx(:,j) = vCF*resi;   % obtain the CF approximants, this takes a long time to run
end

%% Compute the Zolotarev approximants
[pol, ~, ~, ~,zer] = lspf_opt_inf(nZ, Z, t);  % Zolotarev poles (pol) and interpolation pts (zer)
zer = zer(:);
f2 = exp(-t.*zer); % f at interpolation pts
subCmat = 1./(mp(zer)-mp(pol));  % Cauchy matrix for interpolation
coef_Z = subCmat\mp(f2);  % Zolotarev residues in multiprecision
coefZ = double(coef_Z);  % Zolotarev residues in double precision

for k = 1:length(pol)
   V(:,k) = (C-pol(k)*M)\q;
end

err_Zolo = V*coefZ-cfapprx;  % Zolotarev approximation error

%% Compute the RKFIT approximants

D = spdiags(Z,0,N,N); u = ones(N,1);  % surrogate 
xistart = inf(1,nR);   % initial poles at infinity

for j = 1:length(t) F{j} = spdiags(f(:,j),0,N,N); end

k = -1;
param.k = k;                     % subdiagonal approximant
param.maxit = 16; param.tol = 0; % 16 RKFIT iterations
param.real = 1;                  % data is real-valued

[xi, ratfun, misfit, out] = rkfit(F, D, u, xistart, param);  % xi: the set of RKFIT poles
Cmat = 1./(mp(Z)-mp(xi)); 
coef_R = Cmat\mp(f);    % RKFIT residues in multiprecision
coefR = double(coef_R); % RKFIT residues in double precision

for j = 1:length(xi)
    V_R(:,j) = (C-xi(j)*M)\q;
end

err_rkf = V_R*coefR - cfapprx; % RKFIT approximation error

%% Compute the M-norm error for both
for k = 1:length(t)
    ee_z = err_Zolo(:,k);
    Zolo_Merr(k) = sqrt((ee_z.')*M*ee_z);  % Zolotarev M-norm error
    ee_r = err_rkf(:,k);              
    rkf_Merr(k) = sqrt((ee_r')*M*ee_r); % RKFIT M-norm error
end

%% Plot the normalized M-norm error for both
figure(1)
loglog(t, Zolo_Merr/Mnorm,'r-*','LineWidth',2);
hold on
loglog(t, rkf_Merr/Mnorm,'k-s','LineWidth',2);
hold off
title('$\|r_t(M^{-1}K)\mathbf{b}-\exp(-tM^{-1}K)\mathbf{b}\|_M / \|\mathbf{b}\|_M$','Interpreter','latex','FontSize',18);
legend('Zolotarev','RKFIT','Location','northeast','interpreter','latex','FontSize',18);
xlabel('time $t$','FontSize',18,'Interpreter','latex');
ylabel('$M$-norm error','FontSize',18,'Interpreter','latex');
ylim([1e-9 1e-5]);
set(gcf, 'color','w');
set(gca,'FontSize',18, 'TickLabelInterpreter','latex')


%% plot the RHS decay functions (for tmax)
CZ = 1./(mp(Z)-mp(pol));
tmax_err_zolo = CZ*coef_Z(:,end)-mp(f(:,end));
tmax_err_rkf = Cmat*coef_R(:,end)-mp(f(:,end));
loglog(Z,abs(tmax_err_zolo),'r-','LineWidth',2); hold on;
loglog(Z,abs(tmax_err_rkf),'k-','LineWidth',2); hold off;
title('$|r_{t_{\max}}(z)-\exp(-t_{\max}z)|$','interpreter','latex','FontSize',18);
set(gca,'FontSize',18, 'TickLabelInterpreter','latex')
xlabel('$z$','FontSize',18,'interpreter','latex');
ylabel('Errors at $t_{\max}$','FontSize', 18, 'interpreter','latex');
legend('Zolotarev','RKFIT','FontSize',18,'Interpreter','latex');
yticks([1e-12 1e-10 1e-8 1e-6 1e-4 1e-2]);
