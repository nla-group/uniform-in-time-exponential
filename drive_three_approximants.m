% This codes creates the plot in Figure 3.
% It requires the Rational Krylov Toolbox 
% and the Advanpix Multiprecision Toolbox 
% to be in the MATLAB path.

clear all, close all

myaddpath;

mp.Digits(128) % set working precision to 128 digits

t = logspace(-3, 0, 31);   % time
C3 = 4.31;                 % optimal linear factor for concentrated real poles
N = 1e3; Z = logspace(-6,6,N); Z = Z(:);  % discretization 
f = exp(-t.*Z);            % functions to be fitted  
n1 = 33;                   % number of concentrated real poles
n2 = 21;                   % number of distinct real poles
n3 = 18;                   % number of RKFIT complex poles

%% Approximation using distinct real poles
[dist_pols, ~, ~, ~,zers] = lspf_opt_inf(n2, Z, t); 
% dist_pols: distinct real poles; zers: interpolation nodes (both Zolotarev-based)
subCmat = 1./(mp(zers)-mp(dist_pols));    % sub Cauchy matrix for interpolation
subf = exp(-t.*zers);                     % interpolated RHS functions
coef = subCmat\subf;                      % residues
Cmat = 1./(mp(Z)-mp(dist_pols));          % full Cauchy matrix for evaluation
errdist = vecnorm(double(Cmat*coef-f), inf);  % error of using distinct-real-pole approximation

%% Approximation using concentrated real poles and RKFIT complex poles
D = spdiags(Z,0,N,N); u = ones(N,1);
conc_poles = -C3*n1*ones(1,n1); % optimal concentrated real poles
xistart = inf(1,n3);

for j = 1:length(t) F{j} = spdiags(f(:,j),0,N,N); end

k = -1;                          % type (n-1,n)                 
param.k = k;                     % subdiagonal approximant
param.maxit = 0; param.tol = 0;  % 0 RKFIT iterations
param.real = 1;                  % data is real-valued

param2.k = k;                     % subdiagonal approximant
param2.maxit = 16; param2.tol = 0; % 16 RKFIT iterations
param2.real = 1;                  % data is real-valued

[xi1, ratfun1, misfit1, out] = rkfit(F, D, u, conc_poles, param);
% Approximation using concentrated poles

[xi2, ratfun2, misfit2, out] = rkfit(F, D, u, xistart, param2);  
% Approximation using RKFIT refined poles

for j = 1:length(t)
    errconc(j) = norm(ratfun1{j}(D, u) - f(:,j), inf); % error of using concentrated-real-pole approximation
    errrkfit(j) = norm(ratfun2{j}(D, u) - f(:,j), inf);% error of using RKFIT-complex-pole approximation
end

%% plot
figure(1)
loglog(t, errconc, '-o','Color','b', 'LineWidth',2);
hold on
loglog(t, errdist, '-*','Color','r','LineWidth',2);
hold on
loglog(t, errrkfit, '-s','Color','k','LineWidth',2);
hold off
set(gcf, 'color', 'w');
xlabel('time $t$','FontSize',14,'Interpreter','latex')
ylabel('$L^{\infty}$-norm error','Interpreter','latex','FontSize',14);
title('Approximate $\exp(-tz)$ for $z \geq 0$', 'Interpreter','latex','FontWeight','normal','FontSize',14);
legend('concentrated ($33$ real poles)','Zolotarev ($21$ real poles)','RKFIT ($18$ complex poles)','Location','northeast','interpreter','latex','FontSize', 14);
ylim([1e-14, 1e-0]); yticks([1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2 1e-0]);
set(gca,'FontSize',14, 'TickLabelInterpreter','latex')
