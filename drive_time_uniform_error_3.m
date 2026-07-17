% This codes creates the third subplot in Figure 4.
% It requires the Rational Krylov Toolbox 
% and the Advanpix Multiprecision Toolbox 
% to be in the MATLAB path.

clear all, close all

myaddpath;

mp.Digits(128) % set working precision to 128 digits

%% Time ratio tmax/tmin = 10^3
t = logspace(-3, 0, 31);   % time
C3 = 4.31;                 % optimal linear factor for concentrated real poles

N = 1e3; Z = logspace(-6,6,N); Z = Z(:);  % discretization 
f = exp(-t.*Z);            % functions to be fitted  

D = spdiags(Z,0,N,N); u = ones(N,1);
for j = 1:length(t) F{j} = spdiags(f(:,j),0,N,N); end

k = -1;                          % type (n-1,n)                 
param.k = k;                     % subdiagonal approximant
param.maxit = 0; param.tol = 0;  % 0 RKFIT iterations
param.real = 1;                  % data is real-valued

param2.k = k;                     % subdiagonal approximant
param2.maxit = 16; param2.tol = 0; % 16 RKFIT iterations
param2.real = 1;                  % data is real-valued

for n = 1:30                  % degree of approximation 

% Approximation using distinct real poles
[dist_pols, ~, ~, ~,zers] = lspf_opt_inf(n, Z, t); 
% dist_pols: distinct real poles; zers: interpolation nodes (both Zolotarev-based)
subCmat = 1./(mp(zers)-mp(dist_pols));    % sub Cauchy matrix for interpolation
subf = exp(-t.*zers);                     % interpolated RHS functions
coef = subCmat\mp(subf);                      % residues
Cmat = 1./(mp(Z)-mp(dist_pols));          % full Cauchy matrix for evaluation
errdist(n) = max(vecnorm(double(Cmat*coef-f), inf));  % time-uniform error of using distinct-real-pole approximation

% Approximation using concentrated real poles and RKFIT complex poles
conc_poles = -C3*n*ones(1,n); % optimal concentrated real poles
xistart = inf(1,n);

[xi1, ratfun1, misfit1, out] = rkfit(F, D, u, conc_poles, param);
% Approximation using concentrated poles

[xi2, ratfun2, misfit2, out] = rkfit(F, D, u, xistart, param2);  
% Approximation using RKFIT refined poles

for j = 1:length(t)
    err_conc(j) = norm(ratfun1{j}(D, u) - f(:,j), inf); % error of using concentrated-real-pole approximation
    err_rkfit(j) = norm(ratfun2{j}(D, u) - f(:,j), inf);% error of using RKFIT-complex-pole approximation
end

errconc(n) = max(err_conc); % time-uniform error of using concentrated-real-pole approximation
errrkfit(n) = max(err_rkfit); % time-uniform error of using RKFIT-complex-pole approximation

end

figure
semilogy(1:30, errconc, '-o','Color','b', 'LineWidth',2);
hold on
semilogy(1:30, errdist, '-*','Color','r','LineWidth',2);
hold on
semilogy(1:30, errrkfit, '-s','Color','k','LineWidth',2);
hold off
set(gcf, 'color', 'w');
xlabel('degree of approximation','FontSize',18,'Interpreter','latex')
ylabel('time-uniform error','Interpreter','latex','FontSize',18);
title('$T = [10^{-3}, 1]$', 'Interpreter','latex','FontWeight','normal','FontSize',18);
%legend('concentrated poles','Zolotarev poles','RKFIT poles','Location','northeast','interpreter','latex','FontSize', 18);
ylim([1e-8 1e0]);
yticks([1e-8 1e-6 1e-4 1e-2 1e0]);
set(gca,'FontSize',18, 'TickLabelInterpreter','latex')
