% This codes creates the plots in Figure 6.
% It requires the Rational Krylov Toolbox 
% and the Advanpix Multiprecision Toolbox 
% to be in the MATLAB path.

clear all, close all

myaddpath;

t = logspace(-6, -3, 31);  % time

N = 1e3; Z = logspace(-12,12,N); Z = Z(:); % discretization

f = exp(-t.*Z);            % functions to be fitted

nZ = 21;   % degree of Zolotarev-real-pole approximation
nR = 18;   % degree of RKFIT-complex-pole approximation

[pols, ~, ~, ~,zers] = lspf_opt_inf(nZ, Z, t);  % pols: Zolotarev poles; zers: Zolotarev interpolation pts
    
%% Now finding the RKFIT poles

D = spdiags(Z,0,N,N); u = ones(N,1);  % surrogate 
xistart = inf(1,nR);   % initial poles at infinity

for j = 1:length(t) F{j} = spdiags(f(:,j),0,N,N); end

k = -1;
param.k = k;                     % subdiagonal approximant
param.maxit = 16; param.tol = 0; % 16 RKFIT iterations
param.real = 1;                  % data is real-valued

[xi, ratfun, misfit, out] = rkfit(F, D, u, xistart, param);  % xi is the set of RKFIT poles
 
load('femdata181302.mat','C','M','q');   % load FEM matrices C, M and vector q (size 181302)
C = C + tril(C,-1)';
M = M + tril(M,-1)';

%% Timing for Zolotarev linear systems
for k = 1:length(pols)
   % solve each shifted lin sys and record the time
   tStart_z = tic;
   xxvec(:,k) = (C-pols(k)*M)\q;
   times_z(k) = toc(tStart_z); % times_z: time for each Zolotarev system
end

%% Timing for RKFIT linear systems
for k = 1:length(xi)/2
   % solve each lin sys and record the time
   tStart_r = tic;
   yyvec(:,k) = (C-xi(2*k-1)*M)\q;  % exploiting the complex conjugacy
   times_r(k) = toc(tStart_r); % times_r: time for each RKFIT system
end

tZ_sum = sum(times_z);
tR_sum = sum(times_r);

%% Plot the timing for Zolotarev linear systems
figure(1)
plot(1:21, times_z,'r-*','LineWidth',2);
hold on
yline(mean(times_z), 'r--', 'LineWidth', 2);
hold off
set(gcf,'color','w');
xlim([1 21])
ylim([0 10])
xticks([1 3 5 7 9 11 13 15 17 19 21])
set(gca, 'TickLabelInterpreter','latex')
xlabel('linear system index','FontSize',18,'Interpreter','latex')
ylabel('time (s)','FontSize',18,'Interpreter','latex')
title('$21$ Zolotarev linear systems','FontWeight','normal','FontSize',18,'Interpreter','latex');

%% Plot the timing for RKFIT linear systems
figure(2)
plot(1:9, times_r, 'k-s','LineWidth',2);
hold on
yline(mean(times_r), 'k--', 'LineWidth', 2);
hold off
set(gcf,'color','w');
xlim([1 9])
ylim([0 50])
xticks([1 3 5 7 9])
set(gca, 'TickLabelInterpreter','latex')
xlabel('linear system index','FontSize',18,'Interpreter','latex')
ylabel('time (s)','FontSize',18,'Interpreter','latex')
title('9 RKFIT linear systems','FontWeight','normal','FontSize',18,'Interpreter','latex')


