% This codes creates the plots in Figure 1.
% It requires the Rational Krylov Toolbox 
% to be in the MATLAB path.

clear all, close all
unzip('http://guettel.com/rktoolbox/rktoolbox.zip'); 
cd('rktoolbox'); addpath(fullfile(cd)); savepath

t = logspace(-3,0,21);           % time points
C3 = 4.31;                       % optimal linear constant for t, used for unweighted approximation
N = 1e4; Z = logspace(-6,6,N); Z = Z(:);   % discretization
f = exp(-t.*Z);                  % functions to be fitted

%% apply RKFIT to find optimal rational approximation with given poles
D = spdiags(Z,0,N,N); u = ones(N,1);

for j = 1:length(t) F{j} = spdiags(exp(-t(j)*Z),0,N,N); end

k = 0;                           % type (n,n)                 
param.k = k;                     % subdiagonal approximant
param.maxit = 0; param.tol = 0;  % 0 RKFIT iterations
param.real = 1;                  % data is real-valued

%%
for k = 1:3
    n = 12*k;                        % degree of approximation
    pol = -C3*n*ones(1,n);           % optimal concentrated poles for unweighted approximation
    polw = (-n./(sqrt(2)*max(t)))*ones(1,n);   % optimal concentrated poles for weighted approximation
    
    [xi, ratfun, misfit, out] = rkfit(F, D, u, pol, param);
    % Compute optimal approximant with unweighted poles using RKFIT
    [xiw, ratfunw, misfitw, out] = rkfit(F, D, u, polw, param);
    % Compute optimal approximant with weighted poles using RKFIT

    for j = 1:length(t) 
        err(k,j) = norm(ratfun{j}(D,u)-f(:,j), inf);   % error of unweighted approximation
        errw(k,j) = norm(ratfunw{j}(D,u)-f(:,j), inf); % error of weighted approximation
    end

end

%% plot (unweighted approximation)
figure(1)
loglog(t, err(1,:), '-*','LineWidth',2);
hold on
loglog(t, err(2,:), '-*','LineWidth',2);
hold on
loglog(t, err(3,:), '-*','LineWidth',2);
hold off
xlabel('time $t$','Interpreter','latex','FontSize',18)
ylabel('$L^\infty$-norm error','Interpreter','latex','FontSize',18);
title('Unweighted approximation','FontWeight','normal','FontSize',18,'Interpreter','latex');
set(gcf,'Color','w');
set(gca,'FontSize',18,'TickLabelInterpreter','latex')
legend('$n = 12$','$n = 24$','$n = 36$', 'Location','northeast','FontSize',18,'interpreter','latex')

%% plot (weighted approximation)
figure(2)
loglog(t, errw(1,:), '-*','LineWidth',2);
hold on
loglog(t, errw(2,:), '-*','LineWidth',2);
hold on
loglog(t, errw(3,:), '-*','LineWidth',2);
hold off
xlabel('time $t$','Interpreter','latex','FontSize',18)
ylabel('$L^\infty$-norm error','Interpreter','latex','FontSize',18);
title('Weighted approximation','FontWeight','normal','FontSize',18,'Interpreter','latex');
set(gcf,'Color','w');
set(gca, 'FontSize', 18, 'TickLabelInterpreter','latex')
legend('$n = 12$','$n = 24$','$n = 36$', 'Location','northeast','FontSize', 18,'interpreter','latex')

