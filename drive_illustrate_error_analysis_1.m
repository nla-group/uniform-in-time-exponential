% This codes creates the first subplot in Figure 5.
% It requires the Advanpix Multiprecision Toolbox 
% to be in the MATLAB path.

clear all, close all

myaddpath;

mp.Digits(128) % set working precision to 128 digits

N = 1e3;
decades = 3;
t = logspace(-decades,0,1+10*decades)/1;  % discretized t
C = 0.1*(N+1)^2*gallery('tridiag',N);     % tridiagonal C
rng('default')
q = randn(N,1); q = q/norm(q);            % normalized q
exact = [];
for j = 1:length(t)
exact(:,j) = expm(-t(j)*full(C))*q;       % exact solution
end

M = speye(1000);            % M is the identity matrix in this example
e0 = M\q;                   % initial vector e0 = e(0)
Mnorm = sqrt((e0.')*M*e0);  % normalization const
K2M = 1;                    % K2M is the condition number of M
N = 1e3; Z = logspace(-6,6,N); Z = Z(:);  % discretized sample pts

f = exp(-t.*Z);   % functions to be fitted

n = 25; % degree of the shared real-pole approx

[pol, ~, ~, ~,zer] = lspf_opt_inf(n, Z, t);  % pol, zer: Zolotarev poles and interpolation pts

zer = zer(:);
f2 = exp(-t.*zer); % f evaluated at the interpolation pts

subCmat = 1./(mp(zer)-mp(pol));  % sub Cauchy matrix for interpolation

coef = subCmat\mp(f2);  % obtain residues using multi-precision
double_coef = double(coef);  % residues cast into double
Cmat = 1./(mp(Z)-mp(pol));   % full Cauchy matrix for evaluation
scalerr = vecnorm(double(Cmat*coef-mp(f)), inf);  % scalar error bound

%% Compute the ingredients for the floating-point error bound E1, E2, E3
for k = 1:length(pol)
    V(:,k) = (C-pol(k)*M)\q;
    rhs = (C-pol(k)*M)*V(:,k)-q;
    scalres = M\rhs;
    res(k) = sqrt((scalres.')*M*scalres);
    v = V(:,k);
    Mv(k) = sqrt((v.')*M*v);
end
abs_sigma = -pol;
u = 1.11e-16;
gam_n = n*u/(1-n*u);

%% Compute E1, E2, E3
errbnd1 = (res./abs_sigma)*abs(coef);  % E1 
errbnd2 = Mv*abs(coef);
errbnd2 = u*errbnd2;                   % E2 

err3mat = gam_n*abs(V)*abs(double_coef);
for j = 1:length(t)
    ee = err3mat(:,j);
    errbnd3(j) = sqrt((ee.')*M*ee);
end
errbnd3 = sqrt(K2M)*errbnd3;          % E3 

comb_errbnd = errbnd1+errbnd2+errbnd3; % Total floating-point error

acterrmat = V*double_coef-exact;       % actual computed error
for j = 1:length(t)
    acterr = acterrmat(:,j);
    actual_err(j) = sqrt((acterr.')*M*acterr); 
end

%% plot
figure(1)
loglog(t, scalerr*Mnorm, 'b-*', 'LineWidth',2);  % scalar error bound of actual_err
hold on
loglog(t, actual_err, 'r-*','LineWidth',2);  % actual error: || fl(r_t(A)b) - expm(-tA)b) ||
hold on
loglog(t, comb_errbnd, 'k-x', 'LineWidth',2); % total bound for the floating err 
hold on
loglog(t, errbnd1, 'k-s', 'LineWidth',2); % bound1 for the floating err 
hold on
loglog(t, errbnd2, 'k-o', 'LineWidth',2); % bound2 for the floating err 
hold on
loglog(t, errbnd3, 'k-+', 'LineWidth',2); % bound3 for the floating err 
hold off
%legend('scalar err','computed err','combined','bound 1','bound 2','bound 3','Location','northeast','NumColumns',3,'Orientation','horizontal','interpreter','latex','FontSize',18)
set(gcf, 'color','w');
ylim([1e-14 1e-2]);
yticks([1e-14 1e-12 1e-10 1e-8 1e-6 1e-4 1e-2]);
set(gca,'FontSize',18, 'TickLabelInterpreter','latex')

