% This codes creates the plots in Figure 8.
% It requires the Advanpix Multiprecision Toolbox 
% to be in the MATLAB path.

clear all, close all

myaddpath;

mp.Digits(128)                   % set working precision to 128 digits

N = 1e3; Z = logspace(-6,6,N); Z = Z(:);  % discretized sample pts
C = 0.1*(N+1)^2*gallery('tridiag',N);     % C is tridiagonal
M = speye(1000);                          % M is the identity matrix
q = randn(N,1); q = q/norm(q);            % normalized random vector q
e0 = M\q;                                 % initial vector e0 = e(0)
Mnorm = sqrt((e0')*M*e0);                 % normalized const
K2M = condest(M);                         % condition number of M
rng('default')

for decades = 1:4

t = logspace(-decades,0,1+10*decades);    % discretized t
f = exp(-t.*Z);                           % functions to be fitted
exact = [];
for j = 1:length(t)
exact(:,j) = expm(-t(j)*full(C))*q;       % exact solutions
end

for n = 1:60                              % degree of approximation

[pol, ~, vals, ss,zer] = lspf_opt_inf2(n, Z, t);  
% Poles (pol) and interpolation pts (zer) are obtained from an optimizer
% that minimizes the total error bound := floating error bound + scalar
% approximation bound

zer = zer(:);
f2 = exp(-t.*zer);                        % f at the interpolation pts
subCmat = 1./(mp(zer)-mp(pol));           % Cauchy matrix for interpolation
coef = subCmat\mp(f2);                    % obtain residues using multiprecision
double_coef = double(coef);               % residues cast into double
Cmat = 1./(mp(Z)-mp(pol));                % full Cauchy matrix for evaluation
scalerr = vecnorm(double(Cmat*coef-mp(f)), inf);  % scalar approximation bound

%% Compute the ingredients of our floating-point error bound E1, E2, E3
for k = 1:length(pol)
    V(:,k) = (C-pol(k)*M)\q;
    rhs = (C-pol(k)*M)*V(:,k)-q;
    scalres = M\rhs;
    res(k) = sqrt((scalres')*M*scalres);
    v = V(:,k);
    Mv(k) = sqrt((v')*M*v);
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
    errbnd3(j) = sqrt((ee')*M*ee);
end
errbnd3 = sqrt(K2M)*errbnd3;           % E3
comb_errbnd = errbnd1+errbnd2+errbnd3; % total floating-point error bound
total_err_bnd(decades, n) = max(double(scalerr*Mnorm)+double(comb_errbnd));
end

V = [];
rhs = [];
scalres = [];
res = [];
v = [];
Mv = [];

end

%% Plot
figure
semilogy(1:60,total_err_bnd(1,:),'LineWidth',2);
hold on
semilogy(1:60,total_err_bnd(2,:),'LineWidth',2);
hold on
semilogy(1:60,total_err_bnd(3,:),'LineWidth',2);
hold on
semilogy(1:60,total_err_bnd(4,:),'LineWidth',2);
hold off
legend('$t_{\max}/t_{\min} = 10^1$','$t_{\max}/t_{\min} = 10^2$','$t_{\max}/t_{\min} = 10^3$','$t_{\max}/t_{\min} = 10^4$','interpreter','latex','Location','northeast','FontSize',18);
set(gcf,'Color','w')
xlabel('degree $n$','Interpreter','latex','FontSize',18)
ylabel('total error bound','interpreter','latex','FontSize',18)

