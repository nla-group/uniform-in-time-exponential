% This codes creates the plots in Figure 2.
% It requires the Advanpix Multiprecision Toolbox 
% to be in the MATLAB path.

clear all, close all

myaddpath;

t = logspace(-3, 0, 31);                 % time                      
N = 1e3; Z = logspace(-6,6,N); Z = Z(:); % discretization
f = exp(-t.*Z);                          % functions to be fitted
m = 25;                                  % degree
[pols, ~, ~, ~,zers] = lspf_opt_inf(m, Z, t); % pols: Zolotarev-based poles; zers: interpolation nodes

r = @(x) 1;    % construct rational nodal function
for j = 1:m
    r = @(x) r(x) .* (x - zers(j))./(x - pols(j));
end

%% Functions evaluated at the mesh grids
xx = -10:0.01:10;
yy = -10:0.01:10;
[X, Y] = meshgrid(xx, yy);
ZZ = X + 1i*Y;
F1 = exp(-min(t).*ZZ)./r(ZZ);
F2 = exp(-max(t).*ZZ)./r(ZZ);

%% Plots
figure(1)
contour(X, Y, log10(abs(F1)), 30, 'LineWidth',2)
colorbar
title('$\log_{10} \vert \exp(-t_{min}z)/s_n(z) \vert$','Interpreter','latex','FontSize',18)
set(gcf, 'Color','w');
hold on
plot(real(pols(21:25)), imag(pols(21:25)),'ko','MarkerFaceColor','k');
hold on
plot(real(zers(20:25)), imag(zers(20:25)),'ro','MarkerFaceColor','r');
hold off
set(gca,'TickLabelInterpreter','latex','FontSize',18)
cb = colorbar;
cb.TickLabelInterpreter = 'latex';

figure(2)
contour(X, Y, log10(abs(F2)), 30, 'LineWidth',2)
colorbar
title('$\log_{10} \vert \exp(-t_{max}z)/s_n(z) \vert$','Interpreter','latex','FontSize',18)
set(gcf, 'Color','w');
hold on
plot(real(pols(21:25)), imag(pols(21:25)),'ko','MarkerFaceColor','k');
hold on
plot(real(zers(20:25)), imag(zers(20:25)),'ro','MarkerFaceColor','r');
hold off
set(gca,'TickLabelInterpreter','latex','FontSize',18)
cb.TickLabelInterpreter = 'latex';
