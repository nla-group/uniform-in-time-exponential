function [pols, res, inf_err, s, interp] = lspf_opt_inf(m, x, t)
% Poles (pols), residues (res), interpolation nodes (interp)
% for the time-uniform partial fraction approximation of exp(-tx) on a possibly 
% unbounded positive interval x.
% Here, t should be a positive bounded time interval.
% The poles and interpolation nodes are obtained by solving a Zolotarev 
% approximation problem on the condenser [-inf,0] U [c,d]. As the interval 
% [c,d] is not known explicitly (depends on the degree, among other things), 
% it is found by numerical optimization, which minimizes the scalar 
% approximation bound. The residues are obtained by the interpolation condition.
% The discrete time-uniform error (inf_err) and refined pole-interval (s)
% are also returned.

mp.Digits(128) % set working precision to 128 digits

x = x(:); t = t(:).';
E = exp(-x*t); % exact values to be fitted

% initial guess for the finite condenser plate (log10 of boundaries)
s = log10([m/sqrt(2)/max(t) , m/sqrt(2)/min(t)]); % guess informed by Anderssen (1981)
%options = optimset('MaxIter', 1500, 'TolX', 1e-8, 'TolFun', 1e-8);

if 1
    disp('LSPF_OPTIM: starting 2-dim optim')
    disp(fun(s, E,x,m))
    s = fminsearch(@(ss) fun(ss, E,x,m),s);
    disp(fun(s, E,x,m))
    %{
    disp('LSPF_OPTIM: follow up with alternating 1-dim optim')
    for it = 1:3
        s1 = fminsearch(@(ss) fun([ss,s(2)], E,x,m),s(1));
        s2 = fminsearch(@(ss) fun([s1,ss], E,x,m),s(2));
        s = [s1,s2];
        disp(fun(s, E,x,m))
    end

    disp('LSPF_OPTIM: 2-dim optim once more')
    s = fminsearch(@(ss) fun(ss,E,x,m),s);
    disp(fun(s, E,x,m))
    %}
end

% Find poles and interpolation nodes from the refined condenser
% method is the same as in the function fun
[interp,~,~,pols] = dkz09_mapped2(m, -10^s(2), -10^s(1)); 
interp = interp(:);
C_sub = 1./(mp(interp)-mp(pols));
E_sub = exp(-interp*t);
res = C_sub\mp(E_sub);
C = 1./(mp(x)-mp(pols));
R = mp(E) - C*res;
inf_err = max(vecnorm(double(R), inf)); % discrete time-uniform error

function val = fun(s,E,x,m)
    [intp,~,~,pols] = dkz09_mapped2(m, -10^s(2), -10^s(1)); 
    % poles (pols) and interpolation nodes (intp) from a given condenser
    intp = intp(:);
    C_sub = 1./(mp(intp)-mp(pols));  % sub Cauchy matrix for interpolation
    E_sub = exp(-intp*t);
    coeff = C_sub\mp(E_sub);         % residues
    C = 1./(mp(x)-mp(pols));         % full Cauchy matrix for evaluation
    R = mp(E) - C*coeff;
    val = log10(max(vecnorm(double(R), inf))); % log10-based time-uniform error 
end
end
