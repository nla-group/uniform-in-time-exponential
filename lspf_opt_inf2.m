function [pols, res, total_err, s, interp] = lspf_opt_inf2(m, x, t)
% Poles (pols), residues (res), interpolation nodes (interp)
% for the time-uniform partial fraction approximation of exp(-tx) on a possibly 
% unbounded positive interval x.
% Here, t should be a positive bounded time interval.
% The poles are obtained by solving a Zolotarev approximation problem on
% the condenser [-inf,0] U [c,d]. As the interval [c,d] is not known
% explicitly (depends on the degree, among other things), it is found
% by numerical optimization, which minimizes the total error bound :=
% floating-point error bound (E1) + scalar approximation bound
% The residues are obtained by the interpolation condition.
% The total error bound (total_err) and refined pole-interval (s)
% are also returned.
% We note that the only difference to the algorithm lspf_opt_inf.m is that
% we now minimize the total error bound rather than the scalar bound only.

mp.Digits(128) % set working precision to 128 digits

x = x(:); t = t(:).';
E = exp(-x*t); % exact values to be fitted

% initial guess for the finite condenser plate, informed by Andersson (log10 of boundaries)
s_1 = log10([m/sqrt(2)/(1*max(t)) , m/sqrt(2)/(1*min(t))]); 

% initial guess taking care of the inaccuracy at late time points, which
% happens sometimes
s_2 = log10([m/sqrt(2)/(3*max(t)) , m/sqrt(2)/(1*min(t))]); 

if 1
    disp('LSPF_OPTIM: starting 2-dim optim')
    disp(fun(s_1, E,x,m))
    s1 = fminsearch(@(ss) fun(ss, E,x,m),s_1);
    s2 = fminsearch(@(ss) fun(ss, E,x,m),s_2);
    % choose the pole-interval with a better initial point
    if fun(s1, E,x,m)>fun(s2, E,x,m)
        s = s2;
    else
        s = s1;
    end
end

% Find poles and interpolation nodes from the refined condenser
% method is the same as in the function fun
[interp,~,~,pols] = dkz09_mapped2(m, -10^s(2), -10^s(1)); 
interp = interp(:);
C = 1./(mp(x)-mp(pols));
C_sub = 1./(mp(interp)-mp(pols));
E_sub = exp(-interp*t);
res = C_sub\mp(E_sub);
R = mp(E) - C*res;
total_err = max(vecnorm(double(R), inf)+((1e-15*ones(1,m))./abs(pols))*abs(res));
% return the total error bound (total_err)

function val = fun(s,E,x,m)
    [intp,~,~,pols] = dkz09_mapped2(m, -10^s(2), -10^s(1)); 
     % poles (pols) and interpolation nodes (intp) from a given condenser
    intp = intp(:);
    C_sub = 1./(mp(intp)-mp(pols));  % sub Cauchy matrix for interpolation
    E_sub = exp(-intp*t);
    coeff = C_sub\mp(E_sub);         % residues
    C = 1./(mp(x)-mp(pols));         % full Cauchy matrix for evaluation
    R = mp(E) - C*coeff;
    val = log10(max(vecnorm(double(R), inf)+((1e-15*ones(1,m))./abs(pols))*abs(coeff)));
    % log10-based total error, here the residual of each shifted linear
    % system is approximated uniformly by 1e-15, as we solve these systems
    % using direct method with standard double precision.
end
end