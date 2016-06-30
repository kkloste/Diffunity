function [x, approx_kappa, gamma, v] = local_fiedler(L, d, S, kappa, lambda2, verbose)
% [x, approx_kappa, gamma, v] = local_fiedler(L, d, S, kappa, lambda2, verbose)
%
% Constructs the "MOV" vector from [Mahoney, Orecchia, Vishnoi 2012],
% which is a locally-biased version of the Fiedler vector.
%
%   L       -  the edge-node incidence matrix
%   d       -  vector of degrees
%   S       -  seed set,
%   kappa   -  correlation parameter, 0 < kappa < 1
%
%   OPTIONAL:
%   lambda2 -  the smallest positive generalized eigenvalue of the Laplacian
%               (This is the Fiedler eigenvalue)
%   verbose -  set to 1 to output gamma and x'*L*x at each step.
%
%   If lambda2 is omitted, this function will compute it automatically via
%   MATLAB's built in eigs command. If this function is going to be called
%   multiple times, some cost can be saved by pre-computing lambda2 a
%   single time and passing it in as a paramter. Compute via
%     [lams] = eigs(L, spdiags(d, 0, n, n), 2, 'SA');
%     lambda2 = lams(2);
%   where L = D-A, n = size(L,1).

% SETUP PARAMETERS
n = size(L,2);
volG = sum(d);
volS = sum(d(S));
volSbar = volG - volS;
temp_scalar = sqrt( volS*volSbar/volG );
v = -(temp_scalar/volSbar).*ones(n,1);
v(S) = temp_scalar/volS;
% dsq = d.^(1/2);
Dsq = spdiags( d.^(1/2), 0, n, n) ;

vd = v.*d;
sqkappa = sqrt(kappa);

if nargin < 6, verbose = 0; end
if nargin < 5,
    [lams] = eigs(L, spdiags(d, 0, n, n), 2, 'SA');
    lambda2 = lams(2);
end
gamma_left = -volG;
gamma_right = lambda2;
gamma_cur = (gamma_left + gamma_right)/2;
gamma = [];
kappa_tol = 1e-2;
approx_kappa = -1;

% BINARY SEARCH FOR GAMMA
while ( abs(approx_kappa - sqkappa) > kappa_tol || approx_kappa < sqkappa ),
    gamma = gamma_cur;
    x = MOV_for_gamma(L, gamma, d, vd);

    approx_kappa = x'*vd;
    if approx_kappa > sqkappa,
        gamma_left = gamma_cur;
    else
        gamma_right = gamma_cur;
    end
    gamma_cur = (gamma_left + gamma_right)/2;
    if verbose==1,
        fprintf('\ngamma_cur=%f   x^TLx=%f',gamma_cur, x'*L*x);
    end
end
    
    
function y = MOV_for_gamma(L, gamma, d, vd)
n = size(L,1);
y = (L - spdiags( gamma.*d, 0, n, n) )\vd;
y = y/sqrt((d'*(y.*y))); % normalize
