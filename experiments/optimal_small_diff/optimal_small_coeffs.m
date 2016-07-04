% Sweep over walk vectors of different lengths and compare conductance
% to rayleigh quotient
%

clear; clc;
addpath ../../util;
addpath ../../diffusion_codes;

save_dir = './results/';
fname = 'netscience-cc';
load(['../../data/', fname, '.mat']);

A = A|A';
A = A-diag(diag(A));
n = size(A,1);
d = full(sum(A,1))';
D = spdiags(d,0,n,n);
Dinv = spdiags(1./d,0,n,n); 
P = colnormout(A);
L = D - A;
nL = Dinv*(L*Dinv);
stdist = ones(n,1)./n;

NUM_SEEDS = 10;
ParamSize = [1000,1000];

trials = struct;
dummy = randn(n,1);
[~,perm] = sort(dummy);
trials.seeds = perm(1:10); % 10 unique seeds
trials.conds = zeros(NUM_SEEDS, ParamSize );
trials.bound = zeros(NUM_SEEDS, ParamSize );
trials.bsetvol = zeros(NUM_SEEDS, ParamSize );
trials.ParamSize = ParamSize;
trials.coeff1 = linspace(0,1,ParamSize(1));
trials.coeff2 = linspace(0,1,ParamSize(2));

for which_seed = 1:NUM_SEEDS,
	seed = trials.seeds(which_seed);
	s = sparse( seed, 1, 1, n, 1 );
%	[bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,s);

	Pvecs = [ s, P*s, P*(P*s) ];

	for k=1:ParamSize(1),
		
		s = P*s;
		supp = find(s);

		dinvs = Dinv*s;

		[bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,s);
		walk_set.conds(k,which_seed) = bestcond;
		walk_set.vols(k,which_seed) = bestvol;
		walk_set.sizes(k,which_seed) = nnz(bestset);
		walk_set.rayleigh(k,which_seed) =  s'*(nL*s) / ( s'* (dinvs) ) ;
		
		c = full(min(dinvs))^2*walk_set.total_volume;
		b = c/( dinvs'*s );
		a = (1 - b);

		walk_set.bound(k,which_seed) =  sqrt( (2*s'*(nL*s) / ( s'* (dinvs) ))/a  ) ;
		walk_set.supp_vol(k,which_seed) = full(sum( d(supp) ) );
		walk_set.distanceInf(k,which_seed) = norm( Dinv*(s - stdist), 'inf' );
		walk_set.distance1(k,which_seed) = norm( Dinv*(s - stdist), 1 );
		walk_set.minf(k,which_seed) = min(dinvs);
	end
	fprintf('Done with %s  seed %d / %d\n', fname, which_seed, NUM_SEEDS );
end
save( [ save_dir, 'optimal_walk_', fname, '.mat'], 'walk_set' );


