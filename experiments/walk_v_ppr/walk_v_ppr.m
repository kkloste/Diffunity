function walk_v_ppr( varargin )
% Sweep over walk vectors of different lengths and compare conductance
% to rayleigh quotient -- then compare the performance vs the conductance
% obtained for personalized PageRank approximations using the same number
% of terms.
%

p = inputParser;
p.addOptional('fname','netscience-cc');
p.addOptional('alpha',0.9,@isnumeric);
p.parse(varargin{:});

alpha = p.Results.alpha;

clc;
addpath ../../util;
addpath ../../diffusion_codes;

save_dir = './results/';
fname = p.Results.fname
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

MAX_TERMS = 100;
NUM_SEEDS = 10;

walk_set = struct;
dummy = randn(n,1);
[~,perm] = sort(dummy);
walk_set.seeds = perm(1:10); % 10 unique seeds
walk_set.conds = zeros(MAX_TERMS, NUM_SEEDS);
walk_set.vols = zeros(MAX_TERMS, NUM_SEEDS);
walk_set.sizes = zeros(MAX_TERMS, NUM_SEEDS);
walk_set.rayleigh = zeros(MAX_TERMS, NUM_SEEDS);
walk_set.bound = zeros(MAX_TERMS, NUM_SEEDS);
walk_set.supp_vol = zeros(MAX_TERMS, NUM_SEEDS);
walk_set.distance1 = zeros(MAX_TERMS, NUM_SEEDS);
walk_set.distanceInf = zeros(MAX_TERMS, NUM_SEEDS);
walk_set.minf = zeros(MAX_TERMS, NUM_SEEDS);
walk_set.total_volume = nnz(A);

ppr_set = struct;
ppr_set.conds = zeros(MAX_TERMS, NUM_SEEDS);
ppr_set.vols = zeros(MAX_TERMS, NUM_SEEDS);
ppr_set.sizes = zeros(MAX_TERMS, NUM_SEEDS);
ppr_set.bound = zeros(MAX_TERMS, NUM_SEEDS);
ppr_set.supp_vol = zeros(MAX_TERMS, NUM_SEEDS);
ppr_set.minf = zeros(MAX_TERMS, NUM_SEEDS);

for which_seed = 1:NUM_SEEDS,
	seed = walk_set.seeds(which_seed);
	s = sparse( seed, 1, 1, n, 1 );
	[bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,s);

	ppr = s;
	sa = s;

	for k=1:MAX_TERMS,
		s = P*s;
		sa = (P*sa).*alpha;
		ppr = ppr + sa;

		supp = find(s);
		dinvs = Dinv*s;

		[bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,s,'halfvol',false);
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

		% now do PPR
		supp = find(ppr);
		dinvppr = Dinv*ppr;

		[bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,ppr,'halfvol',false);
		ppr_set.conds(k,which_seed) = bestcond;
		ppr_set.vols(k,which_seed) = bestvol;
		ppr_set.sizes(k,which_seed) = nnz(bestset);
		
		c = full(min(dinvppr))^2*walk_set.total_volume;
		b = c/( dinvppr'*ppr );
		a = (1 - b);

		ppr_set.bound(k,which_seed) =  sqrt( (2*ppr'*(nL*ppr) / ( ppr'* (dinvppr) ))/a  ) ;
		ppr_set.supp_vol(k,which_seed) = full(sum( d(supp) ) );

	end
	fprintf('Done with %s  seed %d / %d\n', fname, which_seed, NUM_SEEDS );
end
save( [ save_dir, 'walk_v_ppr', fname, '.mat'], 'walk_set', 'ppr_set' );


