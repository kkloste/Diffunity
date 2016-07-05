function walk_v_ppr( varargin )
% Sweep over walk vectors of different lengths and compare conductance
% to rayleigh quotient -- then compare the performance vs the conductance
% obtained for personalized PageRank approximations using the same number
% of terms.
%

p = inputParser;
p.addOptional('fname','netscience-cc');
p.addOptional('alpha',0.9,@isnumeric);
p.addOptional('hk_t',1,@isnumeric);
p.parse(varargin{:});

hk_t = p.Results.hk_t;
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

hk_set = struct;
hk_set.conds = zeros(MAX_TERMS, NUM_SEEDS);
hk_set.vols = zeros(MAX_TERMS, NUM_SEEDS);
hk_set.sizes = zeros(MAX_TERMS, NUM_SEEDS);
hk_set.bound = zeros(MAX_TERMS, NUM_SEEDS);
hk_set.supp_vol = zeros(MAX_TERMS, NUM_SEEDS);
hk_set.minf = zeros(MAX_TERMS, NUM_SEEDS);

for which_seed = 1:NUM_SEEDS,
	seed = walk_set.seeds(which_seed);
	s = sparse( seed, 1, 1, n, 1 );
	[bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,s);

	temp_ppr = s;
	temp_hk = s;
	sa = s;
	sk = s;

	hk_coeff = 1;
	t = hk_t;

	for k=1:MAX_TERMS,
		s = P*s;
		sa = (P*sa).*alpha;
		sk = (P*sk).*(t/hk_coeff);
		hk_coeff = hk_coeff + 1;
		temp_hk = temp_hk + sk;
		temp_ppr = temp_ppr + sa;

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
		ppr = temp_ppr./sum(temp_ppr);
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

		% now do HK
		hk = temp_ppr./sum(temp_hk);
		supp = find(hk);
		dinvhk = Dinv*hk;

		[bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,hk,'halfvol',false);
		hk_set.conds(k,which_seed) = bestcond;
		hk_set.vols(k,which_seed) = bestvol;
		hk_set.sizes(k,which_seed) = nnz(bestset);
		
		c = full(min(dinvhk))^2*walk_set.total_volume;
		b = c/( dinvhk'*hk );
		a = (1 - b);

		hk_set.bound(k,which_seed) =  sqrt( (2*hk'*(nL*hk) / ( hk'* (dinvhk) ))/a  ) ;
		hk_set.supp_vol(k,which_seed) = full(sum( d(supp) ) );

	end
	fprintf('Done with %s  seed %d / %d\n', fname, which_seed, NUM_SEEDS );
end
save( [ save_dir, 'walk_v_ppr', fname, '.mat'], 'walk_set', 'ppr_set', 'hk_set' );


