function plotting_walk_v_ppr( varargin )
% Sweep over walk vectors of different lengths and compare conductance
% to rayleigh quotient
%
% p.addOptional('fname','netscience-cc');

p = inputParser;
p.addOptional('fname','netscience-cc');
p.parse(varargin{:});

fname = p.Results.fname;

load_dir = './results/';
image_dir = './images/';

load( [load_dir, 'walk_v_ppr', fname, '.mat'] );
addpath ../../util;

% MAX_TERMS = n;
% walk_set:  seed, conds, vols, sizes, rayleigh, bound
%		supp_vol, distance1, distanceInf, minf
% ppr_set:   conds, vols, sizes, bound, supp_vol
NUM_SEEDS = length(graph.seeds);
for which_seed = 1:NUM_SEEDS,
	clf;

% conductances
	plot( walk_set.conds(:, which_seed), '-r' );
	hold all;
	plot( ppr_set.conds(:, which_seed) , '-b' );
	plot( hk_set.conds(:, which_seed), '-g'  );
	plot( lazy_set.conds(:, which_seed) , '-y' );

	ylim([0,1]);

	ind = find( walk_set.vols(:,which_seed) >= graph.volume/2, 1 );
	if numel(ind) ==0, ind = length(walk_set.vols); end
	plot( [ind,ind], [0,1], ':k' );


	title( sprintf( '%s, seed %d', fname, num2str(which_seed) ) );
	xlabel('Walk term');
	legend('walk-cond',  'ppr-cond','hk-cond', 'lazy-cond', 'half vol', 'location','Northeast');
	print(gcf,[ image_dir, 'walk-v-ppr-', fname, '-', num2str(which_seed), 'cond.png'],'-dpng');




% bounds
	clf;
	plot( walk_set.bound(:, which_seed) , '--r' );
	hold all
	plot( ppr_set.bound(:, which_seed) , '--b' );
	plot( hk_set.bound(:, which_seed) , '--g' );
	plot( lazy_set.bound(:, which_seed),  '--y' );

	ylim([0,1]);

	ind = find( walk_set.vols(:,which_seed) >= graph.volume/2, 1 );
	if numel(ind) ==0, ind = length(walk_set.vols); end
	plot( [ind,ind], [0,1], ':k' );


	title( sprintf( '%s, seed %d', fname, num2str(which_seed) ) );
	xlabel('Walk term');
	legend('walk-bound', 'ppr-bound', 'hk-bound', 'lazy-bound', 'half vol', 'location','Northeast');
	print(gcf,[ image_dir, 'walk-v-ppr-', fname, '-', num2str(which_seed), 'bound.png'],'-dpng');

	fprintf('Done plotting %s  seed %d / %d\n', fname, which_seed, NUM_SEEDS );
end

