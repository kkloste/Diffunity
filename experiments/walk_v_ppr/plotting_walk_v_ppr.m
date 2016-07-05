function plotting_walk_v_ppr( varargin )
% Sweep over walk vectors of different lengths and compare conductance
% to rayleigh quotient
%

p = inputParser;
p.addOptional('fname','netscience-cc');
p.parse(varargin{:});

clear; clc; 
fname = p.Results.fname;

load_dir = './results/';
image_dir = './images/';

load( [load_dir, 'walk_v_ppr', fname, '.mat'] );
addpath ../../util;

% MAX_TERMS = n;
% walk_set:  seed, conds, vols, sizes, rayleigh, bound
%		supp_vol, distance1, distanceInf, minf
% ppr_set:   conds, vols, sizes, bound, supp_vol
NUM_SEEDS = length(walk_set.seeds);
for which_seed = 1:NUM_SEEDS,
	clf;
	plot( walk_set.conds(:, which_seed) );
	hold all;
	plot( walk_set.bound(:, which_seed) );

	plot( ppr_set.conds(:, which_seed) );
	hold all;
	plot( ppr_set.bound(:, which_seed) );

	plot( hk_set.conds(:, which_seed) );
	hold all;
	plot( hk_set.bound(:, which_seed) );


	ylim([0,1]);

%	ind = find( walk_set.vols(:,which_seed) >= walk_set.total_volume/2, 1 );
%	if numel(ind) ==0, ind = length(walk_set.vols); end
%	plot( [ind,ind], [0,1] );

%	ind = find( ppr_set.vols(:,which_seed) >= walk_set.total_volume/2, 1 );
%	if numel(ind) ==0, ind = length(ppr_set.vols); end
%	plot( [ind,ind], [0,1] );

	title( sprintf( '%s, seed %d', fname, num2str(which_seed) ) );
	xlabel('Walk term');
	legend('walk-cond', 'walk-bound', 'ppr-cond','ppr-bound', 'hk-cond', 'hk-bound', 'location','Northeast');
	print(gcf,[ image_dir, 'walk-v-ppr-', fname, '-', num2str(which_seed), '.png'],'-dpng');

	fprintf('Done plotting %s  seed %d / %d\n', fname, which_seed, NUM_SEEDS );
end

