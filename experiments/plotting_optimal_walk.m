% Sweep over walk vectors of different lengths and compare conductance
% to rayleigh quotient
%

clear; clc; 
fname = 'netscience-cc';

load_dir = './results/';
image_dir = './images/';

load( [load_dir, 'optimal_walk_', fname, '.mat'] );
addpath ../util;

% MAX_TERMS = n;
% walk_set:  seed, conds, vols, sizes, rayleigh, bound
%		supp_vol, distance1, distanceInf, minf
NUM_SEEDS = length(walk_set.seeds);
for which_seed = 1:NUM_SEEDS,
	clf;
	plot( walk_set.conds(:, which_seed) );
	hold all;
	plot( walk_set.bound(:, which_seed) );
	plot( walk_set.supp_vol(:, which_seed) ./ walk_set.total_volume );
	plot( walk_set.vols(:, which_seed) ./ walk_set.total_volume );
	plot( walk_set.minf(:, which_seed) );
	ylim([0,1]);

	ind = find( walk_set.vols(:,which_seed) >= walk_set.total_volume/2, 1 );
	if numel(ind) ==0, ind = length(walk_set.vols); end
	plot( [ind,ind], [0,1] );

	title( sprintf( '%s, seed %d', fname, num2str(which_seed) ) );
	xlabel('Walk term');
	legend('conductance', 'our bound', 'supp vol', 'vols', 'minf', 'location','Northeast');
	print(gcf,[ image_dir, 'cond-v-bound-', fname, '-', num2str(which_seed), '.png'],'-dpng');

	fprintf('Done plotting %s  seed %d / %d\n', fname, which_seed, NUM_SEEDS );
end

%% PLOT DISTANCE

range = 1:50;
for which_seed = 1:NUM_SEEDS,
	clf;
	plot( walk_set.conds(range, which_seed) );
	hold all;
	plot( walk_set.bound(range, which_seed) );
	plot( log10(walk_set.distance1(range, which_seed) ) );
	plot( log10(walk_set.distanceInf(range, which_seed) ) );

yl = ylim;
	ind = find( walk_set.vols(:,which_seed) >= walk_set.total_volume/2, 1 );
	if numel(ind) ==0, ind = max(range); end
	plot( [ind,ind], yl );


	title( sprintf( '%s, seed %d', fname, num2str(which_seed) ) );
	xlabel('Walk term');
	legend('conductance', 'our bound', 'dist1', 'dist_inf', 'location','Northeast');
	print(gcf,[ image_dir, 'cond-v-dist-', fname, '-', num2str(which_seed), '.png'],'-dpng');

	fprintf('Done plotting %s  seed %d / %d\n', fname, which_seed, NUM_SEEDS );
end
