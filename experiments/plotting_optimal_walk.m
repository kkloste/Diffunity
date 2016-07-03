% Sweep over walk vectors of different lengths and compare conductance
% to rayleigh quotient
%

clear; clc; clf;


fname = 'netscience-cc';

load_dir = './results/';
image_dir = './images/';

load( [load_dir, 'optimal_walk_', fname, '.mat'] );
addpath ../util;

% MAX_TERMS = n;
% walk_set:  seed, conds, vols, sizes, rayleigh, bound
%		supp_vol, distance1, distanceInf, minf
which_seed = 1:length(walk_set.seeds),
	plot( walk_set.conds(:, which_seed) );
	hold on;
	plot( walk_set.bound(:, which_seed) );
	ylim([0,1]);


	title( sprintf( '%s, seed %d', fname, num2str(which_seed) ) );
	xlabel('Walk term');
	legend('conductance', 'our bound', 'location','Northeast');
	print(gcf,[ image_dir, 'cond-v-bound-', fname, '-', num2str(which_seed), '.png'],'-dpng');
end
