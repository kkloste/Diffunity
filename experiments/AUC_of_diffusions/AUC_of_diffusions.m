function AUC_of_diffusions( varargin )
% AUC_of_diffusions( varargin )
%i
% p.addOptional('fname','netscience-cc');

clf;

p = inputParser;
p.addOptional('fname','senate');
p.parse(varargin{:});

fname = p.Results.fname;

load_dir = '../../data/';
image_dir = './images/';
addpath ../../util;

load( [load_dir, fname] );

n = size(A,1);
NUM_COMM = size(C,2);
which_comm = randi(NUM_COMM);
labels =  C(:,which_comm);
community = find(labels);

% Get a random member of that community.
j = community( randi(length(community) ) ) ;
eyej = sparse(j,1,1,n,1);
seed_set = j;

% Prep diffusions, set parameters
addpath ../../diffusion_codes
target_eps = 1e-4;
alpha = 0.99;
hk_t =5;
accuracy = target_eps;
debugflag = false;
AUCs = zeros(30,3);

fprintf( 'Inputs and parameters set, about to start computation\n');

for N_terms=1:30,
	which_alg = 0;
	% RANDOM WALK
	which_alg = which_alg+1;
	vec = reg_power_mex(A,eyej, N_terms);
	[X,Y,T,AUC] = perfcurve(labels,vec,1);
	AUCs(N_terms, which_alg) = AUC;

% [bestset,cond,cut,vol,prvec] = pprgrow_mex(A,j,targetvol,alpha)
% [bestset,cond,cut,vol,y,npushes] = hkgrow_mex(A,set,t,eps,debugflag)
% [bestset,cond,cut,vol,y,npushes] = gendiff_mex(A,seed_set,coeffs,eps,debugflag)

	% PAGERANK
	which_alg = which_alg+1;
	coefficients = alpha.^( 0:1:N_terms );
	
	[vec, bestset] = gendiff_mex1(A, seed_set, coefficients, accuracy, debugflag);
	vec = full(vec);
	[X,Y,T,AUC] = perfcurve(labels,vec,1);	
	AUCs(N_terms, which_alg) = AUC;

	%HEAT KERNEL
	which_alg = which_alg+1;
	coefficients = ones( N_terms+1, 1 );
	for j=1:N_terms, coefficients(j+1) = hk_t*coefficients(j)/j; end

	[vec, bestset] = gendiff_mex1(A, seed_set, coefficients, accuracy, debugflag);
	vec = full(vec);
	[X,Y,T,AUC] = perfcurve(labels,vec,1);
	AUCs(N_terms, which_alg) = AUC;
end

fprintf('Done computing, now to print.\n');

	plot( AUCs(:, 1), '-k' );
	hold all;
	plot( AUCs(:, 2), '-b' );
	plot( AUCs(:, 3), '-r'  );

	ylim([0,1]);
	xlim([1,30]);

	title( sprintf( 'Graph: %s, community: %d', fname, which_comm ) );
	xlabel('Number of walk terms');
	legend('walk',  'ppr','hk', 'location','Southeast');

	set_figure_size( [3,3] );
	print(gcf,[ image_dir, 'diff-AUC-', fname, '-', num2str(which_comm), '.png'],'-dpng');
fprintf('Done printing  %s \n', fname);
