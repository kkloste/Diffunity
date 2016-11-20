function two_term_diffusions( varargin )
% AUC_of_diffusions( varargin )
%i
% p.addOptional('fname','senate');

clf;

p = inputParser;
p.addOptional('fname','senate');
p.parse(varargin{:});

fname = p.Results.fname;

save_dir = './results/';
load_dir = '../../data/';
image_dir = './images/';
addpath ../../util;
addpath ../../diffusion_codes

load( [load_dir, fname] );

n = size(A,1);
NUM_COMM = size(C,2);
which_comm = randi(NUM_COMM);
labels =  C(:,which_comm);
community = find(labels);
vert = community( randi(length(community) ) ) ;
seed_vec = sparse(vert,1,1,n,1);



% [1]  Setup parameters, variables

NUM_TERMS = 2;
NUM_POINTS = 20;
mesh_points = linspace(0,1,NUM_POINTS);
coeffs = zeros(NUM_TERMS, 1);

% [2]  Get partial Krylov matrix
%
[rows, cols, vals] = find(seed_vec);
temp_vec = sparse_degpower_mex(A,seed_vec,-1);
temp_vec = sparse( temp_vec(:,1), 1, temp_vec(:,2), n, 1);
for which_term = 1:(NUM_TERMS-1),
  temp_vec = A*temp_vec;
  temp_vec = sparse_degpower_mex(A,temp_vec,-1);
  temp_vec = sparse( temp_vec(:,1), 1, temp_vec(:,2), n, 1);
  [rowst,colst,valst] = find(temp_vec);
  rows = [rows; rowst];
  cols = [cols; (which_term+1)*ones(size(rowst))];
  vals = [vals; valst];
end
Krylov_matrix = sparse( rows, cols, vals, n, NUM_TERMS);


fprintf('About to begin computing.\n');

% tracking stats:
% (full vs push) x (coord setting) x (AUC, cond, set size)
% diffusion_stats( which_diff, coeff_1, coeff_2, deg_normalize_flag, which_diff)
diffusion_stats = zeros(2,NUM_POINTS,2,3);
diffusion_coeffs = zeros(NUM_POINTS,2);
% [3]  Perform computations
%
dflag = true;
dflag2 = 1;
for which_point=1:NUM_POINTS,
  coeffs(1) = mesh_points(which_point);
  coeffs(2) = 1 - coeffs(1);
  diffusion_coeffs(which_point, :) = coeffs';
  diff_vec = Krylov_matrix*coeffs;



  fprintf('Size of diff_vec = %d %d \n', size(diff_vec,1), size(diff_vec,2) );

  % which_diff = full diffusion
  which_diff = 1;
  deg_scale = 1;
  [bestset,cond] = sweepcut(A,diff_vec,'halfvol',true, 'debugflag', dflag);
  [X,Y,T,AUC] = perfcurve(labels,diff_vec,1);
  diffusion_stats(which_diff, which_point, deg_scale, :) = [cond, AUC, numel(bestset)];


  fprintf('Done point %d  diff %d  deg %d \n', which_point, which_diff, deg_scale);

  fprintf('Size of diff_vec = %d %d \n', size(diff_vec,1), size(diff_vec,2) );

  fprintf('Is diff_vec sparse = %d \n', issparse(diff_vec) );

  deg_scale = 2;
  diff_vec = sparse_degpower_mex(A,diff_vec,-1.0, dflag2);
  diff_vec = sparse( diff_vec(:,1), 1, diff_vec(:,2), n, 1);
  [bestset,cond] = sweepcut(A,diff_vec,'halfvol',true, 'debugflag', dflag);
  [X,Y,T,AUC] = perfcurve(labels,diff_vec,1);
  diffusion_stats(which_diff, which_point, deg_scale, :) = [cond, AUC, numel(bestset)];


  fprintf('Done point %d  diff %d \n', which_point, which_diff);

  [diff_vec,bestset,~,~,~,npushes] = gendiff_mex1(A,vert,coeffs,1e-3);

  % which_diff = push diffusion
  which_diff = 2;
  deg_scale = 2;
  [bestset,cond] = sweepcut(A,diff_vec,'halfvol',true);
  [X,Y,T,AUC] = perfcurve(labels,diff_vec,1);
  diffusion_stats(which_diff, which_point, deg_scale, :) = [cond, AUC, numel(bestset)];

  deg_scale = 1;
  diff_vec = sparse_degpower_mex(A,diff_vec,1);
  diff_vec = sparse( diff_vec(:,1), 1, diff_vec(:,2), n, 1);
  [bestset,cond] = sweepcut(A,diff_vec,'halfvol',true);
  [X,Y,T,AUC] = perfcurve(labels,diff_vec,1);
  diffusion_stats(which_diff, which_point, deg_scale, :) = [cond, AUC, numel(bestset)];
end

save( [save_dir, fname, '-two_term_diff.mat'], 'diffusions_stats','diffusion_coeffs','which_comm','vert','-v7.3');

fprintf('Done computing, now to print.\n');

  Xs = diffusion_coeffs(:,1);
  linstyles = { '-', ':' };
  colors = [ 1, 0.5, 0.5; 0, 0, 1 ];
  for which_stat =1:2,
    clf
    for which_diff = 1:2,
      for deg_scale = 1:2,
          Ys = squeeze(diffusions_stats(which_diff, :, deg_scale, which_stat));
        	plot( [Xs,Ys], 'Color', colors(which_diff,:), 'LineStyle', linstyles{deg_scale} );
        	hold all;
      end
    end
    ylim([0,1]);
    xlim([0,1]);

    title( sprintf( 'Graph: %s, community: %d , vert: %d', fname, which_comm, vert ) );
    xlabel('coeff 1');
    ylabel( stat_string{which_stat} );
    legend('full', 'full-n', 'push','push-n', 'location','Southeast');
    set_figure_size( [3,3] );
    print(gcf,[ image_dir, 'gen-diff-AUC-', stat_string{which_stat}, fname, '-', num2str(which_comm), '-', num2str(vert), '.png'],'-dpng');
  end



fprintf('Done printing  %s \n', fname);
