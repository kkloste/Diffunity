function two_term_diffusions( varargin )
% AUC_of_diffusions( varargin )
%i
% p.addOptional('fname','senate');

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
COEFF_MAX = 0.5;
mesh_points = linspace(0,COEFF_MAX,NUM_POINTS);
coeffs = zeros(NUM_TERMS+1, 1);

% [2]  Get partial Krylov matrix
%
[rows, cols, vals] = find(seed_vec);
temp_vec = sparse_degpower_mex(A,seed_vec,-1);
temp_vec = sparse( temp_vec(:,1), 1, temp_vec(:,2), n, 1);
for which_term = 1:(NUM_TERMS),
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
dflag = false;
dflag2 = 0;
coeffs(1) = COEFF_MAX;
for which_point=1:NUM_POINTS,
  coeffs(2) = mesh_points(which_point);
  coeffs(3) = COEFF_MAX - coeffs(2);
  diffusion_coeffs(which_point, :) = coeffs';
  diff_vec = Krylov_matrix*sparse(coeffs);

  % which_diff = full diffusion
  which_diff = 1;
  deg_scale = 1;
  [bestset,cond] = sweepcut(A,diff_vec,'halfvol',true, 'debugflag', dflag);
  [X,Y,T,AUC] = perfcurve(labels,diff_vec,1);
  diffusion_stats(which_diff, which_point, deg_scale, :) = [cond, AUC, numel(bestset)];

  deg_scale = 2;
  diff_vec = sparse_degpower_mex(A,diff_vec,-1.0, dflag2);
  diff_vec = sparse( diff_vec(:,1), 1, diff_vec(:,2), n, 1);
  [bestset,cond] = sweepcut(A,diff_vec,'halfvol',true, 'debugflag', dflag);
  [X,Y,T,AUC] = perfcurve(labels,diff_vec,1);
  diffusion_stats(which_diff, which_point, deg_scale, :) = [cond, AUC, numel(bestset)];


  [diff_vec,bestset,~,~,~,npushes] = gendiff_mex1(A,vert,coeffs,'accuracy',1e-3);
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

save( [save_dir, fname, '-', num2str(which_comm), '-', num2str(vert), '-two_term_diff.mat'], 'diffusion_stats','diffusion_coeffs','which_comm','vert','-v7.3');

fprintf('Done computing.\n');
