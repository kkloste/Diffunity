% [X,Y,T,AUC] = perfcurve(labels,scores,posclass)
clc

load ../../data/senate;
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
power_k = 7;
hk_t = 20;

% RANDOM WALK
vec = reg_power_mex(A,eyej, power_k);
[X,Y,T,AUC] = perfcurve(labels,vec,1);
AUC

vec = reg_power_mex(A,eyej, power_k+1);
[X,Y,T,AUC] = perfcurve(labels,vec,1);
AUC

% [bestset,cond,cut,vol,prvec] = pprgrow_mex(A,j,targetvol,alpha)
% [bestset,cond,cut,vol,y,npushes] = hkgrow_mex(A,set,t,eps,debugflag)
% [bestset,cond,cut,vol,y,npushes] = gendiff_mex(A,seed_set,coeffs,eps,debugflag)

% PAGERANK
accuracy = target_eps;
N_terms = power_k; % 1+floor(  log(accuracy)/log(alpha) - 1 );
coefficients = alpha.^( 0:1:(N_terms-1) );
coefficients = coefficients.*(1-alpha);

[vec, bestset] = gendiff_mex1(A, seed_set, coefficients, accuracy);
vec = full(vec);
[X,Y,T,AUC] = perfcurve(labels,vec,1);
AUC

%HEAT KERNEL
%[~,~,~,~,vec] = hkgrow_mex(A,j,hk_t,target_eps);
%[X,Y,T,AUC] = perfcurve(labels,vec,1);
%AUC
