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
[~,~,~,~,vec] = pprgrow_mex(A,j,1/target_eps, alpha);
[X,Y,T,AUC] = perfcurve(labels,vec,1);
AUC

% [bestset,cond,cut,vol,y,npushes] = hkgrow_mex(A,set,t,eps,debugflag)
[~,~,~,~,vec] = hkgrow_mex(A,j,hk_t,target_eps);
[X,Y,T,AUC] = perfcurve(labels,vec,1);
AUC
