% FIEDLER CUT
%
% Display Fiedler cut on usps

clc;
addpath ../../util; % for set_figure_size
load ../../data/usps3nn.mat ;
addpath ../../diffusion_codes ; %for sweep cut

A = A|A';
A = A-diag(diag(A));
n = size(A,1);
d = full(sum(A,1))';
D = spdiags(d,0,n,n);
Dsqinv = spdiags(d.^(-1/2),0,n,n); 
P = colnormout(A);
L = D - A;
nL = Dsqinv*(L*Dsqinv);

[V, lam] = eigs(nL,2,'SM');
lam = diag(lam);

[min_lam, ind_min]  = max(lam); % because 0 is min

fiedler = V(:,ind_min);

fprintf( ' min_lam = %f ,   f^T*nL*f = %f \n' , min_lam, fiedler'*nL*fiedler );

[bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,fiedler,'halfvol',true) ;

save( './results/fiedler_usps.mat', 'bestset', 'fiedler', 'bestcond' );

