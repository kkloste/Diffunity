% FIEDLER CUT
%
% Display Fiedler cut on usps

clc;
addpath ../util; % for set_figure_size
load ../data/usps3nn.mat ;
addpath ../diffusion_codes ; %for sweep cut

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


spy(A,'k');
title('sparsity pattern, USPS');
print( './images/usps3nn_spy','-dpng','-r600','-painters');

for j=1:length(bestset),
	vert = bestset(j);
	neighb = A(:,bestset);
	neighb_in = intersect(neighb,bestset);
	scatter( neighb_in, vert*ones(length(neighb_in) ), 'b' );
	scatter( vert*ones(length(neighb_in) ), neighb_in, 'b' );

	neighb_out = setdiff(neighb, neighb_in);
	scatter( neighb_out, vert*ones(length(neighb_out) ), 'r' );
	scatter( vert*ones(length(neighb_out) ), neighb_out, 'r' );

end

title('fiedler cut, USPS')
print( './images/usps3nn_fiedler','-dpng','-r600','-painters');


