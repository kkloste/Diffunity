% Testing sweepcut_mex
% by comparing the output of sweepcut_mex to the output of pprgrow_mex.
% Check that both codes identify the same set of nodes as the best cut,
% and that the two codes agree on conductance, cut, and volume.

load ../data/netscience-cc.mat;

A = spones(A);
assert( norm(A-A',1) == 0 , 'Not symmetric\n.');
A = A - diag(diag(A));

addpath ../diffusion_codes;
addpath ../util;

n = size(A,1);
degrees = full(sum(A,2));

for vert = 1:n,
	[bestset,condu,cut,vol,prvec] = pprgrow_mex(A,vert, nnz(A)/6 ,0.85);
	prvec = sparse( prvec(:,1), 1, prvec(:,2), n, 1);
	[bset, cond1, cut1, vol1, noderank] = sweepcut(A,prvec);
	[cond2 cut2 vol2] = cut_cond(A,bset);
	[cond3 cut3 vol3] = cut_cond(A,bestset);

	jsim = jaccard_similarity( bestset, bset );

	if jsim ~= 1 , sprintf( 'error node = %d \t jac sim = %f ', vert, jsim ) ; end
	if condu - cond1 ~= 0, sprintf( 'error node = %d \t pprgrow - sweepcut = %f .', vert, condu-cond1 ) ; end
	if condu - cond3 ~= 0, sprintf( 'error node = %d \t pprgrow - cut_cond = %f .', vert, condu-cond3 ); end
	if cond1 - cond2 ~= 0, sprintf( 'error node = %d \t sweepcut - cut_cond = %f .', vert, cond1-cond2 ); end

end

fprintf( '\n test_sweepcut done on netscience-cc -- everything worked, unless output above said otherwise. \n.' );
