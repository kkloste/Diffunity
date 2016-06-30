% Testing sweepcut

load ../data/netscience-cc.mat;
addpath ../diffusion_codes;
addpath ../util;

n = size(A,1);

for vert = 1:n,
	[bestset,cond,cut,vol,prvec] = pprgrow_mex(A,vert, nnz(A)/6 ,0.85);
	[bset, cond1, cut1, vol1] = sweepcut(A,prvec);
	assert( numel(setdiff(bestset,bset)) == 0 , fprintf( 'Failure on $d,  cond_pr = %f  v  cond1 = %f 'n\, vert, cond, cond1 ) );
end
