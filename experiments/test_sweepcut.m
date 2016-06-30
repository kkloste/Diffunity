% Testing sweepcut

load ../data/netscience-cc.mat;
addpath ../diffusion_codes;
addpath ../util;

n = size(A,1);
degrees = full(sum(A,2));

for vert = 1:n,
	[bestset,condu,cut,vol,prvec] = pprgrow_mex(A,vert, nnz(A)/6 ,0.85);
	densepr = zeros(n,1);
	densepr( prvec(:,1) ) = prvec(:,2)./degrees( prvec(:,1) ) ;
	[bset, cond1, cut1, vol1] = sweepcut(A,densepr);
	[cond2 cut2 vol2] = cut_cond(A,bset);
%	assert( numel(setdiff(bestset,bset)) == 0 , fprintf( 'Failure on %d,  cond_pr = %f  v  cond1 = %f \n', vert, condu, cond1 ) );
	fprintf( ' %f   \t %f  \t %f \n', cond1-condu, cut1-cut, vol1-vol );
	fprintf( ' %f   \t %f  \t %f \n', cond1-cond2, cut1-cut1, vol1-vol1 );
	fprintf( ' setdiff size %d \n', numel( setdiff( bset, bestset) ) );

end
