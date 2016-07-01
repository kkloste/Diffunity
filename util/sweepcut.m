function [bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,vector,debugflag)
% [bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,vector,debugflag)
%
%	A       -- symmetric, unweighted adjacency matrix
%   vector  -- vector of node scores to perform sweep on. Need not be length n.
%               If vector is length k < n, then only those k nodes will be swept.
%
%   OUTPUTS
%       bestset     -- set of nodes obtaining best conductance found during sweep.
%       bestcond    -- conductance obtained. (bestcut, bestvol same).
%       noderank    -- the ordering of nodes used during the sweep procedure.
%
% Kyle Kloster
% NCSU, 2016

if nargin < 3, debugflag = 0; end

if issparse(vector) == true, 
	[rows, cols, vals] = find(vector);
	[~, perm] = sort( vals, 'descend');
	noderank = rows(perm);
else
	[ vals, noderank ] = sort( vector, 'descend' );
	noderank = noderank( find(vals) );
end


% TO DO: offer variable argument to allow degree normalizing

[bestindex,bestcond,bestcut,bestvol] = sweepcut_mex(A,noderank,debugflag);

bestset = [];
if bestindex >= 1, bestset = noderank(1:bestindex); end


