function [bestset,bestcond,bestcut,bestvol] = sweepcut(A,vector,debugflag)
% [bestset,bestcond,bestcut,bestvol] = sweepcut(A,noderank)
%
% Kyle Kloster
% NCSU, 2016

if nargin < 3, debugflag = 0; end

[ ~, noderank ] = sort( vector, 'descend' );
noderank = noderank( find(noderank) );

% TO DO: offer variable argument to allow degree normalizing

[bestindex,bestcond,bestcut,bestvol] = sweepcut_mex(A,noderank,debugflag);

bestset = [];
if bestindex >= 1, bestset = noderank(1:bestindex); end


