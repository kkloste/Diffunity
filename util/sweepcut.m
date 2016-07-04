function [bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,vector,varargin)
% [bestset,bestcond,bestcut,bestvol,noderank] = sweepcut(A,vector,debugflag)
%
%	INPUTS:
%		A       	-- symmetric, unweighted adjacency matrix
%		vector  	-- vector of node scores to perform sweep on. Need not be length n.
%               	If vector is length k < n, then only those k nodes will be swept.
%
%	OPTIONAL INPUTS:
%		varargin	-- 'debugflag':		set to true to have MEX code output errors messages.
%					-- 'halfvol':		set to true to have MEX code stop sweep process once
%						half the volume of the graph is reached.
%					-- 'inputnodes':	set to true to have input "vector" be an
%						ordering on the nodes instead of a vector of values to be sorted.
%						This option does not have safety checks, so only use if you know
%						exactly what you are doing.
%
%	OUTPUTS:
%       bestset     -- set of nodes obtaining best conductance found during sweep.
%       bestcond    -- conductance obtained (similarly for outputs bestcut, bestvol).
%       noderank    -- the ordering of nodes used during the sweep procedure.
%
% Kyle Kloster
% NCSU, July 2016

p = inputParser;
p.addOptional('debugflag',false,@islogical);
p.addOptional('halfvol',false,@islogical);
p.addOptional('inputnodes',false,@islogical);
p.parse(varargin{:});

debugflag = p.Results.debugflag;
halfvol = 0;
if p.Results.halfvol, halfvol = 1; end

if p.Results.inputnodes,
	noderank = vector;
else,
	if issparse(vector) == true, 
		[rows, cols, vals] = find(vector);
		[~, perm] = sort( vals, 'descend');
		noderank = rows(perm);
	else
		[ vals, noderank ] = sort( vector, 'descend' );
		noderank = noderank( find(vals) );
	end
end

[bestindex,bestcond,bestcut,bestvol] = sweepcut_mex(A,noderank,halfvol,debugflag);

bestset = [];
if bestindex >= 1, bestset = noderank(1:bestindex); end


