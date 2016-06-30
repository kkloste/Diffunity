function [bestset,bestcond,bestcut,bestvol] = gendiff_grow(A,vert,varargin)
% GENDIFF_GROW Grow a cluster around a vertex using user-defined diffusion algorithm
%
% [bestset,cond,cut,vol,y,npushes] = gendiff_mex(A,seed_set,coeffs,eps,debugflag)
%
% Computes q(P)*v where
% v is a sparse, nonnegative seed vector, input as "vert",
% and q(x) is some polynomial whose coefficients are input as "coeffs".
% The algorithm uses various values of accuracy, eps,
%
% ... gendiff_grow(A,verts,'key',value,'key',value) specifies optional argument
%
%    'neighborhood' : [false | true] to use the neighborhood of the given
%    vertex as the seed. The default is false.
%
%    'debug' : [false | true] to enable debugging info. The default is
%    false.
%
%
%
% Kyle Kloster
% Purdue University, 2016

p = inputParser;
p.addOptional('debug',false,@islogical);
p.addOptional('neighborhood',false,@islogical);
p.parse(varargin{:});

debugflag = p.Results.debug;

t_vals = [10 20 40 80];
eps_vals = [1e-4 1e-3 5*1e-3 1e-2];


if p.Results.neighborhood
    neighs = find(A(:,vert));
    vert = union(vert,neighs);
end

bestcond = Inf;
bestset = [];
for ei=1:numel(t_vals)

    if debugflag==1, fprintf('gendiff_grow.m: Called gendiff_mex on set of size=%i with t=%f  ;  eps=%f \n', ...
            numel(vert), t_vals(ei), eps_vals(ei)); end

    [curset cond cut vol hk npushes] = gendiff_mex(A, vert, t_vals(ei), eps_vals(ei), debugflag);

    if debugflag==1, fprintf('gendiff_grow.m: gendiff_mex done on set of size=%i with t=%f  ;  eps=%f \n', ...
            numel(vert), t_vals(ei), eps_vals(ei)); end

    if cond < bestcond
        bestcond = cond;
        bestset = curset;
        bestvol = vol;
        bestcut = cut;
    end
end