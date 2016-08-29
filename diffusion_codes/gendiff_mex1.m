function [diff_vec,bestset,conductance,cut,volume,num_pushes] = gendiff_mex1(A,seed_set,coefficients,accuracy,varargin)
% GENDIFF_MEX1 outputs cluster and diffusion vector computed from 'coefficients' at vert
%
% [diff_vec,bestset,cond,cut,vol,npushes] = gendiff_mex1(A,seed_set,coefficients,accuracy,varargin)
%
% Computes q(P)*v where
% v is a sparse, nonnegative seed vector, input as "vert",
% and q(x) is some polynomial whose coefficients are input as "coeffs".
% "coeffs" should be nonnegative and have sum S satisfying  (1-accuracy) < S <= 1
%
% ... gendiff_grow(A,verts,'key',value,'key',value) specifies optional argument
%
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
p.parse(varargin{:});

debugflag = p.Results.debug;

assert( (1.0 - accuracy < sum(coefficients)), ...
 "\nInput error: accuracy and coefficients must satisfy (1 - accuracy < sum(coefficients) )\n" );
assert( (sum(coefficients) <= 1) , "\nInput error: coefficients must satisfy sum(coefficients) <= 1 \n" );
for j=1:length(coefficients),
	assert(   coefficients(j) >= 0 , "\nInput error: coefficients must satisfy coefficients[j] >= 0 \n" );
end

[bestset conductance cut volume diff_vec num_pushes] = gendiff_mex(A, seed_set, coeffs, accuracy, debugflag);
inds = diff_vec(:,1);
vals = diff_vec(:,2);
diffusion = sparse(inds, 1, vals, size(A,1), 1);


end
