function [diffusion,bestset,conductance,cut,volume,num_pushes] = gendiff_mex1(A,seed_set,coefficients,varargin)
% GENDIFF_MEX1 outputs cluster and diffusion vector computed from 'coefficients' at vert
%
% [diffusion,bestset,conductance,cut,volume,num_pushes] = gendiff_mex1(A,seed_set,coefficients,varargin)
%
% Computes q(P)*v where
% v is a sparse, nonnegative seed vector, input as "vert",
% and q(x) is some polynomial whose coefficients are input as "coeffs".
% "coeffs" should be nonnegative and have sum S > 0
%
% ... gendiff_grow(A,verts,'key',value,'key',value) specifies optional argument
%
%    'acccuracy' : > 0, degree-normalized infinity norm error. Default is 1e-4.
%    'debug' 		 : [false | true] to enable debugging info. The default is
%    false.
%
%
%
% Kyle Kloster
% Purdue University, 2016

p = inputParser;
p.addOptional('debug',false,@islogical);
p.addOptional('accuracy',1e-4,@isscalar);
p.parse(varargin{:});

debugflag = p.Results.debug;
accuracy = p.Results.accuracy;
assert( accuracy > 0 ), sprintf('Input error: accuracy must satisfy accuracy > 0;  accuracy = %f ', accuracy ) );
for j=1:length(coefficients),
	assert(   coefficients(j) >= 0 , 'Input error: coefficients must satisfy coefficients[j] >= 0 ' );
end

coeff_sum = sum(coefficients);
coefficients = coefficients./coeff_sum;

[bestset conductance cut volume diff_vec num_pushes] = gendiff_mex(A, seed_set, coefficients, accuracy, debugflag);
inds = diff_vec(:,1);
vals = diff_vec(:,2);
diffusion = sparse(inds, 1, vals, size(A,1), 1);


end
