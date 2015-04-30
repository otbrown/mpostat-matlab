% TenProd.m
% A function which calculates the tensor product of any number of input
% vectors. Must have at least two inputs. Calculates from left to right.
% Oliver Thomson Brown
% 19/06/2014

function product = TenProd(varargin)

narginchk(2,Inf);
product = 1;
for i = 1 : nargin
    product = kron(product, varargin{i});
end
