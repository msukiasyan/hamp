% BINSEARCHX Computes the zero of univariate function on interval via Binary Search
%   Vectorized version of Binary Search
% USAGE
%   [x,fval] = binsearchx(f,a,b,P1,P2,...);
% INPUTS
%   f         : name of function of form fval=f(x)
%   a,b       : left, right endpoints of interval
%   P1,P2,... : optional additional arguments for f
% OUTPUTS
%   x       : zero of f
%
% Note: Assumes the function is increasing over the interval.

function x = binsearchx(f, a, b, varargin)

tol = 1e-11;

while any(b - a > tol)
      
  x     = (a + b) / 2;
  fx    = feval(f, x, varargin{:});
  i     = fx > 0;
  a     = a .* i + x .* (~i);
  b     = x .* i + b .* (~i);

end

x = (a + b) / 2;
