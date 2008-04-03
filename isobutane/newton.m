function x = newton(fun,x0,y,res,iter)
%NEWTON     Newton iteration.
%
%  NEWTON(FUN,X0) finds a solution to F(X) = 0 by Newton iteration, starting
%  from X = X0. NEWTON expects FUN to be the handle to a function that returns a
%  vector [F(X) DF(X)], where DF is the derivative of F at X.
%
%  NEWTON(FUN,X0,Y,RES,ITER) finds a solution to F(X) = Y. NEWTON iterates until
%  ABS( F(X) - Y ) < RES or more than ITER iterations are done. Default values
%  are RES = 1e-6 and ITER = 1000.
% 
%  The exact X is at a distance of approximately RES/DF from the value returned.
%
%  Example
%    Given
%      function [s ds] = senus(x)
%      s = sin(x); ds = cos(x);
%    should give pi, using a function handle to senus.
%      probablypi = newton(@senus,3.1,0)
%
%  See also the MATLAB-functions function_handle, feval.

if nargin < 5 iter = 1000; end
if nargin < 4 res = 1e-6; end
if nargin < 3 y = 0; end

x = x0;
for i = 1:iter
  [y0 dy] = feval(fun,x);
  if abs(y0-y) < res return
  else
    x = x - (y0-y)/dy;
  end
end

func2str(fun)
x0
y
x
y0
res
iter
error('NEWTON not succesful. Increase RES or ITER. Type help newton.')
