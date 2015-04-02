function [minterval,pinterval] = findinterval(shoot,mguess,pzero)
%FINDINTERVAL Find an interval in which a function changes sign.
%  FINDINTERVAL(FUNCTION,MGUESS,PZERO) returns the interval of the argument
%  to a function in which the function FUNCTION changes sign. MGUESS is a
%  guess for the zero crossing, where F(M) = 0. PZERO is the function value
%  for M = 0. PZERO must be smaller than zero.
%  [MINTERVAL, PINTERVAL] = FINDINTERVAL(FUNCTION,MGUESS,PZERO) returns the
%  interval of argument MINTERVAL and the function values at the interval
%  boundaries PINTERVAL.
%
%  For the GLOBAL variable VERBOSE > 0, diagnostic output is printed.
%
%  See also FINDZERO.

global VERBOSE

% Initialize
mold = 0; pold = pzero;
mnow = mguess;
overshoot = true;
fcount = 0;
while overshoot
  try
    fcount = fcount + 1;
    pnow = shoot(mnow);
    overshoot = false;
  catch err
    % this happens most probably, when the temperature is overshoot (by too high
    % mass flux). overshoot remains true
    if strcmp(err.identifier,'MATLAB:deval:SolOutsideInterval')
      if VERBOSE > 0
	fprintf('%s: Reduce mass flux. Probably too large initial guess.\n', mfilename);
      end
      mnow = mnow/2;
    else
      rethrow(err);
    end
  end
end

% pres ^
%      |            * new value, we double the distance to the intersection
%      |          +   intersection
%      |        +     old value
%      |        | | |
%    __|___________________________________>
%      |         /
% pnow |        *
%      |       /      straight line equation, for m:
%      |      /          m = mnow + (p-pnow)*(mnow-mold)/(pnow-pold)
%      |     /        p = 0:
% pold |    *            m = mnow - pnow*(mnow-mold)/(pnow-pold)
%      |  mold mnow
%      |
%

% mguess is too small, pnow is negative
while pnow < 0
% Here, 2, we double the distance.
  madd = - 2*pnow*(mnow-mold)/(pnow-pold);
  mold = mnow; pold = pnow;
  mnow = mnow + madd;
  try
    fcount = fcount + 1;
    pnow = shoot(mnow);
  catch err
    if strcmp(err.identifier,'MATLAB:deval:SolOutsideInterval') %...
       %  || strcmp(err.identifier,'MATLAB:badCellRef')
      if VERBOSE > 0
	fprintf('%s: Too large mass flux. Reduce added distance.\n', mfilename);
      end
      toolarge = true;
      while toolarge
	mnow = mnow - madd;
	madd = madd/2;
	mnow = mnow + madd;
	try
	  fcount = fcount + 1;
	  pnow = shoot(mnow);
	  toolarge = false;
	catch err
	  if strcmp(err.identifier,'MATLAB:deval:SolOutsideInterval') %...
	     % || strcmp(err.identifier,'MATLAB:badCellRef')
	    if fcount > 50
	      error('More than 50 iterations, decreasing added mass flux!');
	    end
	    continue;
	  else
	    rethrow(err);
	  end
	end % end try-catch
      end % end while
    else
      rethrow(err);
    end
  end % end try-catch
  if fcount > 50, error('More than 50 iterations!'); end
end

minterval = [mold mnow];
pinterval = [pold pnow];

% Now the interval could be [0 mguess]; but we can not run with m = 0.
% Look for an interval [m>0 mguess].
% Simply intersect a straight line between (0,pold) and (mnow,pnow) with
% pres = 0; presumably, pres = shoot(m) is concave
if mold == 0 % if ~mold
  if VERBOSE > 0
    fprintf('%s: Looking for m > 0, after %d iterations.\n', mfilename,fcount);
  end
  % pnow is used as plow, henceforth
  while pnow > 0
    mnow = mnow - pnow*mnow/(pnow-pold);
    pnow = shoot(mnow);
    % only once provide for a convex function shoot(m)
    if pnow > 0
      mnow = mnow/2;
      pnow = shoot(mnow);
    end
    fcount = fcount + 1;
    if fcount > 50, error('More than 50 iterations!'); end
  end
  minterval(1) = mnow;
  pinterval(1) = pnow;
end

if VERBOSE > 0
  fprintf('%s: Found an interval, %.3g - %.3g g/m2s.\n', mfilename,minterval*1e3);
  fprintf('%s: %u function calls.\n', mfilename,fcount);

  if VERBOSE > 2
    % Plot pressure-residuum versus mass flux over the found mass flux interval
    mm = [minterval(1):(minterval(2)-minterval(1))/11:minterval(2) minterval(2)];
    pp = mm;
    for i = 1:length(mm)
      pp(i) = shoot(mm(i));
    end
    fprintf('%12.6g g/m2s,  %12.6g Pa\n',[mm*1e3; pp]);
    figure('Name',upper(mfilename));
    plot(mm,pp,'k*');
  end
end
