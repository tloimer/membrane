function [minterval pinterval] = findinterval(shoot,mguess,pzero);
%FINDINTERVAL Find an interval for the mass flux in which pres changes sign.
%  FINDINTERVAL(FUNCTION,MGUESS,PZERO) returns an interval for M  in which the
%  function FUNCTION changes sign. MGUESS is the first guess, PZERO is the
%  function value for M = 0. PZERO must be smaller than zero. For the GLOBAL
%  variable VERBOSE > 0, diagnostic output is printed.

global VERBOSE

% Initialize
mold = 0; pold = pzero;
mnow = mguess;
pnow = shoot(mguess);
fcount = 1;

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

while pnow < 0
% Here, 2, we double the distance.
  mnew = mnow - 2*pnow*(mnow-mold)/(pnow-pold);
  mold = mnow; pold = pnow;
  mnow = mnew;
  pnow = shoot(mnow);
  fcount = fcount + 1;
  if fcount > 50, error([upper(mfilename) ': More than 50 iterations!']); end
end

minterval = [mold mnow];
pinterval = [pold pnow];

if VERBOSE > 0
  fprintf([upper(mfilename) ': Found an interval, %.3g - %.3g kg/m2s.\n'],...
    mold,mnow);
  fprintf([upper(mfilename) ': %u function calls.\n'],fcount);
end
