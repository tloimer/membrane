function [m,f2ph1ph] = mflow(deltap,T0,zrange)
%MFLOW(DELTAP,T0,ZRANGE) calculates the mass flow rate m for a given
%  pressure difference.
%
%  [M FLOWSTRUCT] = MFLOW(DELTAP,T0,ZRANGE) returns m and the solution
%  flowstruct.

% input check
p0=ps(T0);
if (deltap>p0)
  msg = sprintf('deltap must be smaller than ps(T0) = %g',p0);
  error(msg);
  return
end

p1 = ps(T0) - deltap;
L = abs(zrange(2)-zrange(1));
% mguess basierend on Darcy für Dampf;
%mguess = kappa*deltap./(L*nu(T0,1));
% mguess nach linearer Theorie
mguess = mlambda(T0,deltap,L);

% find an interval for m
% initially, mguess should be too small
pres=m2ph(mguess,T0,p1,zrange);
if (pres>0)
  fac=1.2; % i carefully extend the range!
  mult=1;
else
  fac=1/1.2;
  mult=-1;
end
mb = mguess;
while (mult*pres>0)
  mb = fac*mb;
  pres = m2ph(mb,T0,p1,zrange);
end

% Display: none iter final notify
options =optimset(optimset('fzero'),'TolX',mguess/1000,'Display','none');
m = fzero(@m2ph,[mguess mb],options,T0,p1,zrange);
f2ph1ph=flow2ph1ph(m,T0,1,0,zrange);


%-----------------------------------------------------------------------
function pres = m2ph(m,T0,p1,zrange)
f2ph1ph=flow2ph1ph(m,T0,1,0,zrange);
pres = f2ph1ph(2).p(end) - p1;
% debug:
% sprintf('m: %g,   pres: %g',m,pres)
