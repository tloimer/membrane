function flowstruct = mback(deltap,T0,L)
%MBACK      Backward integration.
%  MBACK(DELTAP,T0,L) calculates the mass flux M for given
%  pressure difference DELTAP.
%
%  [M FLOWSTRUCT] = MFLOW(DELTAP,T0,L) returns m and the solution
%  flowstruct.

% input check
p0=ps(T0);
if (deltap>p0)
  error( sprintf('deltap must be smaller than ps(T0) = %g',p0) );
  return
end

pe = p0 - deltap;
Te = T0 - jt(T0,p0,deltap);
% mguess basierend auf Darcy für Dampf;
%mguess = kappa*deltap./(L*nu(T0,1));
% mguess nach linearer Theorie
mguess = mlin(T0,p0,deltap,L);

% find an interval for m
% initially, mguess should be too small
%warning off;
pres=m2ph(mguess,Te,pe,L,p0);
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
  pres = m2ph(mb,Te,pe,L,p0);
end

% Display: none iter final notify
options =optimset(optimset('fzero'),'TolX',mguess/1000,'Display','none');
m = fzero(@m2ph,[mb/fac mb],options,Te,pe,L,p0);
flowstruct=flowback(m,Te,pe,L);
flowstruct.info.T0 = T0;
flowstruct.info.p0 = p0;
flowstruct.info.dp = deltap;

warning on;
lastwarn;

%-----------------------------------------------------------------------
function pres = m2ph(m,Te,pe,L,p0)
flowstruct=flowback(m,Te,pe,L);
pres = p0 - flowstruct.sol.p0;
% debug:
% sprintf('m: %g,   pres: %g',m,pres)
