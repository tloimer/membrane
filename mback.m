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
warning off;
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
warning backtrace;
flowstruct=flowback(m,Te,pe,L);
flowstruct.info.T0 = T0;
flowstruct.info.p0 = p0;
flowstruct.info.dp = deltap;

%lastwarn;

% add linear theory values
[mdot dedlin tmplin pelin] = mlin(T0,p0,deltap,L);

%flowstruct.lin = struct('m',
kkc = kappa/kappac(T0);
Telin = T0-jt(T0,p0)*deltap;
if ( kkc>1 ) %2ph
  a3 = tmplin;
  x3 = x(T0,a3);
  dfL = 0;
  T4 = [];
else % liq.film, tmplin  = film thickness
  if tmplin==0
    T4 = T0;
  else
    kldf=kl(T0)/tmplin;
    kmde=k(T0,0)/dedlin;
    T4 = ( Telin*kmde - T0*kldf )/( kmde-kldf );
  end
  dfL = -tmplin; % i put the z-position
  x3 = 0;
  a3 = 0;
end
flowstruct.lin = struct('m',mdot,'Te',Telin,'T4',T4,...
  'p6',pelin,'deL',dedlin,'dfL',dfL,'a3',a3,'x3',x3);

%-----------------------------------------------------------------------
function pres = m2ph(m,Te,pe,L,p0)
flowstruct=flowback(m,Te,pe,L);
pres = p0 - flowstruct.sol.p0;
% debug:
% sprintf('m: %g,   pres: %g',m,pres)
