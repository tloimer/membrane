function [m,ded,x,pe] = mlin(T,deltap,L)
%MLIN(T,deltap,L) Mass flow to be expected from linear theory.
%  The mass flow is calculated for initial state T and  ps(T), for the
%  applied pressure drop deltap and membrane thickness L.
%
%  [M,DED,X,PE] = MLIN(...) gives, for initial 2ph-flow, the location of
%  the evaporation front de/L, the pressure at this point and the inital
%  vapor content x.
%
%  [M,DED,DFD,PE] = MLIN(...) gives the location of the evaporation front
%  de/L and the thickness of the liquid film df/L.

n = jt(T,ps(T)).*dpsdT(T);
if ( n>=1 )
  n=1;
end

if (kappa>kappac(T))
  % supercritical permeability

  % determine x from  (1-x) nu'  k' / (nu k) = kc/k
  kkc = kappac(T)./(kappa*nu(T,0).*k(T,0));

  % calculate vapor content x
  x=fzero(@xres,[0 1],optimset('fzero'),T,kkc);

  % calculate m
  nn = 1-n+n.*nu(T,1)./nu(T,x);
  m = nn.*deltap.*kappa./(L.*nu(T,1));
  % de
  ded = n.*nu(T,1)./(nu(T,x).*nn);
  pe = ps(T) - m.*nu(T,x).*ded*L/kappa;

else
% subcritical permeability

  nn = 1-n+n.*nu(T,1)./nu(T,0);
  m = nn.*deltap.*kappa./(L.*nu(T,1));
  % de
  ded = n.*nu(T,1)./(nu(T,0).*nn);
  dfd = n.*nu(T,1).*kl(T).*(kappac(T)/kappa-1)./(nn.*nu(T,0).*k(T,0));
  pe = ps(T) - m.*nu(T,1).*ded*L/kappa;
  % dirty outpt:
  x=dfd;

end

function xr = xres(x,T,kkc)

xr = (1-x) - nu(T,x).*k(T,x).*kkc;
