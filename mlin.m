function [m,ded,a,pe] = mlin(T,p0,deltap,L)
%MLIN       Mass flux from linear theory [kg/m2s].
%  MLIN(T,P0,DELTAP,L) calculates the mass flux for initial state T and
%  p0 for the applied pressure drop DELTAP and membrane thickness L.
%
%  [M,DED,A,PE] = MLIN(...) gives, for initial 2ph-flow, the location of
%  the evaporation front de/L, the pressure at this point and the inital
%  vapor content a.
%
%  [M,DED,DFD,PE] = MLIN(...) gives the location of the evaporation front
%  de/L and the thickness of the liquid film df/L.

%Te = T-jt(T,p0).*deltap;

n = jt(T,p0).*dpsdT(T);
if ( n>=1 )
  n=1;
end

if (kappa>kappac(T))
  % supercritical permeability

  % determine a from  (1-lambda) nu'  k' / (nu k) = kc/k
  kkc = kappac(T)./(kappa*nu(T,0).*k(T,0));

  % calculate vapor content a
  a=fzero(@ares,[0 1],optimset('fzero'),T,kkc);

  % calculate m
  nn = 1-n+n.*nu(T,1)./nu(T,a);
  m = nn.*deltap.*kappa./(L.*nu(T,1));
  % de
  ded = n.*nu(T,1)./(nu(T,a).*nn);
  pe = p0 - m.*nu(T,a).*ded*L/kappa;

else
% subcritical permeability

  nn = 1-n+n.*nu(T,1)./nu(T,0);
  m = nn.*deltap.*kappa./(L.*nu(T,1));
  % de
  ded = n.*nu(T,1)./(nu(T,0).*nn);
  dfd = -n.*nu(T,1).*kl(T).*(kappac(T)/kappa-1)./(nn.*nu(T,0).*k(T,0));
  pe = p0 - m.*nu(T,0).*ded*L/kappa;
  % dirty outpt:
  a=dfd;

end
%%%%
%mb=m;
%dedb=ded;
%ab=a;
%peb=pe;
%
%% backward
%% find the end state
%Tee = fzero(...
%  inline('T-(te+jt(te,ps(te))*deltap)','te','T','deltap'),...
%  Te,optimset('fzero'),T,deltap);
%T=Tee;
%p0=ps(Tee);
%
%% now the same as above;
%n = jt(T,p0).*dpsdT(T);
%if ( n>=1 )
%  n=1;
%end
%
%if (kappa>kappac(T))
%  % supercritical permeability
%
%  % determine a from  (1-lambda) nu'  k' / (nu k) = kc/k
%  kkc = kappac(T)./(kappa*nu(T,0).*k(T,0));
%
%  % calculate vapor content a
%  a=fzero(@ares,[0 1],optimset('fzero'),T,kkc);
%
%  % calculate m
%  nn = 1-n+n.*nu(T,1)./nu(T,a);
%  m = nn.*deltap.*kappa./(L.*nu(T,1));
%  % de
%  ded = n.*nu(T,1)./(nu(T,a).*nn);
%  pe = p0 - m.*nu(T,a).*ded*L/kappa;
%
%else
%% subcritical permeability
%
%  nn = 1-n+n.*nu(T,1)./nu(T,0);
%  m = nn.*deltap.*kappa./(L.*nu(T,1));
%  % de
%  ded = n.*nu(T,1)./(nu(T,0).*nn);
%  dfd = -n.*nu(T,1).*kl(T).*(kappac(T)/kappa-1)./(nn.*nu(T,0).*k(T,0));
%  pe = p0 - m.*nu(T,0).*ded*L/kappa;
%  % dirty outpt:
%  a=dfd;
%
%end
%
%me=m;
%dede=ded;
%ae=a;
%pee=pe;
%
%m=mb;
%ded=dedb;
%a=ab;
%pe=peb;

function ar = ares(a,T,kkc)
ar = (1-xdot(T,a)) - nu(T,a).*k(T,a).*kkc;
% Schneider:
% xr = (1-x) - nu(T,x).*k(T,x).*kkc;
