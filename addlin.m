function flowstruct = addlin(fstruct)
%ADDLIN     Add solution from linear theory.
%  ADDLIN(FS) Add the solution from linear theory to the flowstruct FS.
%  ADDLIN returns a new flowstruct and leaves the old one unchanged.

%we write a new flowstruct
flowstruct = fstruct;

% set to flowstruct values
kapold = kappa;
kappa(flowstruct.info.kap);
%mlin calls xdot that depends on the 2ph-model
if ( ~strcmp(fmodel,flowstruct.info.ph) )
  oldmodel = fmodel;
  phmodel(flowstruct.info.ph);
end

% values given
m = flowstruct.info.m;
L = flowstruct.info.L;
T0 = flowstruct.info.T0;
p0 = flowstruct.info.p0;
deltap = flowstruct.info.dp;

i = length(flowstruct.flow) + 1;
%i = abs(flowstruct.sol.len) + 1;

% from the front
Telin = T0 - jt(T0,p0)*deltap;
[mdot dedlin tmp pelin] = mlin(T0,p0,deltap,L);
%disp('linear theory:');
%disp(sprintf('  forward: k/kc = %.3g; m = %.3gkg/m2s, DeltaT = %.3gK.',...
%disp(sprintf(' linear theory: k/kc = %.3g; m = %.3gkg/m2s, DeltaT = %.3gK.',...
%  kappa/kappac(T0),mdot,T0-Telin));

if (tmp>0) % 2ph-flow
  a3=tmp;
  x3=x(T0,a3);
  q1 = mdot*(1-xdot(T0,a3))*r(T0);
  zlin1 = [0 dedlin]*L;
  zlin2 = [dedlin 1]*L;
  Tlin1 = [T0 Telin];
  Tlin2 = [Telin Telin];
  plin1 = [p0 pelin];
  plin2 = [pelin p0-deltap];
  qlin1 = [q1 q1]; %/(mdot*r(T0));
  qlin2 = [0 0]; %/(mdot*r(T0));
  alin1 = [a3 x3/( (1-x3)/(v(Telin,pelin)*rho(Telin)) + x3 )];
  alin2 = [1 1];
  xlin1 = [x3 x3];
  xlin2 = [1 1];
else % liq. film
  dfd=tmp;
  q1 = mdot*r(T0);
  % kliq*(T0-T2)/dfd = k*(T2-Telin)/ded.  Achtung: dfd neg.
  kldf=kl(T0)/dfd;
  kmde=k(T0,0)/dedlin;
  T2 = ( Telin*kmde - T0*kldf )/( kmde-kldf );
  zlin1 = [dfd 0 dedlin]*L;
  zlin2 = [dedlin 1]*L;
  Tlin1 = [T0 T2 Telin];
  Tlin2 = [Telin Telin];
  plin1 = [p0 p0 pelin];
  plin2 = [pelin p0-deltap];
  qlin1 = [kldf*(T2-T0)/L q1 kmde*(T2-Telin)/L];% /q1;
  qlin2 = [0 0];% /q1;
  alin1 = [0 0 0];
  alin2 = [1 1];
  xlin1 = [0 0 0];
  xlin2 = [1 1];
end

flowstruct.flow(i).z = zlin1;
flowstruct.flow(i).T = Tlin1;
flowstruct.flow(i).p = plin1;
flowstruct.flow(i).q = qlin1;
flowstruct.flow(i).a = alin1;
flowstruct.flow(i).x = xlin1;
flowstruct.flow(i).color = 'y';
flowstruct.flow(i+1).z = zlin2;
flowstruct.flow(i+1).T = Tlin2;
flowstruct.flow(i+1).p = plin2;
flowstruct.flow(i+1).q = qlin2;
flowstruct.flow(i+1).a = alin2;
flowstruct.flow(i+1).x = xlin2;
flowstruct.flow(i+1).color = 'y';

%% backwards
%% find the end temp.
%Telin = fzero(...
%  inline('T0-(te+jt(te,p0-deltap)*deltap)','te','T0','p0','deltap'),...
%  Telin,optimset('fzero'),T0,p0,deltap);
%
%[mdot dedlin tmp pelin] = mlambda(Telin,p0,p0-deltap,deltap,L);
%% now nearly the same as above
%disp(sprintf('  back - : k/kc = %.3g: m = %.3gkg/m2s, DeltaT = %.3gK.',...
%  kappa/kappac(Telin),mdot,T0-Telin));
%
%if (tmp>0) % 2ph-flow
%  alin=tmp;
%  xlin=x(Telin,alin);
%  q1 = mdot*(1-lambda(Telin,alin))*r(Telin);
%  zlin = [0 dedlin 1]*L;
%  Tlin = [T0 Telin Telin];
%  plin = [p0 pelin p0-deltap];
%  qlin = [q1 q1 0];% /(mdot*r(Telin));
%  alin = [alin alin 0];
%  xlin = [xlin xlin 0];
%else
%  dfd=tmp;
%  q1 = mdot*r(Telin);
%  % kliq*(T0-T2)/dfd = k*(T2-Telin)/ded.  Achtung: dfd neg.
%  kldf=kl(Telin)/dfd;
%  kmde=k(Telin,0)/dedlin;
%  T2 = ( Telin*kmde - T0*kldf )/( kmde-kldf );
%  zlin = [dfd 0 dedlin 1]*L;
%  Tlin = [T0 T2 Telin Telin];
%  plin = [p0 p0 pelin p0-deltap];
%  qlin = [kldf*(T2-T0) q1 kmde*(T2-Telin) 0];% /q1;
%  alin = [0 0 0 1];
%  xlin = [0 0 0 1];
%end
%
%i=i+1;
%flowstruct.flow(i).z = zlin;
%flowstruct.flow(i).T = Tlin;
%flowstruct.flow(i).p = plin;
%flowstruct.flow(i).q = qlin;
%flowstruct.flow(i).a = alin;
%flowstruct.flow(i).x = xlin;
%flowstruct.flow(i).color = ':k+';
%
%%end backwards

% reset old values
kappa(kapold);
if ( exist('oldmodel','var') )
  phmodel(oldmodel);
end
