function k = kappa(kin,T)
%KAPPA      Permeability of the membrane [m2].
%  KAPPA without arguments returns the permeability of the membrane.
%
%  KAPPA(VAL) sets the permeability to VAL.
%
%  KAPPA(K,T) sets the permeability to K times the critical
%  permeabiltiy, K*KAPPAC(T).
%
%  Calls KAPPAC.

persistent kstore;

if nargin==0
  % return kappa
  if isempty(kstore)
    % default value
    kstore = 3.8259e-15; % 1.5*kappac
  end
elseif nargin==1
  % dimensional permeability given
  kstore = kin;
elseif nargin==2
  % relative permeability to kappac given
  kstore = kin*kappac(T);
end

k = kstore;

%k = 0.8 * 2.553e-15;
%k = 1.14 * 2.553e-15;
%k = 5 * 2.553e-15;
