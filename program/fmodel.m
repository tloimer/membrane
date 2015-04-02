function f = fmodel(fname)
%FMODEL     Two-phase flow model.
%  FMODEL('plug') returns a struct F that contains functions to calculate
%  the effective properties of a two-phase mixture. At present, only a
%  homogeneous ('plug' or 'homogeneous') flow model is supported. All
%  functions in F can be vectorized.
%
%  F contains the fields
%    F.name                            Identifier of the flow model
%  Functions
%    F.xdot(a,vgas,vliq)               Vapor mass flow fraction
%    F.x(a,vgas,vliq)                  Vapor mass fraction
%    F.nu2ph(a,vgas,vliq,mugas,muliq)  Effective kinematic viscosity [m2/s]
%    F.k2ph(a,epsilon,km,kgas,kliq)    Effective thermal conductivity of the
%                                      mixture-filled membrane [W/mK]
%  Auxiliary functions, for convenience:
%    F.kmgas(epsilon,km,kgas)  Thermal conductivity, gas-filled membrane [W/mK]
%    F.kmliq(epsilon,km,kliq)  Thermal conductivity, liquid-filled mem. [W/mK]

if ~strcmp(fname,'plug') && ~strcmp(fname,'homogeneous')
  error('Only ''plug'' and ''homogeneous'' flow model is supported.');
end
% Struct constructor.
f = struct('name',fname,'xdot',[],'x',[],'nu2ph',[],'k2ph',[],'kmgas',[],...
  'kmliq',[]);

%  Equations see Loimer (2007).
% eq. (20) and above eq. (23)
f.xdot = @(a,vgas,vliq) a ./ ( a + (1-a).*vgas./vliq );
f.x = f.xdot;
% eq. (21)
f.nu2ph = @(a,vgas,vliq,mugas,muliq) ( a.*mugas + (1-a).*muliq ) ...
  ./ ( (1-a)./vliq + a./vgas );
% thermal conductivity model
keff = @(eps,km,kf) (1-eps).*km + eps.*kf;	% Parallel model
%keff = @(eps,km,kf) 1./((1-eps)./km + eps./kf);	% Series model
% eq. (22)
f.k2ph = @(a,epsilon,km,kgas,kliq) keff(epsilon,km,a.*kgas+(1-a).*kliq);

% Auxiliary functions
f.kmgas = keff;
f.kmliq = keff;
