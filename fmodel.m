function f = fmodel(fname)
%FMODEL     Two-phase flow model.
%  F = FMODEL('plug') returns a struct F that contains functions to calculate
%  the effective properties of the two-phase mixture. At present, only a
%  homogeneous ('plug') flow model is supported. All functions in F can be
%  vectorized.
%
%  F returns its name
%    F.name                            Identifier of the flow model.
%  F provides the functions
%    F.xdot(A,VGAS,VLIQ)               Vapor mass flow fraction.
%    F.x(A,VGAS,VLIQ)                  Vapor mass fraction.
%    F.nu2ph(A,VGAS,VLIQ,MUGAS,MULIQ)  Effective kinematic viscosity [m2/s].
%    F.k2ph(A,EPSILON,KM,KGAS,KLIQ)    Effective thermal conductivity of the
%                                      mixture-filled membrane [W/mK].
%  Auxiliary functions, for convenience:
%    F.kmgas(EPSILON,KM,KGAS)  Thermal conductivity, gas-filled membrane [W/mK].
%    F.kmliq(EPSILON,KM,KLIQ)  Thermal conductivity, liquid-filled mem. [W/mK].

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
% eq. (22)
f.k2ph = @(a,epsilon,km,kgas,kliq) (1-epsilon).*km ...
  + epsilon.*(a.*kgas+(1-a).*kliq);

% Auxiliary functions
f.kmgas = @(epsilon,km,kgas) (1-epsilon).*km + epsilon.*kgas;
f.kmliq = @(epsilon,km,kliq) (1-epsilon).*km + epsilon.*kliq;