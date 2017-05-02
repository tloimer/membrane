function f = fmodel(fname,cond)
%FMODEL     Thermal conductivity and two-phase flow model.
%  FMODEL returns a struct F that contains functions to calculate the
%  effective thermal conductivity of a porous medium, and the effective
%  properties of a two-phase mixture. The effective thermal conductivity is
%  computed according to the parallel model. The properties of the two-phase
%  flow are computed according to the homogeneous model.
%
%  All functions in F can be vectorized.
%
%  FMODEL('plug') or FMODEL('homogeneous') is equivalent to FMODEL().
%
%  FMODEL(COND), where COND is one of 'series', 'parallel', 'me1', 'me2' or
%  'EMT', returns the effective thermal conductivity of the fluid-filled
%  membrane according to the series, parallel, Maxwell-Eucken 1,
%  Maxwell-Eucken 2 or Effective Medium Theory model. These correspond to
%  eqs. (3) to (6), and eq. (8) in Carson (2005), respectively.
%
%  FMODEL(FNAME, COND), where FNAME is either 'plug' or 'homogeneous', is
%  equivalent to FMODEL(COND).
%
%  Currently, only a homogeneous two-phase flow model is supported.
%
%  J.K. Carson, S.J. Lovatt, D.J. Tanner, A.C. Cleland: Thermal conductivity
%  bounds for isotropic, porous materials. Int. J. Heat Mass Transfer 48,
%  2150-2158 (2005).
%
%  F contains the fields
%    F.name                            Identifier of the 2ph-flow model
%    F.cond                            Thermal conductivity model
%  Functions
%    F.xdot(a,vgas,vliq)               Vapor mass flow fraction
%    F.x(a,vgas,vliq)                  Vapor mass fraction
%    F.nu2ph(a,vgas,vliq,mugas,muliq)  Effective kinematic viscosity [m2/s]
%    F.k2ph(a,epsilon,km,kgas,kliq)    Effective thermal conductivity of the
%                                      mixture-filled membrane [W/mK]
%  Auxiliary functions, for convenience:
%    F.kmgas(epsilon,km,kgas)  Thermal conductivity, gas-filled membrane [W/mK]
%    F.kmliq(epsilon,km,kliq)  Thermal conductivity, liquid-filled mem. [W/mK]

if nargin < 2
  if nargin == 0
    fname = 'plug';
    cond = 'parallel';
  else % nargin == 1
    if any(strcmp(fname,{'plug', 'homogeneous'}))
      cond = 'parallel';
    else
      cond = fname;
      fname = 'plug';
    end
  end
elseif ~any(strcmp(fname,{'plug', 'homogeneous'}))
  error('Only ''plug'' and ''homogeneous'' flow model is supported.');
end

% Struct constructor.
f = struct('name',fname,'cond',cond,'xdot',[],'x',[],'nu2ph',[],...
    'k2ph',[],'kmgas',[],'kmliq',[]);

% eq. (20) and above eq. (23), Loimer (2007).
f.xdot = @(a,vgas,vliq) a ./ ( a + (1-a).*vgas./vliq );
f.x = f.xdot;
% eq. (21), Loimer (2007).
f.nu2ph = @(a,vgas,vliq,mugas,muliq) ( a.*mugas + (1-a).*muliq ) ...
  ./ ( (1-a)./vliq + a./vgas );
% thermal conductivity model
switch cond
    case 'series'   % Eq. (3) in Carson et al. (2005)
        keff = @(eps,km,kf) 1./((1-eps)./km + eps./kf);
    case 'parallel' % Eq. (4), ibid.
        keff = @(eps,km,kf) (1-eps).*km + eps.*kf;
    case 'me1'      % Eq. (5), ibid.
        % Upper bound for internal porous media, e.g., a track-etched
        % membrane, Vycor glass.
        keff = @(eps,km,kf) km .* (2*km + kf - 2*(km-kf).*eps)...
            ./ (2*km + kf + (km-kf).*eps);
    case 'me2'      % Eq. (6), ibid.
        % Lower bound for external porous media, e.g., compacted
        % particles like sand., ceramics?.
        keff = @(eps,km,kf) kf .* (2*kf + km - 2*(kf-km).*(1-eps))...
            ./ (2*kf + km + (kf-km).*(1-eps));
    case {'emt', 'EMT'} % Eq. (8), ibid.
        % Boundary between external and internal porous media
        keff = @(eps,km,kf) 0.25*((3*eps-1).*kf + (3*(1-eps)-1).*km...
            + sqrt(((3*eps-1).*kf + (3*(1-eps)-1).*km).^2 + 8*km.*kf));
	% To test for equality between eq. (7) and (8), ibid,
	% it is really (1/4)(...), not 1/(4*(...)), try mathematica:
	% Solve[(1 - n2) (k1 - ke)/(k1 + 2*ke) + n2 (k2 - ke)/(k2 + 2*ke) == 0
	%	&& k1 > 0 && k2 > 0 && n2 > 0 && n2 < 1, ke, Reals]
	% (1/4) ((3 n2 - 1) k2 + (3 (1 - n2) - 1) k1 +
	%	sqrt[((3 n2 - 1) k2 + (3 (1 - n2) - 1) k1)^2 + 8 k1 k2]) == %
	% Test the expression under the sqrt() separately.
    otherwise
        error('Effective thermal conductivity model %s not supported',...
            cond);
end

% eq. (22), Loimer (2007).
f.k2ph = @(a,epsilon,km,kgas,kliq) keff(epsilon,km,a.*kgas+(1-a).*kliq);

% Auxiliary functions
f.kmgas = keff;
f.kmliq = keff;
