function resplot(flowstruct,var)
%RESPLOT    Pretty print the flow result.
%  RESPLOT(FLOWSTRUCT) prints the temperature, pressure and mass
%  fraction distribution (for 2ph-flow) or heat flux (for fully
%  condensed flow) across the membrane. RESPLOT also calculates some
%  values from linear theory.
%
%  RESPLOT(FLOWSTRUCT,VAR) prints temperature, pressure and variable
%  VAR distribution, where VAR is one of 'x', 'a' or 'q'.
%
%  Calls MEMPLOT, PHMODEL, MLIN, KAPPA, KAPPAC.

m = flowstruct.info.m;
L = flowstruct.info.L;
ph = flowstruct.info.ph;
% T0, p0, and dp are not set if flowback produced the flowstruct
T0 = flowstruct.info.T0;
p0 = flowstruct.info.p0;
deltap = flowstruct.info.dp;

T0s = flowstruct.sol.T0;
Te = flowstruct.sol.Te;
a3 = flowstruct.sol.a3;
x3 = x(T0s,a3);
%p0s = flowstruct.sol.p0;
%pe = flowstruct.sol.pe;

if isempty(T0)
  if ~isempty(p0) | ~isempty(deltap)
    warning('strange flowstruct: empty T0, but p0 or deltap exist');
  end
  T0 = T0s;
  p0 = flowstruct.sol.p0;
  deltap = p0 - flowstruct.sol.pe;
end

flow = flowstruct.flow;

% nondimensionalize
for i=1:length(flow)
  flow(i).z = flow(i).z/L;
  flow(i).q = flow(i).q./(m*r(flow(i).T));
end

% set to flowstruct values
kapold = kappa;
kappa(flowstruct.info.kap);
% some function  calls, e.g. mlin, depend on the 2ph-model
if ( ~strcmp(fmodel,flowstruct.info.ph) )
  oldmodel = fmodel;
  phmodel(flowstruct.info.ph);
end

[mdot dedlin tmplin pelin] = mlin(T0,p0,deltap,L);

% some diagnostics
kkc = kappa/kappac(T0);
if ( kkc>1 ) %2ph
  tmplmsg = ', x_3 = %.3g';
  tmplin = x(T0s,tmplin);
%  tmplmsg = ', \\alpha_1 = %.3g';
else
  tmplmsg = ', d_f/L = %.3g';
  tmplin = -tmplin;
end
if ( a3>0 ) %2ph
  tmpsmsg = ', x_3 = %.3g';
  tmpsol = x3;
%  tmpsmsg = ', \\alpha_1 = %.3g';
%  tmpsol = a3;
else
  tmpsmsg = ', d_f/L = %.3g';
  tmpsol = -flowstruct.sol.df/L;
end

%begin the plots
subplot(3,1,1);
xl = floor(flowstruct.sol.df*10/L)/10;
memplot(flow,'T');

title(sprintf(...
  ['\\kappa/\\kappa_c = %.3g, \\Deltap = %.3gbar, \\DeltaT = %.3gK, ' ph ...
  ': m = %.3gkg/m^2s, d_e/L = %.3g' tmpsmsg ';\n'...
  '\\Deltap_{err} = %.3gPa, \\DeltaT_{err} = %.3gK, q_0/(mr) = %.3g.'],...
  kkc, deltap/1e5, T0-Te, m, flowstruct.sol.de/L, tmpsol,...
  flowstruct.sol.p0-p0,T0s-T0,...
  (flowstruct.sol.q3/(m*(1-xdot(T0s,a3))*r(T0s))-1) ));
ylabel('T [K]');
xlim([xl 1]);
ylim([floor(Te) ceil(T0s)]);

subplot(3,1,2);
memplot(flow,'p');
ylabel('p [Pa]');
xlim([xl 1]);

subplot(3,1,3);
if nargin==1
  if (a3>0)
    % 2ph flow
    var = 'x';
  else
    % complete condensation
    var = 'q';
  end
end

switch var
  case 'x'
    memplot(flow,'x');
    ylabel('x');
  case 'a'
    memplot(flow,'a');
    ylabel('\alpha');
%    ylim([0 1]);
  case 'q'
    memplot(flow,'q');
    ylabel('q/(m*r(T))');
%  ylim([0.95 1.05]);
  otherwise
    error('unsupported variable VAR given.');
end

xlabel(sprintf( ['z/L\nlinear theory: \\DeltaT = %.3gK, '...
  'm = %.3gkg/m^2s, d_e/L = %.3g' tmplmsg '.'], ...
  jt(T0,p0)*deltap, mdot, dedlin, tmplin ));
xlim([xl 1]);

% reset old values
kappa(kapold);
if ( exist('oldmodel','var') )
  phmodel(oldmodel);
end
